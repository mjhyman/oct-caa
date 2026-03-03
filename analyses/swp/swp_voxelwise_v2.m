function [subsampled, interpolated, interpolated_vs_rm] = ...
    swp_voxelwise_v2(epvs, mask, ves, radius, p, mode)
% Optimized for 1500x2000x150 volumes (~1.8GB per single-precision volume)
dims = size(epvs);
% Immediate conversion to logical to save RAM
epvs = logical(epvs);
mask = logical(mask);
if ~isempty(ves), ves = logical(ves); end

%% 1. Build points and weights
fprintf('Building points and weights...\n')
if strcmpi(mode, 'centroids')
    cc = bwconncomp(epvs, 26);
    points = zeros(cc.NumObjects, 3, 'single');
    weights = zeros(cc.NumObjects, 1, 'single');
    for i = 1:cc.NumObjects
        [x, y, z] = ind2sub(dims, cc.PixelIdxList{i});
        points(i, :) = [mean(x), mean(y), mean(z)];
        weights(i) = numel(x); 
    end
    clear cc;
    is_voxel_mode = false;
else
    epi_inds = find(epvs);
    [ex, ey, ez] = ind2sub(dims, epi_inds);
    points = single([ex, ey, ez]);
    weights = 1; 
    is_voxel_mode = true;
    clear ex ey ez epi_inds;
end
clear epvs; 

%% 2. Define sampling grid
fprintf('Defining sampling grid\n')
gv = {1:4:dims(1), 1:4:dims(2), 1:4:dims(3)};
[xg, yg, zg] = ndgrid(gv{1}, gv{2}, gv{3});
sample_coords = single([xg(:), yg(:), zg(:)]);
clear xg yg zg;

fprintf('Converting sampling grid to valid indices\n')
sample_inds = sub2ind(dims, sample_coords(:,1), sample_coords(:,2), sample_coords(:,3));
is_valid = mask(sample_inds); 
valid_coords = sample_coords(is_valid, :);
clear sample_coords sample_inds is_valid;
clear mask; 

%% 3. KD-Tree Search (Chunked)
fprintf('Performing KD Tree Searcher\n')
tree = KDTreeSearcher(points);
num_valid = size(valid_coords, 1);
data = zeros(num_valid, 1, 'single');
radius = single(radius);
p = single(p);

% Use a slightly smaller chunk size to keep the cell arrays manageable
chunk_size = 5000; 
fprintf('Processing %d grid samples...\n', num_valid);

for start_idx = 1:chunk_size:num_valid
    end_idx = min(start_idx + chunk_size - 1, num_valid);
    chunk_c = valid_coords(start_idx:end_idx, :);
    
    fprintf('\tPerforming range search\n')
    [IdxC, DC] = rangesearch(tree, chunk_c, radius);
    
    % Force every cell to be column vector to prevent dimension mismatch
    IdxC = cellfun(@(x) x(:), IdxC, 'UniformOutput', false);
    DC = cellfun(@(x) x(:), DC, 'UniformOutput', false);
    
    % KNN Fallback for empty regions
    num_n = cellfun(@numel, IdxC);
    empty_idx = find(num_n == 0);
    if ~isempty(empty_idx)
        fprintf('\tKNN search for missing values\n')
        [fIdx, fD] = knnsearch(tree, chunk_c(empty_idx, :), 'K', 1);
        for i = 1:numel(empty_idx)
            ii = empty_idx(i);
            IdxC{ii} = fIdx(i); 
            DC{ii} = fD(i); 
        end
    end
    
    % NEW FIX: Calculate contributions per-voxel to avoid massive 'allD' array
    fprintf('\tCalculating and summing contributions per voxel\n')
    vals = zeros(numel(IdxC), 1, 'single');
    for i = 1:numel(IdxC)
        d = DC{i};
        % Apply physical limit to prevent infinity
        d(d < 0.5) = 0.5; 
        
        if is_voxel_mode
            % Fast path for weights = 1
            vals(i) = sum(1 ./ (d.^p));
        else
            % Multiply each neighbor distance by its component weight
            vals(i) = sum(weights(IdxC{i}) ./ (d.^p));
        end
    end
    
    data(start_idx:end_idx) = vals;
    chunk_end_idx = start_idx + chunk_size;
    fprintf('\tProcessed chunk indices [%i - %i] of %i\n',...
                start_idx, chunk_end_idx, num_valid)
    clear IdxC DC vals d; 
end
clear tree points weights;

%% 4. Interpolation
fprintf('Building subsampled volume\n')
subsampled = zeros(dims, 'single');
v_inds = sub2ind(dims, valid_coords(:,1), valid_coords(:,2), valid_coords(:,3));
subsampled(v_inds) = data;
clear valid_coords data v_inds;

fprintf('Preparing interpolant\n')
Sgrid = subsampled(gv{1}, gv{2}, gv{3});
F = griddedInterpolant({single(gv{1}), single(gv{2}), single(gv{3})}, Sgrid, 'linear', 'nearest');
clear Sgrid; 

% Force garbage collection before the big 3.6GB allocation
java.lang.System.gc(); 

fprintf('Allocating final volumes (~3.6 GB)...\n');
interpolated = zeros(dims, 'single');
interpolated_vs_rm = zeros(dims, 'single');
[X2D, Y2D] = ndgrid(single(1:dims(1)), single(1:dims(2)));

fprintf('Interpolating slices...\n');
for z = 1:dims(3)
    z_coords = single(z) * ones(dims(1), dims(2), 'single');
    slice = F(X2D, Y2D, z_coords);
    
    interpolated(:,:,z) = slice;
    if ~isempty(ves)
        % Apply vessel mask slice-by-slice to keep memory usage low
        slice(ves(:,:,z)) = 0; 
    end
    interpolated_vs_rm(:,:,z) = slice;
    
    if mod(z, 25) == 0, fprintf('\tSlice %d/%d done\n', z, dims(3)); end
end
fprintf('Finished.\n');
end