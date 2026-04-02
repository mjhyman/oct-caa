function [subsampled, interpolated, interpolated_vs_rm] = ...
    planar_size_weighted_proximity(epvs, mask, ves, radius, p, mode)
% Compute SWP every ds-th voxel in 2D space (flip-flopped EPVS->samples)
% with textual progress reporting.
%
% INPUTS:
%   epvs   - logical 2D matrix (1 = EPVS voxel)
%   mask   - logical 2D matrix (1 = include in computation). vessels should
%               be set to 0 in this mask.
%   radius - search radius (voxels)
%   p      - exponent in the weighting function
%   ds     - downsample integer (sample every ds voxels)
%   parpool_flag (boolean) - True: use parallel computing
%
% OUTPUTS:
%   subsampled_volume    - sparse volume with EPVS density at ds-subsampled grid
%   interpolated_volume  - interpolated full-volume EPVS density

% 2D version - operates on single slices (e.g., 1500x2000)
dims = size(epvs);  % Now [rows, cols]

epvs = logical(epvs);
mask = logical(mask);
if ~isempty(ves), ves = logical(ves); end

%% 1. Build points and weights
fprintf('Building points and weights...\n')
if strcmpi(mode, 'centroids')
    cc = bwconncomp(epvs, 8);  % 8-connectivity for 2D (was 26 for 3D)
    points = zeros(cc.NumObjects, 2, 'single');  % 2 cols, not 3
    weights = zeros(cc.NumObjects, 1, 'single');
    for i = 1:cc.NumObjects
        [x, y] = ind2sub(dims, cc.PixelIdxList{i});  % No z
        points(i, :) = [mean(x), mean(y)];
        weights(i) = numel(x);
    end
    clear cc;
    is_voxel_mode = false;
else
    epi_inds = find(epvs);
    [ex, ey] = ind2sub(dims, epi_inds);  % No ez
    points = single([ex, ey]);
    weights = 1;
    is_voxel_mode = true;
    clear ex ey epi_inds;
end
clear epvs;

%% 2. Define sampling grid
fprintf('Defining sampling grid\n')
gv = {1:4:dims(1), 1:4:dims(2)};  % No 3rd dimension
[xg, yg] = ndgrid(gv{1}, gv{2});  % 2D grid
sample_coords = single([xg(:), yg(:)]);
clear xg yg;

fprintf('Converting sampling grid to valid indices\n')
sample_inds = sub2ind(dims, sample_coords(:,1), sample_coords(:,2));  % 2D
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

chunk_size = 5000;
fprintf('Processing %d grid samples...\n', num_valid);

for start_idx = 1:chunk_size:num_valid
    end_idx = min(start_idx + chunk_size - 1, num_valid);
    chunk_c = valid_coords(start_idx:end_idx, :);

    fprintf('\tPerforming range search\n')
    [IdxC, DC] = rangesearch(tree, chunk_c, radius);

    IdxC = cellfun(@(x) x(:), IdxC, 'UniformOutput', false);
    DC   = cellfun(@(x) x(:), DC,   'UniformOutput', false);

    num_n = cellfun(@numel, IdxC);
    empty_idx = find(num_n == 0);
    if ~isempty(empty_idx)
        fprintf('\tKNN search for missing values\n')
        [fIdx, fD] = knnsearch(tree, chunk_c(empty_idx, :), 'K', 1);
        for i = 1:numel(empty_idx)
            ii = empty_idx(i);
            IdxC{ii} = fIdx(i);
            DC{ii}   = fD(i);
        end
    end

    fprintf('\tCalculating and summing contributions per voxel\n')
    vals = zeros(numel(IdxC), 1, 'single');
    for i = 1:numel(IdxC)
        d = DC{i};
        d(d < 0.5) = 0.5;
        if is_voxel_mode
            vals(i) = sum(1 ./ (d .^ p));
        else
            vals(i) = sum(weights(IdxC{i}) ./ (d .^ p));
        end
    end

    data(start_idx:end_idx) = vals;
    fprintf('\tProcessed chunk indices [%i - %i] of %i\n', ...
        start_idx, start_idx + chunk_size, num_valid)
    clear IdxC DC vals d;
end
clear tree points weights;

%% 4. Interpolation
fprintf('Building subsampled volume\n')
subsampled = zeros(dims, 'single');
v_inds = sub2ind(dims, valid_coords(:,1), valid_coords(:,2));  % 2D
subsampled(v_inds) = data;
clear valid_coords data v_inds;

fprintf('Preparing interpolant\n')
Sgrid = subsampled(gv{1}, gv{2});  % 2D indexing
F = griddedInterpolant({single(gv{1}), single(gv{2})}, Sgrid, 'linear', 'nearest');
clear Sgrid;

fprintf('Interpolating...\n');
[X2D, Y2D] = ndgrid(single(1:dims(1)), single(1:dims(2)));
interpolated = F(X2D, Y2D);  % Single call, no loop needed

interpolated_vs_rm = interpolated;
if ~isempty(ves)
    interpolated_vs_rm(ves) = 0;  % Apply mask directly (2D)
end

fprintf('Finished.\n');
end