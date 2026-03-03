function [subsampled_volume, interpolated_volume] =...
    epvs_density_variable_p(epvs, mask, radius, p)
% Memory-optimized EPVS density at 4x subsampled grid + interpolation
% Vectorized neighbor accumulation (Option A) + micro-optimizations.
% Date: February 24, 2026

%% Parameters

% memory available
bytes = get_approx_free_bytes();
if isempty(bytes)
    chunk_size = 100000; % fallback
else
    chunk_size = min(max(2000, floor(bytes/200/8)), 20000);
end

% Size of volume
dims = size(epvs);

%% EPVS components -> use centroid per component (weight = size)
cc = bwconncomp(epvs, 26);
nComp = cc.NumObjects;

if nComp == 0
    warning('No EPVS voxels found.');
    subsampled_volume = sparse([]);
    interpolated_volume = single([]);
    return;
end

centroids = zeros(nComp,3,'single');
weights = zeros(nComp,1,'single');

for i = 1:nComp
    idxs = cc.PixelIdxList{i};
    [x,y,z] = ind2sub(dims, idxs);
    cent = round(mean([x(:), y(:), z(:)],1));
    centroids(i,:) = single(cent);
    weights(i) = single(numel(idxs));
end

%% Define sampling grid: every 4th voxel in 3D
[xg, yg, zg] = ndgrid(1:4:dims(1), 1:4:dims(2), 1:4:dims(3));
grid_coords = [xg(:), yg(:), zg(:)];
sample_inds = sub2ind(dims, grid_coords(:,1), grid_coords(:,2), grid_coords(:,3));
is_valid_sample = mask(sample_inds) & ~epvs(sample_inds);
valid_coords = grid_coords(is_valid_sample, :);
clear grid_coords;

%%% Ensure sampled coordinates are unique
flat_valid = sub2ind(dims, valid_coords(:,1), valid_coords(:,2), valid_coords(:,3));
[~, uniq_map, ~] = unique(flat_valid, 'stable');
valid_coords = valid_coords(uniq_map, :);
num_valid = size(valid_coords,1);

rows = valid_coords(:,1);
cols = valid_coords(:,2);
deps = valid_coords(:,3);
data = zeros(num_valid,1,'single');

%% Build KD-tree on centroids in single precision
% The centroids are used for the weights and distances, rather than using
% all voxels within all EPVS. This accelerates the computation and provides
% a similar metric.
tree = KDTreeSearcher(single(centroids));

%% Precompute cast values for speed
radius = single(radius);
p = single(p);

%% Chunking
chunk_start_indices = 1:chunk_size:num_valid;
nChunks = numel(chunk_start_indices);

fprintf('Computing EPVS density at every 4th voxel (%d sampled)...\n', num_valid);
tic;

for chunk_idx = 1:nChunks
    start_idx = chunk_start_indices(chunk_idx);
    stop_idx = min(start_idx + chunk_size - 1, num_valid);
    chunk_inds = start_idx:stop_idx;
    chunk_coords = single(valid_coords(chunk_inds,:)); 

    % rangesearch once per chunk
    [IdxChunk, DChunk] = rangesearch(tree, chunk_coords, radius);

    % Ensure each cell is a column vector (fast in MATLAB)
    nonemptyIdx = find(~cellfun('isempty', DChunk));
    for k = nonemptyIdx.'              % iterate over row vector for slight speed
        di = DChunk{k};
        if size(di,1) == 1 && size(di,2) > 1
            DChunk{k} = di(:);
        end
        ii = IdxChunk{k};
        if size(ii,1) == 1 && size(ii,2) > 1
            IdxChunk{k} = ii(:);
        end
    end
    
    % Vectorized accumulation
    len = cellfun(@numel, IdxChunk);
    totalNeigh = sum(len);
    
    if totalNeigh == 0
        vals = zeros(numel(chunk_inds),1,'single');
    else
        % Concatenate neighbor indices and distances into column vectors
        allIdx = vertcat(IdxChunk{:});    % (totalNeigh x 1)
        allD   = vertcat(DChunk{:});     % (totalNeigh x 1)
    
        % Sanity check
        if numel(allIdx) ~= totalNeigh || numel(allD) ~= totalNeigh
            error('Concatenation mismatch after transpose fix.');
        end
    
        w_all = weights(allIdx);
        d_all = single(allD);
        contrib = w_all ./ (d_all .^ p);
    
        % Build group vector via cumulative lengths
        cumlen = [0; cumsum(len(:))];
        group = zeros(totalNeigh,1,'uint32');
        for q = 1:numel(len)
            if len(q) > 0
                r1 = cumlen(q) + 1;
                r2 = cumlen(q+1);
                group(r1:r2) = q;
            end
        end
    
        vals = single(splitapply(@sum, contrib, group));
    end
    
    data(chunk_inds) = vals;

    clear IdxChunk DChunk allIdx allD contrib group len vals;

    % Report status to console
    fprintf('Chunk %d/%d complete (%.1f%%)\n', chunk_idx, nChunks,...
             100 * chunk_idx / nChunks);
end

fprintf('Subsampled computation done in %.2f seconds.\n', toc);

%% Build subsampled outputs
flat_idx = sub2ind(dims, rows, cols, deps);
subsampled_volume = zeros(dims, 'single');
subsampled_volume(flat_idx) = data;

%% Prepare coarse grid for interpolation
coarseX = single(1:4:dims(1));
coarseY = single(1:4:dims(2));
coarseZ = single(1:4:dims(3));

Sgrid = subsampled_volume(1:4:end, 1:4:end, 1:4:end);

F = griddedInterpolant({coarseX, coarseY, coarseZ}, Sgrid, 'linear', 'nearest');

interpolated_volume = zeros(dims, 'single');
[X2D, Y2D] = ndgrid(single(1:dims(1)), single(1:dims(2)));
for z = 1:dims(3)
    Z2D = single(z) * ones(size(X2D),'single');
    slice = F(X2D, Y2D, Z2D);
    interpolated_volume(:,:,z) = single(slice);
    if mod(z, max(1,ceil(dims(3)/10))) == 0
        fprintf('Interpolation slice %d/%d\n', z, dims(3));
    end
end

fprintf('Interpolation complete.\n');

end
