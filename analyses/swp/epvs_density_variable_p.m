function [subsampled_volume, interpolated_volume] = epvs_density_variable_p(epvs, mask, radius, p)
% Computes EPVS density at every 4th voxel in 3D space and interpolates full map
%
% INPUTS:
%   epvs   - logical 3D matrix (1 = EPVS voxel)
%   mask   - logical 3D matrix (1 = include in computation)
%   radius - search radius (voxels)
%   p      - exponent in the weighting function
%
% OUTPUTS:
%   subsampled_volume    - sparse volume with EPVS density at 4x subsampled grid
%   interpolated_volume  - interpolated full-volume EPVS density

%% Parameters
epsilon = 1e-3;
chunk_size = 10000;
dims = size(epvs);

%% Efficient: Preallocate for EPVS coords and weights
% Get all EPVS voxel indices
epvs_idx = find(epvs);
n_epvs_vox = numel(epvs_idx);

if n_epvs_vox == 0
    warning('No EPVS voxels found.');
    subsampled_volume = sparse([]);
    interpolated_volume = single([]);
    return;
end

% Preallocate for efficiency
all_EPVS_coords = zeros(n_epvs_vox, 3, 'uint32'); % Use uint32 for larger images
all_EPVS_weights = zeros(n_epvs_vox, 1, 'uint32');

cc = bwconncomp(epvs, 26);

curr_idx = 1;
for i = 1:cc.NumObjects
    idxs = cc.PixelIdxList{i};
    sz = numel(idxs);
    [x, y, z] = ind2sub(dims, idxs);
    n = length(idxs);
    all_EPVS_coords(curr_idx:curr_idx+n-1, :) = [x(:), y(:), z(:)];
    all_EPVS_weights(curr_idx:curr_idx+n-1) = sz;
    curr_idx = curr_idx + n;
end

%% Define sampling grid: every 4th voxel in 3D space
[xg, yg, zg] = ndgrid(1:4:dims(1), 1:4:dims(2), 1:4:dims(3));
grid_coords = [xg(:), yg(:), zg(:)];

% Ensure type consistency for sub2ind; must be double
sample_inds = sub2ind(dims, double(grid_coords(:,1)), double(grid_coords(:,2)), double(grid_coords(:,3)));
valid_idx = mask(sample_inds) & ~epvs(sample_inds);
valid_coords = grid_coords(valid_idx, :);
num_valid = size(valid_coords, 1);

if num_valid == 0
    warning('No valid sampled voxels found.');
    subsampled_volume = sparse([]);
    interpolated_volume = single([]);
    return;
end

% Instead of allocating a full-size array, store sampled grid in sparse format
data = zeros(num_valid, 1, 'single');
rows = double(valid_coords(:,1)); % Must be double for sub2ind
cols = double(valid_coords(:,2));
deps = double(valid_coords(:,3));

%% Build KD-tree (KDTreeSearcher needs double, so cast)
tree = KDTreeSearcher(double(all_EPVS_coords));

%% Chunking
chunk_start_indices = 1:chunk_size:num_valid;
nChunks = numel(chunk_start_indices);

fprintf('Computing EPVS density at every 4th voxel (%d sampled)...\n', num_valid);
tic;

% Compute density at sampled voxels (in chunks)
for chunk_idx = 1:nChunks
    start_idx = chunk_start_indices(chunk_idx);
    stop_idx = min(start_idx + chunk_size - 1, num_valid);

    chunk_inds = start_idx:stop_idx;
    chunk_coords = double(valid_coords(chunk_inds, :)); % KDTree requires double
    
    % rangesearch returns cell arrays Idx and D
    [Idx, D] = rangesearch(tree, chunk_coords, radius);

    for i = 1:numel(Idx) % You had parfor; remove unless you have Parallel Toolbox.
        if isempty(Idx{i}), data(chunk_inds(i)) = 0; continue; end
        d = D{i};
        w = double(all_EPVS_weights(Idx{i}));
        d = d(:);
        w = w(:);
        val = sum(w ./ (d + epsilon).^p); % d+epsilon avoids divide by zero
        data(chunk_inds(i)) = single(val);
    end
    fprintf('Chunk %d/%d complete (%.1f%%)\n', ...
            chunk_idx, nChunks, 100 * chunk_idx / nChunks);
end

fprintf('Subsampled computation done in %.2f seconds.\n', toc);

%% Build sparse 1D vector, then convert to full and reshape to 3D
S = sparse(sub2ind(dims, rows, cols, deps), 1, double(data), prod(dims), 1);
subsampled_volume = single(reshape(full(S), dims));

% Prepare coarsely sampled grid for interpolation
[xg, yg, zg] = ndgrid(1:4:dims(1), 1:4:dims(2), 1:4:dims(3));
Sgrid = single(subsampled_volume(1:4:end, 1:4:end, 1:4:end));

% Build interpolant
F = griddedInterpolant(xg, yg, zg, Sgrid, 'linear', 'nearest');
[xq, yq, zq] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
interpolated_volume = F(xq, yq, zq);
interpolated_volume = single(interpolated_volume);
fprintf('Interpolation complete.\n');
end
