function [subsampled_volume, interpolated_volume] =...
    epvs_density(epvs, mask, radius)
% Computes EPVS density at every 4th voxel in 3D space and interpolates full map
%
% INPUTS:
%   epvs   - logical 3D matrix (1 = EPVS voxel)
%   mask   - logical 3D matrix (1 = include in computation)
%   radius - search radius (voxels)
%
% OUTPUTS:
%   subsampled_volume    - sparse volume with EPVS density at 4x subsampled grid
%   interpolated_volume  - interpolated full-volume EPVS density

%% Parameters
epsilon = 1e-3;
p = 2;
chunk_size = 10000;
dims = size(epvs);

%% Get EPVS voxel coordinates and weights
num_epvs_voxels = sum(cellfun(@numel, CC.PixelIdxList));
all_EPVS_coords = zeros(num_epvs_voxels, 3, 'uint16'); % use uint16 if dims < 2^16
all_EPVS_weights = zeros(num_epvs_voxels, 1, 'single');

ptr = 1;

for i = 1:CC.NumObjects
    idxs = CC.PixelIdxList{i};
    sz = numel(idxs);
    m = length(idxs);
    [x, y, z] = ind2sub(dims, idxs);
    all_EPVS_coords(ptr:ptr+m-1, :) = [x(:), y(:), z(:)];
    all_EPVS_weights(ptr:ptr+m-1) = sz;
    ptr = ptr + m;
end

%% Define sampling grid: every 4th voxel in 3D space
[xg, yg, zg] = ndgrid(1:4:dims(1), 1:4:dims(2), 1:4:dims(3));
grid_coords = [xg(:), yg(:), zg(:)];

% Only keep sampled voxels that are in the mask and not in EPVS
valid_idx = mask(sub2ind(dims, grid_coords(:,1), grid_coords(:,2), grid_coords(:,3))) & ...
            ~epvs(sub2ind(dims, grid_coords(:,1), grid_coords(:,2), grid_coords(:,3)));
valid_coords = grid_coords(valid_idx, :);
num_valid = size(valid_coords, 1);

% Initialize volumes
subsampled_volume = zeros(dims);
interpolated_volume = zeros(dims);

if isempty(all_EPVS_coords) || isempty(valid_coords)
    warning('No EPVS or valid sampled voxels found.');
    return;
end

%% Build KD-tree
tree = KDTreeSearcher(all_EPVS_coords);

%% Chunking
chunk_start_indices = 1:chunk_size:num_valid;
nChunks = numel(chunk_start_indices);

fprintf('Computing EPVS density at every 4th voxel (%d sampled)...\n', ...
        num_valid);
tic;

% Compute density at sampled voxels
chunk_results = cell(nChunks,1);

parfor chunk_idx = 1:nChunks
    start_idx = chunk_start_indices(chunk_idx);
    stop_idx = min(start_idx + chunk_size - 1, num_valid);

    chunk_inds = start_idx:stop_idx;
    chunk_coords = valid_coords(chunk_inds, :); % OK to read

    [Idx, D] = rangesearch(tree, chunk_coords, radius);

    % Store results as [x, y, z, value]
    vals = zeros(numel(chunk_inds),4);
    for i = 1:numel(chunk_inds)
        d = D{i};
        w = all_EPVS_weights(Idx{i});
        if ~isempty(d)
            w = w(:);
            d = d(:);
            val = sum(w ./ (d + epsilon).^p);
        else
            val = 0;
        end
        vals(i,:) = [chunk_coords(i,:) val];
    end
    chunk_results{chunk_idx} = vals;
end

% After parfor, update subsampled_volume
for chunk_idx = 1:nChunks
    vals = chunk_results{chunk_idx};
    inds = sub2ind(dims, vals(:,1), vals(:,2), vals(:,3));
    subsampled_volume(inds) = vals(:,4);
end

fprintf('Subsampled computation done in %.2f seconds.\n', toc);

%% Interpolation to full grid
fprintf('Interpolating to full volume...\n');
[xq, yq, zq] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));

% Interpolate using griddedInterpolant
F = griddedInterpolant(xg, yg, zg,...
    subsampled_volume(1:4:end, 1:4:end, 1:4:end), 'linear', 'nearest');
interpolated_volume = F(xq, yq, zq);

fprintf('Interpolation complete.\n');
end
