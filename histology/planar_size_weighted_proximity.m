function [subsampled_volume, interpolated_volume] =...
    planar_size_weighted_proximity(epvs, mask, radius, p, ds, parpool_flag)
% Compute SWP every 4th voxel in 2D space - interpolates full map
%
% INPUTS:
%   epvs   - logical 2D matrix (1 = EPVS voxel)
%   mask   - logical 2D matrix (1 = include in computation). vessels should
%               be set to 0 in this mask.
%   radius - search radius (voxels)
%   p      - exponent in the weighting function
%   ds     - downsample integer
%   parpool_flag (boolean) - True: use parallel computing
%
% OUTPUTS:
%   subsampled_volume    - sparse volume with EPVS density at 4x subsampled grid
%   interpolated_volume  - interpolated full-volume EPVS density

%% Parameters

% Value in denominator to prevent dividing by zero
epsilon = 1e-3;

try
    user_mem = memory;
    chunk_size = min(max(20000,floor(user_mem.MaxPossibleArrayBytes/100/8)),...
                    200000);
catch
    chunk_size = 20000;
end

% Retrieve dimensions of the EPVS matrix. These are the dimensions 
dims = size(epvs);

%% EPVS voxel coordinates and weights
epvs_idx = find(epvs);
n_epvs_vox = numel(epvs_idx);

if n_epvs_vox == 0
    warning('No EPVS voxels found.');
    subsampled_volume = sparse([]);
    interpolated_volume = single([]);
    return;
end

% Initialize arrays for improved speed
all_EPVS_coords = zeros(n_epvs_vox,2,'single'); 
all_EPVS_weights = zeros(n_epvs_vox,1,'single');

%%% Find EPVS coordinates and weights of each disjoint group
% Find EPVS connected components
cc = bwconncomp(epvs, 26);
% Iterate over each disjoint EPVS
curr_idx = 1;
for i = 1:cc.NumObjects
    idxs = cc.PixelIdxList{i};
    [x, y] = ind2sub(dims, idxs);
    n = numel(idxs);
    all_EPVS_coords(curr_idx:curr_idx+n-1, :) = single([x(:), y(:)]);
    all_EPVS_weights(curr_idx:curr_idx+n-1) = single(n);
    curr_idx = curr_idx + n;
end

%% Define sampling grid: every 4th voxel in 3D

% Initilize grid at every fourth voxel
[xg, yg] = ndgrid(1:ds:dims(1), 1:ds:dims(2)); 
% Retrieve grid coordinates of every fourth voxel
grid_coords = [xg(:), yg(:)];
% Convert subscripts to indices
sample_inds = sub2ind(dims, grid_coords(:,1), grid_coords(:,2));
% Find indices of parenchyma tissue that are NOT EPVS
valid_idx = mask(sample_inds) & ~epvs(sample_inds);
% Find the respective grid coordinates
valid_coords = grid_coords(valid_idx, :); 
% Find the number of valid grid coordinates
num_valid = size(valid_coords, 1);

% If no valid parenchyma, then return warning
if num_valid == 0
    warning('No valid sampled voxels found.');
    subsampled_volume = sparse([]);
    interpolated_volume = single([]);
    return;
end

% Initialize rows, columns, and data matrix
rows = valid_coords(:,1);
cols = valid_coords(:,2);
data = zeros(num_valid,1,'single');

%% Build KD-tree (KDTreeSearcher needs double)
tree = KDTreeSearcher(double(all_EPVS_coords));

%% Chunking
chunk_start_indices = 1:chunk_size:num_valid;
nChunks = numel(chunk_start_indices);

fprintf('Computing EPVS density at every 4th voxel (%d sampled)...\n',...
        num_valid);
tic;

for chunk_idx = 1:nChunks
    start_idx = chunk_start_indices(chunk_idx);
    stop_idx = min(start_idx + chunk_size - 1, num_valid);
    chunk_inds = start_idx:stop_idx;
    chunk_coords = double(valid_coords(chunk_inds,:)); 

    [Idx, D] = rangesearch(tree, chunk_coords, radius);

    vals = zeros(numel(chunk_inds), 1, 'single');
    
    % Parallel loop for density computation
    if parpool_flag
        parfor i = 1:numel(chunk_inds)
            idx = Idx{i};
            d = D{i};
            if isempty(idx)
                vals(i) = 0;
            else
                w = all_EPVS_weights(idx);
                % 'omitnan' safe for rare cases
                vals(i) = sum(w ./ (single(d) + epsilon).^p, "all",'omitnan'); 
            end
        end
    else
        for i = 1:numel(chunk_inds)
            idx = Idx{i};
            d = D{i};
            if isempty(idx)
                vals(i) = 0;
            else
                w = all_EPVS_weights(idx);
                % 'omitnan' safe for rare cases
                vals(i) = sum(w ./ (single(d) + epsilon).^p, "all",'omitnan'); 
            end
        end
    end

    data(chunk_inds) = vals;

    clear Idx D;

    % Print chunk # to console
    fprintf('Chunk %d/%d complete (%.1f%%)\n', chunk_idx, nChunks,...
            100 * chunk_idx / nChunks);
end

fprintf('Subsampled computation done in %.2f seconds.\n', toc);

%% Build sparse 1D vector, then convert to full and reshape to 3D
fprintf('Starting Interpolation\n');

flat_idx = sub2ind(dims, rows, cols);
S = sparse(flat_idx, 1, double(data), prod(dims), 1);
subsampled_volume = single(reshape(full(S), dims));

% Prepare coarsely sampled grid for interpolation
Sgrid = single(subsampled_volume(1:ds:end,1:ds:end));

F = griddedInterpolant(xg, yg, Sgrid, 'linear', 'nearest');
[xq, yq] = ndgrid(1:dims(1), 1:dims(2));
interpolated_volume = single(F(xq, yq));

fprintf('Interpolation complete.\n');
end