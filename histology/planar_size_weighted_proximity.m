function [subsampled_volume, interpolated_volume] = ...
planar_size_weighted_proximity(epvs, mask, radius, p, ds, parpool_flag)
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

%% Parameters
epsilon = 1e-3; % small value to prevent division by zero

% Tuning: memory threshold for dense per-worker local_data (bytes)
MAX_DENSE_BYTES_PER_WORKER = 200e6;

% Tuning: target EPVS per block (helps choose a sensible block count)
TARGET_EPVS_PER_BLOCK = 5000;

% Progress reporting frequency inside each block (report every this many EPVS)
N_EPVS_REPORT = 500;

%% Basic checks and dims
if islogical(epvs)
    epvs = logical(epvs);
end
if islogical(mask)
    mask = logical(mask);
end

dims = size(epvs);

% Locate EPVS voxels
epvs_idx = find(epvs);
n_epvs_vox = numel(epvs_idx);

if n_epvs_vox == 0
    warning('No EPVS voxels found.');
    subsampled_volume = sparse([]);
    interpolated_volume = single([]);
    return;
end

%% Build EPVS coords and weights
all_EPVS_coords = zeros(n_epvs_vox,2,'double');
all_EPVS_weights = zeros(n_epvs_vox,1,'single');

cc = bwconncomp(epvs, 8); % 8-connectivity for 2D
curr_idx = 1;
for k = 1:cc.NumObjects
    idxs = cc.PixelIdxList{k};
    n = numel(idxs);
    [x, y] = ind2sub(dims, idxs);
    all_EPVS_coords(curr_idx:curr_idx+n-1, :) = [x(:), y(:)];
    all_EPVS_weights(curr_idx:curr_idx+n-1) = single(n);
    curr_idx = curr_idx + n;
end

%% Sampling grid (every ds voxels)
[xg_coarse, yg_coarse] = ndgrid(1:ds:dims(1), 1:ds:dims(2));
grid_coords = [xg_coarse(:), yg_coarse(:)];
sample_inds = sub2ind(dims, grid_coords(:,1), grid_coords(:,2));
valid_mask_samples = mask(sample_inds) & ~epvs(sample_inds);
valid_coords = grid_coords(valid_mask_samples, :);
num_valid = size(valid_coords, 1);

if num_valid == 0
    warning('No valid sampled voxels found.');
    subsampled_volume = sparse([]);
    interpolated_volume = single([]);
    return;
end

rows = valid_coords(:,1);
cols = valid_coords(:,2);

%% Build KD-tree on sampled points
tree_samples = KDTreeSearcher(double(valid_coords)); % requires double

%% Decide block and accumulation strategy
bytes_per_single = 4;
dense_bytes = num_valid * bytes_per_single;

numWorkers = 1;
if parpool_flag
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        poolobj = parpool; % start default pool
    end
    numWorkers = poolobj.NumWorkers;
end

nBlocks = max(1, ceil(n_epvs_vox / TARGET_EPVS_PER_BLOCK));
nBlocks = min(nBlocks, max(1, 2*numWorkers));
if parpool_flag
    nBlocks = max(nBlocks, numWorkers);
end

use_dense = (dense_bytes <= MAX_DENSE_BYTES_PER_WORKER);

if use_dense
    fprintf('Using dense per-worker accumulators (%.1f MB per vector)\n', dense_bytes/1e6);
else
    fprintf('Using sparse per-worker accumulators (dense %.1f MB too large)\n', dense_bytes/1e6);
end

%% Partition EPVS indices into blocks
block_edges = round(linspace(1, n_epvs_vox+1, nBlocks+1));
blocks = cell(nBlocks,1);
for b = 1:nBlocks
    blocks{b} = block_edges(b):(block_edges(b+1)-1);
end

%% Main accumulation with text progress
fprintf('Starting accumulation: %d EPVS -> %d samples (blocks=%d, workers=%d)\n', ...
n_epvs_vox, num_valid, nBlocks, numWorkers);
tic;

partials = cell(nBlocks,1);

if parpool_flag
    % Parallel blocks
    parfor b = 1:nBlocks
        idx_block = blocks{b};
        n_block = numel(idx_block);
        last_report = 0;
        if use_dense
            local_data = zeros(num_valid,1,'single');
            for jj = 1:n_block
                j = idx_block(jj);
                ep_coord = all_EPVS_coords(j, :);
                w = all_EPVS_weights(j);
    
                [idxs_cell, dists_cell] = rangesearch(tree_samples, ep_coord, radius);
                idxs = idxs_cell{1};
                dists = dists_cell{1};
    
                if ~isempty(idxs)
                    d_s = single(dists(:));
                    contrib = w ./ (d_s + single(epsilon)).^single(p);
                    local_data(idxs) = local_data(idxs) + contrib;
                end
    
                if mod(jj, N_EPVS_REPORT) == 0 || jj == n_block
                    % Progress message for this worker-block
                    fprintf('Worker block %d: processed %d/%d EPVS in block (global %d/%d)\n', ...
                            b, jj, n_block, idx_block(jj), n_epvs_vox);
                end
            end
            partials{b} = local_data;
        else
            idxs_all = zeros(0,1,'int32');
            contrib_all = zeros(0,1,'single');
            for jj = 1:n_block
                j = idx_block(jj);
                ep_coord = all_EPVS_coords(j, :);
                w = all_EPVS_weights(j);
    
                [idxs_cell, dists_cell] = rangesearch(tree_samples, ep_coord, radius);
                idxs = idxs_cell{1};
                dists = dists_cell{1};
    
                if ~isempty(idxs)
                    d_s = single(dists(:));
                    contrib = w ./ (d_s + single(epsilon)).^single(p);
                    idxs_all = [idxs_all; int32(idxs(:))]; %#ok<AGROW>
                    contrib_all = [contrib_all; single(contrib(:))]; %#ok<AGROW>
                end
    
                if mod(jj, N_EPVS_REPORT) == 0 || jj == n_block
                    fprintf('Worker block %d: processed %d/%d EPVS in block (global %d/%d)\n', ...
                            b, jj, n_block, idx_block(jj), n_epvs_vox);
                end
            end
            if isempty(idxs_all)
                partials{b} = sparse([],[],[],num_valid,1,0);
            else
                partials{b} = sparse(double(idxs_all), ones(numel(idxs_all),1), double(contrib_all), num_valid, 1);
            end
        end
    end
else

        % Serial blocks with progress
    total_processed = 0;
    for b = 1:nBlocks
        idx_block = blocks{b};
        n_block = numel(idx_block);
    
        fprintf('Starting block %d/%d (EPVS %d to %d)\n', b, nBlocks, idx_block(1), idx_block(end));
    
        if use_dense
            local_data = zeros(num_valid,1,'single');
            for jj = 1:n_block
                j = idx_block(jj);
                ep_coord = all_EPVS_coords(j, :);
                w = all_EPVS_weights(j);
    
                [idxs_cell, dists_cell] = rangesearch(tree_samples, ep_coord, radius);
                idxs = idxs_cell{1};
                dists = dists_cell{1};
    
                if ~isempty(idxs)
                    d_s = single(dists(:));
                    contrib = w ./ (d_s + single(epsilon)).^single(p);
                    local_data(idxs) = local_data(idxs) + contrib;
                end
    
                total_processed = total_processed + 1;
                if mod(total_processed, N_EPVS_REPORT) == 0 || jj == n_block
                    fprintf(' Block %d: processed %d/%d EPVS in block (global %d/%d)\n', ...
                            b, jj, n_block, total_processed, n_epvs_vox);
                end
            end
            partials{b} = local_data;
        else
            idxs_all = zeros(0,1,'int32');
            contrib_all = zeros(0,1,'single');
            for jj = 1:n_block
                j = idx_block(jj);
                ep_coord = all_EPVS_coords(j, :);
                w = all_EPVS_weights(j);
    
                [idxs_cell, dists_cell] = rangesearch(tree_samples, ep_coord, radius);
                idxs = idxs_cell{1};
                dists = dists_cell{1};
    
                if ~isempty(idxs)
                    d_s = single(dists(:));
                    contrib = w ./ (d_s + single(epsilon)).^single(p);
                    idxs_all = [idxs_all; int32(idxs(:))]; %#ok<AGROW>
                    contrib_all = [contrib_all; single(contrib(:))]; %#ok<AGROW>
                end
    
                total_processed = total_processed + 1;
                if mod(total_processed, N_EPVS_REPORT) == 0 || jj == n_block
                    fprintf(' Block %d: processed %d/%d EPVS in block (global %d/%d)\n', ...
                            b, jj, n_block, total_processed, n_epvs_vox);
                end
            end
            if isempty(idxs_all)
                partials{b} = sparse([],[],[],num_valid,1,0);
            else
                partials{b} = sparse(double(idxs_all), ones(numel(idxs_all),1), double(contrib_all), num_valid, 1);
            end
        end
        fprintf('Finished block %d/%d\n', b, nBlocks);
    end
end

% Reduce partials into final data vector
fprintf('Reducing partial results into final sampled vector...\n');
if use_dense
    data = zeros(num_valid,1,'single');
    for b = 1:nBlocks
        data = data + single(partials{b});
        partials{b} = [];
    end
else
    data = zeros(num_valid,1,'single');
    for b = 1:nBlocks
        if isa(partials{b},'double') || isa(partials{b},'single') || isa(partials{b},'logical')
            data = data + single(partials{b});
        else
            [ii, ~, vv] = find(partials{b});
            if isempty(ii)
                data(ii) = data(ii) + single(vv);
            end
        end
        partials{b} = [];
    end
end

fprintf('Accumulation and reduction complete (elapsed %.2f s)\n', toc);

%% Build subsampled volume from data
flat_idx = sub2ind(dims, rows, cols);
S = sparse(flat_idx, 1, double(data), prod(dims), 1);
subsampled_volume = single(reshape(full(S), dims));

%% Interpolated full-resolution map
fprintf('Starting interpolation\n');
Sgrid = single(subsampled_volume(1:ds:end, 1:ds:end));
F = griddedInterpolant(xg_coarse, yg_coarse, Sgrid, 'linear', 'nearest');
[xq, yq] = ndgrid(1:dims(1), 1:dims(2));
interpolated_volume = single(F(xq, yq));
fprintf('Interpolation done\n');

end