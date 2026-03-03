function [subsampled, interpolated, interpolated_ves_rm] = ...
    swp_voxelwise(epvs, mask, ves, radius, p, mode, chunk_size_in)
% EPVS density at 4x subsampled grid + interpolation
% INPUTS:
%   epvs (logical): binary EPVS mask
%   mask (logical): binary WM mask with vessels removed
% OUTPUTS:
%   subsampled (matrix): subsampled matrix
%   interpolated (matrix): interpolated SWP matrix
%   interpolated_ves (matrix): interpolated SWP matrix with the vessel
%                               lumens removed
% Modes:
%   'voxels'    - use every EPVS voxel (weight = 1) [default]
%   'centroids' - use connected-component centroids (weight = component size)
%
% Usage:
%   [Ssub, Sint] = epvs_density_variable_p(epvs, mask, radius, p)
%   [Ssub, Sint] = epvs_density_variable_p(..., mode)
%   [Ssub, Sint] = epvs_density_variable_p(..., mode, chunk_size)
%
% Date: February 24, 2026

if nargin < 6 || isempty(mode)
    mode = 'voxels';
end
if nargin < 7
    chunk_size_in = [];
end

dims = size(epvs);

%% Determine chunk_size robustly
if ~isempty(chunk_size_in)
    chunk_size = chunk_size_in;
else
    % Try memory(); if unavailable, try /proc/meminfo, Java, env, fallback
    chunk_size = [];
    try
        user_mem = memory; %#ok<NASGU>
        chunk_size = min(max(2000, floor(user_mem.MaxPossibleArrayBytes/200/8)), 20000);
    catch
        % Try /proc/meminfo (Linux)
        try
            fid = fopen('/proc/meminfo','r');
            if fid > 0
                s = textscan(fid, '%s %s %s', 'Delimiter', ':');
                fclose(fid);
                keys = s{1};
                vals = s{2};
                ii = find(strcmp(keys,'MemAvailable') | strcmp(keys,'MemFree') | strcmp(keys,'MemTotal'), 1);
                if ~isempty(ii)
                    v = vals{ii};
                    tok = textscan(v, '%f %s');
                    kb = tok{1};
                    bytes = single(kb * 1024);
                    chunk_size = min(max(2000, floor(bytes/200/8)), 20000);
                end
            end
        catch
            chunk_size = [];
        end
        if isempty(chunk_size)
            try
                rt = java.lang.Runtime.getRuntime();
                bytes = single(rt.maxMemory());
                chunk_size = min(max(2000, floor(bytes/200/8)), 20000);
            catch
                env_cs = getenv('EPVS_CHUNK_SIZE');
                if ~isempty(env_cs)
                    chunk_size = str2double(env_cs);
                else
                    chunk_size = 20000; % conservative fallback
                end
            end
        end
    end
end

%% Build points and weights depending on mode
switch lower(mode)
    case 'centroids'
        cc = bwconncomp(epvs, 26);
        nComp = cc.NumObjects;
        
        if nComp == 0
            warning('No EPVS voxels found.');
            subsampled = sparse([]);
            interpolated = single([]);
            return;
        end
        
        centroids = zeros(nComp, 3, 'single');
        % Initialize sparse weights
        weights = single(nComp, 1);  % Create a sparse vector for weights
        
        for i = 1:nComp
            idxs = cc.PixelIdxList{i};
            [x, y, z] = ind2sub(dims, idxs);
            cent = round(mean([x(:), y(:), z(:)], 1));
            centroids(i, :) = single(cent);
            weights(i) = single(numel(idxs));  % Assign number of elements as weight
        end
        points = centroids;  % (N x 3)
        
    case 'voxels'
        epi_inds = find(epvs);
        if isempty(epi_inds)
            warning('No EPVS voxels found.');
            subsampled = sparse([]);
            interpolated = single([]);
            return;
        end
        
        [ex, ey, ez] = ind2sub(dims, epi_inds);
        points = single([ex, ey, ez]);              % (S x 3)
        % Initialize sparse weights
        weights = sparse(size(points, 1), 1,'logical'); % Create a sparse vector for weights
        weights(:) = true;  % Assign weight = 1 for all EPVS voxels
end


%% Define sampling grid: every 4th voxel in 3D
[xg, yg, zg] = ndgrid(1:4:dims(1), 1:4:dims(2), 1:4:dims(3));
grid_coords = [xg(:), yg(:), zg(:)];
sample_inds = sub2ind(dims, grid_coords(:,1), grid_coords(:,2), grid_coords(:,3));
is_valid_sample = mask(sample_inds) & ~epvs(sample_inds);
valid_coords = grid_coords(is_valid_sample, :);
clear grid_coords;

% Ensure sampled coordinates are unique
flat_valid = sub2ind(dims, valid_coords(:,1), valid_coords(:,2), valid_coords(:,3));
mask_values = mask(flat_valid);
valid_coords = valid_coords(mask_values, :);
num_valid = size(valid_coords,1);

if num_valid == 0
    warning('No valid sampled voxels found.');
    subsampled = sparse([]);
    interpolated = single([]);
    return;
end

rows = valid_coords(:,1);
cols = valid_coords(:,2);
deps = valid_coords(:,3);
data = zeros(num_valid,1,'single');

%% Build KD-tree on points (single precision)
tree = KDTreeSearcher(single(points));

%% Precompute single casts
radius_s = single(radius);
p_s = single(p);

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
    [IdxChunk, DChunk] = rangesearch(tree, chunk_coords, radius_s);

    % Ensure DChunk and IdxChunk entries are column vectors (reshape only nonempty cells)
    nonemptyIdx = find(~cellfun('isempty', DChunk));
    for k = nonemptyIdx.'  % iterate over row vector for slight speed
        dk = DChunk{k};
        if size(dk,1) == 1 && size(dk,2) > 1
            DChunk{k} = dk(:);
        end
        ik = IdxChunk{k};
        if size(ik,1) == 1 && size(ik,2) > 1
            IdxChunk{k} = ik(:);
        end
    end

    % Vectorized accumulation using vertcat + splitapply
    len = cellfun(@numel, IdxChunk);
    totalNeigh = sum(len);

    if totalNeigh == 0
        vals = zeros(numel(chunk_inds),1,'single');
    else
        %%% Concatenate neighbor indices and distances into column vectors
        % Preallocate
        len = cellfun(@numel, IdxChunk);
        totalNeigh = sum(len);
        allIdx = zeros(totalNeigh,1,'like',IdxChunk{find(len,1)});
        allD   = zeros(totalNeigh,1,'like',DChunk{find(len,1)});
        % Fill in
        pos = 1;
        for i = 1:numel(IdxChunk)
            ni = len(i);
            if ni>0
                allIdx(pos:pos+ni-1) = IdxChunk{i}(:);
                allD(pos:pos+ni-1)   = DChunk{i}(:);
                pos = pos + ni;
            end
        end
        clear DChunk IdxChunk
        % Sanity check
        if numel(allIdx) ~= totalNeigh || numel(allD) ~= totalNeigh
            error('Concatenation mismatch after transpose fix: expected %d but got %d/%d.', totalNeigh, numel(allIdx), numel(allD));
        end

        %%% Compute contributions (single precision)
        % optimize denominator call
        if p == 1
            denom = single(allD);
        elseif p == 2
            denom = single(allD) .* single(allD);
        else
            denom = single(allD) .^ p_s;
        end
        
        %%% Calculate SWP
        % The weights have a value of "1" assigned to each EPVS voxel. This
        % will divide each weight by the respective distance scalar. The
        % subsequent code will sum over all neighborhood weights.
        contrib = single(weights(allIdx)) ./ denom;       % (totalNeigh x 1)
        
        % Calculate values for each
        vals = zeros(numel(len),1,'single');
        pos = 1;
        for q = 1:numel(len)
            nq = len(q);
            if nq > 0
                vals(q) = vals(q) + sum(contrib(pos:pos+nq-1));
                pos = pos + nq;
            end
        end
    end

    data(chunk_inds) = vals;

    clear allIdx allD contrib group len vals;

    % Print to console
    fprintf('Chunk %d/%d complete (%.1f%%)\n', chunk_idx, nChunks, 100 * chunk_idx / nChunks);
end

fprintf('Subsampled computation done in %.2f seconds.\n', toc);

%% Build subsampled outputs
flat_idx = sub2ind(dims, rows, cols, deps);
subsampled = zeros(dims, 'single');
subsampled(flat_idx) = data;

%% Prepare coarse grid for interpolation
% subsample coordinates
coarseX = single(1:4:dims(1));
coarseY = single(1:4:dims(2));
coarseZ = single(1:4:dims(3));

% subb-sampled grid coordinates
Sgrid = subsampled(1:4:end, 1:4:end, 1:4:end);

% Interpolated volume
F = griddedInterpolant({coarseX, coarseY, coarseZ}, Sgrid, 'linear', 'nearest');

interpolated = zeros(dims, 'single');
[X2D, Y2D] = ndgrid(single(1:dims(1)), single(1:dims(2)));
for z = 1:dims(3)
    Z2D = single(z) * ones(size(X2D),'single');
    slice = F(X2D, Y2D, Z2D);
    interpolated(:,:,z) = single(slice);
    if mod(z, max(1,ceil(dims(3)/10))) == 0
        fprintf('Interpolation slice %d/%d\n', z, dims(3));
    end
end

% Remove the vessel lumens from the SWP heatmap
interpolated_ves_rm = interpolated;
interpolated_ves_rm(ves) = 0;

fprintf('Interpolation complete.\n');
end
