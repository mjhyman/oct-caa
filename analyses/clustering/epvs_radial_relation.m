function results = epvs_radial_relation(brain_mask, epvs_mask,...
                                        voxel_size_mm, varargin)
% EPVS_RADIAL_RELATION Computes spatial relationship between large and small EPVS.
%
%   'LargePercentile' (90) - Top X% of volumes treated as "seeds"
%   'RadialBins' ([0 0.5 1 2 4]) - Distance bins in mm for shell analysis
%   'NumPermutations' (500) - Number of random-label shuffles
%   'MinSeeds' (5) - Minimum number of seeds required to run
%   'IncludePermutations' (true) - Run null distribution
%   'Verbose' (true) - Print progress to console
%   'MinVolumeMM3' (0.0005) - Filter out objects smaller than this
%   'Debug' (false) - Show 3D visualization of the largest split object
try
% ---- Parse inputs ----
p = inputParser;
addRequired(p,'brain_mask',@(x)islogical(x) || isnumeric(x));
addRequired(p,'epvs_mask',@(x)islogical(x) || numeric(x));
addRequired(p,'voxel_size_mm',@(x)isnumeric(x) && (isscalar(x)||numel(x)==3));
addParameter(p,'LargePercentile',90,@isnumeric);
addParameter(p,'RadialBins', [0 0.5 1 2 4], @isnumeric);
addParameter(p,'NumPermutations',500,@isnumeric);
addParameter(p,'MinSeeds',5,@isnumeric);
addParameter(p,'IncludePermutations',true,@islogical);
addParameter(p,'Verbose',true,@islogical);
addParameter(p,'MinVolumeMM3', 0.0005, @isnumeric); 
addParameter(p,'Debug', false, @islogical);
parse(p,brain_mask,epvs_mask,voxel_size_mm,varargin{:});
params = p.Results;

if isscalar(params.voxel_size_mm)
    vx = repmat(params.voxel_size_mm,1,3); 
else 
    vx = params.voxel_size_mm(:)'; 
end
voxel_vol = prod(vx);

% ---- Watershed Segmentation ----
if params.Verbose, fprintf('Splitting trunks and branches...\n'); end
D = bwdist(~epvs_mask);
h_threshold = 2 * min(vx); 
D_sloped = imhmax(D, h_threshold);
L = watershed(-D_sloped);
L(~epvs_mask) = 0; 

% ---- Filtering and Centroids ----
stats = regionprops(L, 'Area', 'Centroid', 'PixelIdxList', 'BoundingBox');
all_volumes_mm3 = [stats.Area]' * voxel_vol;
keep_idx = all_volumes_mm3 >= params.MinVolumeMM3;
filtered_stats = stats(keep_idx);
volumes_mm3 = all_volumes_mm3(keep_idx);
M = numel(filtered_stats);
raw_centroids = cat(1, filtered_stats.Centroid);
centroids_mm = (raw_centroids - 0.5) .* vx; 

% ---- Optional Debugging Figure ----
if params.Debug && M > 0
    [~, biggest] = max(volumes_mm3);
    bb = floor(filtered_stats(biggest).BoundingBox);
    pad = 5;
    rowR = max(1, bb(2)-pad):min(size(L,1), bb(2)+bb(5)+pad);
    colR = max(1, bb(1)-pad):min(size(L,2), bb(1)+bb(4)+pad);
    pagR = max(1, bb(3)-pad):min(size(L,3), bb(3)+bb(6)+pad);
    crop = L(rowR, colR, pagR);
    
    figure('Color', 'w', 'Name', 'Watershed Split Debug');
    p_fig = patch(isosurface(crop > 0, 0.5));
    p_fig.FaceColor = 'interp';
    p_fig.FaceVertexCData = crop(isosurface(crop > 0, 0.5)); 
    set(p_fig, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    view(3); camlight; lighting gouraud; grid on;
    title('Largest EPVS: Trunk vs. Branches');
end

%% Define Binning and Seeds
binEdges = params.RadialBins;
nbins = numel(binEdges) - 1;
binMids = (binEdges(1:end-1) + binEdges(2:end)) / 2;
maxr = binEdges(end);

large_thresh = prctile(volumes_mm3, params.LargePercentile);
seed_idx = find(volumes_mm3 >= large_thresh);      
neighbor_idx = find(volumes_mm3 < large_thresh);   
numSeeds = numel(seed_idx);

if params.Verbose
    fprintf('Detected %d EPVS; seeds=%d (thresh=%.4g mm^3)\n', M, numSeeds, large_thresh);
end

if numSeeds < params.MinSeeds
    warning('Insufficient seeds. Output will be NaN.');
    results = struct(); return;
end

%% Prepare voxel coords for NON-EPVS brain voxels
non_epvs_brain_mask = brain_mask & ~epvs_mask;
[idxX, idxY, idxZ] = ind2sub(size(brain_mask), find(non_epvs_brain_mask));
voxel_coords = [(idxX - 0.5) * vx(1), (idxY - 0.5) * vx(2), (idxZ - 0.5) * vx(3)];

% Build KD-tree
useKD = license('test','statistics_toolbox') && exist('createns','file');
if useKD
    if params.Verbose, fprintf('Building KD-Trees...\n'); end
    voxelTree = createns(voxel_coords,'NSMethod','kdtree','Distance','euclidean');
    centTree = createns(centroids_mm,'NSMethod','kdtree','Distance','euclidean');
end

%% Precompute Shell Voxel Counts
seed_shell_vox_counts_all = zeros(M, nbins); 
if params.Verbose, fprintf('Computing shell volumes (voxel counting)...\n'); end

for i = 1:M
    pt = centroids_mm(i,:);
    if useKD
        [~, dists] = rangesearch(voxelTree, pt, maxr);
        dists = dists{1}(:);
        if ~isempty(dists)
            [~,~,binidx] = histcounts(dists, binEdges);
            for b = 1:nbins
                seed_shell_vox_counts_all(i,b) = sum(binidx==b);
            end
        end
    else
        dists_all = sqrt(sum((voxel_coords - pt).^2,2));
        maskv = dists_all <= maxr;
        if any(maskv)
            [~,~,binidx] = histcounts(dists_all(maskv), binEdges);
            for b = 1:nbins
                seed_shell_vox_counts_all(i,b) = sum(binidx==b);
            end
        end
    end
    if params.Verbose && mod(i,50)==0, fprintf('\tObject %i of %i\n',i,M); end
end
seed_shell_vol_all = seed_shell_vox_counts_all * voxel_vol;

%% Precompute Neighbor Centroid Distances
neighbor_list = cell(M,1);
if params.Verbose, fprintf('Computing neighbor distances...\n'); end

if useKD
    for i = 1:M
        [idxs, dists] = rangesearch(centTree, centroids_mm(i,:), maxr);
        idxs = idxs{1}(:); dists = dists{1}(:);
        self = (idxs == i);
        idxs(self) = []; dists(self) = [];
        if ~isempty(idxs), neighbor_list{i} = [idxs, dists]; end
    end
else
    D_mat = squareform(pdist(centroids_mm));
    for i = 1:M
        idxsLocal = find(D_mat(i,:) <= maxr & (1:M) ~= i);
        if ~isempty(idxsLocal)
            neighbor_list{i} = [idxsLocal(:), D_mat(i,idxsLocal)'];
        end
    end
end

%% Observed Density Calculation

seed_epvs_counts = zeros(numSeeds, nbins);
seed_shell_volumes = zeros(numSeeds, nbins);

for s = 1:numSeeds
    idx = seed_idx(s);
    nl = neighbor_list{idx};
    if ~isempty(nl)
        isSmall = ismember(nl(:,1), neighbor_idx);
        if any(isSmall)
            seed_epvs_counts(s,:) = histcounts(nl(isSmall,2), binEdges);
        end
    end
    seed_shell_volumes(s,:) = seed_shell_vol_all(idx,:);
end

seed_epvs_density = seed_epvs_counts ./ seed_shell_volumes;
seed_epvs_density(~isfinite(seed_epvs_density)) = NaN;
observed_mean_density = mean(seed_epvs_density, 1, 'omitmissing');
observed_sem = std(seed_epvs_density, 0, 1, 'omitmissing') ./ sqrt(sum(~isnan(seed_epvs_density),1));

%% Permutations
if params.IncludePermutations && params.NumPermutations > 0
    if params.Verbose, fprintf('Running %d Permutations...\n', params.NumPermutations); end
    nperm = params.NumPermutations;
    perm_mean = NaN(nperm, nbins);
    rng(0); 
    
    for p = 1:nperm
        perm_idx = randperm(M, numSeeds); 
        perm_small_mask = true(M,1);
        perm_small_mask(perm_idx) = false;
        
        per_seed_counts = zeros(numSeeds, nbins);
        per_seed_vols = zeros(numSeeds, nbins);
        
        for k = 1:numSeeds
            objidx = perm_idx(k);
            nl = neighbor_list{objidx};
            if ~isempty(nl)
                keepMask = perm_small_mask(nl(:,1));
                if any(keepMask)
                    per_seed_counts(k,:) = histcounts(nl(keepMask,2), binEdges);
                end
            end
            per_seed_vols(k,:) = seed_shell_vol_all(objidx,:);
        end
        dens = per_seed_counts ./ per_seed_vols;
        dens(~isfinite(dens)) = NaN;
        perm_mean(p,:) = mean(dens,1,'omitmissing');
    end
    
    perm_median = median(perm_mean, 1, 'omitmissing');
    perm_lo = prctile(perm_mean, 2.5, 1);
    perm_hi = prctile(perm_mean, 97.5, 1);
    pvals = (sum(perm_mean >= observed_mean_density, 1) + 1) / (nperm + 1);
else
    perm_mean = []; perm_median = []; perm_lo = []; perm_hi = []; pvals = [];
end

%% Collect outputs
results.observed_mean_density = observed_mean_density;
results.observed_sem = observed_sem;
results.pvals = pvals;
results.binMids = binMids;
results.perm_median = perm_median;
results.params = params;
results.perm_lo = perm_lo;
results.perm_hi = perm_hi;
results.perm_mean = perm_mean;
catch
    pause(0.1)
end
end