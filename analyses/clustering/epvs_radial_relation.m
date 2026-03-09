function results = epvs_radial_relation(brain_mask, epvs_mask,...
                                        voxel_size_mm, varargin)
% EPVS_RADIAL_RELATION Computes spatial relationship between large and small EPVS.
%
%   'LargePercentile'     (90)    - Top X% of volumes treated as "seeds"
%   'RadialBins'    ([0 0.5 1 2 4]) - Distance bins in mm for shell analysis
%   'NumPermutations'     (500)   - Number of random-label shuffles
%   'MinSeeds'            (5)     - Minimum number of seeds required to run
%   'IncludePermutations' (true)  - Run null distribution
%   'Verbose'             (true)  - Print progress to console
%   'MinVolumeMM3'        (0.0005)- Filter out objects smaller than this
%   'Debug'               (false) - Show 3D visualization
%   'WatershedFraction'   (0.35)  - Adaptive split threshold as fraction of local radius
%   'WatershedH'          ([])    - Fixed watershed h in mm; [] = adaptive
%   'WatershedMinRadius'  ([])    - Min radius in mm below which no splitting occurs; [] = auto
%   'SkelSkipN'           (3)     - Sample every Nth skeleton point for rangesearch

try
%% ---- Parse inputs ----
p = inputParser;
addRequired(p,'brain_mask',    @(x)islogical(x)||isnumeric(x));
addRequired(p,'epvs_mask',     @(x)islogical(x)||isnumeric(x));
addRequired(p,'voxel_size_mm', @(x)isnumeric(x)&&(isscalar(x)||numel(x)==3));
addParameter(p,'LargePercentile',    90,            @isnumeric);
addParameter(p,'RadialBins',         [0 0.5 1 2 4], @isnumeric);
addParameter(p,'NumPermutations',    500,           @isnumeric);
addParameter(p,'MinSeeds',           5,             @isnumeric);
addParameter(p,'IncludePermutations',true,          @islogical);
addParameter(p,'Verbose',            true,          @islogical);
addParameter(p,'MinVolumeMM3',       0.0005,        @isnumeric);
addParameter(p,'Debug',              false,         @islogical);
addParameter(p,'WatershedFraction',  0.35,          @isnumeric);
addParameter(p,'WatershedH',         [],            @(x)isempty(x)||isnumeric(x));
addParameter(p,'WatershedMinRadius', [],            @(x)isempty(x)||isnumeric(x));
addParameter(p,'SkelSkipN',          3,             @isnumeric);
parse(p, brain_mask, epvs_mask, voxel_size_mm, varargin{:});
params = p.Results;

if isscalar(params.voxel_size_mm)
    vx = repmat(params.voxel_size_mm, 1, 3);
else
    vx = params.voxel_size_mm(:)';
end
voxel_vol = prod(vx);

if isempty(params.WatershedMinRadius)
    params.WatershedMinRadius = 1.5 * min(vx);
end

%% ---- Watershed Segmentation ----
if params.Verbose, fprintf('Splitting trunks and branches...\n'); end

epvs_mask  = logical(epvs_mask);
brain_mask = logical(brain_mask);
D          = bwdist(~epvs_mask);

if ~isempty(params.WatershedH)
    if params.Verbose
        fprintf('  Watershed mode: fixed h = %.4g mm\n', params.WatershedH);
    end
    h_threshold = params.WatershedH / min(vx);
    D_sloped    = imhmax(D .* epvs_mask, h_threshold);
else
    if params.Verbose
        fprintf('  Watershed mode: adaptive fraction=%.2f, minRadius=%.4g mm\n',...
            params.WatershedFraction, params.WatershedMinRadius);
    end
    frac       = params.WatershedFraction;
    min_radius = params.WatershedMinRadius / min(vx);
    D_masked   = D .* epvs_mask;
    D_sloped   = D_masked;
    local_max  = imregionalmax(D_masked) & epvs_mask;
    lm_labels  = bwlabeln(local_max);
    lm_stats   = regionprops(lm_labels, D_masked, 'MaxIntensity','PixelIdxList');

    n_suppressed = 0;
    for lm = 1:numel(lm_stats)
        local_radius = lm_stats(lm).MaxIntensity;
        if local_radius < min_radius
            continue
        end
        h_local       = frac * local_radius;
        pix           = lm_stats(lm).PixelIdxList;
        D_sloped(pix) = min(D_sloped(pix), local_radius - h_local);
        n_suppressed  = n_suppressed + 1;
    end
    if params.Verbose
        fprintf('  Suppressed %d of %d local maxima (%.1f%%)\n',...
            n_suppressed, numel(lm_stats),...
            100*n_suppressed/max(numel(lm_stats),1));
    end
end

fprintf('   Applying watershed\n')
L = watershed(-D_sloped);
L(~epvs_mask) = 0;

%% ---- Filtering and Centroids ----
fprintf('   Filtering and finding centroids\n')
stats          = regionprops(L,'Area','Centroid','PixelIdxList','BoundingBox');
all_volumes    = [stats.Area]' * voxel_vol;
keep_idx       = all_volumes >= params.MinVolumeMM3;
filtered_stats = stats(keep_idx);
volumes_mm3    = all_volumes(keep_idx);
M              = numel(filtered_stats);

if M == 0
    warning('No EPVS objects survived filtering. Output will be NaN.');
    results = struct(); return;
end

raw_centroids = cat(1, filtered_stats.Centroid);
centroids_mm  = (raw_centroids - 0.5) .* vx;

%% ---- Classify Seeds and Neighbours ----
fprintf('   Classifying seeds and neighbors\n')
binEdges     = params.RadialBins;
nbins        = numel(binEdges) - 1;
binMids      = (binEdges(1:end-1) + binEdges(2:end)) / 2;
maxr         = binEdges(end);

large_thresh = prctile(volumes_mm3, params.LargePercentile);
seed_idx     = find(volumes_mm3 >= large_thresh);
neighbor_idx = find(volumes_mm3 <  large_thresh);
numSeeds     = numel(seed_idx);

if params.Verbose
    fprintf('Detected %d EPVS; seeds=%d (thresh=%.4g mm^3)\n',...
        M, numSeeds, large_thresh);
end

if numSeeds < params.MinSeeds
    warning('Insufficient seeds (%d < %d). Output will be NaN.',...
        numSeeds, params.MinSeeds);
    results = struct(); return;
end

%% ---- Skeletonise Large EPVS ----
if params.Verbose, fprintf('Skeletonising %d large EPVS...\n', numSeeds); end

seed_skel_pts  = cell(numSeeds, 1);
seed_branchpts = cell(numSeeds, 1);
seed_endpoints = cell(numSeeds, 1);
skipN          = max(1, round(params.SkelSkipN));
vol_sz         = size(L);
conn26         = conndef(3,'maximal');

for s = 1:numSeeds
    idx = seed_idx(s);
    pix = filtered_stats(idx).PixelIdxList;

    bb  = floor(filtered_stats(idx).BoundingBox);
    pad = 3;
    r1 = max(1,bb(2)-pad);   r2 = min(vol_sz(1),bb(2)+bb(5)+pad);
    c1 = max(1,bb(1)-pad);   c2 = min(vol_sz(2),bb(1)+bb(4)+pad);
    z1 = max(1,bb(3)-pad);   z2 = min(vol_sz(3),bb(3)+bb(6)+pad);

    obj_crop          = false(vol_sz);
    obj_crop(pix)     = true;
    obj_crop          = obj_crop(r1:r2, c1:c2, z1:z2);
    skel_crop         = bwskel(obj_crop,'MinBranchLength',2);

    nb_count = uint8(imfilter(uint8(skel_crop), double(conn26),'same')) ...
               - uint8(skel_crop);
    bp_crop  = skel_crop & (nb_count >= 3);
    ep_crop  = skel_crop & (nb_count == 1);

    [sr,sc,sz] = ind2sub(size(skel_crop), find(skel_crop));
    sr = sr+r1-1; sc = sc+c1-1; sz = sz+z1-1;

    [br,bc,bz] = ind2sub(size(bp_crop), find(bp_crop));
    br = br+r1-1; bc = bc+c1-1; bz = bz+z1-1;

    [er,ec,ez] = ind2sub(size(ep_crop), find(ep_crop));
    er = er+r1-1; ec = ec+c1-1; ez = ez+z1-1;

    ns   = numel(sr);
    keep = 1:skipN:ns;
    sr   = sr(keep); sc = sc(keep); sz = sz(keep);

    skel_mm = [(sr-0.5)*vx(1), (sc-0.5)*vx(2), (sz-0.5)*vx(3)];
    bp_mm   = [(br-0.5)*vx(1), (bc-0.5)*vx(2), (bz-0.5)*vx(3)];
    ep_mm   = [(er-0.5)*vx(1), (ec-0.5)*vx(2), (ez-0.5)*vx(3)];

    if isempty(skel_mm)
        skel_mm = centroids_mm(idx,:);
    end

    seed_skel_pts{s}  = skel_mm;
    seed_branchpts{s} = bp_mm;
    seed_endpoints{s} = ep_mm;

    if params.Verbose && mod(s,20)==0
        fprintf('  Skeleton %d of %d\n', s, numSeeds);
    end
end

% Build lookup: object index -> seed list index
seed_obj_to_s = containers.Map(seed_idx, num2cell(1:numSeeds));

%% ---- Debug Figure ----
if params.Debug
    display_vol            = zeros(size(L),'single');
    display_vol(epvs_mask) = 1;
    display_vol(L > 0)     = 2;
    kept_mask              = false(size(L));
    for k = 1:numel(filtered_stats)
        kept_mask(filtered_stats(k).PixelIdxList) = true;
    end
    display_vol(kept_mask) = 3;

    iso_removed   = isosurface(display_vol>=1 & display_vol<2, 0.5);
    iso_watershed = isosurface(display_vol>=2 & display_vol<3, 0.5);
    iso_kept      = isosurface(display_vol>=3, 0.5);

    colors = struct(...
        'removed',   [0.85, 0.33, 0.10],...
        'watershed', [0.30, 0.60, 0.90],...
        'kept',      [0.47, 0.87, 0.47]);

    figure('Color','k','Name','EPVS Mask vs Watershed vs Kept');
    ax = axes('Parent',gcf,'Color','k');
    hold(ax,'on');
    legend_handles = [];
    legend_entries = {};

    if ~isempty(iso_removed.vertices)
        p1 = patch(ax,iso_removed);
        p1.FaceColor=colors.removed; p1.EdgeColor='none'; p1.FaceAlpha=0.35;
        legend_handles(end+1) = patch(ax,'XData',[],'YData',[],'ZData',[],...
            'FaceColor',colors.removed,'EdgeColor','none');
        legend_entries{end+1} = 'Removed by watershed';
    end
    if ~isempty(iso_watershed.vertices)
        p2 = patch(ax,iso_watershed);
        p2.FaceColor=colors.watershed; p2.EdgeColor='none'; p2.FaceAlpha=0.50;
        legend_handles(end+1) = patch(ax,'XData',[],'YData',[],'ZData',[],...
            'FaceColor',colors.watershed,'EdgeColor','none');
        legend_entries{end+1} = 'Removed by size filter';
    end
    if ~isempty(iso_kept.vertices)
        p3 = patch(ax,iso_kept);
        p3.FaceColor=colors.kept; p3.EdgeColor='none'; p3.FaceAlpha=0.80;
        legend_handles(end+1) = patch(ax,'XData',[],'YData',[],'ZData',[],...
            'FaceColor',colors.kept,'EdgeColor','none');
        legend_entries{end+1} = 'Kept EPVS';
    end

    all_skel = cell2mat(seed_skel_pts);
    all_bp   = cell2mat(seed_branchpts);
    all_ep   = cell2mat(seed_endpoints);

    if ~isempty(all_skel)
        sv = all_skel ./ vx + 0.5;
        scatter3(ax,sv(:,2),sv(:,1),sv(:,3),2,[1 1 0],'filled');
        legend_handles(end+1) = scatter3(ax,nan,nan,nan,20,[1 1 0],'filled');
        legend_entries{end+1} = 'Skeleton';
    end
    if ~isempty(all_bp)
        bv = all_bp ./ vx + 0.5;
        scatter3(ax,bv(:,2),bv(:,1),bv(:,3),25,[1 0 1],'filled');
        legend_handles(end+1) = scatter3(ax,nan,nan,nan,25,[1 0 1],'filled');
        legend_entries{end+1} = 'Branch points';
    end
    if ~isempty(all_ep)
        ev = all_ep ./ vx + 0.5;
        scatter3(ax,ev(:,2),ev(:,1),ev(:,3),25,[1 0.5 0],'filled');
        legend_handles(end+1) = scatter3(ax,nan,nan,nan,25,[1 0.5 0],'filled');
        legend_entries{end+1} = 'End points';
    end

    view(ax,3); camlight(ax,'headlight'); camlight(ax,'left');
    lighting(ax,'gouraud'); axis(ax,'equal','tight','off');
    legend(ax,legend_handles,legend_entries,...
        'TextColor','w','Color','none','EdgeColor','w','Location','northeast');

    n_orig    = sum(epvs_mask(:));
    n_kept_v  = sum(kept_mask(:));
    n_removed = n_orig - n_kept_v;
    title(ax,sprintf('EPVS: %d original | %d kept | %d removed (%.1f%%)',...
        n_orig,n_kept_v,n_removed,100*n_removed/n_orig),...
        'Color','w','FontSize',11);
end

%% ---- KD-Trees ----
fprintf('   Building KD-Trees\n')
non_epvs_brain_mask = brain_mask & ~epvs_mask;
[idxX,idxY,idxZ]    = ind2sub(size(brain_mask),find(non_epvs_brain_mask));
voxel_coords        = [(idxX-0.5)*vx(1),(idxY-0.5)*vx(2),(idxZ-0.5)*vx(3)];

useKD = license('test','statistics_toolbox') && exist('createns','file');
if useKD
    if params.Verbose, fprintf('Building KD-Trees...\n'); end
    voxelTree = createns(voxel_coords, 'NSMethod','kdtree','Distance','euclidean');
    centTree  = createns(centroids_mm, 'NSMethod','kdtree','Distance','euclidean');
end

%% ---- Shell Volumes ----
% Seeds   -> minimum distance from any skeleton point to each brain voxel
% Neighbours -> centroid distance (they are approximately spherical)
if params.Verbose, fprintf('Computing shell volumes...\n'); end
seed_shell_vol_all = zeros(M, nbins);

for i = 1:M
    isSeed = isKey(seed_obj_to_s, i);

    if isSeed
        s_pos    = seed_obj_to_s(i);
        skel_pts = seed_skel_pts{s_pos};

        if useKD
            [nb_idx, nb_dist] = rangesearch(voxelTree, skel_pts, maxr);
            n_vox    = size(voxel_coords,1);
            min_dist = inf(n_vox,1,'single');
            for sp = 1:numel(nb_idx)
                vi = nb_idx{sp}(:);
                di = single(nb_dist{sp}(:));
                improve          = di < min_dist(vi);
                min_dist(vi(improve)) = di(improve);
            end
        else
            n_vox    = size(voxel_coords,1);
            min_dist = inf(n_vox,1,'single');
            for sp = 1:size(skel_pts,1)
                d        = single(sqrt(sum((voxel_coords - skel_pts(sp,:)).^2,2)));
                min_dist = min(min_dist, d);
            end
        end

        valid = isfinite(min_dist) & min_dist <= maxr;
        if any(valid)
            [~,~,binidx] = histcounts(min_dist(valid), binEdges);
            for b = 1:nbins
                seed_shell_vol_all(i,b) = sum(binidx==b) * voxel_vol;
            end
        end

    else
        pt = centroids_mm(i,:);
        if useKD
            [~,dists] = rangesearch(voxelTree, pt, maxr);
            dists = dists{1}(:);
        else
            dists = sqrt(sum((voxel_coords - pt).^2,2));
            dists = dists(dists <= maxr);
        end
        if ~isempty(dists)
            [~,~,binidx] = histcounts(dists, binEdges);
            for b = 1:nbins
                seed_shell_vol_all(i,b) = sum(binidx==b) * voxel_vol;
            end
        end
    end

    if params.Verbose && mod(i,50)==0
        fprintf('  Shell volumes: object %d of %d\n', i, M);
    end
end

%% ---- Neighbour Distances ----
% Initial pass: centroid-to-centroid (catches all pairs within maxr)
% Then refine any edge involving a seed to use skeleton min-distance
%   - seed -> neighbour : rangesearch from seed skeleton outward
%   - neighbour -> seed : update stored distance to min skel distance

if params.Verbose, fprintf('Computing neighbour distances...\n'); end

neighbor_list = cell(M,1);

% --- Centroid pass (all objects) ---
fprintf('   Finding neighbor distance between centroids\n')
if useKD
    for i = 1:M
        [idxs,dists] = rangesearch(centTree, centroids_mm(i,:), maxr);
        idxs  = idxs{1}(:); dists = dists{1}(:);
        self  = (idxs == i);
        idxs(self) = []; dists(self) = [];
        if ~isempty(idxs)
            neighbor_list{i} = [idxs, dists];
        end
    end
else
    D_mat = squareform(pdist(centroids_mm));
    for i = 1:M
        cols = find(D_mat(i,:) <= maxr & (1:M)~=i);
        if ~isempty(cols)
            neighbor_list{i} = [cols(:), D_mat(i,cols)'];
        end
    end
end

% --- Skeleton refinement: seed -> neighbours (outward pass) ---
if params.Verbose, fprintf('Refining seed neighbour distances via skeleton...\n'); end

for s = 1:numSeeds
    idx      = seed_idx(s);
    skel_pts = seed_skel_pts{s};

    if isempty(skel_pts), continue; end

    % Rangesearch from ALL skeleton points to all centroids
    if useKD
        [nb_idx_all, nb_dist_all] = rangesearch(centTree, skel_pts, maxr);
    else
        nb_idx_all  = cell(size(skel_pts,1),1);
        nb_dist_all = cell(size(skel_pts,1),1);
        for sp = 1:size(skel_pts,1)
            d  = sqrt(sum((centroids_mm - skel_pts(sp,:)).^2,2));
            in = find(d <= maxr);
            nb_idx_all{sp}  = in';
            nb_dist_all{sp} = d(in)';
        end
    end

    % Accumulate minimum distance per neighbour centroid
    min_dist = inf(M,1,'single');
    for sp = 1:numel(nb_idx_all)
        vi = nb_idx_all{sp}(:);
        di = single(nb_dist_all{sp}(:));
        improve           = di < min_dist(vi);
        min_dist(vi(improve)) = di(improve);
    end
    min_dist(idx) = inf;   % remove self

    valid = find(isfinite(min_dist) & min_dist <= maxr);
    if ~isempty(valid)
        neighbor_list{idx} = [valid(:), double(min_dist(valid))];
    else
        neighbor_list{idx} = [];
    end
end

% --- Skeleton refinement: neighbour -> seed (inward pass) ---
% For every non-seed that has a seed in its neighbor_list,
% recompute that distance as min dist from neighbour centroid to seed skeleton.
for nb = 1:M
    nl = neighbor_list{nb};
    if isempty(nl), continue; end

    for row = 1:size(nl,1)
        obj_j = nl(row,1);
        if ~isKey(seed_obj_to_s, obj_j), continue; end

        s_j      = seed_obj_to_s(obj_j);
        skel_pts = seed_skel_pts{s_j};
        if isempty(skel_pts), continue; end

        pt_nb     = centroids_mm(nb,:);
        d_skel    = sqrt(sum((skel_pts - pt_nb).^2,2));
        nl(row,2) = min(d_skel);
    end
    neighbor_list{nb} = nl;
end

%% ---- Observed Density ----
seed_epvs_counts   = zeros(numSeeds, nbins);
seed_shell_volumes = zeros(numSeeds, nbins);

for s = 1:numSeeds
    idx = seed_idx(s);
    nl  = neighbor_list{idx};
    if ~isempty(nl)
        isSmall = ismember(nl(:,1), neighbor_idx);
        if any(isSmall)
            seed_epvs_counts(s,:) = histcounts(nl(isSmall,2), binEdges);
        end
    end
    seed_shell_volumes(s,:) = seed_shell_vol_all(idx,:);
end

seed_epvs_density             = seed_epvs_counts ./ seed_shell_volumes;
seed_epvs_density(~isfinite(seed_epvs_density)) = NaN;
observed_mean_density         = mean(seed_epvs_density, 1, 'omitmissing');
observed_sem                  = std(seed_epvs_density,  0, 1, 'omitmissing') ...
                                ./ sqrt(sum(~isnan(seed_epvs_density),1));

%% ---- Permutations ----
if params.IncludePermutations && params.NumPermutations > 0
    if params.Verbose
        fprintf('Running %d permutations...\n', params.NumPermutations);
    end
    nperm     = params.NumPermutations;
    perm_mean = NaN(nperm, nbins);
    rng(0);

    for perm = 1:nperm
        perm_idx        = randperm(M, numSeeds);
        perm_small_mask = true(M,1);
        perm_small_mask(perm_idx) = false;

        per_seed_counts = zeros(numSeeds, nbins);
        per_seed_vols   = zeros(numSeeds, nbins);

        for k = 1:numSeeds
            objidx = perm_idx(k);
            nl     = neighbor_list{objidx};
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
        perm_mean(perm,:) = mean(dens,1,'omitmissing');

        if params.Verbose && mod(perm,100)==0
            fprintf('  Permutation %d of %d\n', perm, nperm);
        end
    end

    perm_median = median(perm_mean, 1, 'omitmissing');
    perm_lo     = prctile(perm_mean, 2.5,  1);
    perm_hi     = prctile(perm_mean, 97.5, 1);
    pvals       = (sum(perm_mean >= observed_mean_density,1)+1) / (nperm+1);
else
    perm_mean=[]; perm_median=[]; perm_lo=[]; perm_hi=[]; pvals=[];
end

%% ---- Collect outputs ----
results.observed_mean_density = observed_mean_density;
results.observed_sem          = observed_sem;
results.pvals                 = pvals;
results.binMids               = binMids;
results.perm_median           = perm_median;
results.perm_lo               = perm_lo;
results.perm_hi               = perm_hi;
results.perm_mean             = perm_mean;
results.params                = params;
results.numSeeds              = numSeeds;
results.numObjects            = M;
results.large_thresh_mm3      = large_thresh;
results.seed_skel_pts         = seed_skel_pts;
results.seed_branchpts        = seed_branchpts;
results.seed_endpoints        = seed_endpoints;

catch ME
    warning(ME.identifier,...
        'epvs_radial_relation failed: %s\nAt: %s line %d',...
        ME.message, ME.stack(1).file, ME.stack(1).line);
    rethrow(ME);
end
end