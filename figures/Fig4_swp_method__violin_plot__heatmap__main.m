%% EPVS-SWP vs vessel-SWP in clean white matter
% Violins (SWP and object volume), joint densities, binned curves, and
% within-subject block-bootstrap Spearman's rho (voxel-level and on the
% binned curve).
%
% Sections: setup/config -> load caa -> SWP source -> pass 1 (data) ->
% violins -> heatmaps -> bootstrap/curves -> save. Every analysis section
% reads only from the data sections, so each can be run on its own in cell
% mode; the RUN flags skip sections when running the whole file.
%
% Design
%   Main:          control = CAA26, CAA6 ; severe = CAA22, CAA25 ; front + occip
%   Supplemental:  CAA17 (moderate), occipital only. Never pooled into a
%                  main group; reported in its own figures/rows.
%
% Inference
%   SWP is smoothed over a 500 um radius on a 20 um grid, so neighbouring
%   voxels are not independent. All intervals come from a spatial block
%   bootstrap within subject/region (block edge = 2 x SWP radius). Group
%   values weight subjects equally (mean of subject statistics; heatmaps
%   weight each subject's voxels to sum to 1/n_subjects).
%   With n = 2 subjects per group, group intervals describe within-subject
%   spatial uncertainty, NOT between-subject variability. Subject panels
%   are drawn beside every group panel so group claims can be checked.
%
% Two Spearman statistics (both block-bootstrapped)
%   Voxel rho : rank correlation of EPVS-SWP and vessel-SWP over voxels
%               (how tightly the two co-vary voxel by voxel).
%   Curve rho : rank correlation across the n_qbins points of the binned
%               curve (bin-median EPVS-SWP vs bin-median vessel-SWP), i.e.
%               how monotonic the curve is. Bins with < MIN_BIN_VOX voxels
%               (full data) are excluded; >= MIN_CURVE_BINS bins required.
%               Subject curve rho: on that subject's curve, recomputed in
%               every replicate. Group curve rho: on the group-mean curve
%               (the curve drawn), recomputed from the replicate group-mean
%               curves. Difference: severe - control per replicate.
%               With ~20 points the curve rho saturates near +/-1 for any
%               monotonic curve, so it describes shape, not strength.
%
% Outputs
%   <fig_root>/<EXPLORATORY|MANUSCRIPT>_<yyyy-MM-dd_HHmm>/
%     violins_main.*, violin_<swp|volume>_<reg>.*, violins_<run_date>.xlsx
%       (sheets: deltas, counts)
%     heatmap_<reg>[_zoom].* (tiled), heatmap_<reg>_<panel>[_zoom].* (single)
%     heatmap_occip_CAA17_supp[_zoom].{svg,png,fig}
%     curves_by_group.* (tiled), curves_by_group_<reg>.* (single),
%     curves_diff_severe_minus_control.*, curves_occip_CAA17_supp.*
%     spearman_block_bootstrap.csv        voxel rho
%     spearman_curve_block_bootstrap.csv  curve rho
%     object_summary.csv, results.mat
%
% Inputs (reused from the workspace when present, so this script reads the
% same data as the other pipeline scripts):
%   caa     caa.(sub).(reg) structs; only missing subjects are loaded
%   swp_op  swp_op.raw.(sub).(reg) = [mus, ret, vessel SWP, EPVS SWP] over
%           mask_wm & ~seg & ~epvs; only missing sub/regions are read from
%           the interpolated SWP heatmaps, using the same mask
%
% No external plotting dependencies (bin_swp, swp_scatterplot_overlay and
% violin_box_whisker_paired are no longer used).

%% Initialization
clearvars -except caa swp_op
clc; close all;

% Add project tree (two levels above the current directory)
topdir = fileparts(fileparts(pwd));
addpath(genpath(topdir));

% Threads: use allocated SCC slots, not the node's physical cores
nslots = str2double(getenv('NSLOTS'));
if isnan(nslots) || nslots < 1
    nslots = feature('numcores');
end
maxNumCompThreads(nslots);

%% Configuration
%%% Sections to run when executing the whole file
RUN = struct('violins', false, 'heatmaps', false, 'bootstrap', true);

EXPLORATORY = false;
if EXPLORATORY
    run_tag   = 'EXPLORATORY';
    n_boot    = 100;
    n_sub_vox = 2e4;     % voxels drawn per bootstrap replicate
else
    run_tag   = 'MANUSCRIPT';
    n_boot    = 5000;
    n_sub_vox = 1e5;
end
rng(1, 'twister');

%%% Directories
mat_dir  = '/projectnb/npbssmic/ns/CAA/';
swp_dir  = '/projectnb/npbssmic/ns/CAA/swp';
fig_root = '/projectnb/npbssmic/ns/CAA/figures/swp_method_distributions_heatmaps/';
run_date = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmm'));
out_dir  = fullfile(fig_root, sprintf('%s_%s', run_tag, run_date));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
fprintf('Output directory: %s\n', out_dir);

%%% SWP filenames
epvs_fbase = '_swp_voxelwise_radius_500_exp_2_interpolated_heatmap.mat';
ves_fbase  = '_swp_voxelwise_ves_radius_500_exp_2_interpolated_heatmap.mat';

%%% Geometry
vox_um        = 20;                       % isotropic voxel edge (um)
vox_vol_um3   = vox_um^3;
swp_radius_um = 500;
block_vox     = round(2 * swp_radius_um / vox_um);   % bootstrap block edge (voxels)
conn          = 26;
MIN_OBJ_VOX   = 100;                      % object-size floor for counts/volumes
MIN_BLOCKS    = 20;                       % warn if a region has fewer blocks

%%% Subject table (single source of truth for group membership)
S = table( ...
    {'caa26'; 'caa6'; 'caa22'; 'caa25'; 'caa17'}, ...
    {'control'; 'control'; 'severe'; 'severe'; 'moderate'}, ...
    {{'front','occip'}; {'front','occip'}; {'front','occip'}; {'front','occip'}; {'occip'}}, ...
    {'caa26/caa26.mat'; 'caa6/caa6.mat'; 'caa22/caa22.mat'; 'caa25/caa25.mat'; 'caa17/occip/caa17.mat'}, ...
    'VariableNames', {'sub', 'group', 'regions', 'matfile'});

main_groups = {'control', 'severe'};
supp_sub    = 'caa17';
regions     = {'front', 'occip'};

grp_disp = struct('control', 'Control', 'severe', 'Severe CAA', 'moderate', 'Moderate CAA');
% Structure colours (violins): EPVS blue, vessel orange
st_col   = struct('epvs', [0.098 0.392 0.690], ...    % #1964B0
                  'ves',  [0.859 0.345 0.161]);       % #DB5829
% Group colours (curves): kept distinct from the structure colours
grp_col  = struct('control',  [0.25 0.25 0.25], ...
                  'severe',   [0.80 0.15 0.45], ...
                  'moderate', [0.55 0.55 0.55]);
% Subject markers: diamond = 1st subject in a group, square = 2nd
sub_mrk  = struct('caa26', 'd', 'caa6', 's', 'caa22', 'd', 'caa25', 's', 'caa17', 'v');

%%% 2D heatmaps: common edges for control and severe within a region
n_bins_2d  = 150;
MIN_COUNT  = 3;          % raw voxel count below which a bin is drawn empty
CLIM_PCT   = 99.5;       % colour ceiling = this percentile of occupied bins
CLIM_DEC   = 3;          % log colour range below the ceiling (decades)
axis_lims.front = struct('epv', [0 600], 'ves', [0 200]);
axis_lims.occip = axis_lims.front;
zoom_lims.front = struct('epv', [0 50],  'ves', [0 50]);
zoom_lims.occip = zoom_lims.front;

%%% Binned curves
n_qbins           = 20;     % EPVS-SWP quantile bins (equal subject weight)
MIN_BIN_VOX       = 50;     % subject bin needs this many voxels to count
MIN_CURVE_BINS    = 5;      % curve rho needs at least this many valid bins
RESAMPLE_SUBJECTS = false;  % n = 2/group -> only 3 distinct subject sets

%%% Violins
VIOLIN_EQUAL_N  = true;          % SWP pooled violins: subsample each subject to the region's smallest n
SWP_VIOLIN_YLIM = [1e-1, 1e4];   % original units, log axis
MIN_N_VOX       = 50;            % subject marker hollow below this many voxels
MIN_N_OBJ       = 10;            % ... or this many objects
KDE_MAX_N       = 5e4;           % values used for each violin's KDE
CLIFF_MAX_N     = 2e4;           % values per side for unpaired Cliff's delta

%% CAA structs: reuse `caa` from the workspace if present
% Only subjects missing from `caa` are loaded, so a partially populated
% struct from another script is completed rather than replaced.
if ~exist('caa', 'var') || ~isstruct(caa)
    caa = struct();
end
for i = 1:height(S)
    sub = S.sub{i};
    if isfield(caa, sub)
        fprintf('caa.%s found in workspace\n', sub);
    else
        fprintf('Loading %s\n', sub);
        t = load(fullfile(mat_dir, S.matfile{i}), sub);
        caa.(sub) = t.(sub);
        clear t
        fprintf('Finished loading %s\n', sub);
    end
    regs = S.regions{i};
    for j = 1:numel(regs)
        assert(isfield(caa.(sub), regs{j}), 'caa.%s has no region "%s".', sub, regs{j});
    end
end

%% SWP source: reuse `swp_op` from the workspace if present
% swp_op.raw.(sub).(reg) rows are the voxels of mask_wm & ~seg & ~epvs in
% column-major order (as built by the optical-property scripts). Only the
% SWP columns are used. Sub/regions absent from swp_op are read from the
% interpolated heatmaps with the same mask, so both routes give the same rows.
SWP_VES_COL  = 3;
SWP_EPVS_COL = 4;
have_swp_op  = exist('swp_op', 'var') && isstruct(swp_op) && isfield(swp_op, 'raw');
if have_swp_op
    fprintf('swp_op found in workspace; using swp_op.raw columns %d (EPVS) and %d (vessel)\n', ...
            SWP_EPVS_COL, SWP_VES_COL);
else
    fprintf('swp_op not in workspace; importing SWP heatmaps\n');
end

%% Pass 1: clean-WM SWP pairs, bootstrap blocks, object volumes
pairs    = struct();
blocks   = struct();
vols     = struct();
obj_rows = {};

for i = 1:height(S)
    sub  = S.sub{i};
    regs = S.regions{i};
    for j = 1:numel(regs)
        reg = regs{j};
        c   = caa.(sub).(reg);

        % Clean WM, identical to the mask used to build swp_op
        clean = logical(c.mask_wm) & ~logical(c.seg) & ~logical(c.epvs);
        assert(ndims(clean) == 3, '%s/%s: expected a 3D volume.', sub, reg);
        idx = find(clean);

        if have_swp_op && isfield(swp_op.raw, sub) && isfield(swp_op.raw.(sub), reg)
            M = swp_op.raw.(sub).(reg);
            assert(size(M, 1) == numel(idx), ...
                   '%s/%s: swp_op has %d rows but mask_wm & ~seg & ~epvs has %d voxels.', ...
                   sub, reg, size(M, 1), numel(idx));
            P = single(M(:, [SWP_EPVS_COL, SWP_VES_COL]));
            clear M
        else
            if have_swp_op
                warning('%s/%s not in swp_op; reading SWP heatmaps.', sub, reg);
            end
            fprintf('Importing SWP for %s %s\n', sub, reg);
            fbase = fullfile(swp_dir, sub, reg, [sub '_' reg]);
            t = load([fbase epvs_fbase], 'interpolated_volume');
            e = t.interpolated_volume;
            t = load([fbase ves_fbase], 'interpolated_volume');
            v = t.interpolated_volume;
            clear t
            assert(isequal(size(e), size(v), size(clean)), ...
                   '%s/%s: SWP volumes and WM mask differ in size.', sub, reg);
            P = single([e(idx), v(idx)]);
            clear e v
        end

        % Finite rows only; keep voxel indices aligned with P
        ok  = all(isfinite(P), 2);
        idx = idx(ok);
        P   = P(ok, :);

        % Spatial block id for each retained voxel
        [ir, ic, is] = ind2sub(size(clean), idx);
        bsz = ceil(size(clean) / block_vox);
        bid = sub2ind(bsz, ceil(ir / block_vox), ceil(ic / block_vox), ceil(is / block_vox));
        clear ir ic is
        [~, ~, bid] = unique(bid);
        n_blocks = max(bid);
        if n_blocks < MIN_BLOCKS
            warning('%s/%s: only %d bootstrap blocks; intervals will be unstable.', ...
                    sub, reg, n_blocks);
        end

        pairs.(sub).(reg)  = P;
        blocks.(sub).(reg) = uint32(bid);

        % Object volumes (vessels: seg_wm_no_epvs within WM; EPVS mask as stored)
        ves_wm = logical(c.seg_wm_no_epvs) & logical(c.mask_wm);
        n_diff = nnz(xor(logical(c.seg) & logical(c.mask_wm), ves_wm));
        if n_diff > 0
            fprintf('  note: seg_wm_no_epvs (objects) and seg (clean-WM mask) differ in WM by %d voxels\n', n_diff);
        end
        cc_e = bwconncomp(logical(c.epvs), conn);
        cc_v = bwconncomp(ves_wm, conn);
        ne = cellfun(@numel, cc_e.PixelIdxList);
        nv = cellfun(@numel, cc_v.PixelIdxList);
        ve = ne(ne >= MIN_OBJ_VOX)' * vox_vol_um3;
        vv = nv(nv >= MIN_OBJ_VOX)' * vox_vol_um3;
        vols.(sub).(reg) = struct('epvs', ve, 'ves', vv);

        obj_rows(end+1, :) = {sub, S.group{i}, reg, size(P, 1), n_blocks, ...
                              numel(ve), median(ve), numel(vv), median(vv)}; %#ok<SAGROW>
        fprintf('%s/%s: %d clean-WM voxels, %d blocks, %d EPVS, %d vessels (>= %d vox)\n', ...
                sub, reg, size(P, 1), n_blocks, numel(ve), numel(vv), MIN_OBJ_VOX);
    end
end

obj_tbl = cell2table(obj_rows, 'VariableNames', ...
    {'subject', 'group', 'region', 'n_clean_wm_vox', 'n_blocks', ...
     'n_epvs', 'median_epvs_vol_um3', 'n_vessels', 'median_ves_vol_um3'});
writetable(obj_tbl, fullfile(out_dir, 'object_summary.csv'));

%% Violins: EPVS vs vessel SWP (paired voxels) and object volumes
% Needs only the three data sections above; run it on its own in cell
% mode, or set RUN.violins = false to skip it.
%
% Per region: Control EPVS | Control vessel || Severe EPVS | Severe vessel
%   Violin        KDE of the pooled data in log10 space (log axis). For SWP,
%                 each subject is subsampled to the region's smallest
%                 clean-WM count when VIOLIN_EQUAL_N, so subjects contribute
%                 equally to the pooled shape, median and IQR.
%   White circle  pooled median        Black bar  pooled Q1-Q3
%   Markers       subject medians from each subject's FULL data
%                 (diamond = 1st subject in the group, square = 2nd; hollow
%                 if n < MIN_N_VOX voxels / MIN_N_OBJ objects)
%   Colour        EPVS blue, vessel orange
% Cliff's delta (A vs B; > 0 means A larger), descriptive only:
%   within group  EPVS vs vessel (SWP: paired voxels, mean sign(A - B);
%                 volume: unpaired objects)
%   between       severe vs control, for EPVS and for vessel (unpaired)
% Unpaired deltas use at most CLIFF_MAX_N values per side.
if RUN.violins
    measures   = {'swp', 'volume'};
    meas_disp  = struct('swp', 'SWP', 'volume', 'Object volume');
    meas_ylab  = struct('swp', 'Size-weighted proximity', 'volume', 'Object volume (\mum^3)');
    st_names   = {'epvs', 'ves'};
    st_disp    = {'EPVS', 'Vessel'};
    mrk_order  = {'d', 's', '^', 'v'};
    vio        = struct();
    delta_rows = {};
    count_rows = {};

    for r = 1:numel(regions)
        reg     = regions{r};
        mem_all = cellfun(@(g) local_members(S, g, reg), main_groups, 'UniformOutput', false);
        if VIOLIN_EQUAL_N
            n_cap = min(cellfun(@(sb) size(pairs.(sb).(reg), 1), [mem_all{:}]));
        else
            n_cap = Inf;
        end

        for mi = 1:numel(measures)
            meas = measures{mi};
            pd = struct('data', {cell(1, 4)}, 'col', {cell(1, 4)}, ...
                        'subj', {cell(1, 4)}, 'xlab', {cell(1, 4)}, ...
                        'glab', {cellfun(@(g) grp_disp.(g), main_groups, 'UniformOutput', false)});
            for g = 1:numel(main_groups)
                grp = main_groups{g};
                mem = mem_all{g};
                nm  = numel(mem);
                pooled = cell(nm, 2);
                smed   = nan(nm, 2);
                sn     = zeros(nm, 2);
                gc     = zeros(1, 6);    % group sums for the counts sheet
                for k = 1:nm
                    sb = mem{k};
                    if strcmp(meas, 'swp')
                        P      = double(pairs.(sb).(reg));
                        n_full = size(P, 1);
                        posf   = all(P > 0, 2);
                        smed(k, :) = median(P(posf, :), 1);          % full data
                        sn(k, :)   = nnz(posf);
                        if n_full > n_cap
                            P = P(randperm(n_full, n_cap), :);
                        end
                        n_capd = size(P, 1);
                        pos    = all(P > 0, 2);                     % log axis
                        P      = P(pos, :);
                        pooled(k, :) = {P(:, 1), P(:, 2)};
                        crow = [n_full, n_capd, nnz(~pos), size(P, 1), NaN, NaN];
                    else
                        ve = vols.(sb).(reg).epvs;
                        vv = vols.(sb).(reg).ves;
                        smed(k, :)   = [median(ve), median(vv)];
                        sn(k, :)     = [numel(ve), numel(vv)];
                        pooled(k, :) = {ve, vv};
                        crow = [NaN, NaN, NaN, NaN, numel(ve), numel(vv)];
                    end
                    if any(sn(k, :) == 0)
                        fprintf('  %s %s/%s: subject has no data for one structure; no marker drawn\n', ...
                                meas, upper(sb), reg);
                    end
                    count_rows(end+1, :) = [{meas, reg, 'subject', grp, upper(sb)}, num2cell(crow)]; %#ok<SAGROW>
                    gc = gc + crow;
                end
                count_rows(end+1, :) = [{meas, reg, 'group', grp, strjoin(upper(mem), '+')}, num2cell(gc)]; %#ok<SAGROW>

                for si = 1:2
                    v = 2*(g - 1) + si;
                    pd.data{v} = vertcat(pooled{:, si});
                    pd.col{v}  = st_col.(st_names{si});
                    pd.xlab{v} = st_disp{si};   % group name is drawn under each pair
                    pd.subj{v} = struct('name', upper(mem), 'med', num2cell(smed(:, si)'), ...
                                        'n', num2cell(sn(:, si)'), 'mrk', mrk_order(1:nm));
                end
            end
            pd.lab_d = sprintf('%s / %s median', upper(mem_all{1}{1}), upper(mem_all{2}{1}));
            pd.lab_s = sprintf('%s / %s median', upper(mem_all{1}{2}), upper(mem_all{2}{2}));
            vio.(meas).(reg) = pd;

            %%% Cliff's delta
            paired = strcmp(meas, 'swp');
            for g = 1:numel(main_groups)            % within group: EPVS vs vessel
                A = pd.data{2*g - 1};
                B = pd.data{2*g};
                if paired
                    dl = mean(sign(A - B));
                    nu = [numel(A), numel(B)];
                    ty = 'paired';
                else
                    [dl, nu(1), nu(2)] = local_cliff(A, B, CLIFF_MAX_N);
                    ty = 'unpaired';
                end
                delta_rows(end+1, :) = {meas, reg, 'EPVS vs vessel', ty, ...
                    sprintf('%s EPVS', grp_disp.(main_groups{g})), ...
                    sprintf('%s vessel', grp_disp.(main_groups{g})), ...
                    dl, numel(A), numel(B), nu(1), nu(2)}; %#ok<SAGROW>
            end
            for si = 1:2                            % between: severe vs control
                A = pd.data{2 + si};
                B = pd.data{si};
                [dl, na, nb] = local_cliff(A, B, CLIFF_MAX_N);
                delta_rows(end+1, :) = {meas, reg, 'severe vs control', 'unpaired', ...
                    sprintf('%s %s', grp_disp.severe, st_disp{si}), ...
                    sprintf('%s %s', grp_disp.control, st_disp{si}), ...
                    dl, numel(A), numel(B), na, nb}; %#ok<SAGROW>
            end
        end
    end

    %%% Y-limits: fixed for SWP, shared across regions for volume
    vylim.swp = SWP_VIOLIN_YLIM;
    allv = [];
    for r = 1:numel(regions)
        allv = [allv; vertcat(vio.volume.(regions{r}).data{:})]; %#ok<AGROW>
    end
    allv = allv(allv > 0);
    vylim.volume = [10^floor(log10(min(allv))), 10^ceil(log10(max(allv)))];

    vcfg = struct('xc', [1 2 3.4 4.4], 'hw', 0.42, 'med_off', 0.12, 'med_inside', 0.75, ...
                  'kde_max_n', KDE_MAX_N, 'n_grid', 256);
    min_n = struct('swp', MIN_N_VOX, 'volume', MIN_N_OBJ);

    %%% Main figure: rows = measure, columns = region
    fig = figure('Color', 'w', 'Position', [100 100 650*numel(regions) 1100]);
    tl  = tiledlayout(fig, numel(measures), numel(regions), ...
                      'TileSpacing', 'compact', 'Padding', 'compact');
    for mi = 1:numel(measures)
        meas = measures{mi};
        vcfg.min_n = min_n.(meas);
        for r = 1:numel(regions)
            reg = regions{r};
            ax  = nexttile(tl);
            hleg = local_violin_panel(ax, vio.(meas).(reg), vylim.(meas), meas_ylab.(meas), ...
                       sprintf('%s - %s', meas_disp.(meas), reg), vcfg);
        end
    end
    lg = legend(ax, hleg, 'Orientation', 'horizontal', 'Box', 'off', 'AutoUpdate', 'off');
    lg.Layout.Tile = 'south';
    local_export(fig, out_dir, 'violins_main');
    close(fig);

    %%% One figure per panel
    for mi = 1:numel(measures)
        meas = measures{mi};
        vcfg.min_n = min_n.(meas);
        for r = 1:numel(regions)
            reg = regions{r};
            fig = figure('Color', 'w', 'Position', [100 100 650 620]);
            ax  = axes(fig);
            hleg = local_violin_panel(ax, vio.(meas).(reg), vylim.(meas), meas_ylab.(meas), ...
                       sprintf('%s - %s', meas_disp.(meas), reg), vcfg);
            legend(ax, hleg, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
                   'NumColumns', 3, 'Box', 'off', 'AutoUpdate', 'off');
            local_export(fig, out_dir, sprintf('violin_%s_%s', meas, reg));
            close(fig);
        end
    end

    %%% Spreadsheet: deltas + counts
    vio_delta_tbl = cell2table(delta_rows, 'VariableNames', ...
        {'measure', 'region', 'comparison', 'type', 'A', 'B', 'cliffs_delta', ...
         'n_A', 'n_B', 'n_A_used', 'n_B_used'});
    vio_count_tbl = cell2table(count_rows, 'VariableNames', ...
        {'measure', 'region', 'level', 'group', 'name', 'n_vox_clean_wm', ...
         'n_vox_after_cap', 'n_pairs_excluded_nonpos', 'n_vox_plotted', ...
         'n_epvs', 'n_vessels'});
    xlsx_path = fullfile(out_dir, sprintf('violins_%s.xlsx', run_date));
    writetable(vio_delta_tbl, xlsx_path, 'Sheet', 'deltas');
    writetable(vio_count_tbl, xlsx_path, 'Sheet', 'counts');
    disp(vio_delta_tbl);
    disp(vio_count_tbl);
    fprintf('Wrote violin deltas and counts to %s\n', xlsx_path);
end

%% 2D joint densities: group panels over subject panels, per region
if RUN.heatmaps
    cmap = local_density_cmap(256);
    clims_reg = struct();

    for r = 1:numel(regions)
        reg = regions{r};
        for zoom = [false true]
            if zoom
                L = zoom_lims.(reg);  sfx = '_zoom';
            else
                L = axis_lims.(reg);  sfx = '';
            end
            xe = linspace(L.epv(1), L.epv(2), n_bins_2d + 1);
            ye = linspace(L.ves(1), L.ves(2), n_bins_2d + 1);

            % Panels: groups in row 1 (each spanning 2 tiles), subjects in row 2
            panels = struct('ttl', {}, 'tag', {}, 'dens', {}, 'counts', {}, 'tile', {}, 'span', {});
            tile_sub = 5;
            for g = 1:numel(main_groups)
                grp = main_groups{g};
                mem = local_members(S, grp, reg);
                [P, w] = local_pool(pairs, mem, reg);
                [d, c] = local_hist2(P, w, xe, ye);
                panels(end+1) = struct('ttl', sprintf('%s (%s)', grp_disp.(grp), strjoin(upper(mem), ', ')), ...
                                       'tag', grp, 'dens', d, 'counts', c, ...
                                       'tile', 2*g - 1, 'span', [1 2]); %#ok<SAGROW>
                for k = 1:numel(mem)
                    [P, w] = local_pool(pairs, mem(k), reg);
                    [d, c] = local_hist2(P, w, xe, ye);
                    panels(end+1) = struct('ttl', upper(mem{k}), 'tag', upper(mem{k}), ...
                                           'dens', d, 'counts', c, ...
                                           'tile', tile_sub, 'span', [1 1]); %#ok<SAGROW>
                    tile_sub = tile_sub + 1;
                end
            end

            % Shared colour limits per region, set on the full-range figure and
            % reused for the zoom (density units are bin-size independent)
            if ~zoom
                occ = [];
                for p = 1:numel(panels)
                    occ = [occ; panels(p).dens(panels(p).counts >= MIN_COUNT)]; %#ok<AGROW>
                end
                top = prctile(occ, CLIM_PCT);
                clims_reg.(reg) = [top / 10^CLIM_DEC, top];
            end
            lims = clims_reg.(reg);

            fig = figure('Color', 'w', 'Position', [100 100 1600 850]);
            tl  = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
            title(tl, sprintf('%s: vessel SWP vs EPVS SWP (clean WM)%s', reg, ...
                  strrep(sfx, '_', ' ')), 'FontName', 'Arial', 'FontSize', 16);
            for p = 1:numel(panels)
                ax = nexttile(tl, panels(p).tile, panels(p).span);
                local_draw_tile(ax, xe, ye, panels(p).dens, panels(p).counts, ...
                                panels(p).ttl, lims, MIN_COUNT, cmap);
            end
            cb = colorbar(ax);
            cb.Layout.Tile = 'east';
            cb.Label.String = '% of voxels per SWP^2 (subjects equally weighted)';
            local_export(fig, out_dir, sprintf('heatmap_%s%s', reg, sfx));
            close(fig);

            % Same panels as single figures (same edges and colour limits)
            for p = 1:numel(panels)
                fig = figure('Color', 'w', 'Position', [100 100 700 600]);
                ax  = axes(fig);
                local_draw_tile(ax, xe, ye, panels(p).dens, panels(p).counts, ...
                                sprintf('%s - %s%s', panels(p).ttl, reg, strrep(sfx, '_', ' ')), ...
                                lims, MIN_COUNT, cmap);
                cb = colorbar(ax);
                cb.Label.String = '% of voxels per SWP^2';
                local_export(fig, out_dir, sprintf('heatmap_%s_%s%s', reg, panels(p).tag, sfx));
                close(fig);
            end

            % Supplemental CAA17 (occipital), same edges and colour limits
            if any(strcmp(local_members(S, 'moderate', reg), supp_sub))
                [P, w] = local_pool(pairs, {supp_sub}, reg);
                [d, c] = local_hist2(P, w, xe, ye);
                fig = figure('Color', 'w', 'Position', [100 100 700 600]);
                ax  = axes(fig);
                local_draw_tile(ax, xe, ye, d, c, sprintf('%s (moderate) %s%s', ...
                                upper(supp_sub), reg, strrep(sfx, '_', ' ')), lims, MIN_COUNT, cmap);
                cb = colorbar(ax);
                cb.Label.String = '% of voxels per SWP^2';
                local_export(fig, out_dir, sprintf('heatmap_%s_%s_supp%s', reg, upper(supp_sub), sfx));
                close(fig);
            end
        end
    end
end

%% Common EPVS-SWP quantile edges per region (main subjects, equal weight)
if RUN.bootstrap
    qedges = struct();
    qp = linspace(0, 1, n_qbins + 1);
    for r = 1:numel(regions)
        reg = regions{r};
        mem = [local_members(S, 'control', reg), local_members(S, 'severe', reg)];
        [P, w] = local_pool(pairs, mem, reg);
        e = unique(local_wquantile(P(:, 1), w, qp));
        e(1) = -Inf;
        e(end) = Inf;
        qedges.(reg) = e(:);
    end
end

%% Block bootstrap per subject/region: voxel rho, binned curve, curve rho
if RUN.bootstrap
    sub_res = struct();
    for i = 1:height(S)
        sub  = S.sub{i};
        regs = S.regions{i};
        for j = 1:numel(regs)
            reg = regs{j};
            fprintf('Bootstrap %s/%s (%d reps)\n', sub, reg, n_boot);
            P  = pairs.(sub).(reg);
            bb = local_block_index(blocks.(sub).(reg));
            e  = qedges.(reg);
            nb = numel(e) - 1;

            % Point estimates on all voxels. p_corr is corr's p-value, which
            % treats every voxel as independent (voxel-level, pseudoreplicated);
            % at these n it underflows to ~0 for any nonzero rho.
            [rho0, p0] = corr(double(P(:, 1)), double(P(:, 2)), 'Type', 'Spearman');
            [xm, ym, nbin] = local_bin_medians(P, e);
            ym(nbin < MIN_BIN_VOX) = NaN;
            xm(nbin < MIN_BIN_VOX) = NaN;
            [rc0, nc0] = local_curve_rho(xm, ym, MIN_CURVE_BINS);
            fprintf('  %s/%s: voxel rho = %.4f (p = %.3g, n = %d voxels); curve rho = %.4f (%d bins)\n', ...
                    sub, reg, rho0, p0, size(P, 1), rc0, nc0);

            % Replicates: resample blocks, draw n_sub_vox voxels, rebuild the
            % curve, and take Spearman's rho across its bins
            rho_b = nan(n_boot, 1);
            xm_b  = nan(n_boot, nb);
            ym_b  = nan(n_boot, nb);
            for b = 1:n_boot
                Pb = P(local_block_draw(bb, n_sub_vox), :);
                rho_b(b) = corr(double(Pb(:, 1)), double(Pb(:, 2)), 'Type', 'Spearman');
                [xb, yb] = local_bin_medians(Pb, e);
                xm_b(b, :) = xb';
                ym_b(b, :) = yb';
            end
            xm_b(:, nbin < MIN_BIN_VOX) = NaN;      % same bins as the point estimate
            ym_b(:, nbin < MIN_BIN_VOX) = NaN;
            rc_b = nan(n_boot, 1);
            for b = 1:n_boot
                rc_b(b) = local_curve_rho(xm_b(b, :), ym_b(b, :), MIN_CURVE_BINS);
            end

            sub_res.(sub).(reg) = struct('rho', rho0, 'p', p0, 'rho_boot', rho_b, ...
                'x', xm, 'y', ym, 'n_bin', nbin, 'x_boot', xm_b, 'y_boot', ym_b, ...
                'rho_curve', rc0, 'rho_curve_boot', rc_b, 'n_curve_bins', nc0, ...
                'n_vox', size(P, 1), 'n_blocks', numel(bb.cnt));
        end
    end
end

%% Group statistics (mean of subjects) and severe - control differences
if RUN.bootstrap
    rows   = {};
    rows_c = {};
    grp_res = struct();
    ci = @(x) prctile(x, [2.5 97.5]);
    % Two-sided block-bootstrap p against 0; see local_p_boot. When no
    % replicate falls on the far side of 0, only the bound 2/(n_boot+1)
    % can be reported (p_boot_bound = true).

    for r = 1:numel(regions)
        reg = regions{r};

        % Subject rows (main + supplemental)
        for i = 1:height(S)
            sub = S.sub{i};
            if ~isfield(sub_res.(sub), reg), continue; end
            sr = sub_res.(sub).(reg);
            q  = ci(sr.rho_boot);
            [pb, pbnd] = local_p_boot(sr.rho_boot, sprintf('%s/%s', sub, reg));
            if isnan(sr.p)
                warning('%s/%s: corr returned NaN (constant SWP column?).', sub, reg);
            end
            rows(end+1, :) = {'subject', upper(sub), S.group{i}, reg, sr.rho, q(1), q(2), ...
                              sr.p, pb, pbnd, sr.n_vox, sr.n_blocks}; %#ok<SAGROW>

            q = ci(sr.rho_curve_boot);
            [pb, pbnd] = local_p_boot(sr.rho_curve_boot, sprintf('%s/%s curve', sub, reg));
            rows_c(end+1, :) = {'subject', upper(sub), S.group{i}, reg, sr.rho_curve, ...
                                q(1), q(2), pb, pbnd, sr.n_curve_bins, ...
                                nnz(isfinite(sr.rho_curve_boot))}; %#ok<SAGROW>
        end

        % Group rows
        for g = 1:numel(main_groups)
            grp = main_groups{g};
            mem = local_members(S, grp, reg);
            grp_res.(grp).(reg) = local_group(sub_res, mem, reg, n_boot, RESAMPLE_SUBJECTS, MIN_CURVE_BINS);
            gr = grp_res.(grp).(reg);
            q  = ci(gr.rho_boot);
            [pb, pbnd] = local_p_boot(gr.rho_boot, sprintf('%s/%s', grp, reg));
            rows(end+1, :) = {'group', strjoin(upper(mem), '+'), grp, reg, gr.rho, q(1), q(2), ...
                              NaN, pb, pbnd, ...
                              sum(cellfun(@(s) sub_res.(s).(reg).n_vox, mem)), ...
                              sum(cellfun(@(s) sub_res.(s).(reg).n_blocks, mem))}; %#ok<SAGROW>

            q = ci(gr.rho_curve_boot);
            [pb, pbnd] = local_p_boot(gr.rho_curve_boot, sprintf('%s/%s curve', grp, reg));
            rows_c(end+1, :) = {'group', strjoin(upper(mem), '+'), grp, reg, gr.rho_curve, ...
                                q(1), q(2), pb, pbnd, gr.n_curve_bins, ...
                                nnz(isfinite(gr.rho_curve_boot))}; %#ok<SAGROW>
        end

        % Severe - control
        d_rho  = grp_res.severe.(reg).rho      - grp_res.control.(reg).rho;
        d_boot = grp_res.severe.(reg).rho_boot - grp_res.control.(reg).rho_boot;
        q = ci(d_boot);
        [pb, pbnd] = local_p_boot(d_boot, sprintf('severe-control/%s', reg));
        rows(end+1, :) = {'difference', 'severe-control', 'severe-control', reg, d_rho, q(1), q(2), ...
                          NaN, pb, pbnd, NaN, NaN}; %#ok<SAGROW>

        d_rc  = grp_res.severe.(reg).rho_curve      - grp_res.control.(reg).rho_curve;
        d_rcb = grp_res.severe.(reg).rho_curve_boot - grp_res.control.(reg).rho_curve_boot;
        q = ci(d_rcb);
        [pb, pbnd] = local_p_boot(d_rcb, sprintf('severe-control/%s curve', reg));
        rows_c(end+1, :) = {'difference', 'severe-control', 'severe-control', reg, d_rc, ...
                            q(1), q(2), pb, pbnd, NaN, nnz(isfinite(d_rcb))}; %#ok<SAGROW>
    end

    rho_tbl = cell2table(rows, 'VariableNames', ...
        {'level', 'name', 'group', 'region', 'spearman_rho', 'ci_lo', 'ci_hi', ...
         'p_corr', 'p_boot', 'p_boot_bound', 'n_vox', 'n_blocks'});
    % p_corr exists only for subject rows (one corr call on that subject's
    % voxels); group and difference rows are NaN by design.
    writetable(rho_tbl, fullfile(out_dir, 'spearman_block_bootstrap.csv'));

    rho_curve_tbl = cell2table(rows_c, 'VariableNames', ...
        {'level', 'name', 'group', 'region', 'curve_rho', 'ci_lo', 'ci_hi', ...
         'p_boot', 'p_boot_bound', 'n_bins', 'n_boot_valid'});
    writetable(rho_curve_tbl, fullfile(out_dir, 'spearman_curve_block_bootstrap.csv'));

    fprintf('\nSpearman rho, EPVS-SWP vs vessel-SWP (clean WM)\n');
    fprintf('p_corr: voxel-level corr p (subject rows only). p_boot: block-bootstrap, two-sided vs 0.\n');
    fprintf('%-10s %-12s %-14s %-6s %8s %18s %10s %10s\n', ...
            'level', 'name', 'group', 'region', 'rho', '95% CI', 'p_corr', 'p_boot');
    for k = 1:height(rho_tbl)
        if isnan(rho_tbl.p_corr(k))
            pc = '--';
        else
            pc = sprintf('%.3g', rho_tbl.p_corr(k));
        end
        fprintf('%-10s %-12s %-14s %-6s %8.4f   [%6.3f, %6.3f] %10s %10s\n', ...
                rho_tbl.level{k}, rho_tbl.name{k}, rho_tbl.group{k}, rho_tbl.region{k}, ...
                rho_tbl.spearman_rho(k), rho_tbl.ci_lo(k), rho_tbl.ci_hi(k), pc, ...
                local_fmt_pboot(rho_tbl.p_boot(k), rho_tbl.p_boot_bound(k)));
    end

    fprintf('\nCurve Spearman rho (across %d EPVS-SWP bins of the binned curve)\n', n_qbins);
    fprintf('%-10s %-12s %-14s %-6s %8s %18s %10s %6s\n', ...
            'level', 'name', 'group', 'region', 'rho', '95% CI', 'p_boot', 'bins');
    for k = 1:height(rho_curve_tbl)
        fprintf('%-10s %-12s %-14s %-6s %8.4f   [%6.3f, %6.3f] %10s %6s\n', ...
                rho_curve_tbl.level{k}, rho_curve_tbl.name{k}, rho_curve_tbl.group{k}, ...
                rho_curve_tbl.region{k}, rho_curve_tbl.curve_rho(k), ...
                rho_curve_tbl.ci_lo(k), rho_curve_tbl.ci_hi(k), ...
                local_fmt_pboot(rho_curve_tbl.p_boot(k), rho_curve_tbl.p_boot_bound(k)), ...
                num2str(rho_curve_tbl.n_bins(k)));
    end
end

%% Binned curves: groups with bands + subject overlays, front | occip
% Main tiled figure plus one figure per region; each panel is annotated
% with the group curve rho and its block-bootstrap 95% CI.
if RUN.bootstrap
    fig = figure('Color', 'w', 'Position', [100 100 1400 600]);
    tl  = tiledlayout(fig, 1, numel(regions), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Median vessel SWP by EPVS-SWP quantile bin (clean WM)', ...
          'FontName', 'Arial', 'FontSize', 16);
    axs = gobjects(1, numel(regions));
    for r = 1:numel(regions)
        reg = regions{r};
        ax  = nexttile(tl);
        axs(r) = ax;
        [h, lbl] = local_curve_panel(ax, reg, grp_res, sub_res, S, main_groups, ...
                                     grp_col, grp_disp, sub_mrk);
        if r == numel(regions)
            legend(ax, h, lbl, 'Location', 'eastoutside', 'Box', 'off');
        end
    end
    linkaxes(axs, 'y');
    yl_curves = ylim(axs(1));
    local_export(fig, out_dir, 'curves_by_group');
    close(fig);

    for r = 1:numel(regions)
        reg = regions{r};
        fig = figure('Color', 'w', 'Position', [100 100 800 600]);
        ax  = axes(fig);
        [h, lbl] = local_curve_panel(ax, reg, grp_res, sub_res, S, main_groups, ...
                                     grp_col, grp_disp, sub_mrk);
        ylim(ax, yl_curves);                 % same y-range as the tiled figure
        legend(ax, h, lbl, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
               'NumColumns', 3, 'Box', 'off');
        local_export(fig, out_dir, sprintf('curves_by_group_%s', reg));
        close(fig);
    end
end

%% Severe - control difference curves (same edges -> bin-aligned)
if RUN.bootstrap
    fig = figure('Color', 'w', 'Position', [100 100 1200 550]);
    tl  = tiledlayout(fig, 1, numel(regions), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Severe - control: median vessel SWP per EPVS-SWP bin', ...
          'FontName', 'Arial', 'FontSize', 16);
    axs = gobjects(1, numel(regions));
    for r = 1:numel(regions)
        reg = regions{r};
        ax  = nexttile(tl);
        axs(r) = ax;
        hold(ax, 'on');
        gs = grp_res.severe.(reg);
        gc = grp_res.control.(reg);
        xd = mean([gs.x, gc.x], 2);           % bins are shared; x = mean bin position
        dy = gs.y - gc.y;
        q  = prctile(gs.y_boot - gc.y_boot, [2.5 97.5], 1)';
        local_band(ax, xd, q(:, 1), q(:, 2), [0.3 0.3 0.3]);
        plot(ax, xd, dy, 'k-', 'LineWidth', 2);
        yline(ax, 0, '-', 'Color', [0.5 0.5 0.5]);
        local_style_curve_axes(ax, reg, xd);
        ylabel(ax, '\Delta median vessel SWP');
    end
    linkaxes(axs, 'y');
    local_export(fig, out_dir, 'curves_diff_severe_minus_control');
    close(fig);
end

%% Supplemental: CAA17 occipital curve against main occipital groups
if RUN.bootstrap
    if isfield(sub_res, supp_sub) && isfield(sub_res.(supp_sub), 'occip')
        fig = figure('Color', 'w', 'Position', [100 100 800 600]);
        ax  = axes(fig);
        hold(ax, 'on');
        h = gobjects(0);
        lbl = {};
        for g = 1:numel(main_groups)
            grp = main_groups{g};
            gr  = grp_res.(grp).occip;
            h(end+1) = plot(ax, gr.x, gr.y, '-', 'Color', grp_col.(grp), 'LineWidth', 2); %#ok<SAGROW>
            lbl{end+1} = grp_disp.(grp); %#ok<SAGROW>
        end
        sr  = sub_res.(supp_sub).occip;
        col = grp_col.moderate;
        q   = prctile(sr.y_boot, [2.5 97.5], 1)';
        local_band(ax, sr.x, q(:, 1), q(:, 2), col);
        h(end+1) = plot(ax, sr.x, sr.y, '-', 'Color', col, 'LineWidth', 2, ...
                        'Marker', sub_mrk.(supp_sub), 'MarkerFaceColor', 'w');
        qc = prctile(sr.rho_curve_boot, [2.5 97.5]);
        lbl{end+1} = sprintf('%s (moderate), \\rho_{curve} = %.2f [%.2f, %.2f]', ...
                             upper(supp_sub), sr.rho_curve, qc(1), qc(2));
        local_style_curve_axes(ax, 'occip', sr.x);
        legend(ax, h, lbl, 'Location', 'northwest', 'Box', 'off');
        local_export(fig, out_dir, sprintf('curves_occip_%s_supp', upper(supp_sub)));
        close(fig);
    end
end

%% Save results and configuration (whatever sections have run)
cfg = struct('run_tag', run_tag, 'run_date', run_date, 'RUN', RUN, 'n_boot', n_boot, ...
    'n_sub_vox', n_sub_vox, 'block_vox', block_vox, 'vox_um', vox_um, ...
    'MIN_OBJ_VOX', MIN_OBJ_VOX, 'n_bins_2d', n_bins_2d, 'MIN_COUNT', MIN_COUNT, ...
    'CLIM_PCT', CLIM_PCT, 'CLIM_DEC', CLIM_DEC, ...
    'axis_lims', axis_lims, 'zoom_lims', zoom_lims, 'n_qbins', n_qbins, ...
    'MIN_BIN_VOX', MIN_BIN_VOX, 'MIN_CURVE_BINS', MIN_CURVE_BINS, ...
    'RESAMPLE_SUBJECTS', RESAMPLE_SUBJECTS, ...
    'VIOLIN_EQUAL_N', VIOLIN_EQUAL_N, 'SWP_VIOLIN_YLIM', SWP_VIOLIN_YLIM, ...
    'MIN_N_VOX', MIN_N_VOX, 'MIN_N_OBJ', MIN_N_OBJ, 'KDE_MAX_N', KDE_MAX_N, ...
    'CLIFF_MAX_N', CLIFF_MAX_N, 'subjects', S);
save_vars = {'cfg', 'obj_tbl', 'vols'};
if exist('vio_delta_tbl', 'var')
    save_vars = [save_vars, {'vio_delta_tbl', 'vio_count_tbl'}];
end
if exist('clims_reg', 'var')
    cfg.clims = clims_reg;
end
if exist('rho_tbl', 'var')
    save_vars = [save_vars, {'sub_res', 'grp_res', 'rho_tbl', 'rho_curve_tbl', 'qedges'}];
end
save(fullfile(out_dir, 'results.mat'), save_vars{:}, '-v7.3');
fprintf('Done. Results in %s\n', out_dir);


%% ======================================================================
%%  Local functions
%% ======================================================================

function mem = local_members(S, grp, reg)
% Subjects in group GRP that have region REG (1 x k cell).
keep = strcmp(S.group, grp) & cellfun(@(rr) any(strcmp(rr, reg)), S.regions);
mem  = S.sub(keep)';
end

function [P, w] = local_pool(pairs, subs, reg)
% Stack subjects' pairs for one region; weights give each subject 1/n.
Pc = {};
for i = 1:numel(subs)
    if isfield(pairs, subs{i}) && isfield(pairs.(subs{i}), reg)
        Pc{end+1} = pairs.(subs{i}).(reg); %#ok<AGROW>
    end
end
n  = numel(Pc);
assert(n > 0, 'No data for region %s.', reg);
wc = cellfun(@(x) ones(size(x, 1), 1) / (size(x, 1) * n), Pc, 'UniformOutput', false);
P  = vertcat(Pc{:});
w  = vertcat(wc{:});
end

function [dens, counts] = local_hist2(P, w, xe, ye)
% Weighted 2D histogram as density (% of weighted voxels per SWP^2).
% Voxels outside the window still count in the denominator.
ix = discretize(double(P(:, 1)), xe);
iy = discretize(double(P(:, 2)), ye);
ok = ~isnan(ix) & ~isnan(iy);
sz = [numel(xe) - 1, numel(ye) - 1];
counts = accumarray([ix(ok) iy(ok)], 1,     sz);
mass   = accumarray([ix(ok) iy(ok)], w(ok), sz);
dens   = 100 * mass / (diff(xe(1:2)) * diff(ye(1:2)));
end

function q = local_wquantile(x, w, p)
% Weighted quantiles (lower interpolation).
[xs, o] = sort(double(x));
cw = cumsum(w(o));
cw = cw / cw(end);
q  = zeros(size(p));
for k = 1:numel(p)
    q(k) = xs(find(cw >= p(k) - eps, 1, 'first'));
end
end

function [xm, ym, n] = local_bin_medians(P, e)
% Per-bin voxel count and median EPVS / vessel SWP.
nb = numel(e) - 1;
k  = discretize(double(P(:, 1)), e);
ok = ~isnan(k);
n  = accumarray(k(ok), 1, [nb 1]);
xm = accumarray(k(ok), double(P(ok, 1)), [nb 1], @median, NaN);
ym = accumarray(k(ok), double(P(ok, 2)), [nb 1], @median, NaN);
end

function [rho, nv] = local_curve_rho(x, y, min_bins)
% Spearman's rho across the points of a binned curve (bin-median EPVS-SWP
% vs bin-median vessel-SWP), using bins where both are finite. NaN if
% fewer than min_bins valid bins or either coordinate is constant.
x  = double(x(:));
y  = double(y(:));
ok = isfinite(x) & isfinite(y);
nv = nnz(ok);
if nv < min_bins || numel(unique(x(ok))) < 2 || numel(unique(y(ok))) < 2
    rho = NaN;
    return
end
rho = corr(x(ok), y(ok), 'Type', 'Spearman');
end

function bb = local_block_index(bid)
% Rows sorted by block, with block sizes and start positions.
bid = double(bid(:));
[~, bb.order] = sort(bid);
bb.cnt   = accumarray(bid, 1);
bb.start = cumsum([1; bb.cnt(1:end-1)]);
end

function idx = local_block_draw(bb, n_draw)
% Resample blocks with replacement, then draw n_draw voxel rows uniformly
% from the resampled multiset (vectorised).
nb    = numel(bb.cnt);
pick  = randi(nb, nb, 1);
sz    = bb.cnt(pick);
cs    = cumsum(sz);
u     = randi(cs(end), n_draw, 1);
bpos  = discretize(u, [0.5; cs + 0.5]);
prev  = [0; cs(1:end-1)];
off   = u - prev(bpos);
idx   = bb.order(bb.start(pick(bpos)) + off - 1);
end

function gr = local_group(sub_res, mem, reg, n_boot, resample_subj, min_bins)
% Equal-weight group statistics from subject point estimates/replicates.
% Voxel rho: mean of subject rhos. Curve: mean of subject curves (bin-wise);
% curve rho: Spearman's rho of that mean curve, per replicate.
n  = numel(mem);
R0 = cellfun(@(s) sub_res.(s).(reg).rho, mem);
RB = cell2mat(cellfun(@(s) sub_res.(s).(reg).rho_boot, mem, 'UniformOutput', false));
X0 = cell2mat(cellfun(@(s) sub_res.(s).(reg).x, mem, 'UniformOutput', false));
Y0 = cell2mat(cellfun(@(s) sub_res.(s).(reg).y, mem, 'UniformOutput', false));
XB = cellfun(@(s) sub_res.(s).(reg).x_boot, mem, 'UniformOutput', false);
YB = cellfun(@(s) sub_res.(s).(reg).y_boot, mem, 'UniformOutput', false);
XB = cat(3, XB{:});                                   % n_boot x nb x n
YB = cat(3, YB{:});

if resample_subj
    sel = randi(n, n_boot, n);
else
    sel = repmat(1:n, n_boot, 1);
end
rho_boot   = zeros(n_boot, 1);
x_boot     = zeros(n_boot, size(XB, 2));
y_boot     = zeros(n_boot, size(YB, 2));
rhoc_boot  = nan(n_boot, 1);
for b = 1:n_boot
    rho_boot(b)  = mean(RB(b, sel(b, :)));
    x_boot(b, :) = mean(XB(b, :, sel(b, :)), 3);      % NaN if any subject lacks the bin
    y_boot(b, :) = mean(YB(b, :, sel(b, :)), 3);
    rhoc_boot(b) = local_curve_rho(x_boot(b, :), y_boot(b, :), min_bins);
end
xg = mean(X0, 2);
yg = mean(Y0, 2);
[rhoc, nbins] = local_curve_rho(xg, yg, min_bins);
gr = struct('rho', mean(R0), 'rho_boot', rho_boot, ...
            'x', xg, 'y', yg, 'x_boot', x_boot, 'y_boot', y_boot, ...
            'rho_curve', rhoc, 'rho_curve_boot', rhoc_boot, 'n_curve_bins', nbins, ...
            'members', {mem});
end

function [p, is_bound] = local_p_boot(x, label)
% Two-sided bootstrap p for H0: statistic = 0, from replicates x.
% p = 2 * min(#{x <= 0} + 1, #{x >= 0} + 1) / (n + 1), on finite replicates.
% If every replicate is on one side of 0, p equals its floor 2/(n+1) and is
% only an upper bound (is_bound = true): more replicates would lower it.
ok = isfinite(x);
if any(~ok)
    warning('%s: %d of %d bootstrap replicates are NaN/Inf and were dropped.', ...
            label, nnz(~ok), numel(x));
end
x = x(ok);
n = numel(x);
if n == 0
    p = NaN;
    is_bound = false;
    return
end
k_lo = sum(x <= 0);
k_hi = sum(x >= 0);
p = min(1, 2 * (min(k_lo, k_hi) + 1) / (n + 1));
is_bound = min(k_lo, k_hi) == 0;
end

function s = local_fmt_pboot(p, is_bound)
if isnan(p)
    s = '--';
elseif is_bound
    s = sprintf('<%.3g', p);
else
    s = sprintf('%.3g', p);
end
end

function [d, na, nb] = local_cliff(a, b, max_n)
% Unpaired Cliff's delta = P(A > B) - P(A < B), via the rank-sum identity.
% Each side is randomly subsampled to at most max_n values.
a = double(a(:));  a = a(isfinite(a));
b = double(b(:));  b = b(isfinite(b));
if numel(a) > max_n, a = a(randperm(numel(a), max_n)); end
if numel(b) > max_n, b = b(randperm(numel(b), max_n)); end
na = numel(a);
nb = numel(b);
if na == 0 || nb == 0
    d = NaN;
    return
end
rk = tiedrank([a; b]);
U  = sum(rk(1:na)) - na * (na + 1) / 2;
d  = 2 * U / (na * nb) - 1;
end

function [h, lbl] = local_curve_panel(ax, reg, grp_res, sub_res, S, main_groups, ...
                                      grp_col, grp_disp, sub_mrk)
% Group curves with 95% block-bootstrap bands, subject curves overlaid, and
% the group curve rho [95% CI] in the legend labels.
hold(ax, 'on');
h = gobjects(0);
lbl = {};
for g = 1:numel(main_groups)
    grp = main_groups{g};
    gr  = grp_res.(grp).(reg);
    col = grp_col.(grp);
    q   = prctile(gr.y_boot, [2.5 97.5], 1)';
    local_band(ax, gr.x, q(:, 1), q(:, 2), col);
    h(end+1) = plot(ax, gr.x, gr.y, '-', 'Color', col, 'LineWidth', 2.5); %#ok<AGROW>
    qc = prctile(gr.rho_curve_boot, [2.5 97.5]);
    lbl{end+1} = sprintf('%s, \\rho_{curve} = %.2f [%.2f, %.2f]', ...
                         grp_disp.(grp), gr.rho_curve, qc(1), qc(2)); %#ok<AGROW>
    mem = local_members(S, grp, reg);
    for k = 1:numel(mem)
        sr = sub_res.(mem{k}).(reg);
        h(end+1) = plot(ax, sr.x, sr.y, ':', 'Color', col, 'LineWidth', 1, ...
                        'Marker', sub_mrk.(mem{k}), 'MarkerSize', 5, ...
                        'MarkerFaceColor', 'w'); %#ok<AGROW>
        lbl{end+1} = sprintf('%s, \\rho_{curve} = %.2f', upper(mem{k}), sr.rho_curve); %#ok<AGROW>
    end
end
local_style_curve_axes(ax, reg, gr.x);
end

function hleg = local_violin_panel(ax, pd, ylims, ylab, ttl, vc)
% Four violins (log axis): pooled KDE, white pooled median, black Q1-Q3
% bar, subject-median markers kept inside the violin.
hold(ax, 'on');
any_hollow = false;
for v = 1:numel(pd.data)
    x = pd.data{v};
    x = x(isfinite(x) & x > 0);
    xc = vc.xc(v);
    if numel(x) < 2, continue; end
    lx = log10(x);
    if numel(lx) > vc.kde_max_n
        lk = lx(randperm(numel(lx), vc.kde_max_n));
    else
        lk = lx;
    end
    yg = linspace(min(lx), max(lx), vc.n_grid);
    f  = ksdensity(lk, yg);
    w  = vc.hw * f / max(f);
    fill(ax, [xc - w, fliplr(xc + w)], 10.^[yg, fliplr(yg)], pd.col{v}, ...
         'FaceAlpha', 0.55, 'EdgeColor', pd.col{v} * 0.7, 'LineWidth', 0.8, ...
         'HandleVisibility', 'off');
    q = prctile(x, [25 50 75]);
    plot(ax, [xc xc], q([1 3]), 'k-', 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(ax, xc, q(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');

    % Subject medians: spread +/- med_off, limited to med_inside x local half-width
    sj = pd.subj{v};
    ns = numel(sj);
    if ns > 1
        offs = vc.med_off * linspace(-1, 1, ns);
    else
        offs = 0;
    end
    for k = 1:ns
        m = sj(k).med;
        if ~isfinite(m) || m <= 0, continue; end
        wl = interp1(yg, w, log10(m), 'linear', 0);
        dx = sign(offs(k)) * min(abs(offs(k)), vc.med_inside * wl);
        hollow = sj(k).n < vc.min_n;
        any_hollow = any_hollow || hollow;
        if hollow, fc = 'w'; else, fc = 'k'; end
        plot(ax, xc + dx, m, sj(k).mrk, 'MarkerSize', 7, 'MarkerFaceColor', fc, ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end
xline(ax, mean(vc.xc(2:3)), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
set(ax, 'YScale', 'log', 'XLim', [vc.xc(1) - 0.7, vc.xc(end) + 0.7], ...
        'XTick', vc.xc, 'XTickLabel', pd.xlab, 'FontName', 'Arial', 'FontSize', 12, ...
        'Box', 'on', 'LineWidth', 1, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
ylim(ax, ylims);
ylabel(ax, ylab);
title(ax, ttl);

% Group name centred under each EPVS/vessel pair. Single-line tick labels:
% a newline inside a tick label is split into separate labels, which
% shifted every label by one position.
xl = get(ax, 'XLim');
for g = 1:numel(pd.glab)
    xm = mean(vc.xc(2*g - 1:2*g));
    text(ax, (xm - xl(1)) / diff(xl), -0.09, pd.glab{g}, 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
         'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
end

% Legend proxies (NaN data, nothing drawn)
hleg = gobjects(0);
hleg(end+1) = fill(ax, NaN, NaN, pd.col{1}, 'FaceAlpha', 0.55, 'EdgeColor', pd.col{1} * 0.7, ...
                   'DisplayName', 'EPVS');
hleg(end+1) = fill(ax, NaN, NaN, pd.col{2}, 'FaceAlpha', 0.55, 'EdgeColor', pd.col{2} * 0.7, ...
                   'DisplayName', 'Vessel');
hleg(end+1) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
                   'DisplayName', 'Pooled median');
hleg(end+1) = plot(ax, NaN, NaN, 'k-', 'LineWidth', 4, 'DisplayName', 'Pooled IQR');
hleg(end+1) = plot(ax, NaN, NaN, 'kd', 'MarkerFaceColor', 'k', 'DisplayName', pd.lab_d);
hleg(end+1) = plot(ax, NaN, NaN, 'ks', 'MarkerFaceColor', 'k', 'DisplayName', pd.lab_s);
if any_hollow
    hleg(end+1) = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor', 'w', ...
                       'DisplayName', sprintf('subject n < %d', vc.min_n));
end
end

function local_band(ax, x, lo, hi, col)
x = x(:); lo = lo(:); hi = hi(:);
ok = isfinite(x) & isfinite(lo) & isfinite(hi);
if nnz(ok) < 2, return; end
fill(ax, [x(ok); flipud(x(ok))], [lo(ok); flipud(hi(ok))], col, ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function local_style_curve_axes(ax, reg, x)
x = x(isfinite(x));
if ~isempty(x) && min(x) > 0
    set(ax, 'XScale', 'log');
end
set(ax, 'FontName', 'Arial', 'FontSize', 14, 'Box', 'on', 'LineWidth', 1, ...
        'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
xlabel(ax, 'EPVS SWP (bin median)');
ylabel(ax, 'Median vessel SWP');
title(ax, reg);
end

function local_draw_tile(ax, xe, ye, dens, counts, ttl, lims, min_count, cmap)
% Empty / under-threshold bins are transparent, not the bottom colour.
xc = xe(1:end-1) + diff(xe) / 2;
yc = ye(1:end-1) + diff(ye) / 2;
D  = dens;
D(counts < min_count) = NaN;
h  = imagesc(ax, xc, yc, D');                 % rows = y
set(h, 'AlphaData', ~isnan(D'));
set(ax, 'YDir', 'normal', 'Color', 'w', 'ColorScale', 'log', ...
        'Box', 'on', 'Layer', 'top', 'LineWidth', 1, ...
        'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], ...
        'FontName', 'Arial', 'FontSize', 12);
colormap(ax, cmap);
local_set_clim(ax, lims);
xlim(ax, xe([1 end]));
ylim(ax, ye([1 end]));
title(ax, ttl);
xlabel(ax, 'EPVS SWP');
ylabel(ax, 'Vessel SWP');
end

function local_set_clim(ax, lims)
try
    clim(ax, lims);
catch
    caxis(ax, lims); %#ok<CAXIS>
end
end

function cmap = local_density_cmap(n)
% Light-to-dark sequential ramp, monotonically decreasing luminance.
anchors = [0.78 0.86 0.94;
           0.42 0.63 0.83;
           0.10 0.39 0.69;
           0.06 0.24 0.44;
           0.03 0.11 0.22];
x  = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = [interp1(x, anchors(:, 1), xi)', ...
        interp1(x, anchors(:, 2), xi)', ...
        interp1(x, anchors(:, 3), xi)'];
cmap = min(max(cmap, 0), 1);
end

function local_export(fig, out_dir, stem)
f = fullfile(out_dir, stem);
exportgraphics(fig, [f '.svg'], 'ContentType', 'vector');
exportgraphics(fig, [f '.png'], 'Resolution', 600);
savefig(fig, [f '.fig']);
fprintf('  saved %s.{svg,png,fig}\n', f);
end
