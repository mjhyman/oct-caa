%% compare_orient_vessel_groups.m
% Fiber (as-PSOCT optic axis) vs vessel orientation in white matter:
% control (CAA26, CAA6) vs severe CAA (CAA22, CAA25), frontal and occipital.
%
%   RUN_MODE = 'check' : orient/retardance diagnostics + glyph overlays per subject x region
%   RUN_MODE = 'full'  : group analysis, figures, CSVs
%
% Geometry: dims 1-2 = en-face plane, dim 3 = through-plane. The optic axis is an
% in-plane axial angle; vessel axes are projected onto the en-face plane, and vessels
% inclined >= MAX_INCL_DEG are excluded from angle-based metrics.
%
% Vessel axis: 3-D structure tensor of the smoothed vessel mask at vessel-wall voxels
% (smallest-eigenvalue eigenvector; integration scale chosen per anchor by max tube
% coherence). Each perivascular WM voxel is anchored to its nearest vessel-wall voxel.
% Shell distance is measured from the nearest structure (vessel or EPVS) when
% DIST_REF = 'nearest' (matches the annular analysis), or from the vessel wall.
%
% Retardance: optic-axis angle is unreliable at low retardance. Fiber voxels below
% RET_FLOOR (a percentile of pooled control WM retardance, per region, applied to
% every subject) are excluded; alignment is also reported by retardance tertile.
%
% Metrics (per anchor, fibers within PRIMARY_MAX_UM):
%   dtheta  = |axial difference| vessel axis vs fiber circular mean, [0,90] deg
%   excess  = median(shuffled-null dtheta) - median(observed dtheta); > 0 = vessels
%             aligned with local fibers beyond what regional orientation statistics give
%   R       = fiber circular resultant length (local orientation coherence)
%
% Statistics: group value = mean of subject values (subjects weighted equally).
% CIs: blocks (BLOCK_UM cubes) resampled within subject, subjects fixed; the
% severe - control CI is computed on the difference itself. No p-values (n = 2/group);
% CIs describe the measured subjects only. Leave-one-subject-out estimates reported.

%% 0. Config ---------------------------------------------------------------
RUN_MODE    = 'full';         % 'check' | 'full'
EXPLORATORY = true;

GROUPS      = struct('control', {{'CAA26','CAA6'}}, 'severe', {{'CAA22','CAA25'}});  % reference first
GROUP_LABEL = struct('control', 'Control', 'severe', 'Severe');
GROUP_COLOR = struct('control', [0.45 0.45 0.45], 'severe', [0.80 0.20 0.15]);
SUBJ_MARKER = struct('CAA26', 'o', 'CAA6', 's', 'CAA22', '^', 'CAA25', 'd');
REGIONS     = {'front','occip'};
REG_LABEL   = struct('front', 'Frontal', 'occip', 'Occipital');

VOX_UM = 20;

ORIENT_UNITS     = 'deg';
ORIENT_ZERO_AXIS = 'dim2';    % confirm with RUN_MODE = 'check'
ORIENT_CCW       = true;

% vessel axis (structure tensor)
SIGMA_D_UM    = 20;
SIGMA_I_UM    = [40 80 160];
MIN_TUBE_COH  = 0.5;
RAD_SEARCH_UM = 160;          % local radius = max interior distance within this cube half-width

% geometry
DIST_REF        = 'nearest';  % 'nearest' (vessel or EPVS) | 'vessel'
SHELL_UM        = 40;
MAX_DIST_UM     = 400;
MAX_WALL_UM     = 600;        % fiber voxels farther than this from any vessel wall are not anchored
PRIMARY_MAX_UM  = 120;
MAX_INCL_DEG    = 45;
EPVS_ADJ_UM     = 40;         % anchor counted as EPVS-adjacent if an EPVS voxel lies within this distance
MIN_VOX         = 5;          % min fiber voxels per (anchor, shell) and per anchor (primary window)
MIN_ROWS        = 20;         % min anchors for any median/excess estimate

% retardance
RET_FLOOR_PCTL  = 10;         % percentile of pooled control WM retardance (0 = no floor)
RET_WEIGHT      = false;      % retardance-weighted circular means
RET_SAMPLE      = 2e6;        % WM voxels sampled per control subject to set floor/tertiles

BLOCK_UM        = 400;
N_NULL_SHUFFLE  = 10;

if EXPLORATORY
    N_BOOT = 100; MAX_ANCHORS = 50000; TAG = 'EXPLORATORY';
else
    N_BOOT = 500; MAX_ANCHORS = Inf;   TAG = 'MANUSCRIPT';
end

ROOT = '/projectnb/npbssmic/ns/CAA/';
OUT_DIR = fullfile(ROOT, 'figures/orientation/orient_vessel_groups', ...
    sprintf('%s_%s_%s', TAG, RUN_MODE, datestr(now, 'yyyymmdd_HHMMSS')));
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
rng(1);

prm = struct('vox', VOX_UM, 'units', ORIENT_UNITS, 'zeroAxis', ORIENT_ZERO_AXIS, ...
    'ccw', ORIENT_CCW, 'sigmaD', SIGMA_D_UM, 'sigmaI', SIGMA_I_UM, 'minCoh', MIN_TUBE_COH, ...
    'radSearch', RAD_SEARCH_UM, 'distRef', DIST_REF, 'shell', SHELL_UM, 'maxDist', MAX_DIST_UM, ...
    'maxWall', MAX_WALL_UM, 'primaryMax', PRIMARY_MAX_UM, 'maxIncl', MAX_INCL_DEG, ...
    'epvAdj', EPVS_ADJ_UM, 'minVox', MIN_VOX, 'retWeight', RET_WEIGHT, ...
    'block', BLOCK_UM, 'maxPts', MAX_ANCHORS);

GN = fieldnames(GROUPS);
assert(numel(GN) == 2, 'GROUPS must contain exactly two groups (reference first).');
ALL_SUBJ = [GROUPS.(GN{1}), GROUPS.(GN{2})];
nSh = ceil(MAX_DIST_UM / SHELL_UM);
xShell = ((1:nSh)' - 0.5) * SHELL_UM;
if strcmp(DIST_REF, 'nearest'), XL = 'Distance from nearest vessel/EPVS (\mum)';
else, XL = 'Distance from vessel wall (\mum)'; end

%% 1. Data: use workspace caa if present; load only missing subjects -------------
MAT_PATH = struct('caa6',  'caa6/caa6.mat', ...
                  'caa17', 'caa17/occip/caa17.mat', ...
                  'caa22', 'caa22/caa22.mat', ...
                  'caa25', 'caa25/caa25.mat', ...
                  'caa26', 'caa26/caa26.mat');
if evalin('base', 'exist(''caa'',''var'')')
    caa = evalin('base', 'caa');
    fprintf('Using caa from workspace.\n');
else
    caa = struct();
end
for s = 1:numel(ALL_SUBJ)
    key = regexprep(lower(ALL_SUBJ{s}), '[^a-z0-9]', '');
    have = ~isempty(fieldnames(caa)) && ...
           any(strcmp(regexprep(lower(fieldnames(caa)), '[^a-z0-9]', ''), key));
    if ~have
        fprintf('Loading %s\n', upper(key));
        tmp = load(fullfile(ROOT, MAT_PATH.(key)), key);
        caa.(key) = tmp.(key); clear tmp
    end
end

%% 2. Check mode ---------------------------------------------------------------
if strcmp(RUN_MODE, 'check')
    for s = 1:numel(ALL_SUBJ)
        for r = 1:numel(REGIONS)
            run_orient_check(get_vol(caa, ALL_SUBJ{s}, REGIONS{r}), prm, ALL_SUBJ{s}, REGIONS{r}, OUT_DIR);
        end
    end
    fprintf('\nCheck figures in %s\n', OUT_DIR);
    return
end

%% 3. Retardance floor and tertiles from pooled control WM (per region) ---------
retCfg = struct();
for r = 1:numel(REGIONS)
    reg = REGIONS{r}; vals = [];
    for s = 1:numel(GROUPS.(GN{1}))
        V = get_vol(caa, GROUPS.(GN{1}){s}, reg);
        idx = find(logical(V.mask_wm) & ~logical(V.seg) & ~logical(V.epvs));
        if numel(idx) > RET_SAMPLE, idx = idx(randperm(numel(idx), RET_SAMPLE)); end
        v = double(V.ret_full(idx)); vals = [vals; v(isfinite(v))]; %#ok<AGROW>
        clear V idx v
    end
    if RET_FLOOR_PCTL > 0, fl = prctile(vals, RET_FLOOR_PCTL); else, fl = -Inf; end
    tt = prctile(vals(vals >= fl), [100/3 200/3]);
    retCfg.(reg) = struct('floor', fl, 'edges', [fl tt Inf]);
    fprintf('%s retardance: floor %.4g (p%d of control WM), tertile edges %.4g / %.4g\n', ...
        reg, fl, RET_FLOOR_PCTL, tt);
end

%% 4. Per subject x region extraction -------------------------------------------
DAT = struct(); allA = table(); allSH = table(); subjSumm = table(); comp = table();
for r = 1:numel(REGIONS)
    reg = REGIONS{r};
    for s = 1:numel(ALL_SUBJ)
        sub = ALL_SUBJ{s}; grp = group_of(sub, GROUPS);
        fprintf('%s %s (%s) ...\n', sub, reg, grp);
        [A, SH, AN, info] = analyze_volume(get_vol(caa, sub, reg), prm, retCfg.(reg));
        A.subject  = repmat(string(sub), height(A), 1);  A.group  = repmat(string(grp), height(A), 1);
        A.region   = repmat(string(reg), height(A), 1);
        SH.subject = repmat(string(sub), height(SH), 1); SH.group = repmat(string(grp), height(SH), 1);
        SH.region  = repmat(string(reg), height(SH), 1);

        prim = A(~A.inclined, :);
        nul = zeros(height(prim) * N_NULL_SHUFFLE, 1);
        for k = 1:N_NULL_SHUFFLE
            ii = (k-1)*height(prim) + (1:height(prim));
            nul(ii) = axial_diff_deg(prim.thetaF_deg(randperm(height(prim))), prim.thetaV_deg);
        end
        DAT.(reg).(sub) = struct('A', A, 'prim', prim, 'sh', SH(~SH.inclined, :), 'shAll', SH, ...
            'AN', AN, 'null', nul, 'info', info, 'group', grp);
        fprintf('  anchors with fibers %d | valid axis %d | primary %d | inclined %d\n', ...
            info.n_anchor_total, info.n_anchor_valid, height(prim), nnz(A.inclined));

        pv = prim_stat(prim, MIN_ROWS, N_NULL_SHUFFLE);
        incl = A(A.inclined, :);
        subjSumm = [subjSumm; table(string(sub), string(grp), string(reg), ...
            info.n_anchor_total, info.n_anchor_valid, height(prim), height(incl), ...
            numel(unique(prim.block)), pv(1), median(nul), pv(2), cliffs_delta(prim.dtheta_deg, nul), ...
            pv(3), median(incl.dtheta_deg), median(prim.retMean), retCfg.(reg).floor, string(info.units), ...
            'VariableNames', {'subject','group','region','n_anchor_total','n_anchor_valid', ...
            'n_primary','n_inclined','n_blocks','median_dtheta','null_median','alignment_excess', ...
            'cliff_obs_vs_null','median_R','median_dtheta_inclined','median_ret','ret_floor','orient_units'})]; %#ok<AGROW>
        comp = [comp; table(string(sub), string(grp), string(reg), ...
            info.n_anchor_valid / max(1, info.n_anchor_total), median(AN.radius_um), ...
            mean(AN.incl_deg >= MAX_INCL_DEG), median(AN.tube_coh), mean(AN.epvAdj), ...
            'VariableNames', {'subject','group','region','frac_anchor_valid','median_radius_um', ...
            'frac_inclined','median_tube_coh','frac_epvs_adjacent'})]; %#ok<AGROW>
        allA = [allA; A]; allSH = [allSH; SH]; %#ok<AGROW>
    end
end

%% 5. Group contrasts ------------------------------------------------------------
COLS = {'block','dtheta_deg','thetaF_deg','thetaV_deg','R','shell','retBin','epvAdj'};
statPrim  = @(T) prim_stat(T, MIN_ROWS, 3);
statProfD = @(T) grouped_median(T.shell, T.dtheta_deg, nSh, MIN_ROWS);
statProfR = @(T) grouped_median(T.shell, T.R, nSh, MIN_ROWS);
statRet   = @(T) strat_excess(T, T.retBin, 3, MIN_ROWS);
statEpv   = @(T) strat_excess(T, double(T.epvAdj) + 1, 2, MIN_ROWS);

res = struct(); gc = table(); cl = table(); loso = table(); profT = table(); stratT = table();
MNAMES = {'median_dtheta','alignment_excess','median_R'};
for r = 1:numel(REGIONS)
    reg = REGIONS{r}; fprintf('Bootstrapping %s ...\n', reg);
    res.(reg).prim  = group_contrast(DAT.(reg), GROUPS, 'prim',  statPrim,  N_BOOT, COLS);
    res.(reg).profD = group_contrast(DAT.(reg), GROUPS, 'sh',    statProfD, N_BOOT, COLS);
    res.(reg).profR = group_contrast(DAT.(reg), GROUPS, 'shAll', statProfR, N_BOOT, COLS);   % R needs no vessel axis
    res.(reg).ret   = group_contrast(DAT.(reg), GROUPS, 'prim',  statRet,   N_BOOT, COLS);
    res.(reg).epv   = group_contrast(DAT.(reg), GROUPS, 'prim',  statEpv,   N_BOOT, COLS);

    P = res.(reg).prim;
    for m = 1:numel(MNAMES)
        gc = [gc; table(string(reg), string(MNAMES{m}), P.(GN{1}).pt(m), P.(GN{1}).lo(m), P.(GN{1}).hi(m), ...
            P.(GN{2}).pt(m), P.(GN{2}).lo(m), P.(GN{2}).hi(m), P.diff.pt(m), P.diff.lo(m), P.diff.hi(m), ...
            'VariableNames', {'region','metric', [GN{1} '_pt'], [GN{1} '_lo'], [GN{1} '_hi'], ...
            [GN{2} '_pt'], [GN{2} '_lo'], [GN{2} '_hi'], 'diff_pt', 'diff_lo', 'diff_hi'})]; %#ok<AGROW>
    end

    % Cliff's delta on observed dtheta (> 0: severe less aligned); pooled and per subject pair
    a = []; b = [];
    for s = 1:numel(GROUPS.(GN{1})), a = [a; DAT.(reg).(GROUPS.(GN{1}){s}).prim.dtheta_deg]; end %#ok<AGROW>
    for s = 1:numel(GROUPS.(GN{2})), b = [b; DAT.(reg).(GROUPS.(GN{2}){s}).prim.dtheta_deg]; end %#ok<AGROW>
    cl = [cl; table(string(reg), "pooled", "pooled", cliffs_delta(b, a), ...
        'VariableNames', {'region', GN{2}, GN{1}, 'cliff_dtheta'})]; %#ok<AGROW>
    for i = 1:numel(GROUPS.(GN{1}))
        for j = 1:numel(GROUPS.(GN{2}))
            ci_ = GROUPS.(GN{1}){i}; sj = GROUPS.(GN{2}){j};
            cl = [cl; table(string(reg), string(sj), string(ci_), ...
                cliffs_delta(DAT.(reg).(sj).prim.dtheta_deg, DAT.(reg).(ci_).prim.dtheta_deg), ...
                'VariableNames', cl.Properties.VariableNames)]; %#ok<AGROW>
        end
    end

    % leave-one-subject-out difference (point estimates)
    for g = 1:2
        subs = GROUPS.(GN{g});
        for s = 1:numel(subs)
            keep = setdiff(1:numel(subs), s);
            v = P.(GN{g}).subj; other = P.(GN{3-g}).pt;
            gm = mean(v(:, keep), 2, 'omitnan');
            if g == 1, d = other - gm; else, d = gm - other; end
            loso = [loso; table(string(reg), string(subs{s}), d(1), d(2), d(3), ...
                'VariableNames', {'region','left_out','diff_median_dtheta','diff_excess','diff_median_R'})]; %#ok<AGROW>
        end
    end

    profT  = [profT;  prof_table(reg, 'dtheta', res.(reg).profD, xShell, GN); ...
                      prof_table(reg, 'R',      res.(reg).profR, xShell, GN)]; %#ok<AGROW>
    stratT = [stratT; prof_table(reg, 'excess_by_ret_tertile', res.(reg).ret, (1:3)', GN); ...
                      prof_table(reg, 'excess_by_epvs_adjacency', res.(reg).epv, (0:1)', GN)]; %#ok<AGROW>
end

disp(subjSumm); disp(gc); disp(cl); disp(loso); disp(comp);
writetable(subjSumm, fullfile(OUT_DIR, 'subject_summary.csv'));
writetable(gc,       fullfile(OUT_DIR, 'group_contrast.csv'));
writetable(cl,       fullfile(OUT_DIR, 'cliff_dtheta.csv'));
writetable(loso,     fullfile(OUT_DIR, 'leave_one_subject_out.csv'));
writetable(comp,     fullfile(OUT_DIR, 'vessel_composition.csv'));
writetable(profT,    fullfile(OUT_DIR, 'distance_profiles.csv'));
writetable(stratT,   fullfile(OUT_DIR, 'stratified_excess.csv'));
writetable(allA,     fullfile(OUT_DIR, 'anchors_primary.csv'));
writetable(allSH,    fullfile(OUT_DIR, 'anchors_by_shell.csv'));

fid = fopen(fullfile(OUT_DIR, 'README.txt'), 'w');
fprintf(fid, ['%s run %s\nGroups: %s = %s; %s = %s. Contrasts are %s - %s.\n' ...
    'Vessel axis: structure tensor (sigma_d %g um, sigma_i %s um, tube coherence >= %g) at vessel-wall voxels.\n' ...
    'Distance reference: %s; shells %g um to %g um; primary window <= %g um; vessels inclined >= %g deg excluded from angle metrics.\n' ...
    'Retardance floor: p%d of pooled control WM per region (applied to all subjects); weighting = %d.\n' ...
    'dtheta in [0,90] (0 parallel, 90 perpendicular). alignment_excess = median(shuffled) - median(observed); > 0 = aligned beyond regional statistics.\n' ...
    'cliff_dtheta > 0: %s dtheta larger (less aligned) than %s.\n' ...
    'Group value = mean of subject values (equal subject weight). CIs: %g um blocks resampled within subject, subjects fixed; ' ...
    'difference CIs computed on the difference. Pseudoreplicated anchors within n = 2 subjects per group; no p-values; ' ...
    'CIs describe the measured subjects only.\nOrientation convention: zero axis %s, ccw %d, units %s. N_BOOT %d, MAX_ANCHORS %g.\n'], ...
    TAG, datestr(now), GN{1}, strjoin(GROUPS.(GN{1}), ','), GN{2}, strjoin(GROUPS.(GN{2}), ','), GN{2}, GN{1}, ...
    SIGMA_D_UM, mat2str(SIGMA_I_UM), MIN_TUBE_COH, DIST_REF, SHELL_UM, MAX_DIST_UM, PRIMARY_MAX_UM, MAX_INCL_DEG, ...
    RET_FLOOR_PCTL, RET_WEIGHT, GN{2}, GN{1}, BLOCK_UM, ORIENT_ZERO_AXIS, ORIENT_CCW, ORIENT_UNITS, N_BOOT, MAX_ANCHORS);
fclose(fid);

%% 6. Figures --------------------------------------------------------------------
GA = {GROUPS, GROUP_LABEL, GROUP_COLOR, SUBJ_MARKER};
allRad = [];
for r = 1:numel(REGIONS), for s = 1:numel(ALL_SUBJ), allRad = [allRad; DAT.(REGIONS{r}).(ALL_SUBJ{s}).AN.radius_um]; end, end %#ok<AGROW>
radLims = [0, max(VOX_UM, prctile(allRad, 99))];

main = struct('title', {'\Delta\theta: observed vs shuffled', 'Alignment excess', ...
                        '\Delta\theta vs distance', 'Fiber coherence vs distance'}, ...
              'fn', {@(ax, reg) draw_violin_groups(ax, DAT.(reg), GA{:}), ...
                     @(ax, reg) draw_point_contrast(ax, res.(reg).prim, 2, GA{:}, 'Alignment excess (deg)'), ...
                     @(ax, reg) draw_profile(ax, res.(reg).profD, xShell, GA{:}, '\Delta\theta vessel–fiber (deg)', XL, [0 90], PRIMARY_MAX_UM), ...
                     @(ax, reg) draw_profile(ax, res.(reg).profR, xShell, GA{:}, 'Fiber coherence R', XL, [0 1], PRIMARY_MAX_UM)}, ...
              'tag', {'violin','excess','profile_dtheta','profile_R'});
build_figure('main_orient_vessel_groups', REGIONS, REG_LABEL, main, OUT_DIR);

s1 = struct('title', {'\Delta\theta: severe − control', 'R: severe − control'}, ...
            'fn', {@(ax, reg) draw_diff_profile(ax, res.(reg).profD, xShell, '\Delta\Delta\theta (deg)', XL, PRIMARY_MAX_UM), ...
                   @(ax, reg) draw_diff_profile(ax, res.(reg).profR, xShell, '\DeltaR', XL, PRIMARY_MAX_UM)}, ...
            'tag', {'diff_dtheta','diff_R'});
build_figure('supp1_difference_profiles', REGIONS, REG_LABEL, s1, OUT_DIR);

s2 = struct('title', {'Alignment excess by retardance tertile', 'Alignment excess by EPVS adjacency'}, ...
            'fn', {@(ax, reg) draw_strat(ax, res.(reg).ret, {'Low ret','Mid ret','High ret'}, GA{:}, 'Alignment excess (deg)'), ...
                   @(ax, reg) draw_strat(ax, res.(reg).epv, {'No EPVS','EPVS-adjacent'}, GA{:}, 'Alignment excess (deg)')}, ...
            'tag', {'strat_ret','strat_epvs'});
build_figure('supp2_stratified_excess', REGIONS, REG_LABEL, s2, OUT_DIR);

s3 = struct('title', {'Vessel radius', 'Vessel inclination', 'Tube coherence'}, ...
            'fn', {@(ax, reg) draw_comp(ax, DAT.(reg), 'radius_um', radLims, GA{:}, 'Radius (\mum)'), ...
                   @(ax, reg) draw_comp(ax, DAT.(reg), 'incl_deg', [0 90], GA{:}, 'Inclination from en-face plane (deg)'), ...
                   @(ax, reg) draw_comp(ax, DAT.(reg), 'tube_coh', [MIN_TUBE_COH 1], GA{:}, 'Tube coherence')}, ...
            'tag', {'comp_radius','comp_incl','comp_coh'});
build_figure('supp3_vessel_composition', REGIONS, REG_LABEL, s3, OUT_DIR);

fprintf('Done. Output: %s\n', OUT_DIR);

%% ========================= Local functions ==================================

function V = get_vol(caa, sub, reg)
    sf = resolve_field(caa, sub, true);
    rf = resolve_field(caa.(sf), reg, false);
    V = caa.(sf).(rf);
end

function g = group_of(sub, GROUPS)
    gn = fieldnames(GROUPS); g = '';
    for k = 1:numel(gn), if any(strcmp(GROUPS.(gn{k}), sub)), g = gn{k}; return; end, end
    error('Subject %s not in GROUPS', sub);
end

function [A, SH, AN, info] = analyze_volume(V, prm, rc)
    wm = logical(V.mask_wm);
    pad = ceil(max(prm.maxWall, 3 * max(prm.sigmaI)) / prm.vox) + 2;
    [r1, r2, r3] = bbox_ranges(wm, pad);
    wm  = wm(r1,r2,r3);
    seg = logical(V.seg(r1,r2,r3));
    epv = logical(V.epvs(r1,r2,r3));
    [th, info] = decode_orient(V.orient(r1,r2,r3), prm.units, wm);
    th  = single(th);
    ret = single(V.ret_full(r1,r2,r3));

    % fiber voxels, shell distance, and nearest vessel-wall anchor
    [dWall, idxW] = bwdist(seg);
    dWall = dWall * prm.vox;
    if strcmp(prm.distRef, 'nearest'), dRef = bwdist(seg | epv) * prm.vox; else, dRef = dWall; end
    fib = wm & ~seg & ~epv & isfinite(th) & isfinite(ret) & ret >= rc.floor & ...
          dRef > 0 & dRef <= prm.maxDist & dWall <= prm.maxWall;
    clear dWall
    vox = find(fib); clear fib
    nSh = ceil(prm.maxDist / prm.shell);
    shV = min(nSh, ceil(double(dRef(vox)) / prm.shell)); clear dRef
    [aIdx, ~, pid] = unique(double(idxW(vox))); clear idxW
    nA = numel(aIdx);
    info.n_anchor_total = nA;

    % vessel axis, radius, EPVS adjacency at anchors
    vt  = vessel_tensor(seg, aIdx, prm);
    t   = vt.t;
    phi = rad2deg(asin(min(1, abs(t(:,3)))));
    thV = rad2deg(display_to_conv(atan2(-t(:,1), t(:,2)), prm.zeroAxis, prm.ccw));
    okA = isfinite(vt.coh) & vt.coh >= prm.minCoh;
    info.n_anchor_valid = nnz(okA);

    dIn  = bwdist(~seg);
    k    = max(1, round(prm.radSearch / prm.vox));
    rmap = imdilate(dIn, strel('cube', 2*k + 1)); clear dIn
    rad  = double(rmap(aIdx)) * prm.vox; clear rmap
    if any(epv(:))
        dE = bwdist(epv); dEpv = double(dE(aIdx)) * prm.vox; clear dE
    else
        dEpv = inf(nA, 1);
    end
    epvAdj = dEpv <= prm.epvAdj;

    keepA = okA;
    if isfinite(prm.maxPts) && nnz(okA) > prm.maxPts
        fA = find(okA); keepA(:) = false;
        keepA(fA(randperm(numel(fA), prm.maxPts))) = true;
    end
    kv = keepA(pid); vox = vox(kv); pid = pid(kv); shV = shV(kv);
    tv = double(th(vox)); rv = double(ret(vox)); clear th ret
    if prm.retWeight, w = max(rv, 0); else, w = ones(size(rv)); end

    key = (pid - 1) * nSh + shV;
    sz  = [nA * nSh, 1];
    acc = @(x) reshape(accumarray(key, x, sz), nSh, nA)';
    C  = acc(w .* cos(2*tv)); S = acc(w .* sin(2*tv)); W = acc(w);
    N  = acc(ones(size(tv)));  RT = acc(rv);

    [a1, a2, a3] = ind2sub(size(seg), aIdx);
    B   = max(1, round(prm.block / prm.vox));
    blk = sub2ind(ceil(size(seg) / B), ceil(a1 / B), ceil(a2 / B), ceil(a3 / B));
    xyz = [a1 + r1(1) - 1, a2 + r2(1) - 1, a3 + r3(1) - 1];
    incl = phi >= prm.maxIncl;

    % primary window
    kP = round(prm.primaryMax / prm.shell);
    Cp = sum(C(:,1:kP), 2); Sp = sum(S(:,1:kP), 2); Wp = sum(W(:,1:kP), 2);
    Np = sum(N(:,1:kP), 2); RTp = sum(RT(:,1:kP), 2);
    ok  = keepA & Np >= prm.minVox;
    thF = mod(rad2deg(0.5 * atan2(Sp, Cp)), 180);
    Rp  = hypot(Cp, Sp) ./ Wp;
    retM = RTp ./ Np;
    A = table(find(ok), blk(ok), xyz(ok,1), xyz(ok,2), xyz(ok,3), thV(ok), phi(ok), rad(ok), ...
        vt.coh(ok), vt.sigma(ok), dEpv(ok), epvAdj(ok), thF(ok), Rp(ok), Np(ok), retM(ok), ...
        discretize(retM(ok), rc.edges), axial_diff_deg(thF(ok), thV(ok)), incl(ok), ...
        'VariableNames', {'anchor','block','x1','x2','x3','thetaV_deg','incl_deg','radius_um', ...
        'tube_coh','sigmaI_um','dist_epvs_um','epvAdj','thetaF_deg','R','nvox','retMean', ...
        'retBin','dtheta_deg','inclined'});

    % per shell
    [pp, ss] = find(N >= prm.minVox);
    li = sub2ind(size(N), pp, ss);
    thFs = mod(rad2deg(0.5 * atan2(S(li), C(li))), 180);
    retS = RT(li) ./ N(li);
    SH = table(pp, blk(pp), ss, (ss - 0.5) * prm.shell, thV(pp), phi(pp), rad(pp), epvAdj(pp), ...
        thFs, hypot(C(li), S(li)) ./ W(li), N(li), retS, discretize(retS, rc.edges), ...
        axial_diff_deg(thFs, thV(pp)), incl(pp), ...
        'VariableNames', {'anchor','block','shell','dist_um','thetaV_deg','incl_deg','radius_um', ...
        'epvAdj','thetaF_deg','R','nvox','retMean','retBin','dtheta_deg','inclined'});

    % composition: all valid (kept) anchors, independent of fiber coverage
    AN = table(rad(keepA), phi(keepA), vt.coh(keepA), epvAdj(keepA), ...
        'VariableNames', {'radius_um','incl_deg','tube_coh','epvAdj'});
end

function vt = vessel_tensor(seg, aIdx, prm)
    Is = imgaussfilt3(single(seg), prm.sigmaD / prm.vox);
    [g2, g1, g3] = gradient(Is); clear Is          % gradient(): dim2, dim1, dim3
    G = {g1, g2, g3}; clear g1 g2 g3
    pairs = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
    nA = numel(aIdx);
    vt = struct('t', nan(nA, 3), 'coh', -inf(nA, 1), 'sigma', nan(nA, 1));
    for k = 1:numel(prm.sigmaI)
        Jv = zeros(nA, 6);
        for p = 1:6
            Pk = imgaussfilt3(G{pairs(p,1)} .* G{pairs(p,2)}, prm.sigmaI(k) / prm.vox);
            Jv(:, p) = double(Pk(aIdx));
        end
        clear Pk
        [t, coh] = tensor_axis(Jv);
        upd = coh > vt.coh;
        vt.t(upd, :) = t(upd, :); vt.coh(upd) = coh(upd); vt.sigma(upd) = prm.sigmaI(k);
    end
    vt.coh(~isfinite(vt.coh)) = NaN;
end

function [t, coh] = tensor_axis(Jv)
    n = size(Jv, 1);
    J = zeros(3, 3, n);
    J(1,1,:) = Jv(:,1); J(1,2,:) = Jv(:,2); J(1,3,:) = Jv(:,3);
    J(2,1,:) = Jv(:,2); J(2,2,:) = Jv(:,4); J(2,3,:) = Jv(:,5);
    J(3,1,:) = Jv(:,3); J(3,2,:) = Jv(:,5); J(3,3,:) = Jv(:,6);
    if exist('pageeig', 'builtin') || exist('pageeig', 'file')
        [Vp, Lp] = pageeig(J);
        mu = [squeeze(Lp(1,1,:)), squeeze(Lp(2,2,:)), squeeze(Lp(3,3,:))];
        [mu, o] = sort(real(mu), 2, 'ascend');
        base = (0:n-1)' * 9 + (o(:,1) - 1) * 3;
        t = real(Vp(base + [1 2 3]));
    else
        t = nan(n, 3); mu = nan(n, 3);
        for i = 1:n
            [V, L] = eig(J(:,:,i));
            [mi, o] = sort(diag(L), 'ascend');
            mu(i,:) = mi'; t(i,:) = V(:, o(1))';
        end
    end
    coh = (mu(:,2) - mu(:,1)) ./ (mu(:,2) + mu(:,1) + eps);
    coh(all(Jv == 0, 2)) = NaN;
end

function [th, info] = decode_orient(o, unitsOpt, wm)
    o = double(o);
    v = o(wm & isfinite(o));
    info.raw_min = min(v); info.raw_max = max(v);
    if strcmp(unitsOpt, 'auto')
        if max(abs(v)) > 2*pi + 1e-6, info.units = 'deg'; else, info.units = 'rad'; end
    else
        info.units = unitsOpt;
    end
    if strcmp(info.units, 'deg'), th = deg2rad(o); else, th = o; end
    th = mod(th, pi);
end

function th = display_to_conv(thDisp, zeroAxis, ccw)
    th0 = 0; if strcmp(zeroAxis, 'dim1'), th0 = -pi/2; end
    sgn = 1; if ~ccw, sgn = -1; end
    th = mod(sgn * (thDisp - th0), pi);
end

function thDisp = conv_to_display(th, zeroAxis, ccw)
    th0 = 0; if strcmp(zeroAxis, 'dim1'), th0 = -pi/2; end
    sgn = 1; if ~ccw, sgn = -1; end
    thDisp = sgn * th + th0;
end

function d = axial_diff_deg(a, b)
    d = abs(mod(a - b + 90, 180) - 90);
end

function d = cliffs_delta(x, y)
    x = x(:); y = y(:); nx = numel(x); ny = numel(y);
    if nx == 0 || ny == 0, d = NaN; return; end
    rk = tiedrank([x; y]);
    U = sum(rk(1:nx)) - nx * (nx + 1) / 2;
    d = (2 * U - nx * ny) / (nx * ny);
end

% ---- statistics ----

function v = prim_stat(T, minN, K)
    % [median dtheta; alignment excess; median R]
    n = height(T); v = nan(3, 1);
    if n < minN, return; end
    md = median(T.dtheta_deg); nm = zeros(K, 1);
    for k = 1:K, nm(k) = median(axial_diff_deg(T.thetaF_deg(randperm(n)), T.thetaV_deg)); end
    v = [md; mean(nm) - md; median(T.R)];
end

function v = strat_excess(T, bins, nb, minN)
    % alignment excess within each stratum (null shuffled within stratum)
    v = nan(nb, 1);
    for k = 1:nb
        i = bins == k; n = nnz(i);
        if n < minN, continue; end
        thF = T.thetaF_deg(i); thV = T.thetaV_deg(i);
        v(k) = median(axial_diff_deg(thF(randperm(n)), thV)) - median(T.dtheta_deg(i));
    end
end

function m = grouped_median(g, x, nG, minN)
    ok = isfinite(x) & g >= 1 & g <= nG; g = double(g(ok)); x = double(x(ok));
    m = nan(nG, 1);
    if isempty(x), return; end
    [~, o] = sortrows([g x]); xs = x(o);
    cnt = accumarray(g, 1, [nG 1]); st = [0; cumsum(cnt(1:end-1))];
    for k = 1:nG
        c = cnt(k); if c < minN, continue; end
        i0 = st(k);
        if mod(c, 2), m(k) = xs(i0 + (c+1)/2); else, m(k) = 0.5 * (xs(i0 + c/2) + xs(i0 + c/2 + 1)); end
    end
end

function out = group_contrast(Sreg, GROUPS, fld, statFun, nBoot, cols)
    gn = fieldnames(GROUPS);
    for g = 1:numel(gn)
        subs = GROUPS.(gn{g});
        Tc = cell(1, numel(subs));
        for s = 1:numel(subs)
            T = Sreg.(subs{s}).(fld);
            Tc{s} = T(:, intersect(cols, T.Properties.VariableNames, 'stable'));
        end
        [p, reps, sv] = group_boot(Tc, statFun, nBoot);
        out.(gn{g}) = struct('pt', p, 'lo', prctile(reps, 2.5, 2), 'hi', prctile(reps, 97.5, 2), ...
            'reps', reps, 'subj', sv);
        out.(gn{g}).subs = subs;
    end
    d = out.(gn{2}).reps - out.(gn{1}).reps;
    out.diff = struct('pt', out.(gn{2}).pt - out.(gn{1}).pt, ...
        'lo', prctile(d, 2.5, 2), 'hi', prctile(d, 97.5, 2));
end

function [pt, reps, sv] = group_boot(Tc, statFun, nBoot)
    % group statistic = mean over subjects of statFun(subject table); blocks resampled within subject
    nS = numel(Tc);
    v1 = statFun(Tc{1}); sv = nan(numel(v1), nS); sv(:, 1) = v1;
    for s = 2:nS, sv(:, s) = statFun(Tc{s}); end
    pt = mean(sv, 2, 'omitnan');
    cells = cell(1, nS);
    for s = 1:nS
        if height(Tc{s}) == 0, cells{s} = {}; continue; end
        [~, ~, gi] = unique(Tc{s}.block);
        cells{s} = accumarray(gi, (1:height(Tc{s}))', [], @(v) {v});
    end
    reps = nan(numel(pt), nBoot);
    for b = 1:nBoot
        accS = nan(numel(pt), nS);
        for s = 1:nS
            c = cells{s};
            if isempty(c), continue; end
            pick = c(randi(numel(c), numel(c), 1));
            accS(:, s) = statFun(Tc{s}(vertcat(pick{:}), :));
        end
        reps(:, b) = mean(accS, 2, 'omitnan');
    end
end

function T = prof_table(reg, metric, P, x, gn)
    n = numel(x);
    T = table(repmat(string(reg), n, 1), repmat(string(metric), n, 1), x(:), ...
        P.(gn{1}).pt, P.(gn{1}).lo, P.(gn{1}).hi, P.(gn{2}).pt, P.(gn{2}).lo, P.(gn{2}).hi, ...
        P.diff.pt, P.diff.lo, P.diff.hi, ...
        'VariableNames', {'region','metric','x', [gn{1} '_pt'], [gn{1} '_lo'], [gn{1} '_hi'], ...
        [gn{2} '_pt'], [gn{2} '_lo'], [gn{2} '_hi'], 'diff_pt', 'diff_lo', 'diff_hi'});
end

function [r1, r2, r3] = bbox_ranges(m, pad)
    sz = size(m);
    a1 = find(any(any(m, 2), 3)); a2 = find(any(any(m, 1), 3)); a3 = find(any(any(m, 1), 2));
    r1 = max(1, a1(1)-pad):min(sz(1), a1(end)+pad);
    r2 = max(1, a2(1)-pad):min(sz(2), a2(end)+pad);
    r3 = max(1, a3(1)-pad):min(sz(3), a3(end)+pad);
end

function f = resolve_field(s, key, exact)
    fn = fieldnames(s);
    nk = regexprep(lower(key), '[^a-z0-9]', '');
    nf = regexprep(lower(fn), '[^a-z0-9]', '');
    hit = find(strcmp(nf, nk));
    if isempty(hit) && ~exact, hit = find(startsWith(nf, nk)); end
    if numel(hit) ~= 1
        error('Could not uniquely resolve "%s" among: %s', key, strjoin(fn', ', '));
    end
    f = fn{hit};
end

% ---- figures ----

function build_figure(name, regs, REG_LABEL, panels, outDir)
    nR = numel(regs); nP = numel(panels);
    f = figure('Color', 'w', 'Position', [40 40 480*nP 430*nR]);
    tl = tiledlayout(f, nR, nP, 'TileSpacing', 'compact', 'Padding', 'compact');
    axs = gobjects(nR, nP);
    for r = 1:nR
        for p = 1:nP
            axs(r, p) = nexttile(tl);
            panels(p).fn(axs(r, p), regs{r});
            title(axs(r, p), sprintf('%s — %s', REG_LABEL.(regs{r}), panels(p).title));
        end
    end
    yl = cell(1, nP);
    for p = 1:nP, linkaxes(axs(:, p), 'y'); yl{p} = ylim(axs(1, p)); end   % matched y across regions
    exportgraphics(f, fullfile(outDir, [name '.png']), 'Resolution', 300);
    savefig(f, fullfile(outDir, [name '.fig']));
    for r = 1:nR
        for p = 1:nP
            fs = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 580 470]);
            ax = axes(fs);
            panels(p).fn(ax, regs{r}); ylim(ax, yl{p});
            title(ax, sprintf('%s — %s', REG_LABEL.(regs{r}), panels(p).title));
            fn = sprintf('%s_%s_%s', name, regs{r}, panels(p).tag);
            exportgraphics(fs, fullfile(outDir, [fn '.png']), 'Resolution', 300);
            savefig(fs, fullfile(outDir, [fn '.fig'])); close(fs);
        end
    end
end

function c = lighten(c, a)
    c = c + (1 - c) * a;
end

function add_subj_legend(ax, GROUPS, MK)
    gn = fieldnames(GROUPS); subs = {};
    for g = 1:numel(gn), subs = [subs, GROUPS.(gn{g})]; end %#ok<AGROW>
    h = gobjects(numel(subs), 1);
    for s = 1:numel(subs)
        h(s) = plot(ax, nan, nan, MK.(subs{s}), 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
    end
    legend(ax, h, subs, 'Location', 'best', 'Box', 'off');
end

function draw_violin(ax, x, y, col, lims)
    y = y(isfinite(y)); if numel(y) < 5, return; end
    y = min(max(y, lims(1)), lims(2));
    yi = linspace(lims(1), lims(2), 200);
    pd = 1e-6 * max(1, diff(lims));
    f = ksdensity(y, yi, 'Support', [lims(1) - pd, lims(2) + pd], 'BoundaryCorrection', 'reflection');
    f = 0.35 * f / max(f);
    fill(ax, [x - f, fliplr(x + f)], [yi, fliplr(yi)], col, 'FaceAlpha', 0.6, 'EdgeColor', col * 0.7);
    q = prctile(y, [25 50 75]);
    plot(ax, [x x], q([1 3]), 'k-', 'LineWidth', 4);
    plot(ax, x, q(2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
end

function draw_violin_groups(ax, Sreg, GROUPS, GLAB, GC, MK)
    hold(ax, 'on'); box(ax, 'off');
    gn = fieldnames(GROUPS); pos = 0; xs = []; labs = {};
    for g = 1:numel(gn)
        subs = GROUPS.(gn{g}); col = GC.(gn{g}); colN = lighten(col, 0.6);
        obs = []; nul = []; mo = nan(1, numel(subs)); mn = mo;
        for s = 1:numel(subs)
            D = Sreg.(subs{s});
            obs = [obs; D.prim.dtheta_deg]; nul = [nul; D.null]; %#ok<AGROW>
            mo(s) = median(D.prim.dtheta_deg); mn(s) = median(D.null);
        end
        pos = pos + 1; draw_violin(ax, pos, obs, col, [0 90]);
        for s = 1:numel(subs), plot(ax, pos + 0.18, mo(s), MK.(subs{s}), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 7); end
        xs(end+1) = pos; labs{end+1} = GLAB.(gn{g}); %#ok<AGROW>
        pos = pos + 1; draw_violin(ax, pos, nul, colN, [0 90]);
        for s = 1:numel(subs), plot(ax, pos + 0.18, mn(s), MK.(subs{s}), 'MarkerFaceColor', colN, 'MarkerEdgeColor', 'k', 'MarkerSize', 7); end
        xs(end+1) = pos; labs{end+1} = [GLAB.(gn{g}) ' shuffled']; %#ok<AGROW>
        pos = pos + 0.5;
    end
    yline(ax, 45, 'k:');
    set(ax, 'XTick', xs, 'XTickLabel', labs); xlim(ax, [0.4 pos + 0.1]);
    ylim(ax, [0 90]); yticks(ax, 0:15:90); ylabel(ax, '\Delta\theta vessel–fiber (deg)');
    add_subj_legend(ax, GROUPS, MK);
end

function draw_point_contrast(ax, P, idx, GROUPS, GLAB, GC, MK, ylab)
    hold(ax, 'on'); box(ax, 'off');
    gn = fieldnames(GROUPS);
    for g = 1:numel(gn)
        e = P.(gn{g}); col = GC.(gn{g});
        errorbar(ax, g, e.pt(idx), e.pt(idx) - e.lo(idx), e.hi(idx) - e.pt(idx), 'o', 'Color', col, ...
            'MarkerFaceColor', col, 'MarkerSize', 9, 'LineWidth', 1.8, 'CapSize', 10);
        for s = 1:numel(e.subs)
            plot(ax, g + 0.2, e.subj(idx, s), MK.(e.subs{s}), 'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', col, 'MarkerSize', 7, 'LineWidth', 1.2);
        end
    end
    xd = numel(gn) + 1.3;
    errorbar(ax, xd, P.diff.pt(idx), P.diff.pt(idx) - P.diff.lo(idx), P.diff.hi(idx) - P.diff.pt(idx), 'd', ...
        'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 9, 'LineWidth', 1.8, 'CapSize', 10);
    yline(ax, 0, 'k:');
    labs = [cellfun(@(x) GLAB.(x), gn', 'UniformOutput', false), {sprintf('%s − %s', GLAB.(gn{2}), GLAB.(gn{1}))}];
    set(ax, 'XTick', [1:numel(gn), xd], 'XTickLabel', labs); xlim(ax, [0.4 xd + 0.6]);
    ylabel(ax, ylab);
    add_subj_legend(ax, GROUPS, MK);
end

function draw_profile(ax, P, x, GROUPS, GLAB, GC, MK, ylab, xl, yl, pm)
    hold(ax, 'on'); box(ax, 'off');
    gn = fieldnames(GROUPS); h = gobjects(numel(gn), 1);
    for g = 1:numel(gn)
        e = P.(gn{g}); col = GC.(gn{g});
        ok = isfinite(e.pt) & isfinite(e.lo) & isfinite(e.hi);
        fill(ax, [x(ok); flipud(x(ok))], [e.lo(ok); flipud(e.hi(ok))], col, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
        h(g) = plot(ax, x(ok), e.pt(ok), '-', 'Color', col, 'LineWidth', 2);
        for s = 1:numel(e.subs)
            plot(ax, x, e.subj(:, s), ':', 'Color', col, 'Marker', MK.(e.subs{s}), ...
                'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth', 0.8);
        end
    end
    xline(ax, pm, ':', 'Color', [0.4 0.4 0.4]);
    if ~isempty(yl), ylim(ax, yl); end
    xlabel(ax, xl); ylabel(ax, ylab);
    legend(ax, h, cellfun(@(x) GLAB.(x), gn, 'UniformOutput', false), 'Location', 'best', 'Box', 'off');
end

function draw_diff_profile(ax, P, x, ylab, xl, pm)
    hold(ax, 'on'); box(ax, 'off');
    d = P.diff; ok = isfinite(d.pt) & isfinite(d.lo) & isfinite(d.hi);
    fill(ax, [x(ok); flipud(x(ok))], [d.lo(ok); flipud(d.hi(ok))], [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(ax, x(ok), d.pt(ok), 'k-', 'LineWidth', 2);
    yline(ax, 0, 'k:'); xline(ax, pm, ':', 'Color', [0.4 0.4 0.4]);
    xlabel(ax, xl); ylabel(ax, ylab);
end

function draw_strat(ax, P, labels, GROUPS, GLAB, GC, MK, ylab)
    hold(ax, 'on'); box(ax, 'off');
    gn = fieldnames(GROUPS); nb = numel(labels); off = linspace(-0.25, 0, numel(gn));
    h = gobjects(numel(gn) + 1, 1);
    for g = 1:numel(gn)
        e = P.(gn{g}); col = GC.(gn{g}); x = (1:nb)' + off(g);
        h(g) = errorbar(ax, x, e.pt, e.pt - e.lo, e.hi - e.pt, 'o', 'Color', col, 'MarkerFaceColor', col, ...
            'MarkerSize', 8, 'LineWidth', 1.6, 'CapSize', 8, 'LineStyle', 'none');
        for s = 1:numel(e.subs)
            plot(ax, x + 0.07, e.subj(:, s), MK.(e.subs{s}), 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col, ...
                'MarkerSize', 6, 'LineStyle', 'none');
        end
    end
    xd = (1:nb)' + 0.3;
    h(end) = errorbar(ax, xd, P.diff.pt, P.diff.pt - P.diff.lo, P.diff.hi - P.diff.pt, 'd', 'Color', 'k', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.6, 'CapSize', 8, 'LineStyle', 'none');
    yline(ax, 0, 'k:');
    set(ax, 'XTick', 1:nb, 'XTickLabel', labels); xlim(ax, [0.5 nb + 0.6]); ylabel(ax, ylab);
    legend(ax, h, [cellfun(@(x) GLAB.(x), gn', 'UniformOutput', false), ...
        {sprintf('%s − %s', GLAB.(gn{2}), GLAB.(gn{1}))}], 'Location', 'best', 'Box', 'off');
end

function draw_comp(ax, Sreg, var, lims, GROUPS, GLAB, GC, MK, ylab)
    hold(ax, 'on'); box(ax, 'off');
    gn = fieldnames(GROUPS);
    for g = 1:numel(gn)
        subs = GROUPS.(gn{g}); col = GC.(gn{g}); y = []; md = nan(1, numel(subs));
        for s = 1:numel(subs)
            v = Sreg.(subs{s}).AN.(var); y = [y; v]; md(s) = median(v, 'omitnan'); %#ok<AGROW>
        end
        draw_violin(ax, g, y, col, lims);
        for s = 1:numel(subs)
            plot(ax, g + 0.18, md(s), MK.(subs{s}), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
        end
    end
    set(ax, 'XTick', 1:numel(gn), 'XTickLabel', cellfun(@(x) GLAB.(x), gn', 'UniformOutput', false));
    xlim(ax, [0.4 numel(gn) + 0.6]); ylim(ax, lims); ylabel(ax, ylab);
    add_subj_legend(ax, GROUPS, MK);
end

% ---- check mode ----

function run_orient_check(V, prm, sub, reg, outDir)
    wmAll = logical(V.mask_wm);
    [~, z] = max(squeeze(sum(sum(wmAll, 1), 2)));
    zpad = ceil(3 * max(prm.sigmaI) / prm.vox);
    zs = max(1, z - zpad):min(size(wmAll, 3), z + zpad); zc = z - zs(1) + 1;
    wm = wmAll(:,:,zs); clear wmAll
    o  = double(V.orient(:,:,zs)); v = o(wm & isfinite(o));
    [th, info] = decode_orient(o, prm.units, wm); clear o
    rt = double(V.ret_full(:,:,z)); rv = rt(wm(:,:,zc) & isfinite(rt));
    fprintf('\n%s %s (slab around z = %d)\n', sub, reg, z);
    fprintf('  orient range [%.4g, %.4g], p1/p50/p99 = %.4g/%.4g/%.4g, frac == 0: %.3f -> %s\n', ...
        info.raw_min, info.raw_max, prctile(v, [1 50 99]), mean(v == 0), info.units);
    fprintf('  ret_full (slice) p1/p10/p50/p99 = %.4g/%.4g/%.4g/%.4g\n', prctile(rv, [1 10 50 99]));

    segS = logical(V.seg(:,:,zs));
    wall = segS(:,:,zc) & ~imerode(segS(:,:,zc), true(3));
    [w1, w2] = find(wall); w1 = w1(1:4:end); w2 = w2(1:4:end);
    aIdx = sub2ind(size(segS), w1, w2, zc * ones(size(w1)));
    vt = vessel_tensor(segS, aIdx, prm);
    good = vt.coh >= prm.minCoh & abs(vt.t(:,3)) < sind(prm.maxIncl);
    fprintf('  %d wall anchors sampled, %d pass coherence + inclination\n', numel(aIdx), nnz(good));

    f = figure('Color', 'w', 'Position', [50 50 1500 700]);
    tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
    ax = nexttile(tl); histogram(ax, v, 180); title(ax, sprintf('%s %s orient (WM slab)', sub, reg));
    xlabel(ax, sprintf('value (%s)', info.units));
    ax = nexttile(tl, [1 2]);
    if isfield(V, 'mus'), bg = double(V.mus(:,:,z)); else, bg = double(wm(:,:,zc)); end
    imagesc(ax, bg); axis(ax, 'image'); colormap(ax, gray);
    lim = prctile(bg(wm(:,:,zc)), [2 98]); if numel(lim) == 2 && diff(lim) > 0, clim(ax, lim); end
    hold(ax, 'on');
    contour(ax, double(segS(:,:,zc)), [0.5 0.5], 'r', 'LineWidth', 0.5);
    B = 10; [n1, n2] = size(bg); ths = th(:,:,zc); wms = wm(:,:,zc);
    for a = 1:B:n1 - B + 1
        for b = 1:B:n2 - B + 1
            blk = ths(a:a+B-1, b:b+B-1); m = wms(a:a+B-1, b:b+B-1) & isfinite(blk);
            if nnz(m) < 0.5 * B^2, continue; end
            tm = 0.5 * atan2(mean(sin(2*blk(m))), mean(cos(2*blk(m))));
            draw_glyph(ax, b + B/2, a + B/2, conv_to_display(tm, prm.zeroAxis, prm.ccw), 0.8*B, 'c');
        end
    end
    for k = find(good)'
        draw_glyph(ax, w2(k), w1(k), atan2(-vt.t(k,1), vt.t(k,2)), 8, [1 0.55 0.1]);
    end
    title(ax, sprintf('z = %d | cyan: optic axis (zero=%s, ccw=%d) | orange: vessel axis', z, prm.zeroAxis, prm.ccw));
    exportgraphics(f, fullfile(outDir, sprintf('check_%s_%s.png', sub, reg)), 'Resolution', 200);
    savefig(f, fullfile(outDir, sprintf('check_%s_%s.fig', sub, reg)));
end

function draw_glyph(ax, x, y, thDisp, L, col)
    dx = 0.5 * L * cos(thDisp); dy = -0.5 * L * sin(thDisp);
    plot(ax, [x - dx, x + dx], [y - dy, y + dy], '-', 'Color', col, 'LineWidth', 1.2);
end