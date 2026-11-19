%% compare_orient_vessel_controls.m
% Fiber (as-PSOCT optic axis) vs vessel orientation in white matter of the control
% subjects (CAA26, CAA6), frontal and occipital.
%
%   RUN_MODE = 'check' : orient/retardance diagnostics + glyph overlays (run first)
%   RUN_MODE = 'full'  : all analyses, figures, CSVs
%
% Design: 2 subjects x 2 regions. Every estimate is computed within one subject x region
% volume with a spatial block bootstrap; frontal - occipital contrasts are within subject.
% A result is "replicated" when both subjects give the same sign with 95% CIs excluding 0.
% CIs describe these specimens only (no population inference with n = 2).
%
% Geometry: dims 1-2 = en-face plane, dim 3 = through-plane. The optic axis is an in-plane
% axial angle; vessel axes are projected onto the en-face plane, and vessels inclined
% >= MAX_INCL_DEG are excluded from angle-based metrics. Vessel axis = 3-D structure tensor
% of the smoothed vessel mask at vessel-wall voxels (integration scale chosen per anchor by
% max tube coherence). Each perivascular WM voxel is anchored to its nearest wall voxel.
% excess = median(null dtheta) - median(observed dtheta); > 0 = aligned beyond the null.
%
% Analyses
%  1. Alignment vs three nulls (per volume)
%       global  : fiber angles shuffled across all anchors (regional statistics only)
%       block   : fiber angles shuffled among anchors in the same BLOCK_UM block
%       spatial : smoothed fiber field rigidly shifted in-plane by SHIFT_UM relative to the
%                 vessel network (autocorrelation of both fields preserved). Excess vs shift
%                 distance = spatial scale of vessel-fiber coupling; compare with the fiber
%                 orientation correlation length (variogram, analysis 7).
%  2. Frontal - occipital contrasts within subject (alignment, coherence, retardance, decay).
%  3. Composition: occipital anchors reweighted to the frontal radius x inclination
%     distribution and the contrast recomputed; excess by radius and by inclination class.
%  4. Distance decay: excess per shell; exponential fit (lambda) and half-decay distance.
%  5. Depth from the GM/WM boundary: excess, fraction of radial vessels (axis within
%     RADIAL_MAX_DEG of the cortical normal), dtheta of radial vessels, and fiber angle to
%     the in-plane cortical normal. Radial vessels with dtheta near 90 deg in deep WM = the
%     perpendicular vessel-fiber geometry assumed by DTI-ALPS.
%  6. Vessel-axis-free metrics: fiber coherence R and retardance vs distance; near - far.
%  7. Robustness: orientation variogram (correlation length; block-size check), split-half
%     reliability (checkerboard of blocks), sensitivity of the excess to tube-coherence
%     threshold, integration scale, retardance floor and block size.

%% 0. Config ---------------------------------------------------------------
RUN_MODE    = 'full';        % 'check' | 'full'
EXPLORATORY = true;

SUBJECTS    = {'CAA26','CAA6'};
REGIONS     = {'front','occip'};               % contrasts are front - occip
REG_LABEL   = struct('front', 'Frontal', 'occip', 'Occipital');
REG_COLOR   = struct('front', [0.20 0.45 0.75], 'occip', [0.85 0.45 0.10]);
SUBJ_MARKER = {'o','s'};
SUBJ_LS     = {'-','--'};
SUBJ_FACE   = {'k','w'};

VOX_UM = 20;

ORIENT_UNITS     = 'deg';
ORIENT_ZERO_AXIS = 'dim2';    % confirm with RUN_MODE = 'check'
ORIENT_CCW       = true;

% vessel axis (structure tensor)
SIGMA_D_UM    = 20;
SIGMA_I_UM    = [40 80 160];
MIN_TUBE_COH  = 0.5;
RAD_SEARCH_UM = 160;

% geometry
DIST_REF       = 'nearest';   % 'nearest' (vessel or EPVS) | 'vessel'
SHELL_UM       = 40;
MAX_DIST_UM    = 400;
MAX_WALL_UM    = 600;
PRIMARY_MAX_UM = 120;
FAR_MIN_UM     = 280;         % shells starting at >= this distance form the "far" zone
MAX_INCL_DEG   = 45;
EPVS_ADJ_UM    = 40;
MIN_VOX        = 5;
MIN_ROWS       = 20;

% retardance floor (percentile of pooled control WM retardance, per region)
RET_FLOOR_PCTL = 10;
RET_SAMPLE     = 2e6;

% nulls and bootstrap
BLOCK_UM       = 400;
N_NULL_SHUFFLE = 10;          % shuffles for point estimates
N_NULL_BOOT    = 3;           % shuffles inside each bootstrap replicate

% spatial (shifted-field) null
FIELD_SIGMA_UM = 60;          % Gaussian smoothing of the fiber orientation field
FIELD_MIN_W    = 0.2;         % min smoothed fraction of eligible fiber voxels
SHIFT_UM       = [80 160 240 400 800 1600];
N_SHIFT        = 20;          % random in-plane directions per shift distance
SHIFT_REF_UM   = 400;         % shift used in the region contrast

% depth from cortex
DEPTH_EDGES_UM     = [0 500 1000 2000 4000 Inf];
NORMAL_SIGMA_UM    = 100;     % smoothing of the GM distance map before taking its gradient
RADIAL_MAX_DEG     = 30;      % vessel axis within this angle of the cortical normal = radial
NORMAL_MIN_INPLANE = 0.7;     % |in-plane component| of the normal for fiber-normal angles

% strata
RADIUS_EDGES_UM = [0 30 50 80 Inf];
INCL_EDGES_DEG  = [0 15 30 45];

% variogram
VARIO_R_UM = [20 40 60 80 120 160 240 320 480 640 960 1280 1920];
VARIO_N    = 2e5;

% sensitivity (base settings must be included)
SENS_COH      = [0.3 0.5 0.7];
SENS_RET_PCTL = [0 10 25];
SENS_BLOCK_UM = [200 400 800];

if EXPLORATORY
    N_BOOT = 100; MAX_ANCHORS = 50000; TAG = 'EXPLORATORY';
else
    N_BOOT = 500; MAX_ANCHORS = Inf;   TAG = 'MANUSCRIPT';
end

assert(isequal(REGIONS, {'front','occip'}), 'REGIONS must be {''front'',''occip''} (contrasts are front - occip).');
assert(numel(SUBJECTS) == 2, 'This script is written for two subjects.');
SENS_RET_PCTL = sort(SENS_RET_PCTL);
assert(ismember(MIN_TUBE_COH, SENS_COH) && ismember(RET_FLOOR_PCTL, SENS_RET_PCTL) && ...
       ismember(BLOCK_UM, SENS_BLOCK_UM) && ismember(SHIFT_REF_UM, SHIFT_UM), ...
       'Base settings must be included in SENS_* / SHIFT_UM.');

ROOT = '/projectnb/npbssmic/ns/CAA/';
OUT_DIR = fullfile(ROOT, 'figures/orientation/orient_vessel_controls', ...
    sprintf('%s_%s_%s', TAG, RUN_MODE, datestr(now, 'yyyymmdd_HHMMSS')));
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
rng(1);

prm = struct('vox', VOX_UM, 'units', ORIENT_UNITS, 'zeroAxis', ORIENT_ZERO_AXIS, 'ccw', ORIENT_CCW, ...
    'sigmaD', SIGMA_D_UM, 'sigmaI', SIGMA_I_UM, 'minCoh', MIN_TUBE_COH, 'minCohSens', min(SENS_COH), ...
    'radSearch', RAD_SEARCH_UM, 'distRef', DIST_REF, 'shell', SHELL_UM, 'maxDist', MAX_DIST_UM, ...
    'maxWall', MAX_WALL_UM, 'primaryMax', PRIMARY_MAX_UM, 'maxIncl', MAX_INCL_DEG, ...
    'epvAdj', EPVS_ADJ_UM, 'minVox', MIN_VOX, 'block', BLOCK_UM, 'maxPts', MAX_ANCHORS, ...
    'fieldSigma', FIELD_SIGMA_UM, 'fieldMinW', FIELD_MIN_W, 'shift', SHIFT_UM, 'nShift', N_SHIFT, ...
    'normalSigma', NORMAL_SIGMA_UM, 'varioR', VARIO_R_UM, 'varioN', VARIO_N);

nS   = numel(SUBJECTS);
nR   = numel(REGIONS);
nSh  = ceil(MAX_DIST_UM / SHELL_UM);
xShell = ((1:nSh)' - 0.5) * SHELL_UM;
kP   = round(PRIMARY_MAX_UM / SHELL_UM);
kFar = floor(FAR_MIN_UM / SHELL_UM) + 1;
nD   = numel(SHIFT_UM);
iRef = find(SHIFT_UM == SHIFT_REF_UM);
nDep = numel(DEPTH_EDGES_UM) - 1;
nRad = numel(RADIUS_EDGES_UM) - 1;
nInc = numel(INCL_EDGES_DEG) - 1;
if strcmp(DIST_REF, 'nearest'), XL = 'Distance from nearest vessel/EPVS (\mum)';
else, XL = 'Distance from vessel wall (\mum)'; end

VOL = struct('sub', {}, 'reg', {}, 'name', {}, 'col', {}, 'mk', {}, 'ls', {});
for s = 1:nS
    for r = 1:nR
        VOL(end+1) = struct('sub', SUBJECTS{s}, 'reg', REGIONS{r}, ...
            'name', sprintf('%s %s', SUBJECTS{s}, REG_LABEL.(REGIONS{r})), ...
            'col', REG_COLOR.(REGIONS{r}), 'mk', SUBJ_MARKER{s}, 'ls', SUBJ_LS{s}); %#ok<SAGROW>
    end
end
nV = numel(VOL);
vi = @(s, r) (s - 1) * nR + r;

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
for s = 1:nS
    key = regexprep(lower(SUBJECTS{s}), '[^a-z0-9]', '');
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
    for v = 1:nV
        run_orient_check(get_vol(caa, VOL(v).sub, VOL(v).reg), prm, VOL(v).sub, VOL(v).reg, OUT_DIR);
    end
    fprintf(['\nCheck figures in %s\nCyan (optic axis) should follow visible tract texture; ' ...
             'if rotated/mirrored, change ORIENT_ZERO_AXIS / ORIENT_CCW.\n'], OUT_DIR);
    return
end

%% 3. Retardance floors from pooled control WM (per region) ----------------------
retCfg = struct();
for r = 1:nR
    reg = REGIONS{r}; vals = [];
    for s = 1:nS
        Vt = get_vol(caa, SUBJECTS{s}, reg);
        idx = find(logical(Vt.mask_wm) & ~logical(Vt.seg) & ~logical(Vt.epvs));
        if numel(idx) > RET_SAMPLE, idx = idx(randperm(numel(idx), RET_SAMPLE)); end
        rv_ = double(Vt.ret_full(idx)); vals = [vals; rv_(isfinite(rv_))]; %#ok<AGROW>
        clear Vt idx rv_
    end
    fl = -inf(size(SENS_RET_PCTL)); pos = SENS_RET_PCTL > 0;
    fl(pos) = prctile(vals, SENS_RET_PCTL(pos));
    retCfg.(reg) = struct('pctl', SENS_RET_PCTL, 'floors', fl, 'floor', fl(SENS_RET_PCTL == RET_FLOOR_PCTL));
    fprintf('%s retardance floors (p%s): %s; base p%g = %.4g\n', reg, mat2str(SENS_RET_PCTL), ...
        mat2str(fl, 4), RET_FLOOR_PCTL, retCfg.(reg).floor);
end

%% 4. Per-volume extraction ------------------------------------------------------
X = cell(1, nV);
for v = 1:nV
    fprintf('%s ...\n', VOL(v).name);
    X{v} = analyze_volume(get_vol(caa, VOL(v).sub, VOL(v).reg), prm, retCfg.(VOL(v).reg));
    I = X{v}.info;
    fprintf(['  anchors %d | valid axis %d | primary %d | GM present %d | ' ...
             'orientation corr. length %.0f um\n  shift retention (by distance) %s\n'], ...
        I.n_anchor_total, I.n_anchor_valid, nnz(X{v}.A.valid & ~X{v}.A.inclined), I.has_gm, ...
        X{v}.vario.corrlen_um, mat2str(round(I.shift_retained', 2)));
end

%% 5. Per-volume statistics (block bootstrap) ------------------------------------
fPrim  = @(T, K) stat_prim(T, K, MIN_ROWS);
fSpat  = @(T, K) stat_spatial(T, nD, N_SHIFT, MIN_ROWS);
fShell = @(T, K) stat_shell(T, nSh, xShell, K, MIN_ROWS);
fFree  = @(T, K) stat_free(T, nSh, kP, kFar, MIN_ROWS);
fDepth = @(T, K) stat_depth(T, DEPTH_EDGES_UM, RADIAL_MAX_DEG, NORMAL_MIN_INPLANE, K, MIN_ROWS);
fStrat = @(T, K) [strat_excess(T, discretize(T.radius_um, RADIUS_EDGES_UM), nRad, K, MIN_ROWS); ...
                  strat_excess(T, discretize(T.incl_deg,  INCL_EDGES_DEG),  nInc, K, MIN_ROWS)];
fComp  = @(T1, T2, K) stat_comp(T1, T2, K, RADIUS_EDGES_UM, INCL_EDGES_DEG, MIN_ROWS);

PCOLS = {'block','bx','by','bz','x1','x2','x3','dtheta_deg','thetaF_deg','thetaV_deg','R', ...
         'radius_um','incl_deg','dField','dShift'};
DCOLS = {'block','depth_um','radial_deg','inclined','nInplane','fibN_deg','dtheta_deg', ...
         'thetaF_deg','thetaV_deg'};
ST = cell(1, nV);
for v = 1:nV
    fprintf('Bootstrapping %s ...\n', VOL(v).name);
    A = X{v}.A;
    st = struct();
    st.P   = A(A.valid & ~A.inclined, PCOLS);
    st.SHp = X{v}.SH(X{v}.SH.valid & ~X{v}.SH.inclined, :);
    Av     = A(A.valid, DCOLS);
    st.prim  = boot_stat(st.P,     fPrim,  N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    st.spat  = boot_stat(st.P,     fSpat,  N_BOOT, 0, 0);
    st.shell = boot_stat(st.SHp,   fShell, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    st.free  = boot_stat(X{v}.SHa, fFree,  N_BOOT, 0, 0);
    st.depth = boot_stat(Av,       fDepth, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    st.strat = boot_stat(st.P,     fStrat, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    ST{v} = st;
end

% composition-adjusted region contrast (per subject)
CMP = cell(1, nS);
for s = 1:nS
    CMP{s} = boot_pair(ST{vi(s,1)}.P, ST{vi(s,2)}.P, fComp, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
end

%% 6. Tables: volume metrics, replication, region contrasts ----------------------
MET = struct('name', {}, 'src', {}, 'row', {}, 'unit', {}, 'null0', {}, 'contrast', {});
addm = @(M, n, sr, rw, u, z, c) [M, struct('name', n, 'src', sr, 'row', rw, 'unit', u, 'null0', z, 'contrast', c)];
MET = addm(MET, 'median_dtheta',  'prim', 1, 'deg', false, true);
MET = addm(MET, 'excess_global',  'prim', 2, 'deg', true,  true);
MET = addm(MET, 'excess_block',   'prim', 3, 'deg', true,  true);
MET = addm(MET, 'median_R',       'prim', 4, 'R',   false, true);
for j = 1:nD
    MET = addm(MET, sprintf('excess_spatial_%gum', SHIFT_UM(j)), 'spat', j, 'deg', true, j == iRef);
end
MET = addm(MET, 'decay_lambda_um',    'shell', nSh + 1,   'um',  false, true);
MET = addm(MET, 'decay_half_um',      'shell', nSh + 2,   'um',  false, true);
MET = addm(MET, 'R_near',             'free',  2*nSh + 3, 'R',   false, true);
MET = addm(MET, 'R_near_minus_far',   'free',  2*nSh + 1, 'R',   true,  true);
MET = addm(MET, 'ret_near_minus_far', 'free',  2*nSh + 2, 'ret', true,  true);

VMT = table();
for v = 1:nV
    for m = 1:numel(MET)
        B = ST{v}.(MET(m).src); k = MET(m).row;
        VMT = [VMT; table(string(VOL(v).sub), string(VOL(v).reg), string(MET(m).name), string(MET(m).unit), ...
            B.pt(k), B.lo(k), B.hi(k), MET(m).null0 && (B.lo(k) > 0 || B.hi(k) < 0), ...
            'VariableNames', {'subject','region','metric','unit','pt','lo','hi','ci_excludes_0'})]; %#ok<AGROW>
    end
end

REP = table();      % replication of within-volume effects (metrics whose null value is 0)
for r = 1:nR
    for m = find([MET.null0])
        rows = VMT(VMT.region == REGIONS{r} & VMT.metric == MET(m).name, :);
        REP = [REP; table(string(REGIONS{r}), string(MET(m).name), rows.pt(1), rows.pt(2), ...
            all(rows.lo > 0) || all(rows.hi < 0), ...
            'VariableNames', {'region','metric', ['pt_' SUBJECTS{1}], ['pt_' SUBJECTS{2}], 'replicated'})]; %#ok<AGROW>
    end
end

CT = table();       % frontal - occipital within subject
for s = 1:nS
    f = vi(s,1); o = vi(s,2);
    for m = find([MET.contrast])
        Bf = ST{f}.(MET(m).src); Bo = ST{o}.(MET(m).src); k = MET(m).row;
        d = Bf.reps(k,:) - Bo.reps(k,:);
        CT = [CT; table(string(SUBJECTS{s}), string(MET(m).name), string(MET(m).unit), Bf.pt(k), Bo.pt(k), ...
            Bf.pt(k) - Bo.pt(k), prctile(d, 2.5), prctile(d, 97.5), ...
            'VariableNames', {'subject','metric','unit','front_pt','occip_pt','diff_pt','diff_lo','diff_hi'})]; %#ok<AGROW>
    end
    Bc = CMP{s};
    CT = [CT; table(string(SUBJECTS{s}), "excess_global_comp_adjusted", "deg", Bc.pt(1), Bc.pt(3), Bc.pt(5), ...
        prctile(Bc.reps(5,:), 2.5), prctile(Bc.reps(5,:), 97.5), 'VariableNames', CT.Properties.VariableNames)]; %#ok<AGROW>
end
CT.replicated = false(height(CT), 1);
for mname = unique(CT.metric)'
    i = CT.metric == mname;
    CT.replicated(i) = all(CT.diff_lo(i) > 0) || all(CT.diff_hi(i) < 0);
end

COMP = table();
for s = 1:nS
    Bc = CMP{s};
    COMP = [COMP; table(string(SUBJECTS{s}), ["excess_front"; "excess_occip"; "excess_occip_adjusted"; ...
        "diff_raw"; "diff_adjusted"], Bc.pt, Bc.lo, Bc.hi, 'VariableNames', {'subject','quantity','pt','lo','hi'})]; %#ok<AGROW>
end

%% 7. Long-format profiles -------------------------------------------------------
PRO = table();
depX = DEPTH_EDGES_UM(1:end-1)'; radX = RADIUS_EDGES_UM(1:end-1)'; incX = INCL_EDGES_DEG(1:end-1)';
DEP_NAMES = {'excess','frac_radial','dtheta_radial_deg','fiber_normal_angle_deg'};
for v = 1:nV
    PRO = [PRO; long_rows(VOL(v), 'excess_vs_shell_um',  xShell,      ST{v}.shell, 1:nSh)]; %#ok<AGROW>
    PRO = [PRO; long_rows(VOL(v), 'R_vs_shell_um',       xShell,      ST{v}.free,  1:nSh)]; %#ok<AGROW>
    PRO = [PRO; long_rows(VOL(v), 'ret_vs_shell_um',     xShell,      ST{v}.free,  nSh + (1:nSh))]; %#ok<AGROW>
    PRO = [PRO; long_rows(VOL(v), 'excess_vs_shift_um',  SHIFT_UM(:), ST{v}.spat,  1:nD)]; %#ok<AGROW>
    for c = 1:4
        PRO = [PRO; long_rows(VOL(v), [DEP_NAMES{c} '_vs_depth_um'], depX, ST{v}.depth, (c-1)*nDep + (1:nDep))]; %#ok<AGROW>
    end
    PRO = [PRO; long_rows(VOL(v), 'excess_vs_radius_um',  radX, ST{v}.strat, 1:nRad)]; %#ok<AGROW>
    PRO = [PRO; long_rows(VOL(v), 'excess_vs_incl_deg',   incX, ST{v}.strat, nRad + (1:nInc))]; %#ok<AGROW>
end

VAR = table();
for v = 1:nV
    S = X{v}.vario; n = numel(S.r_um);
    VAR = [VAR; table(repmat(string(VOL(v).sub), n, 1), repmat(string(VOL(v).reg), n, 1), S.r_um, S.C, S.rho, ...
        S.n_pairs, repmat(S.C_far, n, 1), repmat(S.corrlen_um, n, 1), ...
        'VariableNames', {'subject','region','r_um','C','rho','n_pairs','C_far','corrlen_um'})]; %#ok<AGROW>
end

%% 8. Robustness: split-half and sensitivity -------------------------------------
SPL = table();
SPL_NAMES = ["median_dtheta"; "excess_global"; "excess_block"; "median_R"; "decay_half_um"];
for v = 1:nV
    P = ST{v}.P;     par  = mod(P.bx + P.by + P.bz, 2) == 0;
    SHp = ST{v}.SHp; parS = mod(SHp.bx + SHp.by + SHp.bz, 2) == 0;
    h0 = fPrim(P(par,:), N_NULL_SHUFFLE);      h1 = fPrim(P(~par,:), N_NULL_SHUFFLE);
    e0 = fShell(SHp(parS,:), N_NULL_SHUFFLE);  e1 = fShell(SHp(~parS,:), N_NULL_SHUFFLE);
    ok = isfinite(e0(1:nSh)) & isfinite(e1(1:nSh));
    rS = NaN; if nnz(ok) >= 3, rS = corr(e0(ok), e1(ok)); end
    n = numel(SPL_NAMES);
    SPL = [SPL; table(repmat(string(VOL(v).sub), n, 1), repmat(string(VOL(v).reg), n, 1), SPL_NAMES, ...
        [h0; e0(nSh+2)], [h1; e1(nSh+2)], repmat(rS, n, 1), ...
        'VariableNames', {'subject','region','metric','half_A','half_B','shell_profile_r'})]; %#ok<AGROW>
end

SEN = table();
SCOLS = {'block','dtheta_deg','thetaF_deg','thetaV_deg','R'};
for v = 1:nV
    fprintf('Sensitivity %s ...\n', VOL(v).name);
    A = X{v}.A; Tl = {}; lab = {}; grp = {};
    for c = SENS_COH
        Tl{end+1} = A(A.tube_coh >= c & ~A.inclined, SCOLS); %#ok<SAGROW>
        lab{end+1} = sprintf('tube coherence >= %.2g', c); grp{end+1} = 'coherence'; %#ok<SAGROW>
    end
    for k = 1:numel(SIGMA_I_UM)
        ok = A.cohS(:,k) >= MIN_TUBE_COH & A.inclS(:,k) < MAX_INCL_DEG;
        T = A(ok, SCOLS); T.thetaV_deg = A.thVs(ok,k); T.dtheta_deg = axial_diff_deg(T.thetaF_deg, T.thetaV_deg);
        Tl{end+1} = T; lab{end+1} = sprintf('sigma_i = %g um only', SIGMA_I_UM(k)); grp{end+1} = 'scale'; %#ok<SAGROW>
    end
    for f = 1:numel(SENS_RET_PCTL)
        ok = A.valid & ~A.inclined & A.nFf(:,f) >= MIN_VOX;
        T = A(ok, SCOLS); T.thetaF_deg = A.thFf(ok,f); T.dtheta_deg = axial_diff_deg(T.thetaF_deg, T.thetaV_deg);
        Tl{end+1} = T; lab{end+1} = sprintf('retardance floor p%g', SENS_RET_PCTL(f)); grp{end+1} = 'retardance'; %#ok<SAGROW>
    end
    for b = SENS_BLOCK_UM
        T = ST{v}.P; Bv = max(1, round(b / VOX_UM));
        [~, ~, g] = unique([ceil(T.x1 / Bv), ceil(T.x2 / Bv), ceil(T.x3 / Bv)], 'rows'); T.block = g;
        Tl{end+1} = T(:, SCOLS); lab{end+1} = sprintf('block %g um', b); grp{end+1} = 'block'; %#ok<SAGROW>
    end
    for q = 1:numel(Tl)
        B = boot_stat(Tl{q}, fPrim, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
        SEN = [SEN; table(repmat(string(VOL(v).sub), 2, 1), repmat(string(VOL(v).reg), 2, 1), ...
            repmat(string(grp{q}), 2, 1), repmat(string(lab{q}), 2, 1), ["excess_global"; "excess_block"], ...
            repmat(height(Tl{q}), 2, 1), B.pt(2:3), B.lo(2:3), B.hi(2:3), ...
            'VariableNames', {'subject','region','group','variant','metric','n_anchors','pt','lo','hi'})]; %#ok<AGROW>
    end
end

%% 9. Write outputs --------------------------------------------------------------
disp(REP); disp(CT); disp(COMP);
writetable(VMT,  fullfile(OUT_DIR, 'volume_metrics.csv'));
writetable(REP,  fullfile(OUT_DIR, 'replication.csv'));
writetable(CT,   fullfile(OUT_DIR, 'region_contrast.csv'));
writetable(COMP, fullfile(OUT_DIR, 'composition_adjusted.csv'));
writetable(PRO,  fullfile(OUT_DIR, 'profiles_long.csv'));
writetable(VAR,  fullfile(OUT_DIR, 'variogram.csv'));
writetable(SPL,  fullfile(OUT_DIR, 'split_half.csv'));
writetable(SEN,  fullfile(OUT_DIR, 'sensitivity.csv'));
AA = table();
for v = 1:nV
    T = removevars(X{v}.A, {'thVs','inclS','cohS','thFf','nFf','dShift'});
    T.subject = repmat(string(VOL(v).sub), height(T), 1); T.region = repmat(string(VOL(v).reg), height(T), 1);
    AA = [AA; T]; %#ok<AGROW>
end
writetable(AA, fullfile(OUT_DIR, 'anchors.csv')); clear AA

fid = fopen(fullfile(OUT_DIR, 'README.txt'), 'w');
fprintf(fid, ['%s run %s\nSubjects: %s; regions: %s. Contrasts are frontal - occipital within subject.\n' ...
    'Vessel axis: structure tensor (sigma_d %g um, sigma_i %s um, tube coherence >= %g) at vessel-wall voxels; ' ...
    'vessels inclined >= %g deg from the en-face plane excluded from angle metrics.\n' ...
    'Distance reference: %s; shells %g um to %g um; primary window <= %g um; far zone >= %g um.\n' ...
    'Retardance floor: p%g of pooled control WM per region.\n' ...
    'dtheta in [0,90] (0 parallel, 90 perpendicular). excess = median(null) - median(observed); > 0 = aligned beyond the null.\n' ...
    'Nulls: global shuffle; shuffle within %g um blocks; rigid in-plane shift of the fiber field (Gaussian sigma %g um) by %s um, ' ...
    '%d directions each (observed and null both use the smoothed field at the vessel-wall anchor).\n' ...
    'Composition adjustment: occipital anchors weighted to the frontal radius x inclination distribution (bins %s um x %s deg).\n' ...
    'Depth: distance to GM (mask & ~mask_wm); radial vessel = axis within %g deg of the cortical normal (gradient of the GM ' ...
    'distance map smoothed at %g um).\n' ...
    'CIs: %g um blocks resampled within volume (N_BOOT %d); contrast CIs from the difference of independent bootstraps; ' ...
    'composition-adjusted CI from a paired bootstrap. Replicated = both subjects same sign with CIs excluding 0. ' ...
    'Anchors are pseudoreplicated within n = 2 subjects; CIs describe these specimens only.\n' ...
    'Orientation convention: zero axis %s, ccw %d, units %s. MAX_ANCHORS %g.\n'], ...
    TAG, datestr(now), strjoin(SUBJECTS, ','), strjoin(REGIONS, ','), SIGMA_D_UM, mat2str(SIGMA_I_UM), ...
    MIN_TUBE_COH, MAX_INCL_DEG, DIST_REF, SHELL_UM, MAX_DIST_UM, PRIMARY_MAX_UM, FAR_MIN_UM, RET_FLOOR_PCTL, ...
    BLOCK_UM, FIELD_SIGMA_UM, mat2str(SHIFT_UM), N_SHIFT, mat2str(RADIUS_EDGES_UM), mat2str(INCL_EDGES_DEG), ...
    RADIAL_MAX_DEG, NORMAL_SIGMA_UM, BLOCK_UM, N_BOOT, ORIENT_ZERO_AXIS, ORIENT_CCW, ORIENT_UNITS, MAX_ANCHORS);
fclose(fid);

%% 10. Figures -------------------------------------------------------------------
names = {VOL.name};

% Fig 1: alignment vs nulls, spatial-shift curve, variogram
[f, tl] = new_fig(1, 3);
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
off = linspace(-0.27, 0.27, nV); h = gobjects(nV, 1);
for v = 1:nV
    B = ST{v}.prim;
    h(v) = plot_pt(ax, [1 2] + off(v), B.pt(2:3), B.lo(2:3), B.hi(2:3), VOL(v).col, VOL(v).mk);
end
yline(ax, 0, 'k:'); set(ax, 'XTick', [1 2], 'XTickLabel', {'Global shuffle','Block shuffle'}); xlim(ax, [0.5 2.5]);
ylabel(ax, 'Alignment excess (deg)'); title(ax, 'Alignment vs non-spatial nulls');
legend(ax, h, names, 'Location', 'best', 'Box', 'off');
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
for v = 1:nV
    B = ST{v}.spat;
    h(v) = plot_ci(ax, SHIFT_UM, B.pt, B.lo, B.hi, VOL(v).col, VOL(v).mk, VOL(v).ls);
    cl = X{v}.vario.corrlen_um;
    if isfinite(cl), xline(ax, cl, ':', 'Color', VOL(v).col, 'LineWidth', 1.2); end
end
set(ax, 'XScale', 'log'); yline(ax, 0, 'k:');
xlabel(ax, 'In-plane shift of fiber field (\mum)'); ylabel(ax, 'Alignment excess (deg)');
title(ax, 'Spatial null (dotted: fiber correlation length)');
legend(ax, h, names, 'Location', 'best', 'Box', 'off');
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
for v = 1:nV
    S = X{v}.vario;
    h(v) = plot(ax, S.r_um, S.rho, 'Color', VOL(v).col, 'LineStyle', VOL(v).ls, 'Marker', VOL(v).mk, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.5);
end
yline(ax, exp(-1), 'k:'); xline(ax, BLOCK_UM, '--', 'Color', [0.5 0.5 0.5]);
set(ax, 'XScale', 'log'); xlabel(ax, 'In-plane separation (\mum)'); ylabel(ax, 'Orientation correlation \rho');
title(ax, 'Fiber orientation variogram (dashed: block size)');
legend(ax, h, names, 'Location', 'best', 'Box', 'off');
save_fig(f, OUT_DIR, 'fig1_alignment_nulls');

% Fig 2: frontal - occipital forest plots
UNITS = {'deg','um','R','ret'};
ULAB  = {'\Delta (deg)', '\Delta (\mum)', '\Delta R', '\Delta retardance'};
[f, tl] = new_fig(1, numel(UNITS));
for u = 1:numel(UNITS)
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
    rows = CT(CT.unit == UNITS{u}, :);
    if isempty(rows), continue; end
    mets = unique(rows.metric, 'stable'); nm = numel(mets); hs = gobjects(nS, 1);
    for s = 1:nS
        Rs = rows(rows.subject == SUBJECTS{s}, :);
        [~, y] = ismember(Rs.metric, mets);
        hs(s) = errorbar(ax, Rs.diff_pt, y + (s - 1.5) * 0.25, Rs.diff_pt - Rs.diff_lo, Rs.diff_hi - Rs.diff_pt, ...
            'horizontal', 'Marker', SUBJ_MARKER{s}, 'Color', 'k', 'MarkerFaceColor', SUBJ_FACE{s}, 'MarkerSize', 7, ...
            'LineStyle', 'none', 'LineWidth', 1.3, 'CapSize', 5);
    end
    xline(ax, 0, 'k:');
    set(ax, 'YTick', 1:nm, 'YTickLabel', strrep(cellstr(mets), '_', ' '), 'YDir', 'reverse'); ylim(ax, [0.4 nm + 0.6]);
    xlabel(ax, [ULAB{u} ', frontal − occipital']);
    legend(ax, hs, SUBJECTS, 'Location', 'best', 'Box', 'off');
end
title(tl, 'Frontal − occipital within subject (95% block-bootstrap CI)');
save_fig(f, OUT_DIR, 'fig2_region_contrast');

% Fig 3: distance decay
[f, tl] = new_fig(1, 3);
for r = 1:nR
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); hs = gobjects(nS, 1);
    for s = 1:nS
        v = vi(s, r); B = ST{v}.shell;
        hs(s) = plot_ci(ax, xShell, B.pt(1:nSh), B.lo(1:nSh), B.hi(1:nSh), VOL(v).col, VOL(v).mk, VOL(v).ls);
        p = decay_fit(xShell, B.pt(1:nSh));
        if all(isfinite(p))
            xx = linspace(xShell(1), xShell(end), 100)';
            plot(ax, xx, p(2) * exp(-xx / p(1)) + p(3), ':', 'Color', VOL(v).col * 0.7, 'LineWidth', 1.3);
        end
    end
    xline(ax, PRIMARY_MAX_UM, ':', 'Color', [0.4 0.4 0.4]); yline(ax, 0, 'k:');
    xlabel(ax, XL); ylabel(ax, 'Alignment excess (deg)');
    title(ax, sprintf('%s: excess vs distance (dotted: exp. fit)', REG_LABEL.(REGIONS{r})));
    legend(ax, hs, SUBJECTS, 'Location', 'best', 'Box', 'off');
end
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
for v = 1:nV
    B = ST{v}.shell;
    plot_pt(ax, v - 0.1, B.pt(nSh+2), B.lo(nSh+2), B.hi(nSh+2), VOL(v).col, VOL(v).mk);
    plot_pt(ax, v + 0.1, B.pt(nSh+1), B.lo(nSh+1), B.hi(nSh+1), VOL(v).col * 0.6, VOL(v).mk);
end
set(ax, 'XTick', 1:nV, 'XTickLabel', names); xtickangle(ax, 30); xlim(ax, [0.5 nV + 0.5]);
ylabel(ax, 'Distance (\mum)'); title(ax, 'Half-decay distance (left) and \lambda (right)');
save_fig(f, OUT_DIR, 'fig3_distance_decay');

% Fig 4: depth from the GM/WM boundary
DEP_Y = {'Alignment excess (deg)', 'Fraction radial vessels', '\Delta\theta of radial vessels (deg)', ...
         'Fiber angle to cortical normal (deg)'};
DEP_REF = [0 NaN 45 45];
[f, tl] = new_fig(1, 4);
for c = 1:4
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
    for v = 1:nV
        B = ST{v}.depth; k = (c-1)*nDep + (1:nDep);
        h(v) = plot_ci(ax, (1:nDep)', B.pt(k), B.lo(k), B.hi(k), VOL(v).col, VOL(v).mk, VOL(v).ls);
    end
    if isfinite(DEP_REF(c)), yline(ax, DEP_REF(c), 'k:'); end
    if c >= 3, ylim(ax, [0 90]); yticks(ax, 0:15:90); end
    if c == 2, ylim(ax, [0 1]); end
    set(ax, 'XTick', 1:nDep, 'XTickLabel', edge_labels(DEPTH_EDGES_UM)); xtickangle(ax, 30);
    xlabel(ax, 'Depth below GM/WM boundary (\mum)'); ylabel(ax, DEP_Y{c});
    legend(ax, h, names, 'Location', 'best', 'Box', 'off');
end
title(tl, sprintf('Depth from cortex (radial = vessel axis within %g deg of cortical normal)', RADIAL_MAX_DEG));
save_fig(f, OUT_DIR, 'fig4_depth');

% Fig 5: vessel-axis-free metrics
[f, tl] = new_fig(1, 2);
FY = {'Fiber coherence R', 'Retardance'};
for c = 1:2
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
    for v = 1:nV
        B = ST{v}.free; k = (c-1)*nSh + (1:nSh);
        h(v) = plot_ci(ax, xShell, B.pt(k), B.lo(k), B.hi(k), VOL(v).col, VOL(v).mk, VOL(v).ls);
    end
    xline(ax, PRIMARY_MAX_UM, ':', 'Color', [0.4 0.4 0.4]); xline(ax, FAR_MIN_UM, ':', 'Color', [0.4 0.4 0.4]);
    xlabel(ax, XL); ylabel(ax, FY{c}); title(ax, [FY{c} ' vs distance (dotted: near / far limits)']);
    legend(ax, h, names, 'Location', 'best', 'Box', 'off');
end
save_fig(f, OUT_DIR, 'fig5_axis_free');

% Fig 6: strata and composition adjustment
[f, tl] = new_fig(1, 3);
SX = {1:nRad, nRad + (1:nInc)}; SL = {edge_labels(RADIUS_EDGES_UM), edge_labels(INCL_EDGES_DEG)};
SXL = {'Vessel radius (\mum)', 'Inclination from en-face plane (deg)'};
for c = 1:2
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
    for v = 1:nV
        B = ST{v}.strat; k = SX{c};
        h(v) = plot_ci(ax, (1:numel(k))', B.pt(k), B.lo(k), B.hi(k), VOL(v).col, VOL(v).mk, VOL(v).ls);
    end
    yline(ax, 0, 'k:'); set(ax, 'XTick', 1:numel(SX{c}), 'XTickLabel', SL{c});
    xlabel(ax, SXL{c}); ylabel(ax, 'Alignment excess (deg)'); title(ax, ['Excess by ' lower(strtok(SXL{c}, '('))]);
    legend(ax, h, names, 'Location', 'best', 'Box', 'off');
end
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
CC = {REG_COLOR.front, REG_COLOR.occip, lighten(REG_COLOR.occip, 0.5)}; hc = gobjects(3, 1);
for s = 1:nS
    Bc = CMP{s};
    for q = 1:3
        hc(q) = plot_pt(ax, s + (q - 2) * 0.2, Bc.pt(q), Bc.lo(q), Bc.hi(q), CC{q}, SUBJ_MARKER{s});
    end
end
yline(ax, 0, 'k:'); set(ax, 'XTick', 1:nS, 'XTickLabel', SUBJECTS); xlim(ax, [0.5 nS + 0.5]);
ylabel(ax, 'Alignment excess (deg)'); title(ax, 'Occipital reweighted to frontal radius x inclination');
legend(ax, hc, {'Frontal', 'Occipital', 'Occipital (adjusted)'}, 'Location', 'best', 'Box', 'off');
save_fig(f, OUT_DIR, 'fig6_strata_composition');

% Fig 7: robustness
[f, tl] = new_fig(1, 2);
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
for v = 1:nV
    Sv = SPL(SPL.subject == VOL(v).sub & SPL.region == VOL(v).reg, :);
    g1 = Sv.metric == "excess_global"; g2 = Sv.metric == "excess_block";
    h(v) = plot(ax, Sv.half_A(g1), Sv.half_B(g1), VOL(v).mk, 'MarkerFaceColor', VOL(v).col, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 9);
    plot(ax, Sv.half_A(g2), Sv.half_B(g2), VOL(v).mk, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', VOL(v).col, ...
        'MarkerSize', 9, 'LineWidth', 1.5);
end
lim = [min([ax.XLim ax.YLim]) max([ax.XLim ax.YLim])]; plot(ax, lim, lim, 'k:'); axis(ax, 'square');
xlabel(ax, 'Excess, half A (deg)'); ylabel(ax, 'Excess, half B (deg)');
title(ax, 'Split-half (filled: global null, open: block null)');
legend(ax, h, names, 'Location', 'best', 'Box', 'off');
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
Sg = SEN(SEN.metric == "excess_global", :);
vars = unique(Sg.variant, 'stable'); nvar = numel(vars); off = linspace(-0.27, 0.27, nV);
for v = 1:nV
    Sv = Sg(Sg.subject == VOL(v).sub & Sg.region == VOL(v).reg, :);
    [~, y] = ismember(Sv.variant, vars);
    h(v) = errorbar(ax, Sv.pt, y + off(v), Sv.pt - Sv.lo, Sv.hi - Sv.pt, 'horizontal', 'Marker', VOL(v).mk, ...
        'Color', VOL(v).col, 'MarkerFaceColor', VOL(v).col, 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 4);
end
xline(ax, 0, 'k:');
set(ax, 'YTick', 1:nvar, 'YTickLabel', cellstr(vars), 'YDir', 'reverse'); ylim(ax, [0.4 nvar + 0.6]);
xlabel(ax, 'Alignment excess, global null (deg)'); title(ax, 'Sensitivity to analysis settings');
legend(ax, h, names, 'Location', 'best', 'Box', 'off');
save_fig(f, OUT_DIR, 'fig7_robustness');

fprintf('Done. Output: %s\n', OUT_DIR);

%% ========================= Local functions ==================================

function V = get_vol(caa, sub, reg)
    sf = resolve_field(caa, sub, true);
    rf = resolve_field(caa.(sf), reg, false);
    V = caa.(sf).(rf);
end

function X = analyze_volume(V, prm, rc)
    wm  = logical(V.mask_wm);
    pad = ceil(max(prm.maxWall, 3 * max(prm.sigmaI)) / prm.vox) + 2;
    [r1, r2, r3] = bbox_ranges(wm, pad);
    wm  = wm(r1,r2,r3);
    seg = logical(V.seg(r1,r2,r3));
    epv = logical(V.epvs(r1,r2,r3));
    if isfield(V, 'mask'), gm = logical(V.mask(r1,r2,r3)) & ~wm & ~seg & ~epv; else, gm = false(size(wm)); end
    [th, info] = decode_orient(V.orient(r1,r2,r3), prm.units, wm);
    th  = single(th);
    ret = single(V.ret_full(r1,r2,r3));
    sz  = size(wm);
    elig0 = wm & ~seg & ~epv & isfinite(th) & isfinite(ret) & ret >= rc.floors(1);   % lowest sensitivity floor
    eligB = elig0 & ret >= rc.floor;                                                % base floor

    % --- anchors: nearest vessel-wall voxel of each perivascular fiber voxel ---
    [dWall, idxW] = bwdist(seg);
    dWall = dWall * prm.vox;
    if strcmp(prm.distRef, 'nearest'), dRef = bwdist(seg | epv) * prm.vox; else, dRef = dWall; end
    fib = elig0 & dRef > 0 & dRef <= prm.maxDist & dWall <= prm.maxWall;
    clear dWall
    vox = find(fib); clear fib
    nSh = ceil(prm.maxDist / prm.shell);
    shV = min(nSh, ceil(double(dRef(vox)) / prm.shell)); clear dRef
    [aIdx, ~, pid] = unique(double(idxW(vox))); clear idxW
    nA = numel(aIdx);
    info.n_anchor_total = nA;

    % --- vessel axis: adaptive scale and each scale separately ---
    vt = vessel_tensor(seg, aIdx, prm);
    [thV, phi] = axis_angles(vt.t, prm);
    nS = numel(prm.sigmaI); thVs = nan(nA, nS); phiS = nan(nA, nS);
    for k = 1:nS, [thVs(:,k), phiS(:,k)] = axis_angles(vt.tS(:,:,k), prm); end
    valid = isfinite(vt.coh) & vt.coh >= prm.minCoh;
    okAx  = isfinite(vt.coh) & vt.coh >= prm.minCohSens;
    incl  = phi >= prm.maxIncl;
    info.n_anchor_valid = nnz(valid);

    % --- radius, EPVS adjacency ---
    dIn  = bwdist(~seg);
    kR   = max(1, round(prm.radSearch / prm.vox));
    rmap = imdilate(dIn, strel('cube', 2*kR + 1)); clear dIn
    rad  = double(rmap(aIdx)) * prm.vox; clear rmap
    if any(epv(:)), dE = bwdist(epv); dEpv = double(dE(aIdx)) * prm.vox; clear dE; else, dEpv = inf(nA, 1); end
    epvAdj = dEpv <= prm.epvAdj;

    % --- anchor subsets: axis-based analyses / axis-free analyses ---
    keepAx  = subsample(okAx, prm.maxPts);
    keepAll = subsample(true(nA, 1), prm.maxPts);
    kv = keepAx(pid) | keepAll(pid);
    vox = vox(kv); pid = pid(kv); shV = shV(kv);
    tv = double(th(vox)); rv = double(ret(vox)); clear vox
    c2 = cos(2*tv); s2 = sin(2*tv); clear tv

    % --- per (anchor, shell) sums at the base retardance floor ---
    inB = rv >= rc.floor;
    key = (pid - 1) * nSh + shV;
    acc = @(x) reshape(accumarray(key(inB), x(inB), [nA*nSh, 1]), nSh, nA)';
    C = acc(c2); S = acc(s2); N = acc(ones(size(c2))); RT = acc(rv);

    % --- primary-window fiber angle at each sensitivity floor ---
    kP = round(prm.primaryMax / prm.shell);
    nF = numel(rc.floors); thFf = nan(nA, nF); nFf = zeros(nA, nF);
    for f = 1:nF
        m  = rv >= rc.floors(f) & shV <= kP;
        cf = accumarray(pid(m), c2(m), [nA 1]); sf = accumarray(pid(m), s2(m), [nA 1]);
        nFf(:,f)  = accumarray(pid(m), 1, [nA 1]);
        thFf(:,f) = mod(rad2deg(0.5 * atan2(sf, cf)), 180);
    end
    clear c2 s2 rv pid shV key inB

    % --- blocks and coordinates ---
    [a1, a2, a3] = ind2sub(sz, aIdx);
    Bv = max(1, round(prm.block / prm.vox));
    bx = ceil(a1 / Bv); by = ceil(a2 / Bv); bz = ceil(a3 / Bv);
    blk = sub2ind(ceil([sz(1) sz(2) size(wm,3)] / Bv), bx, by, bz);
    xyz = [a1 + r1(1) - 1, a2 + r2(1) - 1, a3 + r3(1) - 1];

    % --- primary window ---
    Cp = sum(C(:,1:kP), 2); Sp = sum(S(:,1:kP), 2); Np = sum(N(:,1:kP), 2); RTp = sum(RT(:,1:kP), 2);
    thF  = mod(rad2deg(0.5 * atan2(Sp, Cp)), 180);
    Rp   = hypot(Cp, Sp) ./ Np;
    retM = RTp ./ Np;
    dth  = axial_diff_deg(thF, thV);
    rowA = find(keepAx & Np >= prm.minVox);

    % --- spatial null: smoothed fiber field, observed at anchor vs rigidly shifted ---
    sg  = prm.fieldSigma / prm.vox;
    tmp = zeros(sz, 'single'); tmp(eligB) = cos(2 * th(eligB)); Fc = imgaussfilt3(tmp, sg);
    tmp(:) = 0;                tmp(eligB) = sin(2 * th(eligB)); Fs = imgaussfilt3(tmp, sg); clear tmp
    Fw  = imgaussfilt3(single(eligB), sg);
    ip  = rowA(valid(rowA) & ~incl(rowA));
    dField = nan(nA, 1);
    dField(ip) = axial_diff_deg(field_angle(Fc, Fs, Fw, wm, a1(ip), a2(ip), a3(ip), prm.fieldMinW, false), thV(ip));
    nD = numel(prm.shift); K = prm.nShift;
    [~, ipLoc] = ismember(ip, rowA);                     % shifted angles stored for table rows only
    dShift = nan(numel(rowA), nD * K, 'single'); kept = nan(nD, K);
    for j = 1:nD
        for k = 1:K
            ang = 2 * pi * rand;
            o1 = round(prm.shift(j) / prm.vox * cos(ang)); o2 = round(prm.shift(j) / prm.vox * sin(ang));
            tf = field_angle(Fc, Fs, Fw, wm, a1(ip) + o1, a2(ip) + o2, a3(ip), prm.fieldMinW, true);
            dShift(ipLoc, (j-1)*K + k) = axial_diff_deg(tf, thV(ip));
            kept(j, k) = mean(isfinite(tf));
        end
    end
    clear Fc Fs Fw
    info.shift_retained = mean(kept, 2);

    % --- depth from the GM/WM boundary and cortical normal ---
    depth = nan(nA, 1); radial = nan(nA, 1); thN = nan(nA, 1); nIn = nan(nA, 1);
    info.has_gm = any(gm(:));
    if info.has_gm
        Dg = single(bwdist(gm)) * prm.vox;
        depth = double(Dg(aIdx));
        Dg = imgaussfilt3(Dg, prm.normalSigma / prm.vox);
        n  = grad_at(Dg, a1, a2, a3); clear Dg
        n  = n ./ vecnorm(n, 2, 2);
        dp = abs(sum(vt.t .* n, 2)); dp(dp > 1) = 1;
        radial = acosd(dp);                                  % 0 = radial (along normal)
        nIn = hypot(n(:,1), n(:,2));
        thN = rad2deg(display_to_conv(atan2(-n(:,1), n(:,2)), prm.zeroAxis, prm.ccw));
    end
    fibN = axial_diff_deg(thF, thN);
    clear gm

    % --- fiber orientation variogram ---
    X.vario = orient_variogram(th, eligB, prm.varioR, prm.vox, prm.varioN);
    clear th ret eligB elig0

    % --- tables ---
    r = rowA;
    X.A = table(r, blk(r), bx(r), by(r), bz(r), xyz(r,1), xyz(r,2), xyz(r,3), thV(r), phi(r), incl(r), ...
        vt.coh(r), valid(r), vt.sigma(r), rad(r), dEpv(r), epvAdj(r), thF(r), Rp(r), Np(r), retM(r), dth(r), ...
        depth(r), radial(r), thN(r), nIn(r), fibN(r), dField(r), thVs(r,:), phiS(r,:), vt.cohS(r,:), ...
        thFf(r,:), nFf(r,:), dShift, ...
        'VariableNames', {'anchor','block','bx','by','bz','x1','x2','x3','thetaV_deg','incl_deg','inclined', ...
        'tube_coh','valid','sigmaI_um','radius_um','dist_epvs_um','epvAdj','thetaF_deg','R','nvox','retMean', ...
        'dtheta_deg','depth_um','radial_deg','thetaN_deg','nInplane','fibN_deg','dField','thVs','inclS', ...
        'cohS','thFf','nFf','dShift'});

    Nk = N; Nk(~keepAx, :) = 0;
    [pp, ss] = find(Nk >= prm.minVox); li = sub2ind(size(N), pp, ss);
    thFs = mod(rad2deg(0.5 * atan2(S(li), C(li))), 180);
    X.SH = table(blk(pp), bx(pp), by(pp), bz(pp), ss, (ss - 0.5) * prm.shell, thV(pp), incl(pp), valid(pp), ...
        thFs, hypot(C(li), S(li)) ./ N(li), RT(li) ./ N(li), N(li), axial_diff_deg(thFs, thV(pp)), ...
        'VariableNames', {'block','bx','by','bz','shell','dist_um','thetaV_deg','inclined','valid', ...
        'thetaF_deg','R','retMean','nvox','dtheta_deg'});

    Nk = N; Nk(~keepAll, :) = 0;
    [pp, ss] = find(Nk >= prm.minVox); li = sub2ind(size(N), pp, ss);
    X.SHa = table(blk(pp), ss, (ss - 0.5) * prm.shell, hypot(C(li), S(li)) ./ N(li), RT(li) ./ N(li), N(li), ...
        'VariableNames', {'block','shell','dist_um','R','retMean','nvox'});
    X.info = info;
end

function m = subsample(m, n)
    if isfinite(n) && nnz(m) > n
        f = find(m); m(:) = false; m(f(randperm(numel(f), n))) = true;
    end
end

function [thV, phi] = axis_angles(t, prm)
    phi = rad2deg(asin(min(1, abs(t(:,3)))));
    phi(any(~isfinite(t), 2)) = NaN;
    thV = rad2deg(display_to_conv(atan2(-t(:,1), t(:,2)), prm.zeroAxis, prm.ccw));
end

function a = field_angle(Fc, Fs, Fw, wm, s1, s2, s3, minW, needWM)
    sz = size(wm); n3 = size(wm, 3); a = nan(numel(s1), 1);
    in = s1 >= 1 & s1 <= sz(1) & s2 >= 1 & s2 <= sz(2) & s3 >= 1 & s3 <= n3;
    li = sub2ind([sz(1) sz(2) n3], s1(in), s2(in), s3(in));
    ok = Fw(li) >= minW;
    if needWM, ok = ok & wm(li); end
    v = nan(numel(li), 1);
    v(ok) = mod(rad2deg(0.5 * atan2(double(Fs(li(ok))), double(Fc(li(ok))))), 180);
    a(in) = v;
end

function g = grad_at(D, a1, a2, a3)
    n1 = size(D, 1); n2 = size(D, 2); n3 = size(D, 3);
    cl = @(x, n) min(max(x, 1), n);
    at = @(i, j, k) double(D(sub2ind([n1 n2 n3], cl(i, n1), cl(j, n2), cl(k, n3))));
    g = [(at(a1+1, a2, a3) - at(a1-1, a2, a3)) / 2, ...
         (at(a1, a2+1, a3) - at(a1, a2-1, a3)) / 2, ...
         (at(a1, a2, a3+1) - at(a1, a2, a3-1)) / 2];
end

function S = orient_variogram(th, elig, rUm, vox, nSamp)
    % rho(r) = (C(r) - C_far) / (1 - C_far), C = mean cos(2 dtheta) for in-plane pairs at r
    idx = find(elig); n1 = size(elig, 1); n2 = size(elig, 2); n3 = size(elig, 3);
    if numel(idx) > nSamp, idx = idx(randperm(numel(idx), nSamp)); end
    [s1, s2, s3] = ind2sub([n1 n2 n3], idx); t1 = double(th(idx));
    C = nan(numel(rUm), 1); np = zeros(numel(rUm), 1);
    for k = 1:numel(rUm)
        a = 2 * pi * rand(numel(idx), 1); d = rUm(k) / vox;
        p1 = round(s1 + d * cos(a)); p2 = round(s2 + d * sin(a));
        in = p1 >= 1 & p1 <= n1 & p2 >= 1 & p2 <= n2;
        li = sub2ind([n1 n2 n3], p1(in), p2(in), s3(in));
        ok = elig(li); t0 = t1(in);
        dt = t0(ok) - double(th(li(ok)));
        C(k) = mean(cos(2 * dt)); np(k) = numel(dt);
    end
    Cfar = mean(cos(2 * (t1 - double(th(idx(randperm(numel(idx))))))));
    rho  = (C - Cfar) / (1 - Cfar);
    S = struct('r_um', rUm(:), 'C', C, 'rho', rho, 'n_pairs', np, 'C_far', Cfar, ...
        'corrlen_um', cross_len(rUm(:), rho, exp(-1)));
end

function L = cross_len(r, rho, thr)
    L = NaN; k = find(rho < thr, 1);
    if isempty(k), return; end
    if k == 1, L = r(1); return; end
    L = r(k-1) + (thr - rho(k-1)) * (r(k) - r(k-1)) / (rho(k) - rho(k-1));
end

function vt = vessel_tensor(seg, aIdx, prm)
    Is = imgaussfilt3(single(seg), prm.sigmaD / prm.vox);
    [g2, g1, g3] = gradient(Is); clear Is          % gradient(): dim2, dim1, dim3
    G = {g1, g2, g3}; clear g1 g2 g3
    pairs = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
    nA = numel(aIdx); nS = numel(prm.sigmaI);
    vt = struct('t', nan(nA, 3), 'coh', -inf(nA, 1), 'sigma', nan(nA, 1), ...
                'tS', nan(nA, 3, nS), 'cohS', nan(nA, nS));
    for k = 1:nS
        Jv = zeros(nA, 6);
        for p = 1:6
            Pk = imgaussfilt3(G{pairs(p,1)} .* G{pairs(p,2)}, prm.sigmaI(k) / prm.vox);
            Jv(:, p) = double(Pk(aIdx));
        end
        clear Pk
        [t, coh] = tensor_axis(Jv);
        vt.tS(:,:,k) = t; vt.cohS(:,k) = coh;
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
            [Ve, L] = eig(J(:,:,i));
            [mi, o] = sort(diag(L), 'ascend');
            mu(i,:) = mi'; t(i,:) = Ve(:, o(1))';
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

function [r1, r2, r3] = bbox_ranges(m, pad)
    sz = size(m); if numel(sz) < 3, sz(3) = 1; end
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

% ---- statistics ----

function B = boot_stat(T, fun, nBoot, Kpt, Kb)
    % point estimate on T; CI from resampling spatial blocks (T.block) with replacement
    B.pt = fun(T, Kpt);
    B.reps = nan(numel(B.pt), nBoot);
    if height(T) > 0
        cells = block_cells(T.block);
        for b = 1:nBoot, B.reps(:, b) = fun(resample_blocks(T, cells), Kb); end
    end
    B.lo = prctile(B.reps, 2.5, 2); B.hi = prctile(B.reps, 97.5, 2);
end

function B = boot_pair(T1, T2, fun, nBoot, Kpt, Kb)
    B.pt = fun(T1, T2, Kpt);
    B.reps = nan(numel(B.pt), nBoot);
    if height(T1) > 0 && height(T2) > 0
        c1 = block_cells(T1.block); c2 = block_cells(T2.block);
        for b = 1:nBoot, B.reps(:, b) = fun(resample_blocks(T1, c1), resample_blocks(T2, c2), Kb); end
    end
    B.lo = prctile(B.reps, 2.5, 2); B.hi = prctile(B.reps, 97.5, 2);
end

function cells = block_cells(block)
    [~, ~, gi] = unique(block);
    cells = accumarray(gi, (1:numel(block))', [], @(v) {v});
end

function T2 = resample_blocks(T, cells)
    nU = numel(cells); pick = cells(randi(nU, nU, 1));
    T2 = T(vertcat(pick{:}), :);
    T2.block = repelem((1:nU)', cellfun(@numel, pick));   % each drawn block is its own unit
end

function p = perm_within(g)
    % random permutation that only exchanges elements with equal g
    n = numel(g);
    [~, o1] = sortrows([g(:), (1:n)']);
    [~, o2] = sortrows([g(:), rand(n, 1)]);
    p = zeros(n, 1); p(o1) = o2;
end

function v = stat_prim(T, K, minN)
    % [median dtheta; excess vs global shuffle; excess vs within-block shuffle; median R]
    v = nan(4, 1); n = height(T);
    if n < minN, return; end
    md = median(T.dtheta_deg); [~, ~, g] = unique(T.block);
    ng = zeros(K, 1); nb = zeros(K, 1);
    for k = 1:K
        ng(k) = median(axial_diff_deg(T.thetaF_deg(randperm(n)), T.thetaV_deg));
        nb(k) = median(axial_diff_deg(T.thetaF_deg(perm_within(g)), T.thetaV_deg));
    end
    v = [md; mean(ng) - md; mean(nb) - md; median(T.R)];
end

function v = stat_spatial(T, nD, K, minN)
    % excess of rigidly shifted fiber field over the field at the anchor, per shift distance
    v = nan(nD, 1);
    if height(T) < minN, return; end
    for j = 1:nD
        e = nan(K, 1);
        for k = 1:K
            c = T.dShift(:, (j-1)*K + k);
            m = isfinite(c) & isfinite(T.dField);
            if nnz(m) >= minN, e(k) = median(c(m)) - median(T.dField(m)); end
        end
        v(j) = mean(e, 'omitnan');
    end
end

function v = strat_excess(T, bins, nb, K, minN)
    % excess within each stratum (null shuffled within stratum)
    v = nan(nb, 1);
    for k = 1:nb
        i = bins == k; n = nnz(i);
        if n < minN, continue; end
        thF = T.thetaF_deg(i); thV = T.thetaV_deg(i);
        nm = zeros(K, 1);
        for q = 1:K, nm(q) = median(axial_diff_deg(thF(randperm(n)), thV)); end
        v(k) = mean(nm) - median(T.dtheta_deg(i));
    end
end

function v = stat_shell(T, nSh, x, K, minN)
    % [excess per shell; decay lambda; half-decay distance]
    e = strat_excess(T, T.shell, nSh, K, minN);
    p = decay_fit(x, e);
    v = [e; p(1); half_distance(x, e)];
end

function p = decay_fit(x, e)
    % e(x) = A exp(-x/lambda) + c; returns [lambda A c] (NaN if no decay)
    p = nan(1, 3);
    ok = isfinite(e); x = x(ok); e = e(ok);
    if numel(e) < 4, return; end
    c0 = mean(e(end-2:end)); A0 = e(1) - c0;
    if A0 <= 0, return; end
    f = @(q) sum((q(2) * exp(-x / exp(q(1))) + q(3) - e).^2);
    q = fminsearch(f, [log(100) A0 c0], optimset('Display', 'off'));
    lam = exp(q(1));
    if lam < 5 || lam > 1e4 || q(2) <= 0, return; end
    p = [lam q(2) q(3)];
end

function d = half_distance(x, e)
    % distance where excess falls halfway from the first shell to the mean of the last three
    d = NaN; ok = isfinite(e); x = x(ok); e = e(ok);
    if numel(e) < 4, return; end
    base = mean(e(end-2:end)); top = e(1);
    if top <= base, return; end
    h = base + 0.5 * (top - base);
    k = find(e < h, 1);
    if isempty(k) || k == 1, return; end
    d = x(k-1) + (h - e(k-1)) * (x(k) - x(k-1)) / (e(k) - e(k-1));
end

function v = stat_free(T, nSh, kP, kFar, minN)
    % [median R per shell; median retardance per shell; R near-far; ret near-far; R near]
    v = [grouped_median(T.shell, T.R, nSh, minN); grouped_median(T.shell, T.retMean, nSh, minN); NaN; NaN; NaN];
    near = T.shell <= kP; far = T.shell >= kFar;
    if nnz(near) >= minN && nnz(far) >= minN
        v(2*nSh + 1) = median(T.R(near)) - median(T.R(far));
        v(2*nSh + 2) = median(T.retMean(near)) - median(T.retMean(far));
    end
    if nnz(near) >= minN, v(2*nSh + 3) = median(T.R(near)); end
end

function v = stat_depth(T, edges, radMax, nMin, K, minN)
    % per depth bin: [excess; fraction radial vessels; median dtheta of radial in-plane
    % vessels; median fiber angle to the in-plane cortical normal], stacked by column
    nb = numel(edges) - 1; out = nan(nb, 4);
    b = discretize(T.depth_um, edges);
    prim = ~T.inclined; rad = T.radial_deg <= radMax;
    out(:,1) = strat_excess(T(prim,:), b(prim), nb, K, minN);
    for k = 1:nb
        i = b == k & isfinite(T.radial_deg);
        if nnz(i) >= minN, out(k,2) = mean(rad(i)); end
        i = b == k & prim & rad;
        if nnz(i) >= minN, out(k,3) = median(T.dtheta_deg(i)); end
        i = b == k & T.nInplane >= nMin & isfinite(T.fibN_deg);
        if nnz(i) >= minN, out(k,4) = median(T.fibN_deg(i)); end
    end
    v = out(:);
end

function v = stat_comp(Tf, To, K, rE, iE, minN)
    % [excess front; excess occip; excess occip reweighted to frontal radius x inclination;
    %  front - occip; front - occip adjusted]
    v = nan(5, 1);
    if height(Tf) < minN || height(To) < minN, return; end
    nb = (numel(rE) - 1) * (numel(iE) - 1);
    bf = comp_bin(Tf, rE, iE); bo = comp_bin(To, rE, iE);
    okf = isfinite(bf); oko = isfinite(bo);
    pf = accumarray(bf(okf), 1, [nb 1]) / nnz(okf);
    po = accumarray(bo(oko), 1, [nb 1]) / nnz(oko);
    w = zeros(height(To), 1); w(oko) = pf(bo(oko)) ./ po(bo(oko));
    eF = excess_w(Tf, ones(height(Tf), 1), K);
    eO = excess_w(To, ones(height(To), 1), K);
    eA = excess_w(To, w, K);
    v = [eF; eO; eA; eF - eO; eF - eA];
end

function b = comp_bin(T, rE, iE)
    b = (discretize(T.radius_um, rE) - 1) * (numel(iE) - 1) + discretize(T.incl_deg, iE);
end

function e = excess_w(T, w, K)
    n = height(T); md = wmedian(T.dtheta_deg, w); nm = zeros(K, 1);
    for k = 1:K, nm(k) = wmedian(axial_diff_deg(T.thetaF_deg(randperm(n)), T.thetaV_deg), w); end
    e = mean(nm) - md;
end

function m = wmedian(x, w)
    ok = isfinite(x) & w > 0; x = x(ok); w = w(ok);
    if isempty(x), m = NaN; return; end
    [x, o] = sort(x); cw = cumsum(w(o)); m = x(find(cw >= 0.5 * cw(end), 1));
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

function T = long_rows(VOLv, analysis, x, B, rows)
    n = numel(rows);
    T = table(repmat(string(VOLv.sub), n, 1), repmat(string(VOLv.reg), n, 1), repmat(string(analysis), n, 1), ...
        x(:), B.pt(rows), B.lo(rows), B.hi(rows), ...
        'VariableNames', {'subject','region','analysis','x','pt','lo','hi'});
end

% ---- figures ----

function [f, tl] = new_fig(nr, nc)
    f  = figure('Color', 'w', 'Position', [40 40 480*nc 430*nr]);
    tl = tiledlayout(f, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
end

function save_fig(f, outDir, name)
    exportgraphics(f, fullfile(outDir, [name '.png']), 'Resolution', 300);
    savefig(f, fullfile(outDir, [name '.fig']));
end

function h = plot_ci(ax, x, pt, lo, hi, col, mk, ls)
    x = x(:); pt = pt(:); lo = lo(:); hi = hi(:);
    ok = isfinite(pt) & isfinite(lo) & isfinite(hi);
    if nnz(ok) >= 2
        fill(ax, [x(ok); flipud(x(ok))], [lo(ok); flipud(hi(ok))], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    okp = isfinite(pt);
    h = plot(ax, x(okp), pt(okp), 'Color', col, 'LineStyle', ls, 'Marker', mk, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
end

function h = plot_pt(ax, x, pt, lo, hi, col, mk)
    h = errorbar(ax, x, pt, pt - lo, hi - pt, 'Marker', mk, 'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 7, 'LineWidth', 1.4, 'CapSize', 6, 'LineStyle', 'none');
end

function c = lighten(c, a)
    c = c + (1 - c) * a;
end

function L = edge_labels(e)
    L = cell(1, numel(e) - 1);
    for k = 1:numel(e) - 1
        if isinf(e(k+1)), L{k} = sprintf('>%g', e(k)); else, L{k} = sprintf('%g–%g', e(k), e(k+1)); end
    end
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