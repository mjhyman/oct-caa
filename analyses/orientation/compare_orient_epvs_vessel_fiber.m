%% compare_orient_epvs_vessel.m
% Fiber organisation (as-PSOCT optic axis, retardance) near EPVS vs near vessels in white
% matter of the control subjects (CAA26, CAA6), frontal and occipital.
%
% Run compare_orient_vessel_controls.m with RUN_MODE = 'check' first: ORIENT_ZERO_AXIS /
% ORIENT_CCW below must match the confirmed convention.
%
% EPVS surround vessels, so EPVS proximity is nested within vessel proximity. Every
% comparison below separates the two.
%
% Per fiber voxel (WM, not vessel/EPVS, retardance >= floor):
%   retardance
%   coherence R      |mean e^{i2theta}| within LOCAL_SIGMA_UM (1 = one direction)
%   crossing index X R4 - R^4, R4 = |mean e^{i4theta}|; ~0 for one fiber population,
%                    > 0 for two populations ~90 deg apart
%   orient gradient  in-plane rate of change of fiber angle (deg per 100 um)
%   psi              fiber angle to the in-plane normal of the nearest vessel / EPVS boundary
%                    (0 = radial, pointing at the structure; 90 = tangential, wrapping around)
%   frac tangential (psi > TANG_DEG) and radial (psi < RAD_DEG); random = 1/3, 1/3, 45 deg
%
% Analyses
%  1. Classes: EPVS (EPVS within NEAR_UM), vessel only (vessel within NEAR_UM, no EPVS within
%     EXCL_UM), far (both >= FAR_UM). Contrasts EPVS - vessel, EPVS - far, vessel - far.
%     Fibers within BOUNDARY_SKIP_UM of the boundary are excluded (partial volume); the
%     classes are repeated without the exclusion as a check.
%  2. 2-D maps (distance to vessel x distance to EPVS) and partial effects: EPVS effect at
%     fixed vessel distance (near - far EPVS within each vessel-distance bin, count-weighted
%     average) and the vessel effect at fixed EPVS distance.
%  3. Profiles vs distance in units of local radius and in um, vessel vs EPVS; size-matched
%     comparison (metrics by local radius of the nearest structure).
%  4. Alignment excess (structure-tensor axis vs fibers within PRIMARY_MAX_UM) for vessels
%     without EPVS, vessels with EPVS, and the EPVS axis; vessels with vs without EPVS
%     matched on radius x inclination.
%  5. EPVS objects (26-connected components): elongation axis vs shell fibers (object-level
%     excess), EPVS axis vs contained vessel axis, Spearman correlations of EPVS volume with
%     tangential fraction, retardance, coherence and crossing.
%  6. Placement: fiber features in a shell around EPVS (and around vessels) vs the same shell
%     rigidly shifted in-plane: are EPVS located in disorganised / crossing / high-gradient /
%     low-retardance white matter?
%
% Statistics: voxel metrics are means; anchor and object metrics are medians. CIs: BLOCK_UM
% blocks resampled within volume. Replicated = both subjects same sign with CIs excluding 0.
% Placement: observed vs the distribution of shifted shells. No population inference (n = 2).

%% 0. Config ---------------------------------------------------------------
clearvars -except caa
RUN_STAGE   = 'all';         % 'all' | 'figures' (redraw figures + tables from a saved results.mat)
RESULTS_MAT = '';            % for 'figures': path to results.mat written by a previous 'all' run
EXPLORATORY = true;

SUBJECTS    = {'CAA26','CAA6'};
REGIONS     = {'front','occip'};
REG_LABEL   = struct('front', 'Frontal', 'occip', 'Occipital');
REG_COLOR   = struct('front', [0.20 0.45 0.75], 'occip', [0.85 0.45 0.10]);
SCOL        = struct('vessel', [0.80 0.20 0.15], 'epvs', [0.15 0.45 0.80]);
SUBJ_MARKER = {'o','s'};
SUBJ_LS     = {'-','--'};

VOX_UM = 20;
ORIENT_UNITS     = 'deg';
ORIENT_ZERO_AXIS = 'dim2';
ORIENT_CCW       = true;

RET_FLOOR_PCTL = 10;          % percentile of pooled control WM retardance, per region
RET_SAMPLE     = 2e6;

% proximity classes and profiles
BOUNDARY_SKIP_UM = 40;
NEAR_UM          = 120;
EXCL_UM          = 200;       % "vessel only": no EPVS within this distance
FAR_UM           = 300;
MAX_DIST_UM      = 400;
D2_EDGES_UM      = [0 40 80 120 160 240 320 400];   % 2-D map bins (both axes)
NORM_EDGES       = [0 0.5 1 1.5 2 3 4 6];           % distance / local radius
DIST_EDGES_UM    = 0:40:400;
SIZE_EDGES_UM    = [0 30 50 80 150 Inf];            % local radius of the nearest structure
EPVS_SIZE_MASK   = 'epvs_or_vessel';                % EPVS size/axis from epvs | seg, or 'epvs'
MIN_N_BIN        = 50;                              % min voxels for a binned mean

% local fiber statistics and boundary normals
LOCAL_SIGMA_UM     = 40;
LOCAL_MIN_W        = 0.3;
NORMAL_SIGMA_UM    = 40;
NORMAL_MIN_INPLANE = 0.5;
TANG_DEG           = 60;
RAD_DEG            = 30;

% axes and alignment (as compare_orient_vessel_controls.m)
SIGMA_D_UM      = 20;
SIGMA_I_UM      = [40 80 160];
MIN_TUBE_COH    = 0.5;
RAD_SEARCH_UM   = 160;
MAX_INCL_DEG    = 45;
PRIMARY_MAX_UM  = 120;
EPVS_ADJ_UM     = 40;         % vessel wall with an EPVS voxel within this distance = "with EPVS"
MIN_VOX         = 5;
MIN_ROWS        = 20;
RADIUS_EDGES_UM = [0 30 50 80 Inf];
INCL_EDGES_DEG  = [0 15 30 45];

% EPVS objects
MIN_OBJ_VOX       = 27;
OBJ_SHELL_UM      = [40 160];
MIN_OBJ_SHELL_VOX = 20;
MIN_ELONG         = 1.5;      % sqrt(lambda1/lambda2) needed for an object axis
VESSEL_IN_EPVS_UM = 40;       % vessel voxels within this distance of an EPVS belong to it
MIN_OBJECTS       = 20;

% placement null
PLACE_SHELL_UM = [40 160];
PLACE_SHIFT_UM = [400 2000];
N_PLACE_SHIFT  = 50;

BLOCK_UM       = 400;
N_NULL_SHUFFLE = 10;
N_NULL_BOOT    = 3;

if EXPLORATORY
    N_BOOT = 100; MAX_VOX = 5e6; MAX_ANCHORS = 50000; PLACE_MAX_VOX = 5e5; TAG = 'EXPLORATORY';
else
    N_BOOT = 500; MAX_VOX = 2e7; MAX_ANCHORS = Inf;   PLACE_MAX_VOX = 2e6; TAG = 'MANUSCRIPT';
end

ROOT = '/projectnb/npbssmic/ns/CAA/';
if strcmp(RUN_STAGE, 'figures')
    assert(isfile(RESULTS_MAT), 'Set RESULTS_MAT to the results.mat of a previous run.');
    OUT_DIR = fullfile(fileparts(RESULTS_MAT), ['redraw_' datestr(now, 'yyyymmdd_HHMMSS')]);
else
    OUT_DIR = fullfile(ROOT, 'figures/orientation/orient_epvs_vessel', ...
        sprintf('%s_%s', TAG, datestr(now, 'yyyymmdd_HHMMSS')));
end
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
rng(1);

prm = struct('vox', VOX_UM, 'units', ORIENT_UNITS, 'zeroAxis', ORIENT_ZERO_AXIS, 'ccw', ORIENT_CCW, ...
    'skip', BOUNDARY_SKIP_UM, 'near', NEAR_UM, 'excl', EXCL_UM, 'far', FAR_UM, 'maxDist', MAX_DIST_UM, ...
    'd2Edges', D2_EDGES_UM, 'normEdges', NORM_EDGES, 'distEdges', DIST_EDGES_UM, 'sizeEdges', SIZE_EDGES_UM, ...
    'epvSizeMask', EPVS_SIZE_MASK, 'localSigma', LOCAL_SIGMA_UM, 'localMinW', LOCAL_MIN_W, ...
    'normalSigma', NORMAL_SIGMA_UM, 'normalMinInplane', NORMAL_MIN_INPLANE, 'tangDeg', TANG_DEG, ...
    'radDeg', RAD_DEG, 'sigmaD', SIGMA_D_UM, 'sigmaI', SIGMA_I_UM, 'minCoh', MIN_TUBE_COH, ...
    'radSearch', RAD_SEARCH_UM, 'maxIncl', MAX_INCL_DEG, 'primaryMax', PRIMARY_MAX_UM, ...
    'epvAdj', EPVS_ADJ_UM, 'minVox', MIN_VOX, 'minObjVox', MIN_OBJ_VOX, 'objShell', OBJ_SHELL_UM, ...
    'minObjShell', MIN_OBJ_SHELL_VOX, 'minElong', MIN_ELONG, 'vesInEpvs', VESSEL_IN_EPVS_UM, ...
    'placeShell', PLACE_SHELL_UM, 'placeShift', PLACE_SHIFT_UM, 'nPlace', N_PLACE_SHIFT, ...
    'placeMax', PLACE_MAX_VOX, 'block', BLOCK_UM, 'maxVox', MAX_VOX, 'maxAnchors', MAX_ANCHORS);

MNAME = {'retardance','coherence_R','crossing_index','orient_gradient','frac_tangential','frac_radial','angle_to_normal'};
MLAB  = {'Retardance','Fiber coherence R','Crossing index','Orientation gradient (deg/100 \mum)', ...
         'Fraction tangential','Fraction radial','Fiber angle to boundary normal (deg)'};
MREF  = [NaN NaN 0 NaN 1/3 1/3 45];
nM    = numel(MNAME);
PNAME = {'coherence_R','crossing_index','orient_gradient','retardance'};
PLAB  = {'coherence R','crossing index','gradient (deg/100 \mum)','retardance'};
CLS_NAME = {'EPVS','vessel_only','far'};
CON_NAME = {'EPVS_minus_vessel','EPVS_minus_far','vessel_minus_far'};

nS  = numel(SUBJECTS); nR = numel(REGIONS);
nH  = numel(D2_EDGES_UM) - 1; nbH = nH * nH;
nN  = numel(NORM_EDGES) - 1;
nDi = numel(DIST_EDGES_UM) - 1;
nZ  = numel(SIZE_EDGES_UM) - 1;
nearB = find(D2_EDGES_UM(1:end-1) >= BOUNDARY_SKIP_UM & D2_EDGES_UM(2:end) <= NEAR_UM);
farB  = find(D2_EDGES_UM(1:end-1) >= FAR_UM);
assert(~isempty(nearB) && ~isempty(farB), 'D2_EDGES_UM must contain near and far bins.');

VOL = struct('sub', {}, 'reg', {}, 'name', {}, 'col', {}, 'mk', {}, 'ls', {});
for s = 1:nS
    for r = 1:nR
        VOL(end+1) = struct('sub', SUBJECTS{s}, 'reg', REGIONS{r}, ...
            'name', sprintf('%s %s', SUBJECTS{s}, REG_LABEL.(REGIONS{r})), ...
            'col', REG_COLOR.(REGIONS{r}), 'mk', SUBJ_MARKER{s}, 'ls', SUBJ_LS{s}); %#ok<SAGROW>
    end
end
nV = numel(VOL);
VARIANTS = {'cls', 'clsB'}; VLAB = ["skip_boundary"; "include_boundary"];
AGRP = {'al_v0','vessel_noEPVS'; 'al_v1','vessel_with_EPVS'; 'al_e','EPVS_axis'};
STRN = {'vessel','epvs'};

if strcmp(RUN_STAGE, 'all')
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

%% 2. Retardance floor from pooled control WM (per region) -----------------------
retFloor = struct();
for r = 1:nR
    reg = REGIONS{r}; vals = [];
    for s = 1:nS
        Vt = get_vol(caa, SUBJECTS{s}, reg);
        idx = find(logical(Vt.mask_wm) & ~logical(Vt.seg) & ~logical(Vt.epvs));
        if numel(idx) > RET_SAMPLE, idx = idx(randperm(numel(idx), RET_SAMPLE)); end
        rv_ = double(Vt.ret_full(idx)); vals = [vals; rv_(isfinite(rv_))]; %#ok<AGROW>
        clear Vt idx rv_
    end
    if RET_FLOOR_PCTL > 0, retFloor.(reg) = prctile(vals, RET_FLOOR_PCTL); else, retFloor.(reg) = -Inf; end
    fprintf('%s retardance floor (p%g of control WM): %.4g\n', reg, RET_FLOOR_PCTL, retFloor.(reg));
end

%% 3. Per-volume extraction ------------------------------------------------------
D = cell(1, nV); CNT = table();
for v = 1:nV
    fprintf('%s ...\n', VOL(v).name);
    D{v} = analyze_volume(get_vol(caa, VOL(v).sub, VOL(v).reg), prm, retFloor.(VOL(v).reg));
    I = D{v}.info;
    fprintf(['  fiber voxels sampled %d | EPVS objects %d (axis valid %d) | vessel anchors %d (with EPVS %d) | ' ...
             'EPVS anchors %d\n'], I.n_fiber_vox, height(D{v}.OBJ), nnz(D{v}.OBJ.axisValid), height(D{v}.AV), ...
        nnz(D{v}.AV.epvAdj), height(D{v}.AE));
    if height(D{v}.OBJ) < MIN_OBJECTS
        warning('%s: only %d EPVS objects; object-level results will be unstable.', VOL(v).name, height(D{v}.OBJ));
    end
    CNT = [CNT; table(string(VOL(v).sub), string(VOL(v).reg), I.n_fiber_vox, I.n_epvs_vox, height(D{v}.OBJ), ...
        nnz(D{v}.OBJ.axisValid), nnz(isfinite(D{v}.OBJ.ves_angle_deg)), height(D{v}.AV), nnz(D{v}.AV.epvAdj), ...
        height(D{v}.AE), 'VariableNames', {'subject','region','n_fiber_vox_sampled','n_epvs_vox', ...
        'n_epvs_objects','n_objects_axis_valid','n_objects_with_vessel','n_vessel_anchors', ...
        'n_vessel_anchors_with_epvs','n_epvs_anchors'})]; %#ok<AGROW>
end

%% 4. Statistics -----------------------------------------------------------------
fCls  = @(S, N) stat_classes(S, N, MIN_N_BIN);
fMean = @(S, N) stat_means(S, N, MIN_N_BIN);
fHeat = @(S, N) stat_heat(S, N, nH, nearB, farB, MIN_N_BIN);
fPrim = @(T, K) stat_prim(T, K, MIN_ROWS);
fComp = @(T1, T2, K) stat_comp(T1, T2, K, RADIUS_EDGES_UM, INCL_EDGES_DEG, MIN_ROWS);
fObj  = @(T, K) stat_obj(T, K, MIN_ROWS);

ST = cell(1, nV);
for v = 1:nV
    fprintf('Bootstrapping %s ...\n', VOL(v).name);
    BS = D{v}.BS; B = struct();
    B.cls   = boot_binned(BS.cls,   fCls,  N_BOOT);
    B.clsB  = boot_binned(BS.clsB,  fCls,  N_BOOT);
    B.heat  = boot_binned(BS.heat,  fHeat, N_BOOT);
    for fn = {'normV','normE','distV','distE','sizeV','sizeE'}
        B.(fn{1}) = boot_binned(BS.(fn{1}), fMean, N_BOOT);
    end
    AV = D{v}.AV; AE = D{v}.AE;
    PV0 = AV(AV.valid & ~AV.inclined & ~AV.epvAdj, :);
    PV1 = AV(AV.valid & ~AV.inclined &  AV.epvAdj, :);
    PE  = AE(AE.valid & ~AE.inclined, :);
    B.al_v0 = boot_stat(PV0, fPrim, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    B.al_v1 = boot_stat(PV1, fPrim, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    B.al_e  = boot_stat(PE,  fPrim, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    B.match = boot_pair(PV1, PV0, fComp, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    B.obj   = boot_stat(D{v}.OBJ, fObj, N_BOOT, N_NULL_SHUFFLE, N_NULL_BOOT);
    ST{v} = B;
end
fprintf('Saving results ...\n');
save(fullfile(OUT_DIR, 'results.mat'), 'D', 'ST', 'CNT', 'retFloor', 'VOL', 'prm', '-v7.3');

else   % RUN_STAGE = 'figures': reuse extraction and bootstrap results
    R_ = load(RESULTS_MAT, 'D', 'ST', 'CNT', 'retFloor', 'VOL');
    assert(isequal({R_.VOL.name}, {VOL.name}), ...
        'SUBJECTS/REGIONS differ from the saved run: %s', strjoin({R_.VOL.name}, ', '));
    D = R_.D; ST = R_.ST; CNT = R_.CNT; retFloor = R_.retFloor; clear R_
    fprintf('Loaded results from %s\n', RESULTS_MAT);
end

%% 5. Figures --------------------------------------------------------------------
names = {VOL.name}; off = linspace(-0.27, 0.27, nV);

% Fig 1: classes (with and without boundary exclusion)
FLAB = {'fig1a_classes', 'fig1b_classes_include_boundary'};
for q = 1:2
    [f, tl] = new_fig(2, 4);
    for j = 1:nM
        ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); h = gobjects(nV, 1);
        for v = 1:nV
            B = ST{v}.(VARIANTS{q}); k = (j-1)*3 + (1:3);
            h(v) = plot_pt(ax, (1:3) + off(v), B.pt(k), B.lo(k), B.hi(k), VOL(v).col, VOL(v).mk);
        end
        if isfinite(MREF(j)), yline(ax, MREF(j), 'k:'); end
        set(ax, 'XTick', 1:3, 'XTickLabel', {'EPVS','Vessel only','Far'}); xlim(ax, [0.5 3.5]); ylabel(ax, MLAB{j});
        if j == 1, legend(ax, h, names, 'Location', 'best', 'Box', 'off'); end
    end
    ax = nexttile(tl); axis(ax, 'off');
    text(ax, 0, 0.5, sprintf(['EPVS: EPVS within %g \\mum\nVessel only: vessel within %g \\mum,\n  no EPVS within %g \\mum\n' ...
        'Far: both \\geq %g \\mum\nDotted: random expectation'], NEAR_UM, NEAR_UM, EXCL_UM, FAR_UM));
    title(tl, sprintf('Fiber metrics by proximity class (%s)', strrep(VLAB(q), '_', ' ')));
    save_fig(f, OUT_DIR, FLAB{q});
end

% Fig 2: partial effects
[f, tl] = new_fig(2, 4);
for j = 1:nM
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); h = gobjects(nV, 1);
    for v = 1:nV
        B = ST{v}.heat; k = nbH*nM + [j, nM + j];
        h(v) = plot_pt(ax, [1 2] + off(v), B.pt(k), B.lo(k), B.hi(k), VOL(v).col, VOL(v).mk);
    end
    yline(ax, 0, 'k:'); set(ax, 'XTick', [1 2], 'XTickLabel', {'EPVS | vessel dist.', 'Vessel | EPVS dist.'});
    xlim(ax, [0.5 2.5]); ylabel(ax, ['\Delta ' MLAB{j} ', near − far']);
    if j == 1, legend(ax, h, names, 'Location', 'best', 'Box', 'off'); end
end
ax = nexttile(tl); axis(ax, 'off');
text(ax, 0, 0.5, sprintf(['Near: %g–%g \\mum; far: \\geq %g \\mum.\nEPVS | vessel dist.: effect of EPVS\n' ...
    '  proximity holding vessel distance fixed.'], D2_EDGES_UM(nearB(1)), D2_EDGES_UM(nearB(end)+1), D2_EDGES_UM(farB(1))));
title(tl, 'Partial effects from the 2-D distance maps');
save_fig(f, OUT_DIR, 'fig2_partial_effects');

% Fig 3: 2-D maps
HM = [2 5 1 3];
[f, tl] = new_fig(numel(HM), nV);
for q = 1:numel(HM)
    j = HM(q); Ms = cell(1, nV); allv = [];
    for v = 1:nV
        Ms{v} = reshape(ST{v}.heat.pt((j-1)*nbH + (1:nbH)), nH, nH);   % (EPVS bin, vessel bin)
        allv = [allv; Ms{v}(isfinite(Ms{v}))]; %#ok<AGROW>
    end
    cl = [0 1]; if numel(allv) > 2, cl = prctile(allv, [2 98]); if diff(cl) <= 0, cl = cl + [-1 1] * 1e-3; end, end
    for v = 1:nV
        ax = nexttile(tl);
        imagesc(ax, Ms{v}, 'AlphaData', isfinite(Ms{v})); axis(ax, 'xy'); axis(ax, 'square'); clim(ax, cl);
        set(ax, 'XTick', 1:nH, 'XTickLabel', D2_EDGES_UM(1:end-1), 'YTick', 1:nH, 'YTickLabel', D2_EDGES_UM(1:end-1));
        xlabel(ax, 'Distance to vessel (\mum)'); ylabel(ax, 'Distance to EPVS (\mum)');
        title(ax, sprintf('%s: %s', VOL(v).name, MLAB{j}));
        if v == nV, colorbar(ax); end
    end
end
save_fig(f, OUT_DIR, 'fig3_maps2d');

% Fig 4-6: profiles, vessel vs EPVS
PM = [5 2 1 3];
prof_figure('fig4_profiles_radius_units', ST, VOL, 'normV', 'normE', nN, edge_labels(NORM_EDGES), ...
    'Distance / local radius', PM, MLAB, MREF, REGIONS, REG_LABEL, SCOL, OUT_DIR);
prof_figure('fig5_profiles_um', ST, VOL, 'distV', 'distE', nDi, edge_labels(DIST_EDGES_UM), ...
    'Distance from boundary (\mum)', PM, MLAB, MREF, REGIONS, REG_LABEL, SCOL, OUT_DIR);
prof_figure('fig6_size_matched', ST, VOL, 'sizeV', 'sizeE', nZ, edge_labels(SIZE_EDGES_UM), ...
    'Local radius of nearest structure (\mum)', PM, MLAB, MREF, REGIONS, REG_LABEL, SCOL, OUT_DIR);

% Fig 7: alignment
[f, tl] = new_fig(1, 3);
GL = {'Vessel, no EPVS', 'Vessel with EPVS', 'EPVS axis'}; GC = {SCOL.vessel, lighten(SCOL.vessel, 0.5), SCOL.epvs};
for c = 2:3
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); hg = gobjects(3, 1);
    for v = 1:nV
        for g = 1:3
            B = ST{v}.(AGRP{g,1});
            hg(g) = plot_pt(ax, v + (g - 2) * 0.22, B.pt(c), B.lo(c), B.hi(c), GC{g}, VOL(v).mk);
        end
    end
    yline(ax, 0, 'k:'); set(ax, 'XTick', 1:nV, 'XTickLabel', names); xtickangle(ax, 30); xlim(ax, [0.5 nV + 0.5]);
    ylabel(ax, 'Alignment excess (deg)');
    if c == 2, title(ax, 'Excess vs global shuffle'); else, title(ax, 'Excess vs within-block shuffle'); end
    legend(ax, hg, GL, 'Location', 'best', 'Box', 'off');
end
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); hg = gobjects(2, 1);
for v = 1:nV
    B = ST{v}.match;
    hg(1) = plot_pt(ax, v - 0.12, B.pt(4), B.lo(4), B.hi(4), [0.4 0.4 0.4], VOL(v).mk);
    hg(2) = plot_pt(ax, v + 0.12, B.pt(5), B.lo(5), B.hi(5), [0 0 0], VOL(v).mk);
end
yline(ax, 0, 'k:'); set(ax, 'XTick', 1:nV, 'XTickLabel', names); xtickangle(ax, 30); xlim(ax, [0.5 nV + 0.5]);
ylabel(ax, 'Excess: with EPVS − without (deg)'); title(ax, 'Vessels with vs without EPVS');
legend(ax, hg, {'Raw', 'Matched radius x inclination'}, 'Location', 'best', 'Box', 'off');
save_fig(f, OUT_DIR, 'fig7_alignment');

% Fig 8: EPVS objects
[f, tl] = new_fig(1, 4);
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
for v = 1:nV
    B = ST{v}.obj; plot_pt(ax, v, B.pt(2), B.lo(2), B.hi(2), VOL(v).col, VOL(v).mk);
end
yline(ax, 0, 'k:'); set(ax, 'XTick', 1:nV, 'XTickLabel', names); xtickangle(ax, 30); xlim(ax, [0.5 nV + 0.5]);
ylabel(ax, 'Object-level excess (deg)'); title(ax, 'EPVS elongation axis vs shell fibers');
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); h = gobjects(nV, 1);
for v = 1:nV
    a = sort(D{v}.OBJ.ves_angle_deg(isfinite(D{v}.OBJ.ves_angle_deg)));
    if isempty(a), a = NaN; end
    h(v) = plot(ax, a, (1:numel(a)) / numel(a), 'Color', VOL(v).col, 'LineStyle', VOL(v).ls, 'LineWidth', 1.5);
end
xlim(ax, [0 90]); xlabel(ax, 'EPVS axis vs contained vessel axis (deg)'); ylabel(ax, 'Cumulative fraction of EPVS');
legend(ax, h, names, 'Location', 'southeast', 'Box', 'off');
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
for v = 1:nV
    O = D{v}.OBJ;
    h(v) = scatter(ax, O.vol_um3, O.fracTang, 14, VOL(v).col, VOL(v).mk, 'filled', 'MarkerFaceAlpha', 0.4);
end
set(ax, 'XScale', 'log'); yline(ax, 1/3, 'k:'); ylim(ax, [0 1]);
xlabel(ax, 'EPVS volume (\mum^3)'); ylabel(ax, 'Fraction tangential (shell)');
legend(ax, h, names, 'Location', 'best', 'Box', 'off');
ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
SPN = {'vol ~ tangential', 'vol ~ retardance', 'vol ~ coherence', 'vol ~ crossing', 'elongation ~ \Delta\theta'};
for v = 1:nV
    B = ST{v}.obj; k = 5:9;
    h(v) = errorbar(ax, B.pt(k), (1:5)' + off(v), B.pt(k) - B.lo(k), B.hi(k) - B.pt(k), 'horizontal', ...
        'Marker', VOL(v).mk, 'Color', VOL(v).col, 'MarkerFaceColor', VOL(v).col, 'LineStyle', 'none', 'CapSize', 4);
end
xline(ax, 0, 'k:'); xlim(ax, [-1 1]); set(ax, 'YTick', 1:5, 'YTickLabel', SPN, 'YDir', 'reverse'); ylim(ax, [0.4 5.6]);
xlabel(ax, 'Spearman \rho across EPVS objects');
legend(ax, h, names, 'Location', 'best', 'Box', 'off');
save_fig(f, OUT_DIR, 'fig8_epvs_objects');

% Fig 9: placement
[f, tl] = new_fig(1, numel(PNAME));
for j = 1:numel(PNAME)
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); hs = gobjects(2, 1);
    for v = 1:nV
        for st = 1:2
            P = D{v}.PL.(STRN{st}); x = v + (st - 1.5) * 0.3;
            nu = P.null(:, j); nu = nu(isfinite(nu)); if isempty(nu), continue; end
            mn = mean(nu); q = prctile(nu, [2.5 97.5]);
            plot(ax, [x x], q - mn, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 6);
            hs(st) = plot(ax, x, P.obs(j) - mn, VOL(v).mk, 'MarkerFaceColor', SCOL.(STRN{st}), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineStyle', 'none');
        end
    end
    yline(ax, 0, 'k:'); set(ax, 'XTick', 1:nV, 'XTickLabel', names); xtickangle(ax, 30); xlim(ax, [0.5 nV + 0.5]);
    ylabel(ax, ['Observed − shifted: ' PLAB{j}]);
    if all(isgraphics(hs)), legend(ax, hs, {'Around vessels (no EPVS)', 'Around EPVS'}, 'Location', 'best', 'Box', 'off'); end
end
title(tl, sprintf('Placement: shell %g–%g \\mum vs shifted shells (gray: 95%% of shifts)', PLACE_SHELL_UM));
save_fig(f, OUT_DIR, 'fig9_placement');


%% 6. Tables ---------------------------------------------------------------------
% classes and class contrasts
CLS = table(); CON = table();
for q = 1:2
    for v = 1:nV
        B = ST{v}.(VARIANTS{q});
        for j = 1:nM
            for c = 1:3
                k = (j-1)*3 + c;
                CLS = [CLS; table(string(VOL(v).sub), string(VOL(v).reg), VLAB(q), string(CLS_NAME{c}), ...
                    string(MNAME{j}), B.pt(k), B.lo(k), B.hi(k), ...
                    'VariableNames', {'subject','region','variant','class','metric','pt','lo','hi'})]; %#ok<AGROW>
                k = 3*nM + (c-1)*nM + j;
                CON = [CON; table(string(VOL(v).sub), string(VOL(v).reg), VLAB(q), string(CON_NAME{c}), ...
                    string(MNAME{j}), B.pt(k), B.lo(k), B.hi(k), ...
                    'VariableNames', {'subject','region','variant','contrast','metric','pt','lo','hi'})]; %#ok<AGROW>
            end
        end
    end
end
CON = add_replication(CON, {'region','variant','contrast','metric'});

% partial effects from the 2-D maps
PART = table(); EFF_NAME = ["EPVS_effect_at_fixed_vessel_distance"; "vessel_effect_at_fixed_EPVS_distance"];
for v = 1:nV
    B = ST{v}.heat;
    for e = 1:2
        for j = 1:nM
            k = nbH*nM + (e-1)*nM + j;
            PART = [PART; table(string(VOL(v).sub), string(VOL(v).reg), EFF_NAME(e), string(MNAME{j}), ...
                B.pt(k), B.lo(k), B.hi(k), 'VariableNames', {'subject','region','effect','metric','pt','lo','hi'})]; %#ok<AGROW>
        end
    end
end
PART = add_replication(PART, {'region','effect','metric'});

% 2-D maps (point estimates and voxel counts)
HEAT = table();
for v = 1:nV
    Ncnt = bin_counts(D{v}.BS.heat);
    for j = 1:nM
        M = ST{v}.heat.pt((j-1)*nbH + (1:nbH));
        [iE_, iV_] = ndgrid(1:nH, 1:nH);
        HEAT = [HEAT; table(repmat(string(VOL(v).sub), nbH, 1), repmat(string(VOL(v).reg), nbH, 1), ...
            repmat(string(MNAME{j}), nbH, 1), D2_EDGES_UM(iV_(:))', D2_EDGES_UM(iE_(:))', M, Ncnt(:, j), ...
            'VariableNames', {'subject','region','metric','dist_vessel_lo_um','dist_epvs_lo_um','mean','n_vox'})]; %#ok<AGROW>
    end
end

% profiles
PROF = table();
PF = {'normV','vessel','dist_over_radius',NORM_EDGES; 'normE','EPVS','dist_over_radius',NORM_EDGES; ...
      'distV','vessel','dist_um',DIST_EDGES_UM;      'distE','EPVS','dist_um',DIST_EDGES_UM; ...
      'sizeV','vessel','radius_um',SIZE_EDGES_UM;    'sizeE','EPVS','radius_um',SIZE_EDGES_UM};
for v = 1:nV
    for p = 1:size(PF, 1)
        B = ST{v}.(PF{p,1}); ed = PF{p,4}; nb = numel(ed) - 1;
        for j = 1:nM
            k = (j-1)*nb + (1:nb);
            PROF = [PROF; table(repmat(string(VOL(v).sub), nb, 1), repmat(string(VOL(v).reg), nb, 1), ...
                repmat(string(PF{p,2}), nb, 1), repmat(string(PF{p,3}), nb, 1), ed(1:end-1)', ed(2:end)', ...
                repmat(string(MNAME{j}), nb, 1), B.pt(k), B.lo(k), B.hi(k), ...
                'VariableNames', {'subject','region','structure','axis','bin_lo','bin_hi','metric','pt','lo','hi'})]; %#ok<AGROW>
        end
    end
end

% alignment
ALN = table(); AMET = ["median_dtheta"; "excess_global"; "excess_block"; "median_R"];
for v = 1:nV
    for g = 1:size(AGRP, 1)
        B = ST{v}.(AGRP{g,1});
        ALN = [ALN; table(repmat(string(VOL(v).sub), 4, 1), repmat(string(VOL(v).reg), 4, 1), ...
            repmat(string(AGRP{g,2}), 4, 1), AMET, B.pt, B.lo, B.hi, ...
            'VariableNames', {'subject','region','group','metric','pt','lo','hi'})]; %#ok<AGROW>
    end
    d = ST{v}.al_e.reps - ST{v}.al_v0.reps;
    ALN = [ALN; table(repmat(string(VOL(v).sub), 4, 1), repmat(string(VOL(v).reg), 4, 1), ...
        repmat("EPVS_axis_minus_vessel_noEPVS", 4, 1), AMET, ST{v}.al_e.pt - ST{v}.al_v0.pt, ...
        prctile(d, 2.5, 2), prctile(d, 97.5, 2), 'VariableNames', ALN.Properties.VariableNames)]; %#ok<AGROW>
end
ALN = add_replication(ALN, {'region','group','metric'});

MATCH = table(); MQ = ["excess_vessel_with_EPVS"; "excess_vessel_noEPVS"; "excess_vessel_noEPVS_matched"; ...
                       "diff_raw"; "diff_matched"];
for v = 1:nV
    B = ST{v}.match;
    MATCH = [MATCH; table(repmat(string(VOL(v).sub), 5, 1), repmat(string(VOL(v).reg), 5, 1), MQ, B.pt, B.lo, B.hi, ...
        'VariableNames', {'subject','region','quantity','pt','lo','hi'})]; %#ok<AGROW>
end
MATCH = add_replication(MATCH, {'region','quantity'});

% objects
OSTAT = table(); ONAME = ["median_dtheta_axis_vs_fibers"; "object_excess"; "median_EPVS_vessel_axis_angle"; ...
    "frac_EPVS_vessel_angle_gt30"; "spearman_volume_tangential"; "spearman_volume_retardance"; ...
    "spearman_volume_coherence"; "spearman_volume_crossing"; "spearman_elongation_dtheta"];
OBJ = table();
for v = 1:nV
    B = ST{v}.obj;
    OSTAT = [OSTAT; table(repmat(string(VOL(v).sub), 9, 1), repmat(string(VOL(v).reg), 9, 1), ONAME, B.pt, B.lo, B.hi, ...
        'VariableNames', {'subject','region','stat','pt','lo','hi'})]; %#ok<AGROW>
    O = D{v}.OBJ; O.subject = repmat(string(VOL(v).sub), height(O), 1); O.region = repmat(string(VOL(v).reg), height(O), 1);
    OBJ = [OBJ; O]; %#ok<AGROW>
end

% placement
PLC = table();
for v = 1:nV
    for st = 1:2
        P = D{v}.PL.(STRN{st});
        for j = 1:numel(PNAME)
            nu = P.null(:, j); nu = nu(isfinite(nu));
            if isempty(nu), q = [NaN NaN]; mn = NaN; ge = NaN; else, q = prctile(nu, [2.5 97.5]); mn = mean(nu); ge = mean(nu >= P.obs(j)); end
            PLC = [PLC; table(string(VOL(v).sub), string(VOL(v).reg), string(STRN{st}), string(PNAME{j}), P.n, ...
                P.obs(j), mn, q(1), q(2), P.obs(j) - mn, ge, mean(P.kept, 'omitnan'), ...
                'VariableNames', {'subject','region','structure','feature','n_shell_vox','observed','null_mean', ...
                'null_lo','null_hi','enrichment','frac_null_ge_observed','shift_retention'})]; %#ok<AGROW>
        end
    end
end

disp(CNT); disp(CON(CON.variant == "skip_boundary", :)); disp(PART); disp(ALN); disp(MATCH); disp(OSTAT); disp(PLC);
writetable(CNT,   fullfile(OUT_DIR, 'counts.csv'));
writetable(CLS,   fullfile(OUT_DIR, 'class_metrics.csv'));
writetable(CON,   fullfile(OUT_DIR, 'class_contrasts.csv'));
writetable(PART,  fullfile(OUT_DIR, 'partial_effects.csv'));
writetable(HEAT,  fullfile(OUT_DIR, 'map2d_long.csv'));
writetable(PROF,  fullfile(OUT_DIR, 'profiles_long.csv'));
writetable(ALN,   fullfile(OUT_DIR, 'alignment.csv'));
writetable(MATCH, fullfile(OUT_DIR, 'alignment_matched.csv'));
writetable(OSTAT, fullfile(OUT_DIR, 'epvs_object_stats.csv'));
writetable(OBJ,   fullfile(OUT_DIR, 'epvs_objects.csv'));
writetable(PLC,   fullfile(OUT_DIR, 'placement.csv'));

fid = fopen(fullfile(OUT_DIR, 'README.txt'), 'w');
fprintf(fid, ['%s run %s\nSubjects %s; regions %s; retardance floor p%g of pooled control WM per region.\n' ...
    'Classes: EPVS = EPVS within %g um; vessel only = vessel within %g um and no EPVS within %g um; far = both >= %g um. ' ...
    'Fibers within %g um of the boundary excluded in the main variant.\n' ...
    'Local statistics: Gaussian sigma %g um. Crossing index = R4 - R^4. psi = fiber angle to the in-plane normal of the ' ...
    'nearest boundary (normal from distance map smoothed at %g um; in-plane fraction >= %g); tangential psi > %g, radial psi < %g.\n' ...
    'Partial effects: (mean at %s um from EPVS - mean at >= %g um) within each vessel-distance bin, count-weighted; and vice versa.\n' ...
    'Alignment: structure tensor (sigma_d %g, sigma_i %s um, coherence >= %g), fibers %g-%g um, inclination < %g deg. ' ...
    'excess = median(shuffled) - median(observed); > 0 = aligned.\n' ...
    'EPVS objects: 26-connected, >= %d voxels; axis = PCA major axis (elongation >= %g); shell %s um.\n' ...
    'Placement: shell %s um around structures vs %d in-plane shifts of %s um.\n' ...
    'CIs: %g um blocks resampled within volume (N_BOOT %d). Replicated = both subjects same sign with CIs excluding 0. ' ...
    'n = 2 subjects; CIs describe these specimens only.\n'], ...
    TAG, datestr(now), strjoin(SUBJECTS, ','), strjoin(REGIONS, ','), RET_FLOOR_PCTL, NEAR_UM, NEAR_UM, EXCL_UM, FAR_UM, ...
    BOUNDARY_SKIP_UM, LOCAL_SIGMA_UM, NORMAL_SIGMA_UM, NORMAL_MIN_INPLANE, TANG_DEG, RAD_DEG, ...
    mat2str(D2_EDGES_UM([nearB(1) nearB(end)+1])), D2_EDGES_UM(farB(1)), SIGMA_D_UM, mat2str(SIGMA_I_UM), ...
    MIN_TUBE_COH, BOUNDARY_SKIP_UM, PRIMARY_MAX_UM, MAX_INCL_DEG, MIN_OBJ_VOX, MIN_ELONG, mat2str(OBJ_SHELL_UM), ...
    mat2str(PLACE_SHELL_UM), N_PLACE_SHIFT, mat2str(PLACE_SHIFT_UM), BLOCK_UM, N_BOOT);
fclose(fid);

fprintf('Done. Output: %s\n', OUT_DIR);

%% ========================= Local functions ==================================

function V = get_vol(caa, sub, reg)
    sf = resolve_field(caa, sub, true);
    rf = resolve_field(caa.(sf), reg, false);
    V = caa.(sf).(rf);
end

function D = analyze_volume(V, prm, floorB)
    wm  = logical(V.mask_wm);
    pad = ceil(max([prm.maxDist, prm.excl, 3 * max(prm.sigmaI)]) / prm.vox) + 2;
    [r1, r2, r3] = bbox_ranges(wm, pad);
    wm  = wm(r1,r2,r3);
    seg = logical(V.seg(r1,r2,r3));
    epv = logical(V.epvs(r1,r2,r3));
    if ~any(epv(:)), error('No EPVS voxels in this volume; remove it from SUBJECTS/REGIONS.'); end
    if ~any(seg(:)), error('No vessel voxels in this volume.'); end
    [th, info] = decode_orient(V.orient(r1,r2,r3), prm.units, wm);
    th  = single(th);
    ret = single(V.ret_full(r1,r2,r3));
    sz  = size(wm); if numel(sz) < 3, sz(3) = 1; end
    elig = wm & ~seg & ~epv & isfinite(th) & isfinite(ret) & ret >= floorB;
    Bv  = max(1, round(prm.block / prm.vox)); nBk = ceil(sz / Bv);
    blockOf = @(a, b, c) sub2ind(nBk, min(max(ceil(a / Bv), 1), nBk(1)), min(max(ceil(b / Bv), 1), nBk(2)), ...
                                       min(max(ceil(c / Bv), 1), nBk(3)));
    if strcmp(prm.epvSizeMask, 'epvs_or_vessel'), objE = epv | seg; else, objE = epv; end
    info.n_epvs_vox = nnz(epv);

    % --- distance maps (um), nearest-structure indices, local radius maps ---
    [dV, iV] = bwdist(seg); dV = dV * prm.vox;
    [dE, iE] = bwdist(epv); dE = dE * prm.vox;
    se  = strel('cube', 2 * max(1, round(prm.radSearch / prm.vox)) + 1);
    rVm = imdilate(bwdist(~seg), se) * prm.vox;
    rEm = imdilate(bwdist(~objE), se) * prm.vox;

    % --- local fiber statistics (coherence, crossing, gradient) ---
    sg = prm.localSigma / prm.vox;
    w  = single(elig); t2 = 2 * th; t2(~elig) = 0;
    Fw  = imgaussfilt3(w, sg);
    Fc  = imgaussfilt3(w .* cos(t2), sg); Fs = imgaussfilt3(w .* sin(t2), sg);
    R2m = hypot(Fc, Fs) ./ Fw;
    uc  = Fc ./ hypot(Fc, Fs); us = Fs ./ hypot(Fc, Fs); clear Fc Fs
    F4c = imgaussfilt3(w .* cos(2 * t2), sg); F4s = imgaussfilt3(w .* sin(2 * t2), sg); clear w t2
    Xm  = hypot(F4c, F4s) ./ Fw - R2m.^4; clear F4c F4s
    Gm  = orient_gradient(uc, us, prm.vox); clear uc us
    bad = Fw < prm.localMinW | ~wm; R2m(bad) = NaN; Xm(bad) = NaN; Gm(bad) = NaN; clear Fw bad

    % --- boundary normals: smoothed distance maps ---
    sgN = prm.normalSigma / prm.vox;
    dVs = imgaussfilt3(dV, sgN); dEs = imgaussfilt3(dE, sgN);

    % --- fiber voxel sample ---
    vox = find(elig);
    if numel(vox) > prm.maxVox, vox = sort(vox(randperm(numel(vox), prm.maxVox))); end
    info.n_fiber_vox = numel(vox);
    [s1, s2, s3] = ind2sub(sz, vox);
    dv = double(dV(vox)); de = double(dE(vox));
    rv = double(rVm(iV(vox))); re = double(rEm(iE(vox)));
    thd  = rad2deg(double(th(vox)));
    base = [double(ret(vox)), double(R2m(vox)), double(Xm(vox)), double(Gm(vox))];
    blk  = blockOf(s1, s2, s3);
    psiV = normal_angle(grad_at(dVs, s1, s2, s3), thd, prm); psiV(dv > prm.maxDist) = NaN;
    psiE = normal_angle(grad_at(dEs, s1, s2, s3), thd, prm); psiE(de > prm.maxDist) = NaN;
    clear s1 s2 s3
    allv = true(numel(vox), 1);
    Ymat = @(psi, m) metric_matrix(base(m, :), psi(m), prm.tangDeg, prm.radDeg);   % rows m only

    % --- 1. classes (with and without boundary exclusion) ---
    skips = [prm.skip, 0]; fld = {'cls', 'clsB'};
    for q = 1:2
        sk = skips(q);
        cls = nan(numel(vox), 1);
        cls(dv >= prm.far & de >= prm.far) = 3;
        cls(dv >= sk & dv <= prm.near & de > prm.excl) = 2;
        cls(de >= sk & de <= prm.near) = 1;
        psiC = nan(numel(vox), 1); psiC(cls == 1) = psiE(cls == 1); psiC(cls == 2) = psiV(cls == 2);
        BS.(fld{q}) = binned_sums(blk, cls, 3, Ymat(psiC, allv));
    end

    % --- 2. 2-D maps: bin = (vessel bin - 1) * nH + EPVS bin; psi from the nearer structure ---
    nH = numel(prm.d2Edges) - 1;
    hb = (discretize(dv, prm.d2Edges) - 1) * nH + discretize(de, prm.d2Edges);
    psiN = psiV; useE = de < dv; psiN(useE) = psiE(useE);
    BS.heat = binned_sums(blk, hb, nH * nH, Ymat(psiN, allv));

    % --- 3. profiles (radius units, um) and size strata ---
    vo = de > prm.excl & dv <= prm.maxDist;      % vessel without EPVS nearby
    en = de <= prm.maxDist;
    BS.normV = binned_sums(blk(vo), discretize(dv(vo) ./ rv(vo), prm.normEdges), numel(prm.normEdges) - 1, Ymat(psiV, vo));
    BS.normE = binned_sums(blk(en), discretize(de(en) ./ re(en), prm.normEdges), numel(prm.normEdges) - 1, Ymat(psiE, en));
    BS.distV = binned_sums(blk(vo), discretize(dv(vo), prm.distEdges), numel(prm.distEdges) - 1, Ymat(psiV, vo));
    BS.distE = binned_sums(blk(en), discretize(de(en), prm.distEdges), numel(prm.distEdges) - 1, Ymat(psiE, en));
    nv = dv >= prm.skip & dv <= prm.near & de > prm.excl;
    ne = de >= prm.skip & de <= prm.near;
    BS.sizeV = binned_sums(blk(nv), discretize(rv(nv), prm.sizeEdges), numel(prm.sizeEdges) - 1, Ymat(psiV, nv));
    BS.sizeE = binned_sums(blk(ne), discretize(re(ne), prm.sizeEdges), numel(prm.sizeEdges) - 1, Ymat(psiE, ne));
    D.BS = BS;

    % --- 5. EPVS objects ---
    D.OBJ = obj_table(epv, seg, dE, iE, vox, de, thd, base, psiE, prm, blockOf);
    clear dv de rv re thd base blk psiV psiE vox

    % --- 6. placement ---
    feats = {R2m, Xm, Gm, ret};
    D.PL.vessel = place_test(elig & dV >= prm.placeShell(1) & dV <= prm.placeShell(2) & dE > prm.excl, elig, feats, prm);
    D.PL.epvs   = place_test(elig & dE >= prm.placeShell(1) & dE <= prm.placeShell(2), elig, feats, prm);
    clear feats R2m Xm Gm

    % --- 4. alignment anchors ---
    D.AV = anchor_table(seg,  dV, iV, dVs, rVm, elig, th, ret, prm, blockOf, dE);
    D.AE = anchor_table(objE, dE, iE, dEs, rEm, elig, th, ret, prm, blockOf, []);
    D.info = info;
end

function Y = metric_matrix(base, psi, tangDeg, radDeg)
    t = double(psi > tangDeg); r = double(psi < radDeg);
    t(~isfinite(psi)) = NaN; r(~isfinite(psi)) = NaN;
    Y = [base, t, r, psi];
end

function G = orient_gradient(uc, us, vox)
    % in-plane |d theta / dx| in deg per 100 um from the doubled-angle unit vector (uc, us)
    g2 = zeros(size(uc), 'single');
    for d = 1:2
        g2 = g2 + ((circshift(uc, -1, d) - circshift(uc, 1, d)) / 2).^2 ...
                + ((circshift(us, -1, d) - circshift(us, 1, d)) / 2).^2;
    end
    G = rad2deg(sqrt(g2) / 2) / vox * 100;
end

function psi = normal_angle(g, thd, prm)
    % axial angle (deg) between fiber and the in-plane projection of the boundary normal g
    nrm = vecnorm(g, 2, 2); nIn = hypot(g(:,1), g(:,2)) ./ nrm;
    thN = rad2deg(display_to_conv(atan2(-g(:,1), g(:,2)), prm.zeroAxis, prm.ccw));
    psi = axial_diff_deg(thd, thN);
    psi(~(nIn >= prm.normalMinInplane)) = NaN;
end

function BS = binned_sums(blk, bin, nb, Y)
    % per-block, per-bin sums and counts of each metric column (NaN ignored)
    ok = isfinite(bin) & bin >= 1 & bin <= nb;
    blk = blk(ok); bin = bin(ok); Y = Y(ok, :);
    [~, ~, bi] = unique(blk);
    nBlk = numel(unique(blk)); m = size(Y, 2);
    S = zeros(nBlk, nb, m); N = zeros(nBlk, nb, m);
    for j = 1:m
        if nBlk == 0, break; end
        f = isfinite(Y(:, j));
        if ~any(f), continue; end
        S(:,:,j) = accumarray([bi(f), bin(f)], Y(f, j), [nBlk nb]);
        N(:,:,j) = accumarray([bi(f), bin(f)], 1, [nBlk nb]);
    end
    BS = struct('S', S, 'N', N, 'nb', nb, 'm', m);
end

function Nc = bin_counts(BS)
    Nc = reshape(sum(reshape(BS.N, size(BS.N, 1), BS.nb * BS.m), 1), BS.nb, BS.m);
end

function B = boot_binned(BS, fun, nBoot)
    % statistic of block-summed means; CI by resampling blocks
    nBlk = size(BS.S, 1); nb = BS.nb; m = BS.m;
    if nBlk == 0
        B.pt = fun(zeros(nb, m), zeros(nb, m)); B.reps = nan(numel(B.pt), nBoot);
        B.lo = nan(size(B.pt)); B.hi = nan(size(B.pt)); return
    end
    S2 = reshape(BS.S, nBlk, nb * m); N2 = reshape(BS.N, nBlk, nb * m);
    B.pt = fun(reshape(sum(S2, 1), nb, m), reshape(sum(N2, 1), nb, m));
    B.reps = nan(numel(B.pt), nBoot);
    if nBlk > 0
        for b = 1:nBoot
            c = accumarray(randi(nBlk, nBlk, 1), 1, [nBlk 1])';
            B.reps(:, b) = fun(reshape(c * S2, nb, m), reshape(c * N2, nb, m));
        end
    end
    B.lo = prctile(B.reps, 2.5, 2); B.hi = prctile(B.reps, 97.5, 2);
end

function [M, Nn] = means_from(S, N, minN)
    M = S ./ N; M(N < minN) = NaN;
    Nn = N; Nn(N < minN) = 0;
end

function v = stat_means(S, N, minN)
    M = means_from(S, N, minN); v = M(:);
end

function v = stat_classes(S, N, minN)
    % [means (class x metric); EPVS - vessel; EPVS - far; vessel - far]
    M = means_from(S, N, minN);
    v = [M(:); (M(1,:) - M(2,:))'; (M(1,:) - M(3,:))'; (M(2,:) - M(3,:))'];
end

function v = stat_heat(S, N, nH, nearB, farB, minN)
    % [2-D means; EPVS effect at fixed vessel distance; vessel effect at fixed EPVS distance]
    [M, Nn] = means_from(S, N, minN); m = size(M, 2);
    Mr = reshape(M, nH, nH, m); Nr = reshape(Nn, nH, nH, m);     % (EPVS bin, vessel bin, metric)
    eE = nan(m, 1); eV = nan(m, 1);
    for j = 1:m
        eE(j) = partial_effect(Mr(:,:,j),  Nr(:,:,j),  nearB, farB);
        eV(j) = partial_effect(Mr(:,:,j)', Nr(:,:,j)', nearB, farB);
    end
    v = [M(:); eE; eV];
end

function e = partial_effect(M, N, nearB, farB)
    % near - far along dim 1 within each column, count-weighted average over columns
    ok = isfinite(M); Mz = M; Mz(~ok) = 0; Nz = N; Nz(~ok) = 0;
    nn = sum(Nz(nearB,:), 1); nf = sum(Nz(farB,:), 1);
    mn = sum(Mz(nearB,:) .* Nz(nearB,:), 1) ./ nn;
    mf = sum(Mz(farB,:)  .* Nz(farB,:),  1) ./ nf;
    use = nn > 0 & nf > 0; w = nn + nf;
    if ~any(use), e = NaN; return; end
    e = sum((mn(use) - mf(use)) .* w(use)) / sum(w(use));
end

function O = obj_table(epv, seg, dE, iE, vox, de, thd, base, psiE, prm, blockOf)
    sz = size(epv); if numel(sz) < 3, sz(3) = 1; end
    L = bwlabeln(epv, 26); nO = max(L(:));
    ie = find(epv); [o1, o2, o3] = ind2sub(sz, ie);
    [cnt, mu, C] = obj_cov(double(L(ie)), [o1 o2 o3], nO); clear ie o1 o2 o3
    sv = find(seg & dE <= prm.vesInEpvs); [v1, v2, v3] = ind2sub(sz, sv);
    [cntV, ~, CV] = obj_cov(double(L(iE(sv))), [v1 v2 v3], nO); clear sv v1 v2 v3

    m  = de >= prm.objShell(1) & de <= prm.objShell(2);
    ol = double(L(iE(vox(m))));
    acc  = @(x) accumarray(ol(isfinite(x)), x(isfinite(x)), [nO 1]);
    cntf = @(x) accumarray(ol(isfinite(x)), 1, [nO 1]);
    avg  = @(x) acc(x) ./ cntf(x);
    c2 = cosd(2 * thd(m)); s2 = sind(2 * thd(m));
    Cs = acc(c2); Ss = acc(s2); Ns = cntf(c2);
    b  = base(m, :);
    retM = avg(b(:,1)); R2M = avg(b(:,2)); XM = avg(b(:,3)); GM = avg(b(:,4));
    p  = psiE(m); tg = double(p > prm.tangDeg); tg(~isfinite(p)) = NaN; TG = avg(tg);

    keep = find(cnt >= prm.minObjVox); n = numel(keep);
    thO = nan(n, 1); inclO = nan(n, 1); elong = nan(n, 1); vesAng = nan(n, 1);
    for q = 1:n
        k = keep(q);
        [Ev, Lv] = eig(C(:,:,k)); [lam, o] = sort(real(diag(Lv)), 'descend'); a = real(Ev(:, o(1)));
        elong(q) = sqrt(max(lam(1), 0) / max(lam(2), eps));
        thO(q)   = rad2deg(display_to_conv(atan2(-a(1), a(2)), prm.zeroAxis, prm.ccw));
        inclO(q) = asind(min(1, abs(a(3))));
        if cntV(k) >= prm.minObjVox
            [Ew, Lw] = eig(CV(:,:,k)); [~, ow] = max(real(diag(Lw))); bv = real(Ew(:, ow));
            vesAng(q) = acosd(min(1, abs(a' * bv)));
        end
    end
    nSh = Ns(keep); okF = nSh >= prm.minObjShell;
    thF = mod(0.5 * atan2d(Ss(keep), Cs(keep)), 180); thF(~okF) = NaN;
    Rsh = hypot(Cs(keep), Ss(keep)) ./ nSh;           Rsh(~okF) = NaN;
    fx  = @(x) nan_where(x(keep), ~okF);
    vol = cnt(keep) * prm.vox^3;
    axisValid = okF & elong >= prm.minElong & inclO < prm.maxIncl;
    blk = blockOf(round(mu(keep,1)), round(mu(keep,2)), round(mu(keep,3)));
    O = table(blk, keep, cnt(keep), vol, (3 * vol / (4 * pi)).^(1/3), elong, thO, inclO, nSh, thF, Rsh, ...
        axial_diff_deg(thF, thO), axisValid, vesAng, fx(retM), fx(R2M), fx(XM), fx(GM), fx(TG), ...
        'VariableNames', {'block','object_id','nvox','vol_um3','eq_radius_um','elong','thetaO_deg','inclO_deg', ...
        'n_shell_vox','thetaF_deg','R_shell','dtheta_deg','axisValid','ves_angle_deg','retMean','R2mean', ...
        'Xmean','Gmean','fracTang'});
end

function x = nan_where(x, m)
    x(m) = NaN;
end

function [n, mu, C] = obj_cov(lab, X, nO)
    % per-label voxel count, centroid and coordinate covariance
    ok = lab >= 1; lab = lab(ok); X = X(ok, :);
    n = accumarray(lab, 1, [nO 1]);
    mu = zeros(nO, 3);
    for d = 1:3, mu(:, d) = accumarray(lab, X(:, d), [nO 1]) ./ n; end
    C = zeros(3, 3, nO);
    for a = 1:3
        for b = a:3
            cab = accumarray(lab, X(:, a) .* X(:, b), [nO 1]) ./ n - mu(:, a) .* mu(:, b);
            C(a, b, :) = cab; C(b, a, :) = cab;
        end
    end
end

function P = place_test(shellMask, elig, feats, prm)
    % mean features in the shell vs the same voxel set rigidly shifted in-plane
    sz = size(elig); n3 = size(elig, 3);
    idx = find(shellMask);
    if numel(idx) > prm.placeMax, idx = idx(randperm(numel(idx), prm.placeMax)); end
    nf = numel(feats);
    P = struct('n', numel(idx), 'obs', nan(1, nf), 'null', nan(prm.nPlace, nf), 'kept', nan(prm.nPlace, 1));
    if isempty(idx), return; end
    [q1, q2, q3] = ind2sub([sz(1) sz(2) n3], idx);
    for j = 1:nf, P.obs(j) = mean(double(feats{j}(idx)), 'omitnan'); end
    for k = 1:prm.nPlace
        ang = 2 * pi * rand;
        mag = (prm.placeShift(1) + rand * diff(prm.placeShift)) / prm.vox;
        p1 = round(q1 + mag * cos(ang)); p2 = round(q2 + mag * sin(ang));
        in = p1 >= 1 & p1 <= sz(1) & p2 >= 1 & p2 <= sz(2);
        li = sub2ind([sz(1) sz(2) n3], p1(in), p2(in), q3(in));
        li = li(elig(li));
        P.kept(k) = numel(li) / numel(idx);
        if numel(li) < 100, continue; end
        for j = 1:nf, P.null(k, j) = mean(double(feats{j}(li)), 'omitnan'); end
    end
end

function T = anchor_table(objMask, dMap, iMap, dS, rMap, elig, th, ret, prm, blockOf, dOther)
    % structure-tensor axis at boundary voxels vs fibers within [skip, primaryMax] anchored to them
    sz = size(elig); if numel(sz) < 3, sz(3) = 1; end
    vox = find(elig & dMap >= prm.skip & dMap <= prm.primaryMax);
    [aIdx, ~, pid] = unique(double(iMap(vox)));
    keep = subsample(true(numel(aIdx), 1), prm.maxAnchors); newId = cumsum(keep);
    kv = keep(pid); vox = vox(kv); pid = newId(pid(kv)); aIdx = aIdx(keep); nA = numel(aIdx);
    [v1, v2, v3] = ind2sub(sz, vox);
    thd = rad2deg(double(th(vox)));
    psi = normal_angle(grad_at(dS, v1, v2, v3), thd, prm); clear v1 v2 v3
    fin = isfinite(psi);
    Cs = accumarray(pid, cosd(2 * thd), [nA 1]); Ss = accumarray(pid, sind(2 * thd), [nA 1]);
    N  = accumarray(pid, 1, [nA 1]); RT = accumarray(pid, double(ret(vox)), [nA 1]);
    TG = accumarray(pid(fin), double(psi(fin) > prm.tangDeg), [nA 1]); NT = accumarray(pid(fin), 1, [nA 1]);
    vt = vessel_tensor(objMask, aIdx, prm);
    [thV, phi] = axis_angles(vt.t, prm);
    [a1, a2, a3] = ind2sub(sz, aIdx);
    thF = mod(rad2deg(0.5 * atan2(Ss, Cs)), 180);
    adj = false(nA, 1); if ~isempty(dOther), adj = double(dOther(aIdx)) <= prm.epvAdj; end
    valid = isfinite(vt.coh) & vt.coh >= prm.minCoh;
    T = table(blockOf(a1, a2, a3), thV, phi, phi >= prm.maxIncl, vt.coh, valid, double(rMap(aIdx)), adj, thF, ...
        hypot(Cs, Ss) ./ N, N, RT ./ N, TG ./ NT, axial_diff_deg(thF, thV), ...
        'VariableNames', {'block','thetaV_deg','incl_deg','inclined','tube_coh','valid','radius_um','epvAdj', ...
        'thetaF_deg','R','nvox','retMean','fracTang','dtheta_deg'});
    T = T(N >= prm.minVox, :);
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

function g = grad_at(D, a1, a2, a3)
    n1 = size(D, 1); n2 = size(D, 2); n3 = size(D, 3);
    cl = @(x, n) min(max(x, 1), n);
    at = @(i, j, k) double(D(sub2ind([n1 n2 n3], cl(i, n1), cl(j, n2), cl(k, n3))));
    g = [(at(a1+1, a2, a3) - at(a1-1, a2, a3)) / 2, ...
         (at(a1, a2+1, a3) - at(a1, a2-1, a3)) / 2, ...
         (at(a1, a2, a3+1) - at(a1, a2, a3-1)) / 2];
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

% ---- anchor / object statistics ----

function B = boot_stat(T, fun, nBoot, Kpt, Kb)
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
    T2.block = repelem((1:nU)', cellfun(@numel, pick));
end

function p = perm_within(g)
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

function v = stat_comp(T1, T2, K, rE, iE, minN)
    % [excess T1; excess T2; excess T2 reweighted to T1 radius x inclination; T1 - T2; T1 - T2 matched]
    v = nan(5, 1);
    if height(T1) < minN || height(T2) < minN, return; end
    nb = (numel(rE) - 1) * (numel(iE) - 1);
    b1 = comp_bin(T1, rE, iE); b2 = comp_bin(T2, rE, iE);
    ok1 = isfinite(b1); ok2 = isfinite(b2);
    p1 = accumarray(b1(ok1), 1, [nb 1]) / nnz(ok1);
    p2 = accumarray(b2(ok2), 1, [nb 1]) / nnz(ok2);
    w = zeros(height(T2), 1); w(ok2) = p1(b2(ok2)) ./ p2(b2(ok2));
    e1 = excess_w(T1, ones(height(T1), 1), K);
    e2 = excess_w(T2, ones(height(T2), 1), K);
    eA = excess_w(T2, w, K);
    v = [e1; e2; eA; e1 - e2; e1 - eA];
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

function v = stat_obj(T, K, minN)
    % [median dtheta (EPVS axis vs shell fibers); object-level excess; median EPVS-vessel axis angle;
    %  fraction angle > 30; Spearman volume~tangential, ~retardance, ~coherence, ~crossing; elongation~dtheta]
    v = nan(9, 1);
    a = T.axisValid; n = nnz(a);
    if n >= minN
        d = T.dtheta_deg(a); thF = T.thetaF_deg(a); thO = T.thetaO_deg(a);
        nm = zeros(K, 1);
        for k = 1:K, nm(k) = median(axial_diff_deg(thF(randperm(n)), thO)); end
        v(1) = median(d); v(2) = mean(nm) - median(d);
        v(9) = spear(T.elong(a), d, minN);
    end
    b = isfinite(T.ves_angle_deg);
    if nnz(b) >= minN, v(3) = median(T.ves_angle_deg(b)); v(4) = mean(T.ves_angle_deg(b) > 30); end
    v(5) = spear(T.vol_um3, T.fracTang, minN); v(6) = spear(T.vol_um3, T.retMean, minN);
    v(7) = spear(T.vol_um3, T.R2mean, minN);   v(8) = spear(T.vol_um3, T.Xmean, minN);
end

function r = spear(x, y, minN)
    ok = isfinite(x) & isfinite(y);
    if nnz(ok) < minN, r = NaN; else, r = corr(x(ok), y(ok), 'type', 'Spearman'); end
end

function T = add_replication(T, keys)
    % replicated = every subject's CI excludes 0 on the same side, within each key group
    [g, ~] = findgroups(T(:, keys));
    rep = splitapply(@(lo, hi) all(lo > 0) || all(hi < 0), T.lo, T.hi, g);
    T.replicated = rep(g);
end

% ---- figures ----

function [f, tl] = new_fig(nr, nc)
    f  = figure('Color', 'w', 'Position', [40 40 430*nc 380*nr]);
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
    h = errorbar(ax, x, pt, pt - lo, hi - pt, 'Marker', mk, 'Color', col, 'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'LineWidth', 1.4, 'CapSize', 6, 'LineStyle', 'none');
end

function prof_figure(fname, ST, VOL, fV, fE, nb, xlabs, xl, mets, MLAB, MREF, REGIONS, REG_LABEL, SCOL, outDir)
    nRg = numel(REGIONS); [f, tl] = new_fig(numel(mets), nRg);
    for q = 1:numel(mets)
        j = mets(q); k = (j-1)*nb + (1:nb);
        for r = 1:nRg
            ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off'); h = gobjects(0); lab = {};
            for v = find(strcmp({VOL.reg}, REGIONS{r}))
                B = ST{v}.(fV);
                h(end+1) = plot_ci(ax, (1:nb)', B.pt(k), B.lo(k), B.hi(k), SCOL.vessel, VOL(v).mk, VOL(v).ls); %#ok<AGROW>
                lab{end+1} = [VOL(v).sub ' vessel']; %#ok<AGROW>
                B = ST{v}.(fE);
                h(end+1) = plot_ci(ax, (1:nb)', B.pt(k), B.lo(k), B.hi(k), SCOL.epvs, VOL(v).mk, VOL(v).ls); %#ok<AGROW>
                lab{end+1} = [VOL(v).sub ' EPVS']; %#ok<AGROW>
            end
            if isfinite(MREF(j)), yline(ax, MREF(j), 'k:'); end
            set(ax, 'XTick', 1:nb, 'XTickLabel', xlabs); xtickangle(ax, 30); xlim(ax, [0.5 nb + 0.5]);
            xlabel(ax, xl); ylabel(ax, MLAB{j}); title(ax, REG_LABEL.(REGIONS{r}));
            legend(ax, h, lab, 'Location', 'best', 'Box', 'off');
        end
    end
    save_fig(f, outDir, fname);
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