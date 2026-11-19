%% compare_orient_vessel_controls.m
% Compare as-PSOCT optic-axis (fiber) orientation with local vessel orientation
% in white matter of the control subjects (CAA26, CAA6), frontal and occipital.
%
%   RUN_MODE = 'check' : orient units/convention diagnostics only (run first)
%   RUN_MODE = 'full'  : full analysis, figures, CSVs
%
% Geometry: dims 1-2 = en-face plane, dim 3 = through-plane (serial-section axis).
% The optic axis is an in-plane (axial) angle, so vessel axes are projected onto
% the en-face plane; steeply inclined vessels are analysed as a separate stratum.
%
% Vessel orientation: 3-D structure tensor of the Gaussian-smoothed vessel mask,
% evaluated at vessel-wall voxels. Axis = eigenvector of the smallest eigenvalue
% (direction of least mask-gradient variation). Integration scale is selected per
% wall voxel from SIGMA_I_UM as the one maximising tube coherence
% (mu2 - mu1)/(mu2 + mu1), mu1 <= mu2 <= mu3.
%
% Anchors: each perivascular fiber voxel is assigned to its nearest vessel-wall
% voxel (bwdist index on seg). Unit of observation = wall anchor (pseudoreplicated
% within subject). CIs: bootstrap over spatial blocks (BLOCK_UM cubes) within
% subject; pooled = subjects, then blocks. CIs describe the measured subjects only.

%% 0. Config ---------------------------------------------------------------
RUN_MODE    = 'full';        % 'check' | 'full'
EXPLORATORY = true;

SUBJECTS = {'CAA26','CAA6','CAA22','CAA25'};
REGIONS  = {'front','occip'};
REG_LABEL = struct('front','Frontal','occip','Occipital');

VOX_UM = 20;

% Orientation encoding of caa.(sub).(reg).orient -- confirm with RUN_MODE='check'
ORIENT_UNITS     = 'deg';    % 'auto' | 'deg' | 'rad'
ORIENT_ZERO_AXIS = 'dim2';    % 0 deg along 'dim2' (image columns) or 'dim1' (image rows)
ORIENT_CCW       = true;      % positive = counter-clockwise as displayed by imagesc

% Structure tensor
SIGMA_D_UM   = 20;            % derivative (pre-smoothing) scale
SIGMA_I_UM   = [40 80 160];   % candidate integration scales (per-anchor max coherence)
MIN_TUBE_COH = 0.5;           % (mu2-mu1)/(mu2+mu1) required for a valid axis
RAD_SEARCH_UM = 160;          % local vessel radius = max interior distance within this ball

SHELL_UM        = 40;         % nearest-structure shell width (matches annular analysis)
MAX_DIST_UM     = 400;        % outermost shell
PRIMARY_MAX_UM  = 120;        % fibers within this distance define the primary Δθ
MAX_INCL_DEG    = 45;         % primary analysis: vessels within 45° of the en-face plane
MIN_VOX         = 5;          % min fiber voxels per (anchor, shell) and per anchor (primary)
BLOCK_UM        = 400;        % spatial block size = bootstrap resampling unit
N_NULL_SHUFFLE  = 10;
RADIUS_EDGES_UM = [0 30 50 80 Inf];
COHER_EDGES     = [0 0.25 0.5 0.75 1.0001];

if EXPLORATORY
    N_BOOT = 100; MAX_ANCHORS = 50000; TAG = 'EXPLORATORY';
else
    N_BOOT = 500; MAX_ANCHORS = Inf;   TAG = 'MANUSCRIPT';
end

ROOT = '/projectnb/npbssmic/ns/CAA/';
OUT_DIR = fullfile(ROOT, 'figures/orientation/orient_vessel', ...
    sprintf('%s_%s_%s', TAG, RUN_MODE, datestr(now, 'yyyymmdd_HHMMSS')));
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
rng(1);

prm = struct('vox', VOX_UM, 'units', ORIENT_UNITS, 'zeroAxis', ORIENT_ZERO_AXIS, ...
    'ccw', ORIENT_CCW, 'sigmaD', SIGMA_D_UM, 'sigmaI', SIGMA_I_UM, 'minCoh', MIN_TUBE_COH, ...
    'radSearch', RAD_SEARCH_UM, 'shell', SHELL_UM, 'maxDist', MAX_DIST_UM, ...
    'primaryMax', PRIMARY_MAX_UM, 'maxIncl', MAX_INCL_DEG, 'minVox', MIN_VOX, ...
    'block', BLOCK_UM, 'maxPts', MAX_ANCHORS, 'nBoot', N_BOOT);

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
for s = 1:numel(SUBJECTS)
    key = regexprep(lower(SUBJECTS{s}), '[^a-z0-9]', '');     % 'CAA26' -> 'caa26'
    have = ~isempty(fieldnames(caa)) && ...
           any(strcmp(regexprep(lower(fieldnames(caa)), '[^a-z0-9]', ''), key));
    if ~have
        fprintf('Loading %s\n', upper(key));
        tmp = load(fullfile(ROOT, MAT_PATH.(key)), key);
        caa.(key) = tmp.(key); clear tmp
    end
end

%% 2. Check mode: orient units/convention + vessel-axis glyphs ------------------
if strcmp(RUN_MODE, 'check')
    for s = 1:numel(SUBJECTS)
        for r = 1:numel(REGIONS)
            sf = resolve_field(caa, SUBJECTS{s}, true);
            rf = resolve_field(caa.(sf), REGIONS{r}, false);
            run_orient_check(caa.(sf).(rf), prm, SUBJECTS{s}, REGIONS{r}, OUT_DIR);
        end
    end
    fprintf(['\nCheck figures in %s\nCyan (optic axis) should follow visible tract texture; ' ...
             'if rotated/mirrored, change ORIENT_ZERO_AXIS / ORIENT_CCW.\n' ...
             'Orange (structure-tensor vessel axis) should run along in-plane vessels.\n'], OUT_DIR);
    return
end

%% 3. Per subject x region analysis ------------------------------------------
R = struct(); allPt = table(); allSh = table();
nSh = ceil(MAX_DIST_UM / SHELL_UM);
for r = 1:numel(REGIONS)
    for s = 1:numel(SUBJECTS)
        sub = SUBJECTS{s}; reg = REGIONS{r};
        sf = resolve_field(caa, sub, true);
        rf = resolve_field(caa.(sf), reg, false);
        fprintf('%s %s ...\n', sub, reg);

        [pt, sh, info] = analyze_volume(caa.(sf).(rf), prm);
        pt.subject = repmat(string(sub), height(pt), 1);  pt.region = repmat(string(reg), height(pt), 1);
        sh.subject = repmat(string(sub), height(sh), 1);  sh.region = repmat(string(reg), height(sh), 1);
        fprintf('  anchors: %d valid (%d inclined); median tube coherence %.2f\n', ...
            height(pt), nnz(pt.inclined), median(pt.tube_coh));

        D = struct('pt', pt, 'sh', sh, 'info', info);
        prim = pt(~pt.inclined, :);
        D.obs = prim.dtheta_deg;

        % shuffled null: fiber angles permuted across anchors (preserves the region's
        % bulk fiber anisotropy and the vessel-axis distribution)
        nul = zeros(height(prim) * N_NULL_SHUFFLE, 1);
        for k = 1:N_NULL_SHUFFLE
            ii = (k-1)*height(prim) + (1:height(prim));
            nul(ii) = axial_diff_deg(prim.thetaF_deg(randperm(height(prim))), prim.thetaV_deg);
        end
        D.null = nul;

        shp = sh(~sh.inclined, :);
        D.shellStats = binned_stats(shp.shell, shp.dtheta_deg, shp.block, 0.5:1:nSh+0.5, N_BOOT);
        D.shellStats.x = (D.shellStats.bin - 0.5) * SHELL_UM;
        D.radStats = binned_stats(prim.radius_um, prim.dtheta_deg, prim.block, RADIUS_EDGES_UM, N_BOOT);
        D.cohStats = binned_stats(prim.R, prim.dtheta_deg, prim.block, COHER_EDGES, N_BOOT);

        R.(reg).(sub) = D;
        allPt = [allPt; pt]; allSh = [allSh; sh]; %#ok<AGROW>
    end
end

%% 4. Summary statistics -------------------------------------------------------
summ = table();
for r = 1:numel(REGIONS)
    reg = REGIONS{r};
    for s = 1:numel(SUBJECTS)
        sub = SUBJECTS{s}; D = R.(reg).(sub);
        prim = D.pt(~D.pt.inclined, :); incl = D.pt(D.pt.inclined, :);
        ci = boot_by_unit(D.obs, prim.block, N_BOOT);
        summ = [summ; table(string(sub), string(reg), string(D.info.units), ...
            height(prim), numel(unique(prim.block)), median(D.obs), ci(1), ci(2), ...
            mean(D.obs < 30), mean(D.obs > 60), median(D.null), cliffs_delta(D.obs, D.null), ...
            height(incl), median(incl.dtheta_deg), ...
            'VariableNames', {'subject','region','orient_units','n_anchors','n_blocks', ...
            'median_dtheta','ci_lo','ci_hi','frac_lt30','frac_gt60','null_median', ...
            'cliff_obs_vs_null','n_inclined','median_dtheta_inclined'})]; %#ok<AGROW>
    end
    % pooled controls: hierarchical bootstrap (subjects, then blocks)
    obs = []; nul = []; subjID = []; blkID = [];
    for s = 1:numel(SUBJECTS)
        D = R.(reg).(SUBJECTS{s}); prim = D.pt(~D.pt.inclined, :);
        obs = [obs; D.obs]; nul = [nul; D.null]; %#ok<AGROW>
        subjID = [subjID; s*ones(numel(D.obs),1)]; blkID = [blkID; prim.block]; %#ok<AGROW>
    end
    ci = boot_hier(obs, subjID, blkID, N_BOOT);
    summ = [summ; table("controls_pooled", string(reg), "", numel(obs), NaN, median(obs), ...
        ci(1), ci(2), mean(obs < 30), mean(obs > 60), median(nul), cliffs_delta(obs, nul), ...
        NaN, NaN, 'VariableNames', summ.Properties.VariableNames)]; %#ok<AGROW>
end

% frontal vs occipital (negative = frontal more aligned)
fo = table();
for s = 1:numel(SUBJECTS)
    a = R.front.(SUBJECTS{s}).obs; b = R.occip.(SUBJECTS{s}).obs;
    fo = [fo; table(string(SUBJECTS{s}), median(a), median(b), cliffs_delta(a, b), ...
        'VariableNames', {'subject','median_front','median_occip','cliff_front_vs_occip'})]; %#ok<AGROW>
end
a = [R.front.(SUBJECTS{1}).obs; R.front.(SUBJECTS{2}).obs];
b = [R.occip.(SUBJECTS{1}).obs; R.occip.(SUBJECTS{2}).obs];
fo = [fo; table("controls_pooled", median(a), median(b), cliffs_delta(a, b), ...
    'VariableNames', fo.Properties.VariableNames)];

disp(summ); disp(fo);
writetable(summ,  fullfile(OUT_DIR, 'summary_dtheta.csv'));
writetable(fo,    fullfile(OUT_DIR, 'front_vs_occip.csv'));
writetable(allPt, fullfile(OUT_DIR, 'anchors_primary.csv'));
writetable(allSh, fullfile(OUT_DIR, 'anchors_by_shell.csv'));

fid = fopen(fullfile(OUT_DIR, 'README.txt'), 'w');
fprintf(fid, ['%s run %s\nVessel axis: structure tensor of seg (sigma_d = %g um, sigma_i in %s um, ' ...
    'per-anchor max tube coherence >= %g), evaluated at vessel-wall voxels.\n' ...
    'Δθ = |axial difference| between en-face vessel axis and local optic-axis circular mean, ' ...
    'in [0,90] deg (0 = parallel, 90 = perpendicular; 45 = random).\n' ...
    'Primary: fibers within %d um of vessel wall, vessels inclined < %d deg from en-face plane.\n' ...
    'cliff_obs_vs_null < 0: vessels more aligned with fibers than shuffled pairing.\n' ...
    'Pseudoreplication: anchors nested within %d um spatial blocks within 2 subjects; CIs are ' ...
    'block bootstrap (pooled: subject then block) and do not generalise beyond these subjects.\n' ...
    'Orientation convention: zero axis %s, ccw = %d, units %s.\nN_BOOT = %d, MAX_ANCHORS = %g.\n'], ...
    TAG, datestr(now), SIGMA_D_UM, mat2str(SIGMA_I_UM), MIN_TUBE_COH, PRIMARY_MAX_UM, ...
    MAX_INCL_DEG, BLOCK_UM, ORIENT_ZERO_AXIS, ORIENT_CCW, ORIENT_UNITS, N_BOOT, MAX_ANCHORS);
fclose(fid);

%% 5. Figures: main tiled + one figure per panel -------------------------------
PANELS = {'violin','shell','radius','coherence'};
PTITLE = {'\Delta\theta: observed vs shuffled', '\Delta\theta vs distance from vessel wall', ...
          '\Delta\theta vs vessel radius', '\Delta\theta vs local fiber coherence'};
figM = figure('Color', 'w', 'Position', [40 40 1900 850]);
tl = tiledlayout(figM, numel(REGIONS), numel(PANELS), 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(REGIONS)
    for p = 1:numel(PANELS)
        ttl = sprintf('%s — %s', REG_LABEL.(REGIONS{r}), PTITLE{p});
        ax = nexttile(tl);
        draw_panel(ax, PANELS{p}, R.(REGIONS{r}), SUBJECTS, prm, RADIUS_EDGES_UM, COHER_EDGES);
        title(ax, ttl);

        f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 560 460]);
        ax2 = axes(f);
        draw_panel(ax2, PANELS{p}, R.(REGIONS{r}), SUBJECTS, prm, RADIUS_EDGES_UM, COHER_EDGES);
        title(ax2, ttl);
        fn = sprintf('panel_%s_%s', REGIONS{r}, PANELS{p});
        exportgraphics(f, fullfile(OUT_DIR, [fn '.png']), 'Resolution', 300);
        savefig(f, fullfile(OUT_DIR, [fn '.fig'])); close(f);
    end
end
exportgraphics(figM, fullfile(OUT_DIR, 'main_orient_vessel.png'), 'Resolution', 300);
savefig(figM, fullfile(OUT_DIR, 'main_orient_vessel.fig'));
fprintf('Done. Output: %s\n', OUT_DIR);

%% ========================= Local functions ==================================

function [pt, sh, info] = analyze_volume(V, prm)
    wm = logical(V.mask_wm); seg = logical(V.seg); epv = logical(V.epvs);
    [th, info] = decode_orient(V.orient, prm.units, wm);

    pad = ceil(max(prm.maxDist, 3 * max(prm.sigmaI)) / prm.vox) + 2;
    [r1, r2, r3] = bbox_ranges(wm, pad);
    wm = wm(r1,r2,r3); seg = seg(r1,r2,r3); epv = epv(r1,r2,r3); th = th(r1,r2,r3);

    % fiber voxels and their nearest vessel-wall voxel (= anchor)
    [dWall, idxW] = bwdist(seg);
    dWall = dWall * prm.vox;
    fib = wm & ~seg & ~epv & isfinite(th) & dWall > 0 & dWall <= prm.maxDist;
    vox = find(fib); clear fib
    [aIdx, ~, pid] = unique(double(idxW(vox))); clear idxW
    nA = numel(aIdx);

    % vessel axis at anchors
    vt  = vessel_tensor(seg, aIdx, prm);
    t   = vt.t;
    phi = rad2deg(asin(min(1, abs(t(:,3)))));
    thV = rad2deg(display_to_conv(atan2(-t(:,1), t(:,2)), prm.zeroAxis, prm.ccw));
    okA = isfinite(vt.coh) & vt.coh >= prm.minCoh;

    % local vessel radius: max interior distance within RAD_SEARCH of the wall voxel
    dIn  = bwdist(~seg);
    rmap = imdilate(dIn, strel('sphere', max(1, round(prm.radSearch / prm.vox)))); clear dIn
    rad  = double(rmap(aIdx)) * prm.vox; clear rmap

    % optional anchor subsample (EXPLORATORY)
    keepA = okA;
    if isfinite(prm.maxPts) && nnz(okA) > prm.maxPts
        fA = find(okA); keepA(:) = false;
        keepA(fA(randperm(numel(fA), prm.maxPts))) = true;
    end
    kv = keepA(pid); vox = vox(kv); pid = pid(kv);

    nSh = ceil(prm.maxDist / prm.shell);
    shell = min(nSh, ceil(double(dWall(vox)) / prm.shell)); clear dWall
    key = (pid - 1) * nSh + shell;
    sz  = [nA * nSh, 1];
    C = reshape(accumarray(key, cos(2*th(vox)), sz), nSh, nA)';
    S = reshape(accumarray(key, sin(2*th(vox)), sz), nSh, nA)';
    N = reshape(accumarray(key, 1, sz),               nSh, nA)';

    [a1, a2, a3] = ind2sub(size(seg), aIdx);
    B   = max(1, round(prm.block / prm.vox));
    blk = sub2ind(ceil(size(seg) / B), ceil(a1 / B), ceil(a2 / B), ceil(a3 / B));
    xyz = [a1 + r1(1) - 1, a2 + r2(1) - 1, a3 + r3(1) - 1];   % original voxel coordinates
    incl = phi >= prm.maxIncl;

    % primary: fibers pooled over shells <= primaryMax
    kP = round(prm.primaryMax / prm.shell);
    Cp = sum(C(:,1:kP), 2); Sp = sum(S(:,1:kP), 2); Np = sum(N(:,1:kP), 2);
    ok = Np >= prm.minVox;
    thF = mod(rad2deg(0.5 * atan2(Sp, Cp)), 180);
    pt = table(find(ok), blk(ok), xyz(ok,1), xyz(ok,2), xyz(ok,3), thV(ok), phi(ok), ...
        rad(ok), vt.coh(ok), vt.sigma(ok), thF(ok), hypot(Cp(ok), Sp(ok)) ./ Np(ok), Np(ok), ...
        axial_diff_deg(thF(ok), thV(ok)), incl(ok), ...
        'VariableNames', {'anchor','block','x1','x2','x3','thetaV_deg','incl_deg', ...
        'radius_um','tube_coh','sigmaI_um','thetaF_deg','R','nvox','dtheta_deg','inclined'});

    % per shell
    [pp, ss] = find(N >= prm.minVox);
    li = sub2ind(size(N), pp, ss);
    thFs = mod(rad2deg(0.5 * atan2(S(li), C(li))), 180);
    sh = table(pp, blk(pp), ss, (ss - 0.5) * prm.shell, thV(pp), phi(pp), rad(pp), ...
        thFs, hypot(C(li), S(li)) ./ N(li), N(li), axial_diff_deg(thFs, thV(pp)), incl(pp), ...
        'VariableNames', {'anchor','block','shell','dist_um','thetaV_deg','incl_deg', ...
        'radius_um','thetaF_deg','R','nvox','dtheta_deg','inclined'});
end

function vt = vessel_tensor(seg, aIdx, prm)
    % Structure tensor J = G_sigmaI * (grad I grad I^T), I = G_sigmaD * seg,
    % sampled at aIdx for each integration scale; keep the scale with max coherence.
    Is = imgaussfilt3(single(seg), prm.sigmaD / prm.vox);
    [g2, g1, g3] = gradient(Is); clear Is          % gradient(): dim2 first, then dim1, dim3
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
    % Jv: n x 6 [J11 J12 J13 J22 J23 J33]. Returns smallest-eigenvalue eigenvector
    % and tube coherence (mu2 - mu1)/(mu2 + mu1), mu ascending.
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
    th = mod(th, pi);                               % axial, [0, pi)
end

function th = display_to_conv(thDisp, zeroAxis, ccw)
    % thDisp: angle CCW from +dim2 as displayed by imagesc (dim1 points down)
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
    d = abs(mod(a - b + 90, 180) - 90);             % [0, 90]
end

function d = cliffs_delta(x, y)
    x = x(:); y = y(:); nx = numel(x); ny = numel(y);
    rk = tiedrank([x; y]);
    U = sum(rk(1:nx)) - nx * (nx + 1) / 2;
    d = (2 * U - nx * ny) / (nx * ny);              % >0: x tends larger
end

function ci = boot_by_unit(x, unit, nBoot)
    [~, ~, gi] = unique(unit);
    cells = accumarray(gi, (1:numel(x))', [], @(v) {v});
    nU = numel(cells); med = nan(nBoot, 1);
    for b = 1:nBoot
        pick = cells(randi(nU, nU, 1));
        med(b) = median(x(vertcat(pick{:})));
    end
    ci = prctile(med, [2.5 97.5]);
end

function ci = boot_hier(x, subj, unit, nBoot)
    us = unique(subj); nS = numel(us); med = nan(nBoot, 1);
    uCells = cell(nS, 1);
    for s = 1:nS
        ii = find(subj == us(s));
        [~, ~, gi] = unique(unit(ii));
        uCells{s} = accumarray(gi, ii, [], @(v) {v});
    end
    for b = 1:nBoot
        idx = [];
        for s = randi(nS, nS, 1)'
            c = uCells{s}; pick = c(randi(numel(c), numel(c), 1));
            idx = [idx; vertcat(pick{:})]; %#ok<AGROW>
        end
        med(b) = median(x(idx));
    end
    ci = prctile(med, [2.5 97.5]);
end

function T = binned_stats(x, y, unit, edges, nBoot)
    b = discretize(x, edges); nb = numel(edges) - 1;
    T = table((1:nb)', nan(nb,1), nan(nb,1), nan(nb,1), zeros(nb,1), ...
        'VariableNames', {'bin','median','lo','hi','n'});
    for k = 1:nb
        i = b == k;
        if nnz(i) < 20, continue; end
        ci = boot_by_unit(y(i), unit(i), nBoot);
        T.median(k) = median(y(i)); T.lo(k) = ci(1); T.hi(k) = ci(2); T.n(k) = nnz(i);
    end
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

function draw_panel(ax, kind, RR, subs, prm, radEdges, cohEdges)
    hold(ax, 'on'); box(ax, 'off');
    C_OBS = [0.90 0.50 0.10]; C_NULL = [0.60 0.60 0.60];
    MK = {'o','s','^','d'}; LS = {'-','--',':','-.'};
    nulAll = [];
    for s = 1:numel(subs), nulAll = [nulAll; RR.(subs{s}).null]; end %#ok<AGROW>
    switch kind
        case 'violin'
            obs = []; mo = nan(1, numel(subs)); mn = mo;
            for s = 1:numel(subs)
                D = RR.(subs{s}); obs = [obs; D.obs]; %#ok<AGROW>
                mo(s) = median(D.obs); mn(s) = median(D.null);
            end
            draw_violin(ax, 1, obs,    C_OBS,  mo, MK);
            draw_violin(ax, 2, nulAll, C_NULL, mn, MK);
            set(ax, 'XTick', [1 2], 'XTickLabel', {'Observed','Shuffled'}); xlim(ax, [0.4 2.6]);
            h = gobjects(numel(subs), 1);
            for s = 1:numel(subs)
                h(s) = plot(ax, nan, nan, MK{s}, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
            end
            legend(ax, h, subs, 'Location', 'northeast', 'Box', 'off');
        otherwise
            switch kind
                case 'shell',     fld = 'shellStats'; xl = 'Distance from vessel wall (\mum)';
                case 'radius',    fld = 'radStats';   xl = 'Vessel radius (\mum)';
                case 'coherence', fld = 'cohStats';   xl = 'Fiber coherence R';
            end
            h = gobjects(numel(subs), 1);
            for s = 1:numel(subs)
                T = RR.(subs{s}).(fld);
                if strcmp(kind, 'shell'), x = T.x; else, x = T.bin; end
                ok = isfinite(T.median);
                fill(ax, [x(ok); flipud(x(ok))], [T.lo(ok); flipud(T.hi(ok))], C_OBS, ...
                    'FaceAlpha', 0.15, 'EdgeColor', 'none');
                h(s) = plot(ax, x(ok), T.median(ok), 'Color', C_OBS, 'LineStyle', LS{s}, ...
                    'Marker', MK{s}, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
            end
            if strcmp(kind, 'radius'), set(ax, 'XTick', 1:numel(radEdges)-1, 'XTickLabel', edge_labels(radEdges)); end
            if strcmp(kind, 'coherence'), set(ax, 'XTick', 1:numel(cohEdges)-1, 'XTickLabel', edge_labels(min(cohEdges, 1))); end
            if strcmp(kind, 'shell'), xline(ax, prm.primaryMax, ':', 'Color', [0.4 0.4 0.4]); end
            yline(ax, median(nulAll), '--', 'Color', C_NULL, 'LineWidth', 1.2);
            legend(ax, h, subs, 'Location', 'northeast', 'Box', 'off');
            xlabel(ax, xl);
    end
    yline(ax, 45, 'k:');
    ylim(ax, [0 90]); yticks(ax, 0:15:90); ylabel(ax, '\Delta\theta vessel–fiber (deg)');
end

function L = edge_labels(e)
    L = cell(1, numel(e) - 1);
    for k = 1:numel(e) - 1
        if isinf(e(k+1)), L{k} = sprintf('>%g', e(k)); else, L{k} = sprintf('%g–%g', e(k), e(k+1)); end
    end
end

function draw_violin(ax, x, y, col, subjMed, mk)
    y = y(isfinite(y));
    yi = linspace(0, 90, 200);
    f = ksdensity(y, yi, 'Support', [-1e-3 90 + 1e-3], 'BoundaryCorrection', 'reflection');
    f = 0.35 * f / max(f);
    fill(ax, [x - f, fliplr(x + f)], [yi, fliplr(yi)], col, 'FaceAlpha', 0.6, 'EdgeColor', col * 0.7);
    q = prctile(y, [25 50 75]);
    plot(ax, [x x], q([1 3]), 'k-', 'LineWidth', 4);
    plot(ax, x, q(2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
    for s = 1:numel(subjMed)
        plot(ax, x + 0.15, subjMed(s), mk{s}, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
    end
end

function run_orient_check(V, prm, sub, reg, outDir)
    wm = logical(V.mask_wm);
    o  = double(V.orient);
    v  = o(wm & isfinite(o));
    p  = prctile(v, [1 50 99]);
    [th, info] = decode_orient(o, prm.units, wm);
    fprintf('\n%s %s: orient %s, size %s\n', sub, reg, class(V.orient), mat2str(size(V.orient)));
    fprintf('  WM voxels with finite orient: %d / %d\n', numel(v), nnz(wm));
    fprintf('  range [%.4g, %.4g]; p1/p50/p99 = %.4g / %.4g / %.4g; frac == 0: %.3f\n', ...
        info.raw_min, info.raw_max, p, mean(v == 0));
    fprintf('  -> interpreted as %s\n', info.units);

    % slab around the slice with most WM (thick enough for the largest integration scale)
    [~, z] = max(squeeze(sum(sum(wm, 1), 2)));
    zpad = ceil(3 * max(prm.sigmaI) / prm.vox);
    zs = max(1, z - zpad):min(size(wm, 3), z + zpad); zc = z - zs(1) + 1;
    segS = logical(V.seg(:,:,zs));
    wall = segS(:,:,zc) & ~imerode(segS(:,:,zc), true(3));
    [w1, w2] = find(wall); w1 = w1(1:4:end); w2 = w2(1:4:end);
    aIdx = sub2ind(size(segS), w1, w2, zc * ones(size(w1)));
    vt = vessel_tensor(segS, aIdx, prm);
    good = vt.coh >= prm.minCoh & abs(vt.t(:,3)) < sind(prm.maxIncl);
    fprintf('  slice %d: %d wall anchors sampled, %d pass coherence + inclination\n', ...
        z, numel(aIdx), nnz(good));

    f = figure('Color', 'w', 'Position', [50 50 1500 700]);
    tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
    ax = nexttile(tl); histogram(ax, v, 180); title(ax, sprintf('%s %s raw orient (WM)', sub, reg));
    xlabel(ax, sprintf('value (interpreted: %s)', info.units));

    ax = nexttile(tl, [1 2]);
    if isfield(V, 'mus'), bg = double(V.mus(:,:,z)); else, bg = double(wm(:,:,z)); end
    imagesc(ax, bg); axis(ax, 'image'); colormap(ax, gray);
    lim = prctile(bg(wm(:,:,z)), [2 98]); if numel(lim) == 2 && diff(lim) > 0, clim(ax, lim); end
    hold(ax, 'on');
    contour(ax, double(V.seg(:,:,z)), [0.5 0.5], 'r', 'LineWidth', 0.5);

    % fiber glyphs: block circular means of orient
    B = 10; [n1, n2] = size(bg); ths = th(:,:,z); wms = wm(:,:,z);
    for a = 1:B:n1 - B + 1
        for b = 1:B:n2 - B + 1
            blk = ths(a:a+B-1, b:b+B-1); m = wms(a:a+B-1, b:b+B-1) & isfinite(blk);
            if nnz(m) < 0.5 * B^2, continue; end
            tm = 0.5 * atan2(mean(sin(2*blk(m))), mean(cos(2*blk(m))));
            draw_glyph(ax, b + B/2, a + B/2, conv_to_display(tm, prm.zeroAxis, prm.ccw), 0.8*B, 'c');
        end
    end
    % vessel glyphs: structure-tensor axis at wall voxels (convention-independent)
    for k = find(good)'
        draw_glyph(ax, w2(k), w1(k), atan2(-vt.t(k,1), vt.t(k,2)), 8, [1 0.55 0.1]);
    end
    title(ax, sprintf('z = %d | cyan: optic axis (zero=%s, ccw=%d) | orange: vessel axis (tensor)', ...
        z, prm.zeroAxis, prm.ccw));
    exportgraphics(f, fullfile(outDir, sprintf('check_%s_%s.png', sub, reg)), 'Resolution', 200);
    savefig(f, fullfile(outDir, sprintf('check_%s_%s.fig', sub, reg)));
end

function draw_glyph(ax, x, y, thDisp, L, col)
    dx = 0.5 * L * cos(thDisp); dy = -0.5 * L * sin(thDisp);   % image rows increase downward
    plot(ax, [x - dx, x + dx], [y - dy, y + dy], '-', 'Color', col, 'LineWidth', 1.2);
end