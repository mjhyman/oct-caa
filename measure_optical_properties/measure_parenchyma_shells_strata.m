%% Measure parenchymal optical properties around EPVS and non-pathological vessels
%{
Each parenchymal voxel is assigned to its single nearest structure (WM
vessel, EPVS, or other vessel) using a 3D Euclidean distance transform, then
binned into 40 um shells.

Changes from the previous version
  1. The white-matter baseline (median mus and retardance) is computed here,
     from the same WM parenchymal voxels the shells are drawn from. The
     external median file and add_percent_change.m are no longer used.
  2. Every parenchymal voxel also carries
       - its depth from the gray-white interface (0.5 mm bands, >2 mm = deep)
       - the LOCAL caliber of the structure it is assigned to (skeleton-based,
         so a branching network is not treated as one size)
  3. Three outputs:
       shells    one row = object x shell               (same layout as before,
                 plus mean depth / mean local caliber as covariates)
       strat     one row = object x shell x depth band x caliber bin
       baseline  one row = subject x region: WM medians, shell coverage,
                 largest-component fractions, caliber bin edges
     object_id is unchanged, so the hierarchical bootstrap still resamples
     objects. Stratified rows carry n_vox so they can be recombined.

Subjects: CAA6, CAA26 (control); CAA17 (moderate, occipital only);
          CAA22, CAA25 (severe).

Requires on path: load_local_params, save_mri.
(Vessel exclusion near EPVS is a local function here, so the result does not
depend on which version of exclude_epvs_ves.m is on the path.)
Requires Image Processing Toolbox (bwdist, bwlabeln, bwskel, imreconstruct)
and Statistics and Machine Learning Toolbox (quantile).

Depth map: fold_depth_mm() reproduces derive_gm() + the interface/depth
construction of network_base() in epvs_network_fold_depth.m (INTERFACE_GM =
'cortex', MIN_COMP_VOX = 27, native 20 um grid), reading seg, epvs, mask and
mask_wm from the original struct, so the bands are the ones in Fig. 2C.
DEPTH_CERTAINTY is not applied (it is off in the network defaults).
%}

% Keep the loaded data: the shared struct 'caa' and any per-subject structs
% (caa6, caa17, ...) stay in the workspace; everything else is cleared.
clearvars -except -regexp ^caa\d*$
clc; close all;
% Paths resolved from this script's location, not the current folder
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end   % e.g. code run via Evaluate Selection
parentDir = fileparts(scriptDir);             % .../oct-caa
fsDir     = fullfile(parentDir, 'freesurfer');
assert(isfolder(fsDir), 'freesurfer folder not found: %s', fsDir);
addpath(genpath(fsDir));                      % freesurfer + all subfolders
addpath(parentDir);

%% Configuration
cfg.exploratory = true;        % tags the output file name
cfg.res_um    = 20;            % isotropic voxel size (um)
cfg.edges_um  = 0:40:520;      % shells [0,40), [40,80), ..., [480,520]
cfg.n_min_obj = 50;            % min voxels for a vessel or EPVS object
cfg.n_min_bin = 25;            % min parenchymal voxels per row (flag only)
cfg.epvs_buffer_vox = 1;       % exclude vessel components within N voxels of an EPVS
cfg.restrict_wm = true;        % measure white-matter parenchyma only
cfg.save_qc   = true;          % write QC NIfTIs
cfg.regions   = {'front', 'occip'};   % struct fields treated as regions
cfg.data_dir  = '/projectnb/npbssmic/ns/CAA/';                    % input structs
cfg.out_base  = '/projectnb/npbssmic/ns/CAA/parenchyma_shells/';  % each run gets its own subfolder

% Depth map inputs, read from the original struct (as in epvs_network_fold_depth.m)
cfg.tissue_field = 'mask';         % TISSUE_FIELD
cfg.wm_field     = 'mask_wm';
cfg.interface_gm = 'cortex';       % INTERFACE_GM: 'cortex' | 'all'
cfg.depth_min_comp_vox = 27;       % MIN_COMP_VOX used when building the tree for the interface

% Raw-value bounds (flag only: mus_ok, ret_ok). Not applied to the baseline.
cfg.mus_upper = 25;
cfg.ret_upper = 45;

% Depth from the gray-white interface (mm): DEPTH_EDGES_MM of the network
% script. WM deeper than the last edge = 'deep'.
cfg.depth_edges_mm = [0 0.5 1 1.5 2];

% Local-caliber bins. Default: quantile bins set within each tissue volume,
% separately for vessels and EPVS. To use fixed edges shared across volumes,
% give the interior edges in um, e.g. cfg.cal_edges_um.EPVS = [150 300];
cfg.n_cal_bins   = 3;
cfg.cal_edges_um = struct('vessel', [], 'EPVS', []);

% file: path relative to data_dir (CAA17 is stored under occip/)
cfg.subjects = struct( ...
    'id',       {'caa6',            'caa26',             'caa17',                 'caa22',             'caa25'}, ...
    'severity', {'control',         'control',           'moderate',              'severe',            'severe'}, ...
    'file',     {'caa6/caa6.mat',   'caa26/caa26.mat',   'caa17/occip/caa17.mat', 'caa22/caa22.mat',   'caa25/caa25.mat'});

% Output: <out_base>/<yyyyMMdd_HHmmss>_<EXPLORATORY|MANUSCRIPT>/
stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
tag   = "MANUSCRIPT";
if cfg.exploratory, tag = "EXPLORATORY"; end
out_dir = char(fullfile(cfg.out_base, stamp + "_" + tag));   % char: save_mri expects char paths
assert(~isfolder(out_dir), 'Run folder already exists: %s', out_dir);
mkdir(out_dir);
cfg.out_dir = out_dir;
fout = fullfile(out_dir, "parenchyma_shells");
fprintf('Writing outputs to %s\n', out_dir);

%% Load all subjects up front (skipped for subjects already in the workspace)
% Data live in caa.(subject id), e.g. caa.caa6.front. A subject is taken, in
% order, from: caa.(id) if present; a workspace variable named <id> (e.g.
% caa6); otherwise its .mat file. Nothing here is cleared afterwards.
if ~exist('caa', 'var') || ~isstruct(caa), caa = struct(); end
for s = 1:numel(cfg.subjects)
    sub = cfg.subjects(s).id;
    if isfield(caa, sub)
        fprintf('%s: found in caa struct; skipping load\n', upper(sub));
    elseif exist(sub, 'var') && isstruct(eval(sub))
        caa.(sub) = eval(sub);              % shared copy, no extra memory
        fprintf('%s: found in workspace; skipping load\n', upper(sub));
    else
        fpath = fullfile(cfg.data_dir, cfg.subjects(s).file);
        fprintf('%s: loading %s\n', upper(sub), fpath);
        tmp = load(fpath);
        assert(isfield(tmp, sub), '%s: variable "%s" not found in %s', upper(sub), sub, fpath);
        caa.(sub) = tmp.(sub);
        clear tmp
    end
end

%% Main loop
allT = {};  allS = {};  allB = {};
for s = 1:numel(cfg.subjects)
    sub = cfg.subjects(s).id;
    sev = cfg.subjects(s).severity;
    fprintf('Processing %s (%s)\n', sub, sev);
    S = caa.(sub);                          % reference to the loaded struct (read only)

    % Only treat known region fields as regions
    regions = intersect(fieldnames(S), cfg.regions, 'stable');
    if isempty(regions)
        warning('%s: no region fields found. Fields present: %s', ...
                sub, strjoin(fieldnames(S), ', '));
        continue
    end

    for r = 1:numel(regions)
        loc = regions{r};
        fprintf('  Processing %s %s\n', sub, loc);
        [mus, ret, seg, seg_wm, ~, epvs] = load_local_params(S, loc);
        wm = S.(loc).mask_wm;

        % Depth from the gray-white interface, built from the original struct
        % with the same rules as epvs_network_fold_depth.m
        depth = fold_depth_mm(S.(loc), cfg);
        assert(isequal(size(depth), size(mus)), ...
               ['%s %s: depth map size differs from the volumes returned by ' ...
                'load_local_params (was the volume cropped or reoriented there?)'], sub, loc);

        [T, Ts, B, qc] = measure_shells(mus, ret, seg, seg_wm, epvs, wm, depth, cfg);
        clear mus ret seg seg_wm epvs wm depth

        if cfg.save_qc && ~isempty(qc)
            qdir = fullfile(out_dir, 'shell_qc', [sub '_' loc]);
            if ~exist(qdir, 'dir'), mkdir(qdir); end
            vox_mm = cfg.res_um * [1 1 1] / 1000;
            save_mri(qc.shell,       fullfile(qdir, 'shell_bin.nii'),   vox_mm, 'uchar', 0);
            save_mri(qc.owner_class, fullfile(qdir, 'owner_class.nii'), vox_mm, 'uchar', 0);
            save_mri(qc.depth_band,  fullfile(qdir, 'depth_band.nii'),  vox_mm, 'uchar', 0);
            save_mri(qc.cal_bin,     fullfile(qdir, 'caliber_bin.nii'), vox_mm, 'uchar', 0);
        end
        clear qc
        if isempty(T)
            warning('%s %s: no rows produced', sub, loc);
            continue
        end

        T  = add_ids(T,  sub, sev, loc);
        Ts = add_ids(Ts, sub, sev, loc);
        B.subject  = string(sub);
        B.severity = string(sev);
        B.region   = string(loc);
        allT{end+1} = T;                          %#ok<SAGROW>
        allS{end+1} = Ts;                         %#ok<SAGROW>
        allB{end+1} = struct2table(B);            %#ok<SAGROW>
        fprintf('  Finished %s %s: %d object x shell rows, %d stratified rows\n', ...
                sub, loc, height(T), height(Ts));
    end
end
clear S    % alias only; caa stays in the workspace

if isempty(allT)
    error('No data produced for any subject/region.');
end
shells   = vertcat(allT{:});
strat    = vertcat(allS{:});
baseline = vertcat(allB{:});
baseline = movevars(baseline, {'subject','severity','region'}, 'Before', 1);

writetable(shells,   fout + "_shells.csv");
writetable(strat,    fout + "_strat.csv");
writetable(baseline, fout + "_baseline.csv");
save(fout + ".mat", 'shells', 'strat', 'baseline', 'cfg', '-v7.3');
fprintf('Saved %s_{shells,strat,baseline}.csv and .mat\n', fout);

% Record the configuration and the script that produced this run
fid = fopen(fullfile(out_dir, 'config.json'), 'w');
fprintf(fid, '%s\n', jsonencode(cfg, 'PrettyPrint', true));
fclose(fid);
src = [mfilename('fullpath') '.m'];
if exist(src, 'file')
    copyfile(src, fullfile(out_dir, 'script_copy.m'));
else
    warning('Script path unavailable (run as selection?); script copy not saved.');
end
disp(baseline(:, {'subject','region','mus_median','ret_median','mus_mean', ...
                  'shell_coverage','epvs_largest_frac','ves_largest_frac'}));

%% ------------------------------------------------------------------------
function T = add_ids(T, sub, sev, loc)
n = height(T);
T.subject   = repmat(string(sub), n, 1);
T.severity  = repmat(string(sev), n, 1);
T.region    = repmat(string(loc), n, 1);
T.object_id = T.subject + "_" + T.region + "_" + T.class + "_" + string(T.obj);
T = movevars(T, {'subject','severity','region','object_id'}, 'Before', 1);
end

%% ------------------------------------------------------------------------
function [T, Ts, B, qc] = measure_shells(mus, ret, seg, seg_wm, epvs, wm, depth, cfg)
% Assign each parenchymal voxel to its nearest structure and 40 um shell,
% then average mus and retardance per structure x shell, and per
% structure x shell x depth band x local-caliber bin.
%
% Label map:  1..nV        reported WM vessels (not overlapping EPVS)
%             nV+1..nV+nE  EPVS
%             nV+nE+1      all other vessels: they claim nearby parenchyma
%                          but are not reported

T = table(); Ts = table(); B = struct(); qc = [];
res = cfg.res_um;
if isempty(epvs), epvs = false(size(seg)); end
assert(isequal(size(mus), size(ret), size(seg), size(seg_wm), size(epvs), size(wm), size(depth)), ...
       'Input volumes have mismatched sizes');
wm = logical(wm);

seg  = keep_large(logical(seg),  cfg.n_min_obj);
epvs = keep_large(logical(epvs), cfg.n_min_obj);

% Reported vessels: size filter applied AFTER WM masking and EPVS exclusion
ves = logical(seg_wm) & seg;
n_raw = nnz(ves);
fprintf('    vessel/EPVS shared voxels: %d\n', nnz(ves & epvs));
n_exc = 0;
if any(epvs(:))
    [ves, n_exc] = drop_vessels_near_epvs(ves, epvs, cfg.epvs_buffer_vox);
    fprintf('    excluded %d vessel components (%.1f%% of WM vessel voxels)\n', ...
            n_exc, 100 * (1 - nnz(ves) / max(n_raw, 1)));
end
B.ves_excluded_n_comp   = n_exc;                          % whole components dropped
B.ves_excluded_vox_frac = 1 - nnz(ves) / max(n_raw, 1);   % share of WM vessel voxels dropped
ves   = keep_large(logical(ves), cfg.n_min_obj);
ves   = ves & ~epvs;                 % EPVS takes precedence on any residual overlap
other = seg & ~ves & ~epvs;

[Lv, nV] = bwlabeln(ves,  26);  Lv = uint32(Lv);
[Le, nE] = bwlabeln(epvs, 26);  Le = uint32(Le);
if nV + nE == 0, return; end
nR = nV + nE;                       % reported objects
nL = nR + 1;

L = zeros(size(seg), 'uint32');
L(other) = nL;
L(ves)   = Lv(ves);
L(epvs)  = Le(epvs) + nV;
struct_mask = L > 0;

% ---- Local caliber of each structure voxel (um) --------------------------
cal = zeros(size(seg), 'single');
if nV, cal(ves)  = local_caliber_um(ves,  res); end
if nE, cal(epvs) = local_caliber_um(epvs, res); end

% ---- Object descriptors --------------------------------------------------
obj_vox = accumarray(double(L(struct_mask)), 1, [nL 1]);
din = zeros(nL, 1);
obj_cal_med = nan(nL, 1);
obj_depth   = nan(nL, 1);
if nV
    Dv = bwdist(~ves);
    din(1:nV) = accumarray(double(Lv(ves)), double(Dv(ves)), [nV 1], @max);
    clear Dv
    obj_cal_med(1:nV) = accumarray(double(Lv(ves)), double(cal(ves)), [nV 1], @median);
    obj_depth(1:nV)   = accumarray(double(Lv(ves)), double(depth(ves)), [nV 1], ...
                                   @(x) mean(x, 'omitnan'));
end
if nE
    De = bwdist(~epvs);
    din(nV + (1:nE)) = accumarray(double(Le(epvs)), double(De(epvs)), [nE 1], @max);
    clear De
    obj_cal_med(nV + (1:nE)) = accumarray(double(Le(epvs)), double(cal(epvs)), [nE 1], @median);
    obj_depth(nV + (1:nE))   = accumarray(double(Le(epvs)), double(depth(epvs)), [nE 1], ...
                                          @(x) mean(x, 'omitnan'));
end
obj_diam_um = 2 * din * res;        % max inscribed diameter of the whole component

% Caliber values used to set the bin edges (structure voxels, per class)
cal_ves_vals  = cal(ves);
cal_epvs_vals = cal(epvs);
clear Lv Le ves epvs other seg seg_wm

% ---- Parenchyma only -----------------------------------------------------
mus(struct_mask) = NaN;
ret(struct_mask) = NaN;

% Restrict measured parenchyma to white matter. Distances are still
% computed through all tissue; only voxels outside WM are dropped.
if cfg.restrict_wm
    mus(~wm) = NaN;
    ret(~wm) = NaN;
end
clear wm

% 3D distance (voxel centre to nearest structure voxel centre) + owner
[D, IDX] = bwdist(struct_mask);
D         = D * res;
owner     = L(IDX);
cal_owner = cal(IDX);               % local caliber of the nearest structure voxel
clear IDX cal L struct_mask

nB   = numel(cfg.edges_um) - 1;
bin  = uint8(fillmissing(discretize(D, cfg.edges_um), 'constant', 0));
base = D > 0 & ~isnan(mus) & ~isnan(ret);    % all WM parenchymal voxels
ok   = base & bin > 0;                        % those inside a shell
clear D

% ---- Baseline: WM medians from the voxels measured here ------------------
B.mus_median = double(median(mus(base)));
B.ret_median = double(median(ret(base)));
B.mus_mean   = mean(double(mus(base)));
B.ret_mean   = mean(double(ret(base)));
B.mus_median_in_shells = double(median(mus(ok)));
B.ret_median_in_shells = double(median(ret(ok)));
B.n_wm_parenchyma_vox  = nnz(base);
B.shell_coverage       = nnz(ok) / max(nnz(base), 1);   % fraction of WM parenchyma within the outermost shell
is_v = ok & owner >= 1  & owner <= nV;
is_e = ok & owner >  nV & owner <= nR;
B.frac_owned_vessel = nnz(is_v) / max(nnz(base), 1);
B.frac_owned_epvs   = nnz(is_e) / max(nnz(base), 1);
B.frac_owned_other  = (nnz(ok) - nnz(is_v) - nnz(is_e)) / max(nnz(base), 1);
B.n_vessel = nV;
B.n_epvs   = nE;
B.ves_largest_frac  = max([obj_vox(1:nV); 0]) / max(sum(obj_vox(1:nV)), 1);
B.epvs_largest_frac = max([obj_vox(nV + (1:nE)); 0]) / max(sum(obj_vox(nV + (1:nE))), 1);
fprintf('    WM median: mus = %.3f, ret = %.3f | shell coverage = %.1f%%\n', ...
        B.mus_median, B.ret_median, 100 * B.shell_coverage);
fprintf('    largest component holds %.1f%% of EPVS voxels, %.1f%% of vessel voxels\n', ...
        100 * B.epvs_largest_frac, 100 * B.ves_largest_frac);

% ---- Depth band and caliber bin of every parenchymal voxel ---------------
% Bands as in the network script: discretize on the finite edges, and
% depth > last edge = deep (index nD). 0 = undefined (no interface found).
dedges = cfg.depth_edges_mm(:)';
nD = numel(dedges);                      % (numel - 1) bands + deep
db = discretize(depth, dedges);
db(depth > dedges(end)) = nD;
db = uint8(fillmissing(db, 'constant', 0));

cb = zeros(size(bin), 'uint8');
[cb(is_v), ev] = caliber_bins(cal_ves_vals,  cal_owner(is_v), cfg.n_cal_bins, cfg.cal_edges_um.vessel);
[cb(is_e), ee] = caliber_bins(cal_epvs_vals, cal_owner(is_e), cfg.n_cal_bins, cfg.cal_edges_um.EPVS);
clear cal_ves_vals cal_epvs_vals
nC = max(1, cfg.n_cal_bins);
if ~isempty(cfg.cal_edges_um.vessel), nC = max(nC, numel(cfg.cal_edges_um.vessel) + 1); end
if ~isempty(cfg.cal_edges_um.EPVS),   nC = max(nC, numel(cfg.cal_edges_um.EPVS) + 1);   end
B.cal_edges_vessel_um = string(mat2str(double(ev), 5));
B.cal_edges_epvs_um   = string(mat2str(double(ee), 5));

% QC volumes
oc = zeros(nL, 1, 'uint8');
oc(1:nV) = 1;  oc(nV + (1:nE)) = 2;  oc(nL) = 3;
qc.shell       = zeros(size(bin), 'uint8');   qc.shell(ok)       = bin(ok);
qc.owner_class = zeros(size(bin), 'uint8');   qc.owner_class(ok) = oc(owner(ok));
qc.depth_band  = zeros(size(bin), 'uint8');   qc.depth_band(ok)  = db(ok);
qc.cal_bin     = zeros(size(bin), 'uint8');   qc.cal_bin(ok)     = cb(ok);

% ---- Pull the reported voxels into vectors -------------------------------
rep = is_v | is_e;
clear is_v is_e base ok
o  = double(owner(rep));
b  = double(bin(rep));
d  = double(db(rep));           % 0..nD
c  = double(cb(rep));           % 0..nC (0 only if a caliber could not be binned)
m  = double(mus(rep));
rt = double(ret(rep));
dp = double(depth(rep));
ca = double(cal_owner(rep));
clear owner bin db cb mus ret depth cal_owner rep

depth_lab = strings(1, nD);                % same labels as the network strata CSV
for i = 1:nD - 1
    depth_lab(i) = string(sprintf('%g-%g mm', dedges(i), dedges(i + 1)));
end
depth_lab(nD) = "deep";
depth_lab = ["undefined", depth_lab];      % index d + 1
e = cfg.edges_um(:);

% ---- Table 1: object x shell (same layout as before) ---------------------
sz2  = [nR, nB];
lin2 = sub2ind(sz2, o, b);
A    = accum_stats(lin2, prod(sz2), m, rt, dp, ca);
k    = find(A.n > 0);
[obj, bb] = ind2sub(sz2, k);
cls = repmat("vessel", numel(obj), 1);
cls(obj > nV) = "EPVS";
T = table(obj, cls, e(bb), e(bb + 1), A.n(k), A.mus(k), A.ret(k), ...
          A.depth(k), A.cal(k), ...
          obj_vox(obj), obj_diam_um(obj), obj_cal_med(obj), obj_depth(obj), ...
          A.n(k) >= cfg.n_min_bin, ...
    'VariableNames', {'obj','class','d_inner_um','d_outer_um','n_vox','mus','ret', ...
                      'depth_mm','local_cal_um', ...
                      'obj_vox','obj_diam_um','obj_cal_median_um','obj_depth_mm', ...
                      'enough_voxels'});
T = add_pct(T, B, cfg);

% ---- Table 2: object x shell x depth band x caliber bin ------------------
sz4  = [nR, nB, nD + 1, nC + 1];
lin4 = sub2ind(sz4, o, b, d + 1, c + 1);
A    = accum_stats(lin4, prod(sz4), m, rt, dp, ca);
k    = find(A.n > 0);
[obj, bb, dd, cc] = ind2sub(sz4, k);
cls = repmat("vessel", numel(obj), 1);
cls(obj > nV) = "EPVS";
depth_band = reshape(depth_lab(dd), [], 1);
Ts = table(obj, cls, e(bb), e(bb + 1), depth_band, dd - 1, cc - 1, ...
           A.n(k), A.mus(k), A.ret(k), A.depth(k), A.cal(k), ...
           obj_vox(obj), obj_diam_um(obj), obj_cal_med(obj), obj_depth(obj), ...
           A.n(k) >= cfg.n_min_bin, ...
    'VariableNames', {'obj','class','d_inner_um','d_outer_um','depth_band','depth_idx','cal_bin', ...
                      'n_vox','mus','ret','depth_mm','local_cal_um', ...
                      'obj_vox','obj_diam_um','obj_cal_median_um','obj_depth_mm', ...
                      'enough_voxels'});
Ts = add_pct(Ts, B, cfg);
end

%% ------------------------------------------------------------------------
function T = add_pct(T, B, cfg)
% Percentage change from the WM median computed in this script, plus flags
T.mus_pct = 100 * (T.mus - B.mus_median) ./ B.mus_median;
T.ret_pct = 100 * (T.ret - B.ret_median) ./ B.ret_median;
T.mus_ok  = T.mus < cfg.mus_upper;
T.ret_ok  = T.ret < cfg.ret_upper;
end

%% ------------------------------------------------------------------------
function A = accum_stats(lin, n, m, rt, dp, ca)
% Voxel count and means of mus, retardance, depth and local caliber per cell
A.n   = accumarray(lin, 1,  [n 1]);
A.mus = accumarray(lin, m,  [n 1]) ./ A.n;
A.ret = accumarray(lin, rt, [n 1]) ./ A.n;
A.cal = accumarray(lin, ca, [n 1]) ./ A.n;
v = ~isnan(dp);                                  % depth is undefined where no interface was found
nd = accumarray(lin(v), 1, [n 1]);
A.depth = accumarray(lin(v), dp(v), [n 1]) ./ nd;   % NaN where nd == 0
end

%% ------------------------------------------------------------------------
function [b, edges] = caliber_bins(struct_vals, query_vals, n_bins, fixed_edges)
% Bin the local caliber of the owning structure voxel.
%   struct_vals : local caliber of every structure voxel of this class (sets the edges)
%   query_vals  : local caliber attached to each parenchymal voxel of this class
% Returns the bin index (uint8, 0 = could not be binned) and the edges used.
b = zeros(size(query_vals), 'uint8');
edges = [-Inf Inf];
if isempty(struct_vals) || isempty(query_vals)
    return
end
if ~isempty(fixed_edges)
    edges = [-Inf, sort(fixed_edges(:)'), Inf];
else
    q = quantile(double(struct_vals), (1:n_bins - 1) / n_bins);
    % Caliber takes few distinct values in small vessels, so quantile edges
    % can coincide; duplicates are dropped and fewer bins result.
    edges = unique([-Inf, q(:)', Inf]);
end
b = uint8(fillmissing(discretize(double(query_vals), edges), 'constant', 0));
end

%% ------------------------------------------------------------------------
function c = local_caliber_um(mask, res)
% Local diameter (um) at every voxel of MASK, returned as a vector in the
% order of find(mask).
% Each voxel takes the inscribed-sphere diameter at its nearest skeleton
% voxel, so a long branching component gets a caliber that varies along its
% length instead of one value for the whole component.
% Approximation: the nearest skeleton voxel is taken over the whole mask, so
% where two structures nearly touch a voxel can pick up its neighbour's
% skeleton.
Din = bwdist(~mask);                 % voxels; 1 at the surface layer
sk  = bwskel(mask);
if any(sk(:))
    [~, idx] = bwdist(sk);
    r = Din(idx(mask));
else
    r = Din(mask);
end
c = (2 * single(r) - 1) * res;       % a w-voxel-wide structure gives w * res
end

%% ------------------------------------------------------------------------
function D = fold_depth_mm(V, cfg)
% Depth (mm) of every voxel from the gray-white interface, constructed as in
% epvs_network_fold_depth.m (derive_gm + network_base), on the native grid:
%   GM            = mask & ~mask_wm
%   interface GM  = GM components connected to background outside the block
%                   (internal holes filled first) or to an array face;
%                   GM enclosed by WM stays GM but is not cortex
%   tree          = (seg & ~epvs) U epvs, components < MIN_COMP_VOX removed
%   WM compartment= mask_wm plus tree voxels outside both masks that are
%                   nearer WM than GM
%   interface     = WM-compartment voxels 26-adjacent to interface GM that is
%                   not tree
%   depth         = Euclidean distance to the interface
% NaN everywhere if no interface is found.
nb26 = true(3, 3, 3);
need = {cfg.tissue_field, cfg.wm_field, 'seg', 'epvs'};
miss = need(~isfield(V, need));
assert(isempty(miss), 'fold_depth_mm: missing field(s): %s', strjoin(miss, ', '));
tis = logical(V.(cfg.tissue_field));
wm  = logical(V.(cfg.wm_field));
seg = logical(V.seg);
epv = logical(V.epvs);
if isempty(epv), epv = false(size(seg)); end

% derive_gm
gm = tis & ~wm;
if strcmp(cfg.interface_gm, 'all')
    gmi = gm;
else
    filled = imfill(tis | wm, 'holes');
    seed = gm & imdilate(~filled, nb26);           % touches background outside the block
    seed([1 end], :, :) = seed([1 end], :, :) | gm([1 end], :, :);   % sectioned cortex
    seed(:, [1 end], :) = seed(:, [1 end], :) | gm(:, [1 end], :);
    seed(:, :, [1 end]) = seed(:, :, [1 end]) | gm(:, :, [1 end]);
    gmi = imreconstruct(seed, gm, 26);
    clear seed filled
end
clear tis

% network_base: tree and WM compartment
ves  = seg & ~epv;                                  % EPVS precedence
epv  = bwareaopen(epv, cfg.depth_min_comp_vox, 26);
ves  = bwareaopen(ves, cfg.depth_min_comp_vox, 26);
tree = ves | epv;
clear ves epv seg
lab1 = wm;
un = tree & ~wm & ~gm;                              % tree voxels outside both masks
if any(un(:))
    dW = bwdist(wm); dG = bwdist(gm);
    lab1(un & dW <= dG) = true;
    clear dW dG
end
clear un wm gm

iface = lab1 & imdilate(gmi & ~tree, nb26);
if ~any(iface(:))
    warning('fold_depth_mm: no gray-white interface found; depth left undefined');
    D = nan(size(lab1), 'single');
    return
end
D = bwdist(iface) * (cfg.res_um / 1000);
end

%% ------------------------------------------------------------------------
function [ves_out, n_exc] = drop_vessels_near_epvs(ves, epvs, buffer_vox)
% Remove every 26-connected vessel component that shares a voxel with an EPVS
% or lies within buffer_vox voxels (Euclidean) of one. buffer_vox = 0 gives
% the overlap-only rule of the original exclude_epvs_ves.m.
% NOTE: the whole component is removed. Where the vessel segmentation is one
% connected tree, a single EPVS contact removes the entire tree.
near = epvs;
if buffer_vox > 0, near = bwdist(epvs) <= buffer_vox; end
[L, n] = bwlabeln(ves, 26);
hit = unique(L(near & ves));
hit(hit == 0) = [];
keep = true(n + 1, 1);
keep(hit + 1) = false;
ves_out = ves & keep(L + 1);
n_exc = numel(hit);
end

%% ------------------------------------------------------------------------
function seg_parsed = keep_large(seg, nmin)
% Keep 26-connected components with at least nmin voxels
CC = bwconncomp(seg, 26);
sizes = cellfun(@numel, CC.PixelIdxList);
seg_parsed = false(size(seg));
seg_parsed(vertcat(CC.PixelIdxList{sizes >= nmin})) = true;
end