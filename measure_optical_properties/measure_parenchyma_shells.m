%% Measure parenchymal optical properties around EPVS and non-pathological vessels
%{
Each parenchymal voxel is assigned to its single nearest structure (WM
vessel, EPVS, or other vessel) using a 3D Euclidean distance transform, then
binned into 40 µm shells. One output row = one structure x one shell, with
an object_id that stays the same across shells (supports (1 | object_id)).

Measures: scattering coefficient (mus) and retardance, in white-matter
parenchyma only (cfg.restrict_wm), as raw values and as percentage change
from the subject/region WM median (mus_pct, ret_pct).
Subjects: CAA6, CAA26 (control); CAA17 (moderate, occipital only);
          CAA22, CAA25 (severe).

Requires on path: load_local_params, exclude_epvs_ves, add_percent_change,
save_mri.
%}

clear; clc; close all;
% Paths resolved from this script's location, not the current folder
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end   % e.g. code run via Evaluate Selection
parentDir = fileparts(scriptDir);             % .../oct-caa
fsDir     = fullfile(parentDir, 'freesurfer');
assert(isfolder(fsDir), 'freesurfer folder not found: %s', fsDir);
addpath(genpath(fsDir));                      % freesurfer + all subfolders
addpath(parentDir);

%% Configuration
cfg.res_um    = 20;            % isotropic voxel size (µm)
cfg.edges_um  = 0:40:520;      % shells [0,40), [40,80), ..., [480,520]
cfg.n_min_obj = 50;            % min voxels for a vessel or EPVS object
cfg.n_min_bin = 25;            % min parenchymal voxels per object x shell (flag only)
cfg.epvs_buffer_vox = 1;       % exclude vessel components within N voxels of an EPVS
cfg.restrict_wm = true;        % measure white-matter parenchyma only
cfg.save_qc   = true;          % write shell / owner-class NIfTIs for visual QC
cfg.regions   = {'front', 'occip'};   % struct fields treated as regions
cfg.data_dir  = '/projectnb/npbssmic/ns/CAA/';

% Percentage change from the subject/region white-matter median
cfg.median_file = fullfile(cfg.data_dir, 'median_white_matter_values_14-Nov-2025.mat');
cfg.mus_upper   = 25;          % raw mus bound (flag only: mus_ok)
cfg.ret_upper   = 45;          % raw retardance bound (flag only: ret_ok)
assert(isfile(cfg.median_file), 'Median file not found: %s', cfg.median_file);

% file: path relative to data_dir (CAA17 is stored under occip/)
cfg.subjects = struct( ...
    'id',       {'caa6',            'caa26',             'caa17',                 'caa22',             'caa25'}, ...
    'severity', {'control',         'control',           'moderate',              'severe',            'severe'}, ...
    'file',     {'caa6/caa6.mat',   'caa26/caa26.mat',   'caa17/occip/caa17.mat', 'caa22/caa22.mat',   'caa25/caa25.mat'});

stamp = string(datetime('now', 'Format', 'yyyy-MM-dd'));
fout  = fullfile(cfg.data_dir, "parenchyma_shells_" + stamp);

%% Main loop
allT = {};
for s = 1:numel(cfg.subjects)
    sub = cfg.subjects(s).id;
    sev = cfg.subjects(s).severity;
    fprintf('Loading %s (%s)\n', sub, sev);
    S = load(fullfile(cfg.data_dir, cfg.subjects(s).file));
    S = S.(sub);

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
        [mus, ret, seg, seg_wm, ~, epvs, wm] = load_local_params(S, loc);

        [T, qc] = measure_shells(mus, ret, seg, seg_wm, epvs, wm, cfg);
        clear mus ret seg seg_wm epvs wm

        if cfg.save_qc && ~isempty(qc)
            qdir = fullfile(cfg.data_dir, sub, loc, 'shell_qc');
            if ~exist(qdir, 'dir'), mkdir(qdir); end
            vox_mm = cfg.res_um * [1 1 1] / 1000;
            save_mri(qc.shell,       fullfile(qdir, 'shell_bin.nii'),   vox_mm, 'uchar', 0);
            save_mri(qc.owner_class, fullfile(qdir, 'owner_class.nii'), vox_mm, 'uchar', 0);
        end
        if isempty(T)
            warning('%s %s: no rows produced', sub, loc);
            continue
        end

        n = height(T);
        T.subject   = repmat(string(sub), n, 1);
        T.severity  = repmat(string(sev), n, 1);
        T.region    = repmat(string(loc), n, 1);
        T.object_id = T.subject + "_" + T.region + "_" + T.class + "_" + string(T.obj);
        allT{end+1} = T; %#ok<SAGROW>
        fprintf('  Finished %s %s: %d object x shell rows\n', sub, loc, n);
    end
    clear S
end

if isempty(allT)
    error('No data produced for any subject/region.');
end
shells = vertcat(allT{:});

% Add percentage change from WM median (requires add_percent_change.m)
fprintf('Adding percentage change from %s\n', cfg.median_file);
M = load(cfg.median_file, 'parench_median');
shells = add_percent_change(shells, M.parench_median, cfg.mus_upper, cfg.ret_upper);
shells = movevars(shells, {'subject','severity','region','object_id'}, 'Before', 1);
writetable(shells, fout + ".csv");
save(fout + ".mat", 'shells', 'cfg', '-v7.3');
fprintf('Saved %s.csv / .mat\n', fout);

%% ------------------------------------------------------------------------
function [T, qc] = measure_shells(mus, ret, seg, seg_wm, epvs, wm, cfg)
% Assign each parenchymal voxel to its nearest structure and 40 µm shell,
% then average mus and retardance per structure x shell.
%
% Label map:  1..nV        reported WM vessels (not overlapping EPVS)
%             nV+1..nV+nE  EPVS
%             nV+nE+1      all other vessels: they claim nearby parenchyma
%                          but are not reported

T = table(); qc = [];
res = cfg.res_um;
if isempty(epvs), epvs = false(size(seg)); end
assert(isequal(size(mus), size(ret), size(seg), size(seg_wm), size(epvs), size(wm)), ...
       'Input volumes have mismatched sizes');

seg  = keep_large(logical(seg),  cfg.n_min_obj);
epvs = keep_large(logical(epvs), cfg.n_min_obj);

% Reported vessels: size filter applied AFTER WM masking and EPVS exclusion
ves = logical(seg_wm) & seg;
n_raw = nnz(ves);
fprintf('    vessel/EPVS shared voxels: %d\n', nnz(ves & epvs));
if any(epvs(:))
    [ves, n_exc] = exclude_epvs_ves(ves, epvs, cfg.epvs_buffer_vox);
    fprintf('    excluded %d vessel components (%.1f%% of WM vessel voxels)\n', ...
            n_exc, 100 * (1 - nnz(ves) / max(n_raw, 1)));
end
ves   = keep_large(logical(ves), cfg.n_min_obj);
other = seg & ~ves & ~epvs;

[Lv, nV] = bwlabeln(ves,  26);  Lv = uint32(Lv);
[Le, nE] = bwlabeln(epvs, 26);  Le = uint32(Le);
if nV + nE == 0, return; end
nL = nV + nE + 1;

L = zeros(size(seg), 'uint32');
L(other) = nL;
L(ves)   = Lv(ves);
L(epvs)  = Le(epvs) + nV;           % EPVS takes precedence on overlap
struct_mask = L > 0;

% Parenchyma only
mus(struct_mask) = NaN;
ret(struct_mask) = NaN;

% Restrict measured parenchyma to white matter. Distances are still
% computed through all tissue; only voxels outside WM are dropped.
% This also drops background voxels if mus/ret are 0 (not NaN) there.
if cfg.restrict_wm
    mus(~wm) = NaN;
    ret(~wm) = NaN;
end

% 3D distance (voxel centre to nearest structure voxel centre) + owner
[D, IDX] = bwdist(struct_mask);
D     = D * res;
owner = L(IDX);
clear IDX

nB  = numel(cfg.edges_um) - 1;
bin = uint8(fillmissing(discretize(D, cfg.edges_um), 'constant', 0));
ok  = D > 0 & bin > 0 & ~isnan(mus) & ~isnan(ret);
clear D

% Means per object x shell
sz  = [nL, nB];
idx = [double(owner(ok)), double(bin(ok))];
nvx = accumarray(idx, 1, sz);
mmu = accumarray(idx, double(mus(ok)), sz) ./ nvx;
mrt = accumarray(idx, double(ret(ok)), sz) ./ nvx;

% Object descriptors (for size-matching or a diameter covariate)
obj_vox = accumarray(double(L(struct_mask)), 1, [nL 1]);
din = zeros(nL, 1);
if nV
    Dv = bwdist(~ves);
    din(1:nV) = accumarray(double(Lv(ves)), double(Dv(ves)), [nV 1], @max);
    clear Dv
end
if nE
    De = bwdist(~epvs);
    din(nV + (1:nE)) = accumarray(double(Le(epvs)), double(De(epvs)), [nE 1], @max);
    clear De
end
obj_diam_um = 2 * din * res;   % max inscribed diameter (approximate)

% Long-format table: reported objects only
nR = nV + nE;
[obj, b] = ndgrid(1:nR, 1:nB);
obj = obj(:);  b = b(:);
lin  = sub2ind(sz, obj, b);
keep = nvx(lin) > 0;
obj = obj(keep);  b = b(keep);  lin = lin(keep);

cls = repmat("vessel", numel(obj), 1);
cls(obj > nV) = "EPVS";
e = cfg.edges_um(:);

T = table(obj, cls, e(b), e(b + 1), nvx(lin), mmu(lin), mrt(lin), ...
          obj_vox(obj), obj_diam_um(obj), nvx(lin) >= cfg.n_min_bin, ...
    'VariableNames', {'obj','class','d_inner_um','d_outer_um','n_vox', ...
                      'mus','ret','obj_vox','obj_diam_um','enough_voxels'});

% QC volumes: shell index, and which class owns each parenchymal voxel
oc = zeros(nL, 1, 'uint8');
oc(1:nV) = 1;  oc(nV + (1:nE)) = 2;  oc(nL) = 3;
qc.shell       = zeros(size(seg), 'uint8');
qc.shell(ok)   = bin(ok);
qc.owner_class = zeros(size(seg), 'uint8');
qc.owner_class(ok) = oc(owner(ok));
end

%% ------------------------------------------------------------------------
function seg_parsed = keep_large(seg, nmin)
% Keep 26-connected components with at least nmin voxels
CC = bwconncomp(seg, 26);
sizes = cellfun(@numel, CC.PixelIdxList);
seg_parsed = false(size(seg));
seg_parsed(vertcat(CC.PixelIdxList{sizes >= nmin})) = true;
end