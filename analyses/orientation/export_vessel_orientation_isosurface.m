%% vessel_orientation_isosurface.m
% 3-D isosurfaces of white-matter vessels and EPVS in the control subjects (CAA26, CAA6),
% coloured by the direction of their structure-tensor axis.
%
% Architecture (one tensor, two volumes):
%   1. Combined mask = (seg | epvs) & white matter. The manual EPVS annotation contains the vessel
%      lumen and the enlarged perivascular space, and the vessel segmentation partly includes EPVS,
%      so the union is the most reliable geometry: its surface is the interface between the
%      vessel/fluid compartment and the white matter that the fibers touch.
%   2. The structure tensor, skeleton and watershed are computed once, on the combined mask. The
%      tensor uses mask gradients, which exist only at the union's outer surface, so EPVS voxels
%      get the EPVS envelope's axis and vessel voxels away from EPVS the vessel wall's axis, with no
%      artificial ends where vessels enter or leave an EPVS.
%   3. Class labels: the manual annotation is authoritative. EPVS = epvs & WM; vessel = combined &
%      ~epvs (seg voxels inside an EPVS annotation are EPVS). Every output is split by class.
%
%   Visualisation 1: tensor at every voxel (integration scale chosen per voxel by max tube coherence,
%                    as in compare_orient_vessel_*.m).
%   Visualisation 2: marker-controlled watershed branches (markers = skeleton pieces between branch
%                    points), each split at the EPVS boundary into a vessel piece and an EPVS piece.
%                    One tensor per piece (BRANCH_MODE = 'branch'): gradients of the combined mask
%                    summed over the piece's own wall voxels (within WALL_VOX of the surface), or a
%                    voxelwise tensor confined to each piece ('local').
%
% Junctions: within about one integration scale of a vessel/EPVS junction the axis blends both
% shapes. dJ = distance from each voxel to the nearest voxel of the other class (um) is written
% out (NIfTI, branch table) so results can be checked with junction regions excluded.
%
% Direction colour: RGB = |t . dim1|, |t . dim2|, |t . dim3| for the unit axis t (dim3 =
% through-plane), as in DTI direction maps; antiparallel axes share a colour. Axes with tube
% coherence < MIN_TUBE_COH are LOWCOH_COLOR.
%
% Outputs per subject x region, in OUT_DIR/<sub>_<reg>/ (summary, README and an unrotated colour key
% in OUT_DIR); <st> = vessel | epvs. PNG renders and <sub>_<reg>_direction_color_key.png are rotated
% so the gyrus is on top (RENDER_ROTATE_UPRIGHT); PLY and NIfTI keep the native voxel frame:
%   <sub>_<reg>_<st>_dir_voxel.{ply,tif,nii.gz}, _view*.png    visualisation 1
%   <sub>_<reg>_<st>_dir_branch.{ply,tif,nii.gz}, _view*.png   visualisation 2
%       PLY: binary, per-vertex RGB, um, x = dim1, y = dim2, z = dim3, voxel (1,1,1) at 0
%       TIF: RGB stack, full volume, slices along dim3, rot90 + scale bar as export_orient_rgb_tif.m
%       NIfTI: uint8 RGB in 3 frames, full volume, native grid (no rot90), voxel size in mm;
%              freeview -v <file>.nii.gz:rgb=1
%   <sub>_<reg>_pieces.csv            one row per piece (branch x class), with class and parent branch
%   <sub>_<reg>_piece_labels.nii.gz   int32 piece ids (0 = background)
%   <sub>_<reg>_class.nii.gz          uint8: 1 = vessel, 2 = EPVS
%   <sub>_<reg>_junction_dist.nii.gz  float32 dJ (um); -1 where the other class is absent
%   <sub>_<reg>_axis_voxel.nii.gz     float32, 3 frames: signed unit axis (dim1, dim2, dim3) of the
%                                     voxelwise tensor; 0 outside vessels/EPVS (Freeview :vector=1)
%   <sub>_<reg>_tensor_voxel.mat      struct TV: TV.vessel and TV.epvs (one row per voxel: full-volume
%                                     index and subscripts, signed axis, en-face angle, inclination,
%                                     tube coherence, chosen sigma_i, junction distance, piece id),
%                                     TV.meta (settings)
%   <sub>_<reg>_pieces.mat            piece label volume (cropped), class, parent, offset
%   direction_color_key.png, summary.csv, README.txt

clearvars -except caa*
close all

%% 0. Config ---------------------------------------------------------------
SUBJECTS = {'CAA26','CAA6'};
REGIONS  = {'front','occip'};
VOX_UM   = 20;

SCOPE       = 'wm';           % combined mask restricted to: 'wm' (mask_wm & mask) | 'tissue' | 'all'
MIN_OBJ_VOX = 27;             % combined-mask fragments smaller than this are removed (26-connected)
JUNCTION_UM = 160;            % voxels within this distance of the other class count as "near junction"

% orientation convention of D.orient (only for the en-face angle column of the piece table)
ORIENT_ZERO_AXIS = 'dim2';
ORIENT_CCW       = true;

% structure tensor (as in the analysis scripts)
SIGMA_D_UM   = 20;
SIGMA_I_UM   = [40 80 160];   % add a larger scale if wide EPVS interiors come out gray
MIN_TUBE_COH = 0.5;

% direction colour
LOWCOH_COLOR = [0.55 0.55 0.55];
RGB_GAMMA    = 1;             % < 1 brightens mixed directions

% branch separation (marker-controlled watershed on the combined mask)
SKEL_MIN_BRANCH_VOX = 3;      % bwskel MinBranchLength: spurs shorter than this are pruned
BRANCHPT_RADIUS_VOX = 1;      % cube of this radius removed around each skeleton branch point
MIN_MARKER_VOX      = 3;      % shortest skeleton piece kept as a branch marker
WALL_VOX            = 2;      % voxels within this distance of the surface carry the gradient
BRANCH_MODE         = 'local';    % 'local' (voxelwise, own piece's wall only) | 'branch' (one axis per piece)

% speed
TENSOR_DOWNSAMPLE = true;     % integration scales >= 4 voxels are filtered on a block-averaged grid
                              % (2x at 80 um, 4x at 160 um); axes differ by < 0.1 deg from full resolution

% meshes and renders
ISO_SMOOTH_VOX   = 1;         % Gaussian sigma (voxels) on the mask before the isosurface
ISO_LEVEL        = 0.5;
MESH_DOWNSAMPLE  = 1;         % PLY mesh resolution: 1 = full, 2 = half, ...
RENDER_MAX_FACES = 1.5e6;     % PNG renders use a coarser mesh above this face count
RENDER_VIEWS     = {[0 -90], [-35 30], [0 0]};   % en-face, oblique, side
RENDER_ROTATE_UPRIGHT = true; % renders rotated in-plane so the gyrus is on top (same rot90 as the TIFs)
WRITE_PLY    = true;
WRITE_PNG    = true;
WRITE_LABELS = true;          % pieces.mat
WRITE_NII    = true;          % compressed NIfTI for Freeview
WRITE_TENSOR = true;          % voxelwise axes: axis_voxel.nii.gz + tensor_voxel.mat
WRITE_VOXEL_VIS = true;       % visualisation 1 meshes / stacks (its axes are always computed for tensor_voxel.mat)

% RGB TIF stacks (full volume, slices along dim3)
WRITE_TIF      = true;
TIF_BACKGROUND = 'tissue';                % 'tissue' (GM / WM in dark grays) | 'black'
TIF_GM_COLOR   = uint8([70 70 70]);
TIF_WM_COLOR   = uint8([35 35 35]);
ROTATE_UPRIGHT = true;                    % same per-volume rot90 as export_orient_rgb_tif.m
BAR_UM         = 5000;                    % 0 = no scale bar
BAR_THICK_PX   = 6;
BAR_MARGIN_PX  = 20;
BAR_COLOR      = uint8([255 255 255]);
COMPRESSION    = Tiff.Compression.LZW;

ROOT    = '/projectnb/npbssmic/ns/CAA/';
OUT_DIR = fullfile(ROOT, 'figures/orientation/vessel_isosurface', datestr(now, 'yyyymmdd_HHMMSS'));
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

prm = struct('vox', VOX_UM, 'sigmaD', SIGMA_D_UM, 'sigmaI', SIGMA_I_UM, 'minCoh', MIN_TUBE_COH, ...
    'lowCol', LOWCOH_COLOR, 'gamma', RGB_GAMMA, 'skelMinBranch', SKEL_MIN_BRANCH_VOX, ...
    'bpRadius', BRANCHPT_RADIUS_VOX, 'minMarker', MIN_MARKER_VOX, 'wallVox', WALL_VOX, ...
    'zeroAxis', ORIENT_ZERO_AXIS, 'ccw', ORIENT_CCW, 'junction', JUNCTION_UM, 'ds', TENSOR_DOWNSAMPLE);
tp = struct('vox', VOX_UM, 'comp', COMPRESSION, 'gmCol', TIF_GM_COLOR, 'wmCol', TIF_WM_COLOR, ...
    'barPx', round(BAR_UM / VOX_UM), 'barThick', BAR_THICK_PX, 'barMargin', BAR_MARGIN_PX, 'barCol', BAR_COLOR);
CLASS_NAME  = {'vessel', 'epvs'};
CLASS_LABEL = {'vessel', 'EPVS'};

%% 1. Data: use caa if present; load only missing control subjects --------------
MAT_PATH = struct('caa6',  'caa6/caa6.mat', ...
                  'caa17', 'caa17/occip/caa17.mat', ...
                  'caa22', 'caa22/caa22.mat', ...
                  'caa25', 'caa25/caa25.mat', ...
                  'caa26', 'caa26/caa26.mat');
if exist('caa', 'var') && isstruct(caa)
    fprintf('Using caa from workspace.\n');
else
    caa = struct();
end
for s = 1:numel(SUBJECTS)
    key = regexprep(lower(SUBJECTS{s}), '[^a-z0-9]', '');
    have = ~isempty(fieldnames(caa)) && ...
           any(strcmp(regexprep(lower(fieldnames(caa)), '[^a-z0-9]', ''), key));
    if ~have
        if exist(key, 'var') && isstruct(eval(key))          % e.g. a separate caa26 variable
            fprintf('Using %s from workspace.\n', key); caa.(key) = eval(key);
        else
            fprintf('Loading %s\n', upper(key));
            tmp = load(fullfile(ROOT, MAT_PATH.(key)), key);
            caa.(key) = tmp.(key); clear tmp
        end
    end
end

write_color_key(fullfile(OUT_DIR, 'direction_color_key.png'), prm, 0);   % unrotated (PLY / NIfTI frame)

%% 2. Per subject x region -------------------------------------------------------
SUMMARY = table();
pad = ceil(3 * max(SIGMA_I_UM) / VOX_UM) + 2;
for s = 1:numel(SUBJECTS)
    for r = 1:numel(REGIONS)
        sub = SUBJECTS{s}; reg = REGIONS{r}; base = sprintf('%s_%s', sub, reg);
        VOL_DIR = fullfile(OUT_DIR, base);                  % one subfolder per subject x region
        if ~exist(VOL_DIR, 'dir'), mkdir(VOL_DIR); end
        fprintf('\n%s %s\n', sub, reg); tV = tic;

        % ---- combined white-matter mask and class labels ----
        Vd = get_vol(caa, sub, reg);
        szFull = size(Vd.seg); if numel(szFull) < 3, szFull(3) = 1; end
        scope = scope_mask(Vd, SCOPE, szFull);
        Mt = []; WMt = [];
        if WRITE_TIF && strcmp(TIF_BACKGROUND, 'tissue')
            if isfield(Vd, 'mask'), Mt = logical(Vd.mask); else, Mt = true(szFull); end
            WMt = logical(Vd.mask_wm) & Mt;
        end
        ves = logical(Vd.seg) & scope;
        if isfield(Vd, 'epvs')
            epv = logical(Vd.epvs) & scope;
        else
            warning('%s: no epvs field; all structures are vessel.', base); epv = false(szFull);
        end
        clear Vd scope
        fprintf('  vessel-segmentation voxels inside the EPVS annotation (counted as EPVS): %d\n', nnz(ves & epv));
        comb = bwareaopen(ves | epv, MIN_OBJ_VOX, 26); clear ves
        if ~any(comb(:)), warning('%s: no vessel or EPVS voxels in scope, skipped.', base); continue; end
        k_rot = 0; if ROTATE_UPRIGHT, k_rot = rot_k(sub, reg); end
        k_ren = 0; if RENDER_ROTATE_UPRIGHT, k_ren = rot_k(sub, reg); end
        if WRITE_PNG, write_color_key(fullfile(VOL_DIR, [base '_direction_color_key.png']), prm, k_ren); end

        [r1, r2, r3] = bbox_ranges(comb, pad);
        combC = comb(r1, r2, r3); epvC = epv(r1, r2, r3) & combC; clear comb epv
        vesC  = combC & ~epvC;
        off   = [r1(1) r2(1) r3(1)] - 1;
        szC   = size(combC); if numel(szC) < 3, szC(3) = 1; end
        vIdx  = find(combC); nVx = numel(vIdx);
        cls   = 1 + double(epvC(vIdx));                     % 1 = vessel, 2 = EPVS
        rowOf = zeros(szC, 'uint32'); rowOf(vIdx) = 1:nVx;
        [~, nIdx] = bwdist(combC);                          % nearest combined voxel (wall assignment)
        fprintf('  combined mask: %d voxels (%d vessel, %d EPVS), crop %s\n', nVx, nnz(cls == 1), ...
            nnz(cls == 2), mat2str(szC));

        % ---- junction distance: each voxel to the nearest voxel of the other class (um) ----
        dJ = inf(nVx, 1);
        if any(epvC(:)) && any(vesC(:))
            dE = bwdist(epvC); dJ(cls == 1) = double(dE(vIdx(cls == 1))) * VOX_UM; clear dE
            dV = bwdist(vesC); dJ(cls == 2) = double(dV(vIdx(cls == 2))) * VOX_UM; clear dV
        end
        clear epvC vesC

        % ---- Visualisation 1: voxelwise tensor on the combined mask ----
        fprintf('  [%5.0f s] voxelwise structure tensor (combined mask) ...\n', toc(tV));
        vt = vessel_tensor(combC, vIdx, prm);
        axT = single(vt.t); axCoh = single(vt.coh); axSig = single(vt.sigma);
        col1 = single(dir_rgb(vt.t, vt.coh, prm));
        okV  = vt.coh >= MIN_TUBE_COH;
        clear vt

        % ---- Visualisation 2: watershed on the combined mask, pieces = branch x class ----
        fprintf('  [%5.0f s] watershed branch separation (combined mask) ...\n', toc(tV));
        [L0, nB0, sk] = branch_watershed(combC, prm);
        [~, ~, pid] = unique((double(L0(vIdx)) - 1) * 2 + cls);
        nP = max(pid);
        pClass  = accumarray(pid, cls, [nP 1], @max);
        pParent = accumarray(pid, double(L0(vIdx)), [nP 1], @max);
        clear L0
        L = zeros(szC, 'uint32'); L(vIdx) = pid;
        fprintf('  %d branches -> %d pieces (%d vessel, %d EPVS)\n', nB0, nP, nnz(pClass == 1), nnz(pClass == 2));
        Bt = branch_tensor(combC, L, nP, nIdx, prm);
        if strcmp(BRANCH_MODE, 'local')
            fprintf('  [%5.0f s] local tensor within each piece ...\n', toc(tV));
            vl = branch_local_tensor(combC, L, nP, nIdx, prm, vIdx);
            axB = single(vl.t); cohB = single(vl.coh); clear vl
        else
            axB = single(Bt.t(pid, :)); cohB = single(Bt.coh(pid));
        end
        col2 = single(dir_rgb(axB, cohB, prm));
        fprintf('  [%5.0f s] tensors done\n', toc(tV));
        clear nIdx

        T = branch_table(L, nP, sk, Bt, off, prm, pClass, pParent, pid, dJ, vIdx, CLASS_NAME);
        writetable(T, fullfile(VOL_DIR, [base '_pieces.csv']));
        clear sk
        if WRITE_NII
            write_nifti_sparse(fullfile(VOL_DIR, [base '_piece_labels.nii.gz']), szFull, off, szC, vIdx, pid, 'int32', VOX_UM);
            write_nifti_sparse(fullfile(VOL_DIR, [base '_class.nii.gz']), szFull, off, szC, vIdx, cls, 'uint8', VOX_UM);
            dJw = dJ; dJw(~isfinite(dJw)) = -1;
            write_nifti_sparse(fullfile(VOL_DIR, [base '_junction_dist.nii.gz']), szFull, off, szC, vIdx, dJw, 'single', VOX_UM);
            clear dJw
        end
        if WRITE_LABELS
            save(fullfile(VOL_DIR, [base '_pieces.mat']), 'L', 'pClass', 'pParent', 'off', 'szFull', 'VOX_UM', '-v7.3');
        end
        clear L
        if WRITE_TENSOR
            axW = axT; axW(~isfinite(axW)) = 0;
            write_nifti_sparse(fullfile(VOL_DIR, [base '_axis_voxel.nii.gz']), szFull, off, szC, vIdx, axW, 'single', VOX_UM);
            clear axW
            TV = tensor_struct(vIdx, szC, off, szFull, cls, axT, axCoh, axSig, axB, cohB, BRANCH_MODE, dJ, pid, prm, CLASS_NAME);
            TV.meta.subject = sub; TV.meta.region = reg; TV.meta.scope = SCOPE;
            save(fullfile(VOL_DIR, [base '_tensor_voxel.mat']), 'TV', '-v7.3');
            fprintf('  wrote %s\n', fullfile(VOL_DIR, [base '_tensor_voxel.mat']));
            clear TV
        end
        clear axT axCoh axSig axB cohB

        % ---- split by class: meshes, stacks, renders ----
        for c = 1:2
            st = CLASS_NAME{c}; stl = CLASS_LABEL{c}; tag = sprintf('%s_%s', base, st);
            rows = find(cls == c);
            if isempty(rows), warning('%s: no %s voxels, skipped.', base, stl); continue; end
            fprintf(' %s: %d voxels\n', stl, numel(rows));
            vIdxC = vIdx(rows);
            segMC = false(szC); segMC(vIdxC) = true;
            [~, nIdxC] = bwdist(segMC);

            Ism = smooth_mask(segMC, ISO_SMOOTH_VOX); clear segMC
            [F, P] = mask_mesh(Ism, ISO_LEVEL, MESH_DOWNSAMPLE);
            Fr = F; Pr = P;
            if WRITE_PNG && size(F, 1) > RENDER_MAX_FACES
                dsR = MESH_DOWNSAMPLE * ceil(sqrt(size(F, 1) / RENDER_MAX_FACES));
                [Fr, Pr] = mask_mesh(Ism, ISO_LEVEL, dsR);
            end
            clear Ism
            vRow  = vertex_voxel(P,  nIdxC, rowOf);            % rows into the combined voxel list
            vRowR = vertex_voxel(Pr, nIdxC, rowOf); clear nIdxC
            Pum = (P + off - 1) * VOX_UM;
            fprintf('  mesh: %d faces (render mesh %d)\n', size(F, 1), size(Fr, 1));

            COL = {col1, col2}; SFX = {'_dir_voxel', '_dir_branch'};
            TTL = {'voxelwise structure tensor', sprintf('per watershed piece (%s)', BRANCH_MODE)};
            for v = 1 + ~WRITE_VOXEL_VIS : 2
                out = fullfile(VOL_DIR, [tag SFX{v}]); Cv = COL{v};
                if WRITE_PLY, write_ply([out '.ply'], Pum, F, Cv(vRow, :)); end
                if WRITE_TIF, write_rgb_stack([out '.tif'], szFull, off, szC, vIdxC, Cv(rows, :), Mt, WMt, k_rot, tp); end
                if WRITE_NII, write_nifti_rgb([out '.nii.gz'], szFull, off, szC, vIdxC, Cv(rows, :), VOX_UM); end
                if WRITE_PNG
                    render_views(Fr, Pr, Cv(vRowR, :), out, RENDER_VIEWS, VOX_UM, k_ren, szC, ...
                        sprintf('%s %s: %s axis, %s (combined vessel + EPVS tensor)', sub, reg, stl, TTL{v}));
                end
            end

            Tc = T(T.class == st, :);
            SUMMARY = [SUMMARY; table(string(sub), string(reg), string(st), numel(rows), mean(okV(rows)), ...
                mean(dJ(rows) <= JUNCTION_UM), height(Tc), median(Tc.nvox), median(Tc.skel_length_um), ...
                mean(Tc.valid_axis), size(F, 1), ...
                'VariableNames', {'subject','region','structure','n_vox','frac_voxel_axis_valid', ...
                'frac_vox_near_junction','n_pieces','median_piece_vox','median_piece_length_um', ...
                'frac_piece_axis_valid','n_faces'})]; %#ok<AGROW>
            clear F P Fr Pr Pum vRow vRowR vIdxC rows
        end
        clear combC rowOf vIdx cls dJ col1 col2 okV pid T Mt WMt
    end
end

disp(SUMMARY);
writetable(SUMMARY, fullfile(OUT_DIR, 'summary.csv'));
fid = fopen(fullfile(OUT_DIR, 'README.txt'), 'w');
fprintf(fid, ['Combined mask = (seg | epvs) & %s, fragments < %d voxels removed. Structure tensor (sigma_d %g um, ' ...
    'sigma_i %s um), skeleton and watershed computed once on the combined mask; classes from the manual EPVS ' ...
    'annotation (EPVS = epvs; vessel = combined & ~epvs).\n' ...
    'Axis direction colour: R = |t.dim1|, G = |t.dim2|, B = |t.dim3| (dim3 = through-plane); gray = tube ' ...
    'coherence < %g.\nPieces: marker-controlled watershed branches (markers = bwskel, MinBranchLength %d, minus ' ...
    'branch points dilated by %d voxel(s), pieces >= %d voxels) split at the EPVS boundary. Piece tensor = sum ' ...
    'of combined-mask gradient outer products over the piece''s wall voxels (<= %d voxels from the surface); ' ...
    'mode %s.\nJunction distance dJ = distance to the nearest voxel of the other class; near junction = ' ...
    'dJ <= %g um.\nPLY: binary, per-vertex RGB, um, x = dim1, y = dim2, z = dim3, voxel (1,1,1) at the origin. ' ...
    'TIF: full volume, slices along dim3, rot90 upright = %d, background %s. NIfTI: native grid, voxel (0,0,0) ' ...
    'at the origin, mm. pieces.mat: L occupies off(1)+(1:size(L,1)), off(2)+(1:size(L,2)), ' ...
    'off(3)+(1:size(L,3)) of a volume of size szFull; pClass 1 = vessel, 2 = EPVS.\n'], ...
    SCOPE, MIN_OBJ_VOX, SIGMA_D_UM, mat2str(SIGMA_I_UM), MIN_TUBE_COH, SKEL_MIN_BRANCH_VOX, ...
    BRANCHPT_RADIUS_VOX, MIN_MARKER_VOX, WALL_VOX, BRANCH_MODE, JUNCTION_UM, ROTATE_UPRIGHT, TIF_BACKGROUND);
fclose(fid);
fprintf('\nDone. Output: %s\n', OUT_DIR);

%% ========================= Local functions ==================================

function m = scope_mask(Vd, sc, szFull)
    switch sc
        case 'wm',     m = logical(Vd.mask_wm); if isfield(Vd, 'mask'), m = m & logical(Vd.mask); end
        case 'tissue', if isfield(Vd, 'mask'), m = logical(Vd.mask); else, m = true(szFull); end
        otherwise,     m = true(szFull);
    end
end

function V = get_vol(caa, sub, reg)
    sf = resolve_field(caa, sub, true);
    rf = resolve_field(caa.(sf), reg, false);
    V = caa.(sf).(rf);
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

function [r1, r2, r3] = bbox_ranges(m, pad)
    sz = size(m); if numel(sz) < 3, sz(3) = 1; end
    a1 = find(any(any(m, 2), 3)); a2 = find(any(any(m, 1), 3)); a3 = find(any(any(m, 1), 2));
    r1 = max(1, a1(1)-pad):min(sz(1), a1(end)+pad);
    r2 = max(1, a2(1)-pad):min(sz(2), a2(end)+pad);
    r3 = max(1, a3(1)-pad):min(sz(3), a3(end)+pad);
end

% ---- direction colour ----

function C = dir_rgb(t, coh, prm)
    % RGB = |unit axis| components (dim1, dim2, dim3); low-coherence axes in prm.lowCol
    C = abs(t) ./ max(vecnorm(t, 2, 2), eps);
    if prm.gamma ~= 1, C = C .^ prm.gamma; end
    bad = ~(coh >= prm.minCoh) | any(~isfinite(t), 2);
    C(bad, :) = repmat(prm.lowCol, nnz(bad), 1);
end

function write_color_key(fname, prm, k)
    % colour sphere drawn in the same display frame as the renders (rot90 k in-plane)
    [x, y, z] = sphere(200);                          % x = dim1, y = dim2, z = dim3 component
    C = reshape(dir_rgb([x(:) y(:) z(:)], ones(numel(x), 1), prm), [size(x) 3]);
    RC = [x(:) y(:)];                                 % [row col] components, rotated about the centre like rot90
    for s = 1:abs(k)
        if k > 0, RC = [-RC(:, 2), RC(:, 1)]; else, RC = [RC(:, 2), -RC(:, 1)]; end
    end
    rr = reshape(RC(:, 1), size(x)); cc = reshape(RC(:, 2), size(x));
    f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 900 800]); ax = axes(f);
    surf(ax, cc, rr, z, C, 'EdgeColor', 'none'); set(ax, 'YDir', 'reverse');
    axis(ax, 'equal'); view(ax, [-35 30]); lighting(ax, 'gouraud'); camlight(ax, 'headlight');
    if mod(k, 2) == 0
        xl = 'horizontal = dim2 (G)'; yl = 'vertical = dim1 (R)';
    else
        xl = 'horizontal = dim1 (R)'; yl = 'vertical = dim2 (G)';
    end
    xlabel(ax, xl); ylabel(ax, yl); zlabel(ax, 'dim3: through-plane (B)');
    title(ax, sprintf('Axis direction: RGB = |dim1|, |dim2|, |dim3| component (display rot90 k = %d)', k));
    exportgraphics(f, fname, 'Resolution', 200); close(f);
end

% ---- meshes ----

function I = smooth_mask(seg, sm)
    if sm > 0, I = imgaussfilt3(single(seg), sm); else, I = single(seg); end
end

function [F, P] = mask_mesh(Ism, level, ds)
    % isosurface of the smoothed mask; P = vertex positions in voxel units (dim1, dim2, dim3)
    I = Ism;
    if ds > 1, I = imresize3(I, max(1, ceil(size(Ism) / ds)), 'linear'); end
    szD = size(I); if numel(szD) < 3, szD(3) = 1; end
    I = padarray(I, [1 1 1], 0);                      % closes surfaces at the crop edge
    [F, Q] = isosurface(I, level);                    % Q = [x y z] = [dim2 dim1 dim3], padded
    clear I
    P = Q(:, [2 1 3]) - 1;
    szO = size(Ism); if numel(szO) < 3, szO(3) = 1; end
    if ds > 1, P = (P - 0.5) .* (szO ./ szD) + 0.5; end
    F = F(:, [1 3 2]);                                % keep outward winding after the x/y swap
end

function row = vertex_voxel(P, nIdx, rowOf)
    % row (into the vessel-voxel list) of the vessel voxel nearest to each vertex
    sz = size(nIdx); if numel(sz) < 3, sz(3) = 1; end
    s  = min(max(round(P), 1), sz);
    li = sub2ind(sz, s(:,1), s(:,2), s(:,3));
    row = double(rowOf(nIdx(li)));
end

function write_ply(fname, P, F, C)
    % binary little-endian PLY with per-vertex RGB
    nVt = size(P, 1); nF = size(F, 1);
    fid = fopen(fname, 'w'); assert(fid > 0, 'Cannot write %s', fname);
    fprintf(fid, ['ply\nformat binary_little_endian 1.0\n' ...
        'comment vessel axis direction; units um; x = dim1, y = dim2, z = dim3\n' ...
        'element vertex %d\nproperty float x\nproperty float y\nproperty float z\n' ...
        'property uchar red\nproperty uchar green\nproperty uchar blue\n' ...
        'element face %d\nproperty list uchar int vertex_indices\nend_header\n'], nVt, nF);
    xyz = single(P)';
    fwrite(fid, [reshape(typecast(xyz(:), 'uint8'), 12, nVt); uint8(round(255 * min(max(C, 0), 1)))'], 'uint8');
    fi = int32(F - 1)';
    fwrite(fid, [repmat(uint8(3), 1, nF); reshape(typecast(fi(:), 'uint8'), 12, nF)], 'uint8');
    fclose(fid);
end

function render_views(F, P, C, base, views, vox, k, szC, ttl)
    % plot X = image columns, Y = image rows (reversed), Z = dim3, after rotating the in-plane
    % coordinates exactly like rot90(slice, k): the en-face view matches the rotated TIF stacks
    RC = rot_rc(P(:, 1:2), k, szC(1:2));
    f = figure('Visible', 'off', 'Color', 'k', 'Position', [50 50 1400 1100]);
    ax = axes(f, 'Color', 'k');
    patch(ax, 'Faces', F, 'Vertices', [RC(:, 2), RC(:, 1), P(:, 3)] * vox, 'FaceVertexCData', C, ...
        'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'BackFaceLighting', 'lit', ...
        'AmbientStrength', 0.45, 'DiffuseStrength', 0.7, 'SpecularStrength', 0.1);
    set(ax, 'YDir', 'reverse'); axis(ax, 'equal'); axis(ax, 'tight'); axis(ax, 'off');
    if k ~= 0, ttl = sprintf('%s | rot90 k = %d', ttl, k); end
    title(ax, ttl, 'Color', 'w', 'Interpreter', 'none');
    lt = camlight(ax, 'headlight');
    for i = 1:numel(views)
        view(ax, views{i}); camlight(lt, 'headlight');
        exportgraphics(f, sprintf('%s_view%d.png', base, i), 'Resolution', 200, 'BackgroundColor', 'k');
    end
    close(f);
end

function RC = rot_rc(RC, k, sz)
    % [row col] of a point after rot90(image of size sz, k)
    for s = 1:abs(k)
        if k > 0, RC = [sz(2) + 1 - RC(:, 2), RC(:, 1)];
        else,     RC = [RC(:, 2), sz(1) + 1 - RC(:, 1)];
        end
        sz = sz([2 1]);
    end
end

% ---- structure tensor ----

function vt = vessel_tensor(seg, aIdx, prm)
    % voxelwise tensor of the smoothed mask at aIdx; scale chosen per voxel by max tube coherence
    G = mask_gradient(seg, prm);
    A = idx2sub3(size(seg), aIdx);
    pairs = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
    nA = numel(aIdx);
    vt = struct('t', nan(nA, 3), 'coh', -inf(nA, 1), 'sigma', nan(nA, 1));
    for k = 1:numel(prm.sigmaI)
        Jv = zeros(nA, 6);
        for p = 1:6
            Jv(:, p) = smooth_sample(G{pairs(p,1)} .* G{pairs(p,2)}, prm.sigmaI(k) / prm.vox, aIdx, A, prm.ds);
        end
        [t, coh] = tensor_axis(Jv);
        upd = coh > vt.coh;
        vt.t(upd, :) = t(upd, :); vt.coh(upd) = coh(upd); vt.sigma(upd) = prm.sigmaI(k);
    end
    vt.coh(~isfinite(vt.coh)) = NaN;
end

function A = idx2sub3(sz, idx)
    if numel(sz) < 3, sz(3) = 1; end
    [a1, a2, a3] = ind2sub(sz, idx); A = [a1 a2 a3];
end

function v = smooth_sample(X, sv, idx, A, useDS)
    % Gaussian (sigma sv voxels) of X sampled at idx / subscripts A. For sv >= 4 voxels (useDS), X is
    % block-averaged by f = 2^floor(log2(sv/2)), filtered with the matching sigma and interpolated.
    f = 1; if useDS, f = 2 ^ floor(log2(max(sv, 1) / 2)); end
    if f <= 1 || any(size(X) < 2 * f)
        Y = imgaussfilt3(X, sv); v = double(Y(idx)); return
    end
    B = imgaussfilt3(block_mean(X, f), sqrt(max(sv^2 - (f^2 - 1) / 12, 0.25)) / f);
    u = (A - 0.5) / f + 0.5;                        % voxel subscript -> block coordinate
    F = griddedInterpolant(B, 'linear', 'nearest');
    v = double(F(u(:,1), u(:,2), u(:,3)));
end

function B = block_mean(A, f)
    sz = size(A); if numel(sz) < 3, sz(3) = 1; end
    n = ceil(sz / f);
    if any(n * f ~= sz), A(n(1)*f, n(2)*f, n(3)*f) = 0; end
    A = reshape(A, f, n(1), f, n(2), f, n(3));
    B = reshape(sum(sum(sum(A, 1), 3), 5), n) / f^3;
end

function [t, coh] = tensor_axis(Jv)
    % Jv: n x 6 [J11 J12 J13 J22 J23 J33]; smallest-eigenvalue eigenvector and tube coherence
    n = size(Jv, 1);
    J = zeros(3, 3, n);
    J(1,1,:) = Jv(:,1); J(1,2,:) = Jv(:,2); J(1,3,:) = Jv(:,3);
    J(2,1,:) = Jv(:,2); J(2,2,:) = Jv(:,4); J(2,3,:) = Jv(:,5);
    J(3,1,:) = Jv(:,3); J(3,2,:) = Jv(:,5); J(3,3,:) = Jv(:,6);
    if exist('pageeig', 'builtin') || exist('pageeig', 'file')
        [Vp, Lp] = pageeig(J);
        mu = [squeeze(Lp(1,1,:)), squeeze(Lp(2,2,:)), squeeze(Lp(3,3,:))];
        if n == 1, mu = mu(:)'; end
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

function G = mask_gradient(seg, prm)
    Is = imgaussfilt3(single(seg), prm.sigmaD / prm.vox);
    [g2, g1, g3] = gradient(Is);
    G = {g1, g2, g3};
end

% ---- branches ----

function [L, nB, sk] = branch_watershed(seg, prm)
    % marker-controlled watershed on -(distance to wall); markers = skeleton minus branch points
    sk = bwskel(seg, 'MinBranchLength', prm.skelMinBranch);
    bp = bwmorph3(sk, 'branchpoints');
    r  = prm.bpRadius;
    if r > 0, bp = imdilate(bp, true(2*r + 1, 2*r + 1, 2*r + 1)); end
    M  = bwareaopen(sk & ~bp, prm.minMarker, 26); clear bp
    D  = bwdist(~seg);
    % every connected vessel component needs a marker: add its deepest voxel where missing
    CC = bwlabeln(seg, 26); nCC = max(CC(:));
    has = false(nCC, 1); has(unique(CC(M))) = true;
    if ~all(has)
        idx = find(seg); [~, o] = sort(D(idx), 'descend'); idx = idx(o);
        [cc, first] = unique(CC(idx), 'stable');
        M(idx(first(~has(cc)))) = true;
    end
    clear CC
    W = imimposemin(-D, M); clear D M
    W(~seg) = max(W(isfinite(W))) + 1;
    Lw = watershed(W, 26); clear W
    Lw(~seg) = 0;
    gap = seg & Lw == 0;                              % watershed lines inside vessels
    if any(gap(:))
        [~, ni] = bwdist(Lw > 0);
        Lw(gap) = Lw(ni(gap)); clear ni
    end
    [~, ~, j] = unique(double(Lw(seg))); clear Lw
    L = zeros(size(seg), 'uint32'); L(seg) = j;
    nB = max([j; 0]);
end

function [wIdx, wLab] = wall_voxels(seg, L, nIdx, wv)
    % voxels within wv of the vessel surface (both sides), labelled by the nearest vessel voxel's branch
    wall = (~seg & bwdist(seg) <= wv) | (seg & bwdist(~seg) <= wv);
    wIdx = find(wall); clear wall
    src = double(wIdx); outside = ~seg(wIdx);
    src(outside) = double(nIdx(wIdx(outside)));
    wLab = double(L(src));
    ok = wLab > 0; wIdx = wIdx(ok); wLab = wLab(ok);
end

function B = branch_tensor(seg, L, nB, nIdx, prm)
    % one structure tensor per branch: gradient outer products summed over the branch's wall voxels
    [wIdx, wLab] = wall_voxels(seg, L, nIdx, prm.wallVox);
    G = mask_gradient(seg, prm);
    pairs = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
    Jv = zeros(nB, 6);
    for p = 1:6
        Jv(:, p) = accumarray(wLab, double(G{pairs(p,1)}(wIdx)) .* double(G{pairs(p,2)}(wIdx)), [nB 1]);
    end
    [B.t, B.coh] = tensor_axis(Jv);
    B.nWall = accumarray(wLab, 1, [nB 1]);
end

function vt = branch_local_tensor(seg, L, nB, nIdx, prm, vIdx)
    % voxelwise tensor using only the wall gradients of the voxel's own branch
    sz = size(seg); if numel(sz) < 3, sz(3) = 1; end
    [wIdx, wLab] = wall_voxels(seg, L, nIdx, prm.wallVox);
    G = mask_gradient(seg, prm);
    wallLab = zeros(sz, 'uint32'); wallLab(wIdx) = wLab; clear wIdx wLab
    labV = double(L(vIdx)); nV = numel(vIdx);
    rowsB = accumarray(labV, (1:nV)', [nB 1], @(x) {x});
    vt = struct('t', nan(nV, 3), 'coh', -inf(nV, 1));
    pad = ceil(3 * max(prm.sigmaI) / prm.vox) + 1;
    pairs = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
    for k = 1:nB
        rows = rowsB{k}; if isempty(rows), continue; end
        [a1, a2, a3] = ind2sub(sz, vIdx(rows));
        lo = max([min(a1) min(a2) min(a3)] - pad, 1); hi = min([max(a1) max(a2) max(a3)] + pad, sz);
        rg = {lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)};
        w  = single(wallLab(rg{:}) == k);
        gc = {G{1}(rg{:}), G{2}(rg{:}), G{3}(rg{:})};   % w restricts the products to this branch's wall
        Al = [a1 - lo(1) + 1, a2 - lo(2) + 1, a3 - lo(3) + 1];
        ci = sub2ind(hi - lo + 1, Al(:,1), Al(:,2), Al(:,3));
        Pp = cell(1, 6);
        for p = 1:6, Pp{p} = w .* gc{pairs(p,1)} .* gc{pairs(p,2)}; end   % products once, all scales
        for q = 1:numel(prm.sigmaI)
            Jv = zeros(numel(rows), 6);
            for p = 1:6
                Jv(:, p) = smooth_sample(Pp{p}, prm.sigmaI(q) / prm.vox, ci, Al, prm.ds);
            end
            [t, coh] = tensor_axis(Jv);
            upd = coh > vt.coh(rows);
            vt.t(rows(upd), :) = t(upd, :); vt.coh(rows(upd)) = coh(upd);
        end
    end
    vt.coh(~isfinite(vt.coh)) = NaN;
end

function T = branch_table(L, nP, sk, Bt, off, prm, pClass, pParent, pid, dJ, vIdx, names)
    % one row per piece (watershed branch x class)
    [a1, a2, a3] = ind2sub(size(L), vIdx);
    n = accumarray(pid, 1, [nP 1]);
    c = [accumarray(pid, a1, [nP 1]), accumarray(pid, a2, [nP 1]), accumarray(pid, a3, [nP 1])] ./ n;
    c = (c + off - 1) * prm.vox;
    skl  = accumarray(double(L(sk & L > 0)), 1, [nP 1]) * prm.vox;
    minJ = accumarray(pid, dJ, [nP 1], @min);
    nearJ = accumarray(pid, double(dJ <= prm.junction), [nP 1]) ./ n;
    t = Bt.t;
    thE  = rad2deg(display_to_conv(atan2(-t(:,1), t(:,2)), prm.zeroAxis, prm.ccw));
    incl = asind(min(1, abs(t(:,3)))); incl(any(~isfinite(t), 2)) = NaN;
    cn = names(pClass); cn = string(cn(:));
    T = table((1:nP)', pParent, cn, n, n * prm.vox^3, skl, c(:,1), c(:,2), c(:,3), t(:,1), t(:,2), t(:,3), ...
        thE, incl, Bt.coh, Bt.coh >= prm.minCoh, Bt.nWall, minJ, nearJ, ...
        'VariableNames', {'piece_id','parent_branch_id','class','nvox','vol_um3','skel_length_um', ...
        'centroid_dim1_um','centroid_dim2_um','centroid_dim3_um','t_dim1','t_dim2','t_dim3', ...
        'enface_angle_deg','incl_deg','tube_coh','valid_axis','n_wall_vox','min_junction_dist_um', ...
        'frac_vox_near_junction'});
end

function TV = tensor_struct(vIdx, szC, off, szFull, cls, axT, axCoh, axSig, axB, cohB, mode, dJ, pid, prm, names)
    % voxelwise tensor results split by class; subscripts and indices refer to the full volume
    [a1, a2, a3] = ind2sub(szC, vIdx);
    sub3 = int32([a1 + off(1), a2 + off(2), a3 + off(3)]); clear a1 a2 a3
    idx  = uint64(sub2ind(szFull, double(sub3(:,1)), double(sub3(:,2)), double(sub3(:,3))));
    thE  = single(rad2deg(display_to_conv(atan2(-double(axT(:,1)), double(axT(:,2))), prm.zeroAxis, prm.ccw)));
    incl = single(asind(min(1, abs(double(axT(:,3))))));
    bad  = any(~isfinite(axT), 2); thE(bad) = NaN; incl(bad) = NaN;
    thB  = single(rad2deg(display_to_conv(atan2(-double(axB(:,1)), double(axB(:,2))), prm.zeroAxis, prm.ccw)));
    inclB = single(asind(min(1, abs(double(axB(:,3))))));
    bad  = any(~isfinite(axB), 2); thB(bad) = NaN; inclB(bad) = NaN;
    for c = 1:numel(names)
        m = cls == c;
        TV.(names{c}) = struct('idx', idx(m), 'sub', sub3(m, :), 'axis', axT(m, :), ...
            'enface_angle_deg', thE(m), 'incl_deg', incl(m), 'tube_coh', axCoh(m), ...
            'valid', axCoh(m) >= prm.minCoh, 'sigma_i_um', axSig(m), ...
            'junction_dist_um', single(dJ(m)), 'piece_id', uint32(pid(m)), ...
            'axis_piece', axB(m, :), 'enface_angle_piece_deg', thB(m), 'incl_piece_deg', inclB(m), ...
            'tube_coh_piece', cohB(m), 'valid_piece', cohB(m) >= prm.minCoh);
    end
    TV.meta = struct('size', szFull, 'vox_um', prm.vox, 'sigma_d_um', prm.sigmaD, 'sigma_i_um', prm.sigmaI, ...
        'min_tube_coh', prm.minCoh, 'orient_zero_axis', prm.zeroAxis, 'orient_ccw', prm.ccw, 'piece_mode', mode, ...
        'notes', ['axis = smallest-eigenvalue eigenvector of the structure tensor of the combined ' ...
        '(seg | epvs) & WM mask, columns (dim1, dim2, dim3), sign arbitrary; enface_angle_deg in the ' ...
        'D.orient convention; incl_deg = angle from the en-face plane; junction_dist_um = distance to ' ...
        'the other class (Inf if absent); piece_id matches <sub>_<reg>_pieces.csv. *_piece fields = visualisation 2 ' ...
        '(piece_mode local: voxelwise axis from the voxel''s own piece wall only; branch: the piece axis).']);
end

% ---- NIfTI (written directly, then gzipped; no toolbox needed) ----

function write_nifti_rgb(fname, szFull, off, szC, vIdx, C, voxUm)
    % uint8 [dim1 dim2 dim3 3]: frame c = colour channel c of the shown voxels, 0 elsewhere
    [a1, a2, a3] = ind2sub(szC, vIdx);
    li2 = sub2ind(szFull(1:2), a1 + off(1), a2 + off(2)); a3 = a3 + off(3); clear a1 a2
    C8  = uint8(round(255 * min(max(C, 0), 1)));
    cnt = accumarray(a3, 1, [szFull(3) 1]); st = [0; cumsum(cnt)];
    raw = fname(1:end-3);
    fid = nifti_open(raw, [szFull 3], 2, 8, voxUm / 1000);
    for c = 1:3
        for z = 1:szFull(3)
            sl = zeros(szFull(1), szFull(2), 'uint8');
            ii = st(z)+1 : st(z+1); sl(li2(ii)) = C8(ii, c);
            fwrite(fid, sl, 'uint8');
        end
    end
    fclose(fid); gzip_replace(raw);
    fprintf('  wrote %s\n', fname);
end

function write_nifti_sparse(fname, szFull, off, szC, vIdx, vals, cls, voxUm)
    % values at the cropped voxels vIdx placed into a full-size volume; cls: 'uint8' | 'int32' | 'single'
    switch cls
        case 'uint8',  code = [2 8];   prec = 'uint8';
        case 'int32',  code = [8 32];  prec = 'int32';
        case 'single', code = [16 32]; prec = 'float32';
    end
    [a1, a2, a3] = ind2sub(szC, vIdx);
    li2 = sub2ind(szFull(1:2), a1 + off(1), a2 + off(2)); a3 = a3 + off(3); clear a1 a2
    cnt = accumarray(a3, 1, [szFull(3) 1]); st = [0; cumsum(cnt)];
    vals = cast(vals, cls);
    raw = fname(1:end-3);
    nFr = size(vals, 2); dims = szFull; if nFr > 1, dims = [szFull nFr]; end   % columns = frames
    fid = nifti_open(raw, dims, code(1), code(2), voxUm / 1000);
    for f = 1:nFr
        for z = 1:szFull(3)
            sl = zeros(szFull(1), szFull(2), cls);
            ii = st(z)+1 : st(z+1); sl(li2(ii)) = vals(ii, f);
            fwrite(fid, sl, prec);
        end
    end
    fclose(fid); gzip_replace(raw);
    fprintf('  wrote %s\n', fname);
end

function fid = nifti_open(raw, dims, dtype, bitpix, voxMm)
    % NIfTI-1 single file: 348-byte header + 4-byte extension flag; sform/qform = diag(voxMm),
    % voxel (0,0,0) at the origin (consistent with the PLY coordinates); units mm
    fid = fopen(raw, 'w', 'ieee-le'); assert(fid > 0, 'Cannot write %s', raw);
    dim = ones(1, 8); dim(1) = numel(dims); dim(2:numel(dims)+1) = dims;
    fwrite(fid, 348, 'int32');
    fwrite(fid, zeros(1, 28), 'uint8');                 % data_type, db_name
    fwrite(fid, 0, 'int32'); fwrite(fid, 0, 'int16');   % extents, session_error
    fwrite(fid, [double('r') 0], 'uint8');              % regular, dim_info
    fwrite(fid, dim, 'int16');
    fwrite(fid, [0 0 0], 'float32');                    % intent_p1..3
    fwrite(fid, [0 dtype bitpix 0], 'int16');           % intent_code, datatype, bitpix, slice_start
    fwrite(fid, [1 voxMm voxMm voxMm 1 1 1 1], 'float32');   % pixdim (qfac = 1)
    fwrite(fid, [352 1 0], 'float32');                  % vox_offset, scl_slope, scl_inter
    fwrite(fid, 0, 'int16'); fwrite(fid, [0 2], 'uint8');    % slice_end, slice_code, xyzt_units = mm
    fwrite(fid, [0 0 0 0], 'float32');                  % cal_max, cal_min, slice_duration, toffset
    fwrite(fid, [0 0], 'int32');                        % glmax, glmin
    d = zeros(1, 80); t = double('as-PSOCT structure-tensor axis'); d(1:numel(t)) = t;
    fwrite(fid, d, 'uint8'); fwrite(fid, zeros(1, 24), 'uint8');   % descrip, aux_file
    fwrite(fid, [1 1], 'int16');                        % qform_code, sform_code (scanner)
    fwrite(fid, zeros(1, 6), 'float32');                % quatern_b/c/d, qoffset_x/y/z
    fwrite(fid, [voxMm 0 0 0, 0 voxMm 0 0, 0 0 voxMm 0], 'float32');   % srow_x/y/z
    fwrite(fid, zeros(1, 16), 'uint8');                 % intent_name
    fwrite(fid, [double('n+1') 0], 'uint8');            % magic
    fwrite(fid, [0 0 0 0], 'uint8');                    % no extensions
end

function gzip_replace(raw)
    % pigz (parallel gzip) when available, else MATLAB gzip
    persistent usePigz
    if isempty(usePigz), [st, ~] = system('command -v pigz'); usePigz = st == 0; end
    if usePigz && system(sprintf('pigz -f "%s"', raw)) == 0, return; end
    gzip(raw); delete(raw);
end

% ---- TIF stacks ----

function write_rgb_stack(fname, szFull, off, szC, vIdx, C, M, WM, k_rot, tp)
    % full-size RGB stack: coloured voxels on black or GM/WM background, slices along dim3
    [a1, a2, a3] = ind2sub(szC, vIdx);
    a1 = a1 + off(1); a2 = a2 + off(2); a3 = a3 + off(3);   % vIdx sorted => a3 ascending
    li2 = sub2ind(szFull(1:2), a1, a2); clear a1 a2
    C8  = uint8(round(255 * min(max(C, 0), 1)));
    cnt = accumarray(a3, 1, [szFull(3) 1]); st = [0; cumsum(cnt)];
    nz = szFull(3); n12 = szFull(1) * szFull(2);
    if 3 * prod(szFull) > 3.5e9, mode = 'w8'; else, mode = 'w'; end
    t = Tiff(fname, mode); cl = onCleanup(@() close(t));
    for z = 1:nz
        img = zeros(szFull(1), szFull(2), 3, 'uint8');
        if ~isempty(M)
            img = paint(img, M(:,:,z) & ~WM(:,:,z), tp.gmCol);
            img = paint(img, WM(:,:,z), tp.wmCol);
        end
        ii = st(z)+1 : st(z+1);
        for c = 1:3, img(li2(ii) + (c-1) * n12) = C8(ii, c); end
        if k_rot ~= 0, img = rot90(img, k_rot); end
        if tp.barPx > 0, img = burn_scale_bar(img, tp.barPx, tp.barThick, tp.barMargin, tp.barCol); end
        write_slice(t, img, z, nz, tp.vox, tp.comp);
    end
    clear cl
    fprintf('  wrote %s\n', fname);
end

function img = paint(img, m, color)
    for c = 1:3
        ch = img(:,:,c); ch(m) = color(c); img(:,:,c) = ch;
    end
end

function img = burn_scale_bar(img, len_px, thick_px, margin_px, color)
    [h, w, ~] = size(img);
    r2 = h - margin_px;  r1 = max(1, r2 - thick_px + 1);
    c2 = w - margin_px;  c1 = max(1, c2 - len_px + 1);
    if r2 < 1 || c2 < 1, return; end
    for c = 1:3, img(r1:r2, c1:c2, c) = color(c); end
end

function k = rot_k(sub, reg)
    % rot90 k making the dim-3 slices upright (top = superior), as in export_orient_rgb_tif.m
    s = lower(strrep(sub, '_', '')); isOcc = contains(lower(reg), 'occ'); k = 0;
    if isOcc && (strcmp(s, 'caa6') || strcmp(s, 'caa22')), k = -1;
    elseif isOcc && strcmp(s, 'caa26'), k = 1; end
end

function write_slice(t, img, z, nz, vox, comp)
    tag.ImageLength         = size(img,1);
    tag.ImageWidth          = size(img,2);
    tag.Photometric         = Tiff.Photometric.RGB;
    tag.BitsPerSample       = 8;
    tag.SamplesPerPixel     = 3;
    tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tag.Compression         = comp;
    tag.XResolution         = 1 / vox;
    tag.YResolution         = 1 / vox;
    tag.ResolutionUnit      = Tiff.ResolutionUnit.None;
    tag.Software            = 'MATLAB';
    if z == 1
        tag.ImageDescription = sprintf( ...
            'ImageJ=1.53t\nimages=%d\nslices=%d\nunit=micron\nspacing=%g\nloop=false\n', nz, nz, vox);
    end
    setTag(t, tag);
    write(t, img);
    if z < nz, writeDirectory(t); end
end

function th = display_to_conv(thDisp, zeroAxis, ccw)
    th0 = 0; if strcmp(zeroAxis, 'dim1'), th0 = -pi/2; end
    sgn = 1; if ~ccw, sgn = -1; end
    th = mod(sgn * (thDisp - th0), pi);
end