%% export_orient_rgb_tif.m
% Export caa.(sub).(reg).orient_rgb (4-D double) as masked RGB TIFF stacks
% with a burned-in scale bar. 20 um isotropic voxels; slices along dim 3.
% Background removed with caa.(sub).(reg).mask (overall tissue mask).

%% ---------------- Config ----------------
VOX_UM        = 20;                 % isotropic voxel size (um)
BAR_UM        = 5000;               % scale bar length (um)
BAR_THICK_PX  = 6;                  % scale bar thickness (px)
BAR_MARGIN_PX = 20;                 % offset from bottom-right corner (px)
BAR_COLOR     = uint8([255 255 255]);
GM_COLOR      = uint8([128 128 128]); % gray matter fill
OUTLINE_COLOR = uint8([255 255 255]); % GM outline
OUTLINE_PX    = 1;                  % outline thickness (px), grown into GM only
ROTATE_UPRIGHT = true;              % apply per-volume rotation so top = superior
COMPRESSION   = Tiff.Compression.LZW;

OUT_ROOT = '/projectnb/npbssmic/ns/CAA/figures/orientation/';
CAA_FILE = '';                      % set if caa is not already in the workspace

%% ---------------- Load ----------------
if ~exist('caa','var')
    if isempty(CAA_FILE)
        error('caa struct not in workspace and CAA_FILE not set.');
    end
    S = load(CAA_FILE, 'caa'); caa = S.caa; clear S
end

out_dir = fullfile(OUT_ROOT, datestr(now, 'yyyymmdd_HHMMSS'));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

bar_len_px = round(BAR_UM / VOX_UM);

%% ---------------- Export ----------------
subs = fieldnames(caa);
for i = 1:numel(subs)
    sub  = subs{i};
    fprintf('\nStarting subject %s',sub)
    regs = fieldnames(caa.(sub));
    for j = 1:numel(regs)
        reg = regs{j};
        D = caa.(sub).(reg);
        if ~isstruct(D) || ~isfield(D, 'orient_rgb'), continue; end
        if ~isfield(D, 'mask') || ~isfield(D, 'mask_wm')
            warning('%s.%s: missing mask or mask_wm, skipped.', sub, reg); continue
        end

        V = D.orient_rgb;
        % Channel-last: [X Y Z 3]
        if size(V,4) == 3
            % already channel-last
        elseif size(V,3) == 3
            V = permute(V, [1 2 4 3]);
        else
            warning('%s.%s: no size-3 channel dim in orient_rgb, skipped.', sub, reg); continue
        end

        M  = logical(D.mask);
        WM = logical(D.mask_wm);
        sz = size(V);
        if ~isequal(size(M), sz(1:3)) || ~isequal(size(WM), sz(1:3))
            warning('%s.%s: mask %s / mask_wm %s vs orient_rgb %s, skipped.', sub, reg, ...
                mat2str(size(M)), mat2str(size(WM)), mat2str(sz(1:3))); continue
        end
        WM = WM & M;          % WM restricted to tissue
        GM = M & ~WM;         % gray matter = mask minus mask_wm

        % Volume-wide scaling to uint8 (same mapping for every slice)
        vmax = max(V(:));
        if vmax <= 1, scl = 255; else, scl = 255 / vmax; end

        k_rot = 0;
        if ROTATE_UPRIGHT, k_rot = rot_k(sub, reg); end

        nz = sz(3);
        fname = fullfile(out_dir, sprintf('%s_%s_orient_rgb_WM_%dum_bar%dum.tif', ...
            sub, reg, VOX_UM, BAR_UM));

        % BigTIFF if the uncompressed stack exceeds ~3.5 GB
        if prod(sz) > 3.5e9, mode = 'w8'; else, mode = 'w'; end
        t = Tiff(fname, mode);
        cleanup = onCleanup(@() close(t));

        for z = 1:nz
            img = squeeze(V(:,:,z,:));              % X x Y x 3
            img(~isfinite(img)) = 0;
            img = uint8(min(max(img * scl, 0), 255));

            wm = WM(:,:,z);
            gm = GM(:,:,z);
            ol = bwperim(gm, 8);                    % GM boundary (GM/WM and GM/background)
            if OUTLINE_PX > 1
                ol = imdilate(ol, strel('disk', OUTLINE_PX - 1)) & gm;  % thicken inward, never over WM
            end
            img = paint(img, ~(wm | gm), uint8([0 0 0]));  % background -> black
            img = paint(img, gm, GM_COLOR);                % GM -> gray
            img = paint(img, ol, OUTLINE_COLOR);           % GM outline -> white
            if k_rot ~= 0, img = rot90(img, k_rot); end
            img = burn_scale_bar(img, bar_len_px, BAR_THICK_PX, BAR_MARGIN_PX, BAR_COLOR);

            tag.ImageLength         = size(img,1);
            tag.ImageWidth          = size(img,2);
            tag.Photometric         = Tiff.Photometric.RGB;
            tag.BitsPerSample       = 8;
            tag.SamplesPerPixel     = 3;
            tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tag.Compression         = COMPRESSION;
            tag.XResolution         = 1 / VOX_UM;   % px per um (ImageJ convention)
            tag.YResolution         = 1 / VOX_UM;
            tag.ResolutionUnit      = Tiff.ResolutionUnit.None;
            tag.Software            = 'MATLAB';
            if z == 1
                tag.ImageDescription = sprintf( ...
                    'ImageJ=1.53t\nimages=%d\nslices=%d\nunit=micron\nspacing=%g\nloop=false\n', ...
                    nz, nz, VOX_UM);
            elseif isfield(tag, 'ImageDescription')
                tag = rmfield(tag, 'ImageDescription');
            end
            setTag(t, tag);
            write(t, img);
            if z < nz, writeDirectory(t); end
        end
        clear cleanup   % closes t
        fprintf('Wrote %s (%d slices, rot90 k=%d)\n', fname, nz, k_rot);
    end
end
fprintf('Done. Output: %s\n', out_dir);

%% ---------------- Local functions ----------------
function img = paint(img, m, color)
    for c = 1:3
        ch = img(:,:,c);
        ch(m) = color(c);
        img(:,:,c) = ch;
    end
end

function img = burn_scale_bar(img, len_px, thick_px, margin_px, color)
    [h, w, ~] = size(img);
    r2 = h - margin_px;            r1 = max(1, r2 - thick_px + 1);
    c2 = w - margin_px;            c1 = max(1, c2 - len_px + 1);
    if r2 < 1 || c2 < 1, return; end
    for c = 1:3
        img(r1:r2, c1:c2, c) = color(c);
    end
end

function k = rot_k(sub, reg)
    % rot90 k needed to make the dim-3 projection upright (top = superior).
    % Appears 90 deg CCW -> rotate CW (k = -1); appears 90 deg CW -> rotate CCW (k = 1).
    s = lower(strrep(sub, '_', ''));
    isOcc = contains(lower(reg), 'occ');
    k = 0;
    if isOcc && (strcmp(s, 'caa6') || strcmp(s, 'caa22'))
        k = -1;
    elseif isOcc && strcmp(s, 'caa26')
        k = 1;
    end
end