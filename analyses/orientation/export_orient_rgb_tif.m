%% export_orient_rgb_tif.m
% Export caa.(sub).(reg).orient_rgb (4-D) as masked RGB TIFF stacks with a burned-in
% scale bar, plus vessel-axis orientation stacks. 20 um isotropic voxels; slices along
% dim 3. Background removed with caa.(sub).(reg).mask (overall tissue mask).
%
% Vessel axis: 3-D structure tensor of the smoothed vessel mask (seg), same settings as
% compare_orient_vessel_*.m, evaluated at every vessel voxel. Axis = smallest-eigenvalue
% eigenvector; integration scale chosen per voxel by max tube coherence. The en-face
% (dims 1-2) projection of the axis is converted to the D.orient convention and coloured
% with the same hue wheel as orient_rgb (fitted from orient vs orient_rgb when
% CALIBRATE_HUE = true), so a vessel parallel to its local fibers has the fibers' colour.
%   hue colour    : valid axis, inclined < MAX_INCL_DEG from the en-face plane
%   INCL_COLOR    : valid axis, inclined >= MAX_INCL_DEG (en-face angle not meaningful)
%   LOWCOH_COLOR  : tube coherence < MIN_TUBE_COH (no reliable axis)
%
% Outputs per subject x region:
%   <sub>_<reg>_orient_rgb_WM_*.tif           fiber orientation (unchanged)
%   <sub>_<reg>_vessel_axis_rgb_*.tif         vessels coloured by axis; GM gray + outline
%   <sub>_<reg>_orient_vessel_overlay_*.tif   orient_rgb with vessels painted in axis colour
%   <sub>_<reg>_colorwheel.png                hue key for that volume (already rotated to match);
%                                             '_vessel_only' suffix if the orient_rgb hue fit failed
%
% Note: rot90 (ROTATE_UPRIGHT) rotates the image but not the colour code; the exported
% colour wheels are rotated with their stacks, so each wheel is read directly on its image.

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

SUBJECTS = {'CAA6','CAA17','CAA22','CAA25','CAA26'};   % e.g. {'CAA26'} to test on one subject

WRITE_ORIENT  = true;               % fiber orientation stack (original output)
WRITE_VESSEL  = true;               % vessel-axis orientation stack
WRITE_OVERLAY = true;               % orient_rgb with vessels painted in vessel-axis colour

% Convention of D.orient (must match compare_orient_vessel_*.m)
ORIENT_UNITS     = 'deg';           % 'deg' | 'rad'
ORIENT_ZERO_AXIS = 'dim2';
ORIENT_CCW       = true;
CALIBRATE_HUE    = true;            % fit hue = f(orient) from orient_rgb so colours match
CAL_N            = 2e5;             % WM voxels sampled for the fit
CAL_MIN_R        = 0.9;             % below this the fit is rejected (non-HSV colormap)

% Vessel axis (same settings as the analysis scripts)
SIGMA_D_UM     = 20;
SIGMA_I_UM     = [40 80 160];
MIN_TUBE_COH   = 0.5;
MAX_INCL_DEG   = 45;
INCL_COLOR     = uint8([220 220 220]);
LOWCOH_COLOR   = uint8([70 70 70]);
VES_EDGE_COLOR = uint8([0 0 0]);    % 1-px vessel outline in the overlay stack

% Colour wheel (PNG, transparent background), one per volume: same hue fit and rot90 as its stacks
WRITE_WHEEL   = true;
WHEEL_PX      = 512;                % image size (px)
WHEEL_INNER   = 0.35;               % inner radius as fraction of outer (0 = full disk)

ROOT     = '/projectnb/npbssmic/ns/CAA/';
OUT_ROOT = fullfile(ROOT, 'figures/orientation/');

%% ---------------- Load: use workspace caa if present; load only missing subjects ----------------
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
    key = regexprep(lower(SUBJECTS{s}), '[^a-z0-9]', '');
    have = ~isempty(fieldnames(caa)) && ...
           any(strcmp(regexprep(lower(fieldnames(caa)), '[^a-z0-9]', ''), key));
    if ~have
        fprintf('Loading %s\n', upper(key));
        tmp = load(fullfile(ROOT, MAT_PATH.(key)), key);
        caa.(key) = tmp.(key); clear tmp
    end
end

out_dir = fullfile(OUT_ROOT, datestr(now, 'yyyymmdd_HHMMSS'));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

bar_len_px = round(BAR_UM / VOX_UM);
fin = @(img, k) finish_slice(img, k, bar_len_px, BAR_THICK_PX, BAR_MARGIN_PX, BAR_COLOR);

prm = struct('vox', VOX_UM, 'units', ORIENT_UNITS, 'zeroAxis', ORIENT_ZERO_AXIS, ...
    'ccw', ORIENT_CCW, 'sigmaD', SIGMA_D_UM, 'sigmaI', SIGMA_I_UM, ...
    'minCoh', MIN_TUBE_COH, 'maxIncl', MAX_INCL_DEG);

%% ---------------- Export ----------------
for i = 1:numel(SUBJECTS)
    sub = resolve_field(caa, SUBJECTS{i});   % actual field name in caa (e.g. 'caa26')
    fprintf('\nStarting subject %s\n', sub);
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
        if vmax <= 1, scl = 255; vnorm = 1; else, scl = 255 / vmax; vnorm = vmax; end

        % ---------- hue wheel: orient angle -> orient_rgb hue (fitted per volume) ----------
        doVes = WRITE_VESSEL || WRITE_OVERLAY;
        if doVes && (~isfield(D, 'seg') || ~isequal(size(D.seg), sz(1:3)))
            warning('%s.%s: seg missing or size mismatch, vessel stacks skipped.', sub, reg);
            doVes = false;
        end
        hueMap = [1 0];   % [sign offset]: hue = mod(sign*theta/180 + offset, 1)
        hueOK  = false;   % true = hueMap reproduces orient_rgb, so the wheel applies to fibers too
        if CALIBRATE_HUE && (doVes || WRITE_WHEEL)
            if isfield(D, 'orient') && isequal(size(D.orient), sz(1:3))
                [hm, Rcal] = calibrate_hue(V, D.orient, WM, ORIENT_UNITS, vnorm, CAL_N);
                fprintf('  hue fit: hue = mod(%+d*theta/180 + %.3f, 1), R = %.3f\n', hm(1), hm(2), Rcal);
                if Rcal >= CAL_MIN_R
                    hueMap = hm; hueOK = true;
                else
                    warning(['%s.%s: orient_rgb hue is not a linear function of orient (R = %.2f); ' ...
                        'using plain HSV wheel. Vessel and fiber colours are NOT directly comparable.'], sub, reg, Rcal);
                end
            else
                warning('%s.%s: no matching D.orient for hue fit; using plain HSV wheel.', sub, reg);
            end
        end

        % ---------- vessel-axis colours (computed once per volume) ----------
        if doVes
            fprintf('  vessel axis (structure tensor) ...\n');
            [vIdx, vRGB, vinfo] = vessel_colors(logical(D.seg) & M, prm, hueMap, INCL_COLOR, LOWCOH_COLOR);
            fprintf('  %d vessel voxels: %.1f%% in-plane axis, %.1f%% inclined, %.1f%% no valid axis\n', ...
                vinfo.n, 100*vinfo.fracGood, 100*vinfo.fracIncl, 100*vinfo.fracLow);
            n12 = sz(1) * sz(2);
            vz  = ceil(vIdx / n12);                    % find() order => already sorted by slice
            cnt = accumarray(vz, 1, [sz(3) 1]);
            st  = [0; cumsum(cnt)];
        end

        k_rot = 0;
        if ROTATE_UPRIGHT, k_rot = rot_k(sub, reg); end

        if WRITE_WHEEL
            if hueOK, wtag = 'colorwheel'; else, wtag = 'colorwheel_vessel_only'; end
            fw = fullfile(out_dir, sprintf('%s_%s_%s.png', sub, reg, wtag));
            write_colorwheel(fw, hueMap, ORIENT_ZERO_AXIS, ORIENT_CCW, k_rot, WHEEL_PX, WHEEL_INNER);
            fprintf('  wheel: %s\n', fw);
        end

        nz = sz(3); nbytes = prod(sz);
        sfx = sprintf('%dum_bar%dum', VOX_UM, BAR_UM);
        tO = []; tV = []; tX = []; written = {};
        if WRITE_ORIENT
            f = fullfile(out_dir, sprintf('%s_%s_orient_rgb_WM_%s.tif', sub, reg, sfx));
            tO = open_stack(f, nbytes); written{end+1} = f; %#ok<SAGROW>
        end
        if doVes && WRITE_VESSEL
            f = fullfile(out_dir, sprintf('%s_%s_vessel_axis_rgb_%s.tif', sub, reg, sfx));
            tV = open_stack(f, nbytes); written{end+1} = f; %#ok<SAGROW>
        end
        if doVes && WRITE_OVERLAY
            f = fullfile(out_dir, sprintf('%s_%s_orient_vessel_overlay_%s.tif', sub, reg, sfx));
            tX = open_stack(f, nbytes); written{end+1} = f; %#ok<SAGROW>
        end
        if isempty(written), continue; end
        cleanup = onCleanup(@() close_stacks({tO, tV, tX}));

        for z = 1:nz
            wm = WM(:,:,z);
            gm = GM(:,:,z);
            ol = bwperim(gm, 8);                    % GM boundary (GM/WM and GM/background)
            if OUTLINE_PX > 1
                ol = imdilate(ol, strel('disk', OUTLINE_PX - 1)) & gm;  % thicken inward, never over WM
            end
            bg = ~(wm | gm);

            % fiber orientation image (needed for the orient and overlay stacks)
            if ~isempty(tO) || ~isempty(tX)
                img = squeeze(V(:,:,z,:));          % X x Y x 3
                img(~isfinite(img)) = 0;
                img = uint8(min(max(img * scl, 0), 255));
                img = paint(img, bg, uint8([0 0 0]));  % background -> black
                img = paint(img, gm, GM_COLOR);        % GM -> gray
                img = paint(img, ol, OUTLINE_COLOR);   % GM outline -> white
                if ~isempty(tO), write_slice(tO, fin(img, k_rot), z, nz, VOX_UM, COMPRESSION); end
            end

            if doVes
                ii = st(z)+1 : st(z+1);
                li = vIdx(ii) - (z-1) * n12;        % in-slice linear indices
                if ~isempty(tV)
                    imv = zeros(sz(1), sz(2), 3, 'uint8');
                    imv = paint(imv, gm, GM_COLOR);
                    imv = paint(imv, ol, OUTLINE_COLOR);
                    imv = paint_idx(imv, li, vRGB(ii,:));
                    write_slice(tV, fin(imv, k_rot), z, nz, VOX_UM, COMPRESSION);
                end
                if ~isempty(tX)
                    ves = false(sz(1), sz(2)); ves(li) = true;
                    edge = imdilate(ves, true(3)) & ~ves & ~bg;
                    imx = paint_idx(img, li, vRGB(ii,:));
                    imx = paint(imx, edge, VES_EDGE_COLOR);
                    write_slice(tX, fin(imx, k_rot), z, nz, VOX_UM, COMPRESSION);
                end
            end
        end
        clear cleanup   % closes all stacks
        fprintf('Wrote (%d slices, rot90 k=%d):\n', nz, k_rot);
        fprintf('  %s\n', written{:});
    end
end
fprintf('Done. Output: %s\n', out_dir);

%% ---------------- Local functions ----------------
function f = resolve_field(s, key)
    fn = fieldnames(s);
    hit = find(strcmp(regexprep(lower(fn), '[^a-z0-9]', ''), regexprep(lower(key), '[^a-z0-9]', '')));
    if numel(hit) ~= 1
        error('Could not uniquely resolve "%s" among: %s', key, strjoin(fn', ', '));
    end
    f = fn{hit};
end

function img = paint(img, m, color)
    for c = 1:3
        ch = img(:,:,c);
        ch(m) = color(c);
        img(:,:,c) = ch;
    end
end

function img = paint_idx(img, li, rgb)
    n = size(img,1) * size(img,2);
    for c = 1:3, img(li + (c-1)*n) = rgb(:,c); end
end

function img = finish_slice(img, k, len_px, thick_px, margin_px, color)
    if k ~= 0, img = rot90(img, k); end
    img = burn_scale_bar(img, len_px, thick_px, margin_px, color);
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

% ---- TIFF writing ----

function t = open_stack(fname, nbytes)
    % BigTIFF if the uncompressed stack exceeds ~3.5 GB
    if nbytes > 3.5e9, mode = 'w8'; else, mode = 'w'; end
    t = Tiff(fname, mode);
end

function write_slice(t, img, z, nz, vox, comp)
    tag.ImageLength         = size(img,1);
    tag.ImageWidth          = size(img,2);
    tag.Photometric         = Tiff.Photometric.RGB;
    tag.BitsPerSample       = 8;
    tag.SamplesPerPixel     = 3;
    tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tag.Compression         = comp;
    tag.XResolution         = 1 / vox;   % px per um (ImageJ convention)
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

function close_stacks(ts)
    for k = 1:numel(ts)
        if ~isempty(ts{k}), try, close(ts{k}); catch, end, end
    end
end

% ---- colour mapping ----

function [map, R] = calibrate_hue(V, orient, WM, units, vnorm, nSamp)
    % Fit hue(orient_rgb) = mod(sgn * theta/180 + h0, 1) over WM voxels.
    % R = circular resultant of the residual (1 = exact linear hue wheel).
    n = numel(WM);
    o = double(orient);
    idx = find(WM & isfinite(o));
    if numel(idx) > nSamp, idx = idx(randperm(numel(idx), nSamp)); end
    rgb = double([V(idx), V(idx + n), V(idx + 2*n)]) / vnorm;
    rgb(~isfinite(rgb)) = 0;
    hsv = rgb2hsv(min(max(rgb, 0), 1));
    keep = hsv(:,2) > 0.2 & hsv(:,3) > 0.05;   % hue undefined for gray/dark voxels
    h  = hsv(keep, 1);
    th = o(idx(keep)); if strcmp(units, 'rad'), th = rad2deg(th); end
    u  = mod(th, 180) / 180;
    R = -1; map = [1 0];
    for sgn = [1 -1]
        zc = mean(exp(2i * pi * (h - sgn * u)));
        if abs(zc) > R, R = abs(zc); map = [sgn, mod(angle(zc) / (2*pi), 1)]; end
    end
end

function write_colorwheel(fname, hueMap, zeroAxis, ccw, k_rot, npx, inner)
    % Each pixel is coloured by the hue of a structure running through the wheel centre
    % at that pixel's displayed angle (axial: opposite sides share a colour). Built in the
    % unrotated display frame (imagesc, rows down) and rotated with the same rot90 as the stacks.
    c = (npx + 1) / 2;
    [x, y] = meshgrid(1:npx, 1:npx);
    dx = x - c; dy = c - y;                         % image rows increase downward
    r  = hypot(dx, dy) / (npx / 2 - 1);
    thDisp = atan2(dy, dx);                         % CCW from +dim2 as displayed
    th  = rad2deg(display_to_conv(thDisp, zeroAxis, ccw));
    hue = mod(hueMap(1) * th / 180 + hueMap(2), 1);
    rgb = hsv2rgb(cat(3, hue, ones(npx), ones(npx)));
    alpha = double(r <= 1 & r >= inner);
    if k_rot ~= 0, rgb = rot90(rgb, k_rot); alpha = rot90(alpha, k_rot); end
    imwrite(rgb, fname, 'Alpha', alpha);
end

function [vIdx, rgb, info] = vessel_colors(seg, prm, hueMap, inclColor, lowColor)
    sz = size(seg);
    vIdx = find(seg);
    nV = numel(vIdx);
    rgb = zeros(nV, 3, 'uint8');
    info = struct('n', nV, 'fracGood', 0, 'fracIncl', 0, 'fracLow', 0);
    if nV == 0, return; end

    % crop to vessel bounding box (+ filter support) to save memory
    pad = ceil(3 * max(prm.sigmaI) / prm.vox) + 2;
    [r1, r2, r3] = bbox_ranges(seg, pad);
    segC = seg(r1, r2, r3);
    [a1, a2, a3] = ind2sub(sz, vIdx);
    cIdx = sub2ind(size(segC), a1 - r1(1) + 1, a2 - r2(1) + 1, a3 - r3(1) + 1);
    clear a1 a2 a3

    vt  = vessel_tensor(segC, cIdx, prm);
    t   = vt.t;
    phi = rad2deg(asin(min(1, abs(t(:,3)))));
    th  = rad2deg(display_to_conv(atan2(-t(:,1), t(:,2)), prm.zeroAxis, prm.ccw));  % orient convention
    hue = mod(hueMap(1) * th / 180 + hueMap(2), 1);
    col = uint8(round(255 * hsv2rgb([hue, ones(nV, 2)])));

    okAx = isfinite(vt.coh) & vt.coh >= prm.minCoh;
    incl = okAx & phi >= prm.maxIncl;
    good = okAx & ~incl;
    rgb(good, :) = col(good, :);
    rgb(incl, :) = repmat(inclColor, nnz(incl), 1);
    rgb(~okAx, :) = repmat(lowColor, nnz(~okAx), 1);
    info.fracGood = mean(good); info.fracIncl = mean(incl); info.fracLow = mean(~okAx);
end

% ---- vessel axis (identical to compare_orient_vessel_*.m) ----

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

function th = display_to_conv(thDisp, zeroAxis, ccw)
    th0 = 0; if strcmp(zeroAxis, 'dim1'), th0 = -pi/2; end
    sgn = 1; if ~ccw, sgn = -1; end
    th = mod(sgn * (thDisp - th0), pi);
end

function [r1, r2, r3] = bbox_ranges(m, pad)
    sz = size(m);
    a1 = find(any(any(m, 2), 3)); a2 = find(any(any(m, 1), 3)); a3 = find(any(any(m, 1), 2));
    r1 = max(1, a1(1)-pad):min(sz(1), a1(end)+pad);
    r2 = max(1, a2(1)-pad):min(sz(2), a2(end)+pad);
    r3 = max(1, a3(1)-pad):min(sz(3), a3(end)+pad);
end