function T = measure_histo_shells(hist, stain_name, cfg)
% MEASURE_HISTO_SHELLS  Nearest-structure shell measurements on 2D histology.
%
% Same approach as measure_parenchyma_shells.m (OCT), in 2D:
%   - every parenchymal pixel is assigned to its single nearest structure
%     (EPVS, reported vessel, or other vessel) with a Euclidean distance
%     transform, then binned into shells (cfg.edges_um)
%   - one output row = one structure x one shell, with a stable object id
%
% Normalisation follows measure_epvs_and_vessel_variable.m:
%   histogram matching to the first section of this stain (within mask),
%   then z-scoring within the section.
%
% INPUTS:
%   hist       (struct array): from import_stain; fields baseName, image,
%                              epvs, ves, mask
%   stain_name (char/string):  e.g. 'LHE'
%   cfg        (struct):
%       pix_um                 pixel size (µm)
%       edges_um               shell edges (µm), e.g. 0:40:520
%       n_min_obj_px           min pixels per EPVS/vessel annotation (0 = keep all)
%       n_min_bin_px           min parenchymal pixels per object x shell (flag)
%       epvs_buffer_px         vessels within this many pixels of an EPVS are
%                              not reported as normal vessels (0 = overlap only)
%       hist_match             true = histogram-match to the first section
%       zscore_parenchyma_only true = z-score statistics from parenchyma only
%                              (mask minus EPVS/vessel pixels); false = whole mask
%                              (as in the previous ring method)
%       save_qc                write shell / owner-class TIFFs
%       qc_dir                 folder for QC TIFFs
% OUTPUT:
%   T (table): one row per object x shell (metadata added by the caller)

pix = cfg.pix_um;
nB  = numel(cfg.edges_um) - 1;
e   = cfg.edges_um(:);

ref_im   = im2single(single_ch(hist(1).image));
ref_mask = logical(single_ch(hist(1).mask));

if cfg.save_qc && ~exist(cfg.qc_dir, 'dir'), mkdir(cfg.qc_dir); end

allT = {};
for ii = 1:numel(hist)
    name = string(hist(ii).baseName);
    fprintf('  %s section %d of %d: %s\n', stain_name, ii, numel(hist), name);

    stain = im2single(single_ch(hist(ii).image));
    epvs  = logical(single_ch(hist(ii).epvs));
    ves   = logical(single_ch(hist(ii).ves));
    mask  = logical(single_ch(hist(ii).mask));
    assert(isequal(size(stain), size(epvs), size(ves), size(mask)), ...
           '%s: image and masks differ in size', name);

    %%% Annotations
    epvs = keep_large2(epvs, cfg.n_min_obj_px);
    ves  = keep_large2(ves,  cfg.n_min_obj_px);
    [ves_rep, n_exc] = exclude_touching(ves, epvs, cfg.epvs_buffer_px);
    other = ves & ~ves_rep & ~epvs;

    [Lv, nV] = bwlabel(ves_rep, 8);
    [Le, nE] = bwlabel(epvs,    8);
    fprintf('    EPVS: %d | vessels: %d reported, %d excluded near EPVS\n', nE, nV, n_exc);
    if nV + nE == 0
        warning('%s: no EPVS or vessel annotations, skipped', name);
        continue
    end
    nL = nV + nE + 1;
    L  = zeros(size(stain), 'uint32');
    L(other)   = nL;
    L(ves_rep) = uint32(Lv(ves_rep));
    L(epvs)    = uint32(Le(epvs)) + nV;     % EPVS takes precedence on overlap
    struct_mask = L > 0;
    parench     = mask & ~struct_mask;

    %%% Normalisation (as in measure_epvs_and_vessel_variable)
    stain_m = stain;
    if cfg.hist_match && ii > 1
        stain_m(mask) = imhistmatch(stain(mask), ref_im(ref_mask));
    end
    if cfg.zscore_parenchyma_only
        norm_px = parench;
    else
        norm_px = mask;
    end
    mu = mean(stain_m(norm_px), 'omitnan');
    sd = std(stain_m(norm_px),  'omitnan');
    z  = (stain_m - mu) ./ sd;

    %%% Nearest structure + distance (µm)
    [D, IDX] = bwdist(struct_mask);
    D     = D * pix;
    owner = L(IDX);
    clear IDX
    bin = discretize(D, cfg.edges_um);
    bin(isnan(bin)) = 0;
    ok  = parench & bin > 0 & ~isnan(z);

    sz  = [nL, nB];
    idx = [double(owner(ok)), bin(ok)];
    npx = accumarray(idx, 1, sz);
    m_raw = accumarray(idx, double(stain(ok)),   sz) ./ npx;
    m_hm  = accumarray(idx, double(stain_m(ok)), sz) ./ npx;
    m_z   = accumarray(idx, double(z(ok)),       sz) ./ npx;

    %%% Object descriptors
    area_px = accumarray(double(L(struct_mask)), 1, [nL 1]);
    din = zeros(nL, 1);
    if nV
        Dv = bwdist(~ves_rep);
        din(1:nV) = accumarray(Lv(ves_rep), double(Dv(ves_rep)), [nV 1], @max);
    end
    if nE
        De = bwdist(~epvs);
        din(nV + (1:nE)) = accumarray(Le(epvs), double(De(epvs)), [nE 1], @max);
    end

    %%% Long table, reported objects only
    nR = nV + nE;
    [obj, b] = ndgrid(1:nR, 1:nB);
    obj = obj(:); b = b(:);
    lin  = sub2ind(sz, obj, b);
    keep = npx(lin) > 0;
    obj = obj(keep); b = b(keep); lin = lin(keep);
    cls = repmat("vessel", numel(obj), 1);
    cls(obj > nV) = "EPVS";

    Ti = table(repmat(name, numel(obj), 1), obj, cls, e(b), e(b + 1), npx(lin), ...
               m_raw(lin), m_hm(lin), m_z(lin), ...
               area_px(obj) * pix^2, 2 * din(obj) * pix, npx(lin) >= cfg.n_min_bin_px, ...
        'VariableNames', {'section','obj','class','d_inner_um','d_outer_um','n_px', ...
                          'raw','hmatch','z','obj_area_um2','obj_diam_um','enough_pixels'});
    allT{end+1} = Ti; %#ok<AGROW>

    %%% QC images: shell index (1..nB) and owner class (1 vessel, 2 EPVS, 3 other)
    if cfg.save_qc
        oc = zeros(nL, 1, 'uint8');
        oc(1:nV) = 1; oc(nV + (1:nE)) = 2; oc(nL) = 3;
        shell_img = zeros(size(stain), 'uint8'); shell_img(ok) = bin(ok);
        class_img = zeros(size(stain), 'uint8'); class_img(ok) = oc(owner(ok));
        imwrite(shell_img, fullfile(cfg.qc_dir, name + "_shell_bin.tif"));
        imwrite(class_img, fullfile(cfg.qc_dir, name + "_owner_class.tif"));
    end
end

if isempty(allT)
    T = table();
else
    T = vertcat(allT{:});
end
end

%% ------------------------------------------------------------------------
function im = single_ch(im)
% Some Fiji exports have 3 identical channels: keep the first
if ~ismatrix(im), im = im(:, :, 1); end
end

function bw = keep_large2(bw, nmin)
if nmin <= 0 || ~any(bw(:)), return; end
bw = bwareaopen(bw, nmin, 8);
end

function [ves_keep, n_excluded] = exclude_touching(ves, epvs, buffer_px)
% Drop vessel components that overlap or lie within buffer_px of an EPVS
ves_keep = ves; n_excluded = 0;
if ~any(epvs(:)) || ~any(ves(:)), return; end
if buffer_px > 0
    w = 2 * buffer_px + 1;
    epvs = imdilate(epvs, true(w, w));
end
cc  = bwconncomp(ves, 8);
Lc  = labelmatrix(cc);
hit = unique(Lc(epvs & ves));
hit(hit == 0) = [];
ves_keep   = ves & ~ismember(Lc, hit);
n_excluded = numel(hit);
end