%% Subfigures for histopathology workflow (LHE and GFAP)
%{
Exported per stain (PNG + SVG), prefixed with stain name:
- <stain>_stain_deconv
- <stain>_histogram_matched
- <stain>_histogram_matched_stain_z_score
- <stain>_histogram_matched_stain_z_score_overlay_donuts_zoom
Also exported (LHE only, unchanged):
- stain_ref
%}

%% Top-level settings
clearvars -except caa swp_struct lhe cd68 gfap
clc; close all;
currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(genpath(parentFolder));

%%% Directories (SCC)
hdir   = '/projectnb/npbssmic/ns/CAA/histology/';
figdir = '/projectnb/npbssmic/ns/CAA/figures/histology_workflow/';

%%% Load the stain and segmentation structs (only if not in workspace)
if ~exist('lhe','var')
    lhe = load(fullfile(hdir,'lhe_rings_15-Oct-2025.mat'));
    lhe = lhe.lhe;
end
if ~exist('gfap','var')
    gfap = load(fullfile(hdir,'gfap_rings_21-Nov-2025.mat'));
    gfap = gfap.gfap;
end

%%% Shared figure properties
slen  = 2000;   % scale bar length (um), full-slide reference
sth   = 25;     % scale bar thickness (line width units)
fsize = 40;     % font size

%%% Overlay colours
col_ves  = [0.89,   0.26,   0.20  ];  % vermillion (vessel)
col_epvs = [0.098,  0.3922, 0.6902];  % dark blue  (EPVS)
alpha    = 0.5;

%% Per-stain configuration
% Crop bounds are in the native (un-rotated) image coordinates.
% transpose = true  -> crop, then flipud(img')   (LHE convention)
% transpose = false -> crop, then flipud(img)    (GFAP convention)
cfg = struct( ...
    'name',      {'lhe',            'gfap'}, ...
    'idx',       {2,                7}, ...
    'x1',        {2274,             1850}, ...
    'x2',        {2274+300,         2600}, ...
    'y1',        {1585,             8980}, ...
    'y2',        {1585+300,         9720}, ...
    'transpose', {true,             false}, ...
    'pix',       {1.3067,           1.3067}, ... % um/pix
    'slen_zoom', {100,              250}, ...    % scale bar (um), zoomed panels
    'smin',      {0,                0}, ...      % deconvolved stain limits
    'smax',      {0.5,              0.5}, ...
    'rmin',      {0,                0}, ...      % histogram-matched limits
    'rmax',      {1.0,              1.0}, ...
    'zmin',      {-2,               -2}, ...     % z-score limits
    'zmax',      {2,                2});
stains = {lhe, gfap};

%% Generate and export figures for each stain
for k = 2:numel(cfg)
    c    = cfg(k);
    s    = stains{k}(c.idx);
    mask = s.mask;
    tag  = upper(c.name);
    rad_sm = ceil(40/c.pix);   % donut width: 40 um

    %%% Stain deconvolved
    img = s.stain; img(~mask) = -4;
    img = crop_orient(img, c);
    plot_panel(img, [c.smin c.smax], 'Stain Intensity (unitless)', ...
        sprintf('%s Stain', tag), fsize);
    sbar(c.slen_zoom, sth, c.pix, img)
    export_png_svg(gca, fullfile(figdir, [c.name '_stain_deconv']));

    %%% Histogram matched
    img = s.stain_matched; img(~mask) = -4;
    img = crop_orient(img, c);
    plot_panel(img, [c.rmin c.rmax], 'Stain Intensity (unitless)', ...
        sprintf('%s Stain Histogram Matched', tag), fsize);
    sbar(c.slen_zoom, sth, c.pix, img)
    export_png_svg(gca, fullfile(figdir, [c.name '_histogram_matched']));

    %%% z-score
    zimg = s.z_stain; zimg(~mask) = -4;
    zimg = crop_orient(zimg, c);
    plot_panel(zimg, [c.zmin c.zmax], 'z-score', ...
        sprintf('%s Z-score Stain', tag), fsize);
    sbar(c.slen_zoom, sth, c.pix, zimg)
    export_png_svg(gca, fullfile(figdir, [c.name '_histogram_matched_stain_z_score']));

    %%% z-score with donut overlays
    % Dilate on the full image, then crop, to avoid edge truncation
    epvs_annot = crop_orient(donut(s.epvs, rad_sm), c);
    ves_annot  = crop_orient(donut(s.ves,  rad_sm), c);

    plot_panel(zimg, [c.zmin c.zmax], 'z-score', '', fsize);
    hold on;
    add_overlay(col_epvs, epvs_annot, alpha);
    add_overlay(col_ves,  ves_annot,  alpha);
    hold off;
    sbar(c.slen_zoom, sth, c.pix, zimg)
    export_png_svg(gca, fullfile(figdir, ...
        [c.name '_histogram_matched_stain_z_score_overlay_donuts_zoom']));
end

%% Reference (CAA 17 occipital, LHE)
ref = lhe(1).stain;
ref(~lhe(1).mask) = -4;
plot_panel(ref, [cfg(1).rmin cfg(1).rmax], 'Stain Intensity (unitless)', ...
    'Reference', fsize);
sbar(slen, sth, cfg(1).pix, ref)
saveas(gca, fullfile(figdir, 'stain_ref'), 'png');

%% Local functions
function img = crop_orient(img, c)
% Crop in native coordinates, then rotate so gyri are at top
img = img(c.y1:c.y2, c.x1:c.x2);
if c.transpose
    img = flipud(img');
else
    img = flipud(img);
end
end

function ring = donut(bw, r)
% Annulus of width r (pixels) outside the binary object
bw   = logical(bw);
ring = xor(bw, imdilate(bw, strel('disk', r)));
end

function plot_panel(img, lims, clabel, ttl, fsize)
figure('Position', [500 500 1000 1000]);
imagesc(img); colormap('gray'); clim(lims);
cb = colorbar('southoutside'); cb.Label.String = clabel;
if ~isempty(ttl), title(ttl); end
set(gca, 'XTick', [], 'YTick', [], 'FontSize', fsize);
end

function add_overlay(rgb, annot, alpha)
ov = repmat(reshape(rgb, 1, 1, 3), [size(annot) 1]);
h  = image(ov);
set(h, 'AlphaData', double(annot) * alpha);
end

function export_png_svg(ax, fout)
exportgraphics(ax, [fout '.png'], 'Resolution', 600);
exportgraphics(ax, [fout '.svg'], 'ContentType', 'vector');
end