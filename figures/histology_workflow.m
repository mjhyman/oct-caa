%% Subfigures for histopathology workflow
%{
- large stain
- stain subset
- stain deconvolved
- histogram matched
- z-score
- z-score with donut annotations
%}

%% Top-level settings
clearvars -except caa6 caa17 caa22 caa25 caa26 swp_struct lhe cd68 gfap
clc; close all;
% Get the current folder
currentFolder = pwd;
% Move one directory up
parentFolder = fileparts(currentFolder);
% Add the parent folder and all its subfolders to the MATLAB search path
addpath(genpath(parentFolder));

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/figures/histology_workflow/';

%%% Load the stain and segmentation structs
if ~exist('lhe','var')
    lhe = load(fullfile(hdir,'lhe_rings_15-Oct-2025.mat')); lhe = lhe.lhe;
end

%%% measurement radius (units = pixels)
pix = 1.3067;               % um/pix
rad_sm = ceil(40/pix);      % 39.2 um ~= 30 pixels
rad_lg = ceil(100/pix);     % 39.2 um ~= 30 pixels
radii_sm = rad_sm : rad_sm : rad_sm*12;
radii_lg = rad_lg : rad_lg : rad_lg*5;

%%% color limits
% stain limits
smin = 0;
smax = 0.5;
% Reference limits
rmin = 0;
rmax = 0.6;
% z-score limits
zmin = -2;
zmax = 2;

%%% Scale bar properties
slen = 2000; % scale bar length in microns
slen_zoom = 100; % scale bar length (um) for zoomed subfigure
sth = 25;   % scale bar thickness in line width units

%%% Font size for color bar
fsize = 40; % Set font size for figures

%%% Define bounds for cropping
x1 = 2274; x2 = x1+300;
y1 = 1585; y2 = y1+300;

%% Debugging Figure
% stain
zstain = lhe(2).z_stain;
mask = lhe(2).mask;
ves =  lhe(2).ves;
epvs = lhe(2).epvs;
zstain(~mask) = -4;

% Flip the axes to place gyri at top
zstain = flipud(zstain');
ves = flipud(ves');
epvs = flipud(epvs');

%%% Dilations of EPVS and Vessel
% EPVS
se1 = strel('disk',0);
se2 = strel('disk',rad_sm);
inner_e = imdilate(epvs,se1);
outter_e = imdilate(epvs,se2);
epvs_annot = xor(inner_e, outter_e);
epvs_annot = epvs_annot(y1:y2,x1:x2);
epvs_annot = flipud(epvs_annot');
% Vessel
se1_v = strel('disk',0);
se2_v = strel('disk',rad_sm);
inner_v = imdilate(ves,se1_v);
outter_v = imdilate(ves,se2_v);
ves_annot = xor(inner_v, outter_v);
ves_annot = ves_annot(y1:y2,x1:x2);
ves_annot = flipud(ves_annot');

% 1. Create Vermillion overlay for vessel [0.89, 0.26, 0.20]
overlay_ves = cat(3, ones(size(zstain)) * 0.89, ... 
                      ones(size(zstain)) * 0.26, ... 
                      ones(size(zstain)) * 0.20);

% 2. Create Dark Blue overlay for EPVS [0.68, 0.85, 0.90]
% Dark Blue 0.0980    0.3922    0.6902
overlay_epvs = cat(3, ones(size(zstain)) * 0.098, ... 
                     ones(size(zstain)) * 0.3922, ... 
                     ones(size(zstain)) * 0.6902);

%%% Plot z-stain and then overlay donuts
figure('Position', [500 500 1000 1000]);
imagesc(zstain); hold on;

% Plot EPVS Overlay (Vermillion)
h1 = imagesc(overlay_epvs); colormap('gray');
set(h1, 'AlphaData', epvs_annot * 0.5); % 50% transparency [cite: 5]

% Plot Vessel Overlay (Light Blue)
h2 = imagesc(overlay_ves); colormap('gray');
set(h2, 'AlphaData', ves_annot * 0.5); % 50% transparency [cite: 5]

set(gca, 'XTick', []); set(gca, 'YTick', []);
c = colorbar; clim([zmin, zmax]);
c.Label.String = 'z-score';
set(gca, 'FontSize', fsize);

%% Create figure of stain deconvolved
% stain
stain = lhe(2).stain;
mask = lhe(2).mask;
stain(~mask) = -4;
% Apply crop and rotate
stain = stain(y1:y2,x1:x2);
stain = flipud(stain');
% Create figure
figure('Position', [500 500 1000 1000]);
imagesc(stain); colormap('gray'); clim([smin, smax]); c = colorbar;
c.Label.String = 'Stain Intensity (unitless)';
% Title, remove tick marks, colorbar
title('Stain')
set(gca,'XTick',[]); set(gca,'YTick',[])
set(gca,'FontSize',fsize)
% Add scale bar (microns)
sbar(slen_zoom, sth, pix, stain)
% Save figure
fout = fullfile(figdir, 'stain_deconv');
png_fout = strcat(fout,'.png');
svg_fout = strcat(fout,'.svg');
exportgraphics(gca,png_fout,'ContentType','vector','Resolution',600);
exportgraphics(gca,svg_fout,'ContentType','vector','Resolution',600);

%% Figure of reference
% The reference is CAA 17 occipital
ref = lhe(1).stain;
mask = lhe(1).mask;
ref(~mask) = -4;
% Create figure for the reference stain
figure('Position', [500 500 1000 1000]);
imagesc(ref);
% clim([lmin, lmax]);
title('Reference Stain - CAA 17 Occipital');
c = colorbar;
c.Label.String = 'Stain Intensity (unitless)';
clim([rmin,rmax]);
% Title, remove tick marks, colorbar
title('Reference')
set(gca,'XTick',[]); set(gca,'YTick',[])
set(gca,'FontSize',fsize)
% Add scale bar (microns)
sbar(slen, sth, pix, ref)
% Save figure
fout = fullfile(figdir, 'stain_ref');
saveas(gca,fout,'png');

%% Figure of histogram matched
% stain 
stain = lhe(2).stain_matched;
mask = lhe(2).mask;
stain(~mask) = -4;
% Apply crop and rotate
stain = stain(y1:y2,x1:x2);
stain = flipud(stain');
% Create figure
figure('Position', [500 500 1000 1000]);
imagesc(stain); colormap('gray'); clim([rmin, rmax]); c = colorbar;
c.Label.String = 'Stain Intensity (unitless)';
% Title, remove tick marks, colorbar
title('Stain Histogram Matched')
set(gca,'XTick',[]); set(gca,'YTick',[])
set(gca,'FontSize',fsize)
% Add scale bar (microns)
sbar(slen_zoom, sth, pix, stain)
% Save figure
fout = fullfile(figdir, 'histogram_matched');
png_fout = strcat(fout,'.png');
svg_fout = strcat(fout,'.svg');
exportgraphics(gca,png_fout,'ContentType','vector','Resolution',600);
exportgraphics(gca,svg_fout,'ContentType','vector','Resolution',600);

%% Figure of z-score
% z stain 
zstain = lhe(2).z_stain;
mask = lhe(2).mask;
zstain(~mask) = -4;
% Apply crop and rotate
zstain = zstain(y1:y2,x1:x2);
zstain = flipud(zstain');
% Create figure
figure('Position', [500 500 1000 1000]);
imagesc(zstain); colormap('gray'); clim([zmin, zmax]); c = colorbar;
c.Label.String = 'z-score';
% Title, remove tick marks, colorbar
title('Z-score Stain')
set(gca,'XTick',[]); set(gca,'YTick',[])
set(gca,'FontSize',fsize)
% Add scale bar (microns)
sbar(slen_zoom, sth, pix, zstain)
% Save figure
fout = fullfile(figdir, 'histogram_matched_stain_z_score');
png_fout = strcat(fout,'.png');
svg_fout = strcat(fout,'.svg');
exportgraphics(gca,png_fout,'ContentType','vector','Resolution',600);
exportgraphics(gca,svg_fout,'ContentType','vector','Resolution',600);

%% Create figure overlaying z-score with donuts

%%% Create EPVS and Vessel annotations
% EPVS
epvs = lhe(2).epvs;
se1 = strel('disk',0);
se2 = strel('disk',rad_sm);
inner_e = imdilate(epvs,se1);
outter_e = imdilate(epvs,se2);
epvs_annot = xor(inner_e, outter_e);
epvs_annot = epvs_annot(y1:y2,x1:x2);
epvs_annot = flipud(epvs_annot');
% Vessel
ves = lhe(2).ves;
se1_v = strel('disk',0);
se2_v = strel('disk',rad_sm);
inner_v = imdilate(ves,se1_v);
outter_v = imdilate(ves,se2_v);
ves_annot = xor(inner_v, outter_v);
ves_annot = ves_annot(y1:y2,x1:x2);
ves_annot = flipud(ves_annot');

% 1. Create Vermillion overlay for vessel [0.89, 0.26, 0.20]
overlay_ves = cat(3, ones(size(zstain)) * 0.89, ... 
                      ones(size(zstain)) * 0.26, ... 
                      ones(size(zstain)) * 0.20);

% 2. Create Dark Blue overlay for EPVS [0.68, 0.85, 0.90]
% Dark Blue 0.0980    0.3922    0.6902
overlay_epvs = cat(3, ones(size(zstain)) * 0.098, ... 
                     ones(size(zstain)) * 0.3922, ... 
                     ones(size(zstain)) * 0.6902);

%%% Plot z-stain and then overlay donuts
figure('Position', [500 500 1000 1000]);
imagesc(zstain); colormap('gray'); hold on;

% Plot EPVS Overlay (Vermillion)
h1 = imagesc(overlay_epvs);
set(h1, 'AlphaData', epvs_annot * 0.5); % 50% transparency [cite: 5]

% Plot Vessel Overlay (Light Blue)
h2 = imagesc(overlay_ves);
set(h2, 'AlphaData', ves_annot * 0.5); % 50% transparency [cite: 5]

set(gca, 'XTick', []); set(gca, 'YTick', []);
c = colorbar; clim([zmin, zmax]);
c.Label.String = 'z-score';
set(gca, 'FontSize', fsize);

%%% Take EPVS annotation and create overlay
% epvs = lhe(2).epvs;
% se1 = strel('disk',0);
% se2 = strel('disk',rad_sm);
% inner = imdilate(epvs,se1);
% outter = imdilate(epvs,se2);
% epvs_annot = xor(inner, outter);
% % Overlay EPVS onto z_stain
% overlay = cat(3, ones(size(zstain)), zeros(size(zstain)), zeros(size(zstain)));
% 
% %%% Plot z-stain and then overlay donuts
% % z-stain
% figure('Position', [500 500 1000 1000]);
% imagesc(zstain); hold on;
% % Overlay with annotation
% h = imagesc(overlay);
% set(h,'AlphaData',epvs_annot * 0.5);
% set(gca,'XTick',[]); set(gca,'YTick',[])
% c = colorbar; clim([zmin,zmax]);
% c.Label.String = 'z-score';
% set(gca,'FontSize',fsize)

%%% Add scale bar (microns)
sbar(slen_zoom, sth, pix, zstain)

%%% Save figure
fout = fullfile(figdir, 'histogram_matched_stain_z_score_overlay_donuts_zoom');
png_fout = strcat(fout,'.png');
svg_fout = strcat(fout,'.svg');
exportgraphics(gca,png_fout,'ContentType','vector','Resolution',600);
exportgraphics(gca,svg_fout,'ContentType','vector','Resolution',600);