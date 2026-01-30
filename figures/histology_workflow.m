%% Import downsampled stains
%{
These are the stains that our collaborator performed. The grayscale images
can be directly used to compute the stain density surrounding the EPVS
because there is no counterstain. The deconvolved images had a
counterstain and therefore had to be deconvolved to extract the target
stain.

The purpose of this script is to create the images for the histology
workflow

%}

%% Top-level settings
clearvars -except caa6 caa17 caa22 caa25 caa26 swp_struct lhe
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
sth = 25;   % scale bar thickness in line width units

%%% Font size for color bar
fsize = 40; % Set font size for figures

%% Create figure of stain deconvolved
% stain
stain = lhe(2).stain;
mask = lhe(2).mask;
stain(~mask) = -4;
% Create figure
figure('Position', [100 100 1500 1000]);
imagesc(stain); clim([smin, smax]); c = colorbar;
c.Label.String = 'Stain Density (a.u.)';
% Title, remove tick marks, colorbar
title('Stain')
set(gca,'XTick',[]); set(gca,'YTick',[])
set(gca,'FontSize',fsize)
% Add scale bar (microns)
sbar(slen, sth, pix, stain)
% Save figure
fout = fullfile(figdir, 'stain_deconv');
saveas(gca,fout,'png');

%% Figure of reference
% The reference is CAA 17 occipital
ref = lhe(1).stain;
mask = lhe(1).mask;
ref(~mask) = -4;
% Create figure for the reference stain
figure('Position', [100 100 1500 1000]);
imagesc(ref);
% clim([lmin, lmax]);
title('Reference Stain - CAA 17 Occipital');
c = colorbar;
c.Label.String = 'Stain Density (a.u.)';
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
% Create figure
figure('Position', [100 100 1500 1000]);
imagesc(stain); clim([rmin, rmax]); c = colorbar;
c.Label.String = 'Stain Density (a.u.)';
% Title, remove tick marks, colorbar
title('Stain Histogram Matched')
set(gca,'XTick',[]); set(gca,'YTick',[])
set(gca,'FontSize',fsize)
% Add scale bar (microns)
sbar(slen, sth, pix, stain)
% Save figure
fout = fullfile(figdir, 'histogram_matched');
saveas(gca,fout,'png');

%% Figure of z-score
% z stain 
zstain = lhe(2).z_stain;
mask = lhe(2).mask;
zstain(~mask) = -4;
% Create figure
figure('Position', [100 100 1500 1000]);
imagesc(zstain); clim([zmin, zmax]); c = colorbar;
c.Label.String = 'z-score';
% Title, remove tick marks, colorbar
title('Z-score Stain')
set(gca,'XTick',[]); set(gca,'YTick',[])
set(gca,'FontSize',fsize)
% Add scale bar (microns)
sbar(slen, sth, pix, zstain)
% Save figure
fout = fullfile(figdir, 'histogram_matched_stain_z_score');
saveas(gca,fout,'png');

%% Create figure overlaying z-score with donuts

%%% Take EPVS annotation and create overlay
epvs = lhe(2).epvs;
se1 = strel('disk',0);
se2 = strel('disk',rad_sm);
inner = imdilate(epvs,se1);
outter = imdilate(epvs,se2);
annot = xor(inner, outter);
% Overlay EPVS onto z_stain
overlay = cat(3, ones(size(zstain)), zeros(size(zstain)), zeros(size(zstain)));

%%% Plot z-stain and then overlay donuts
% z-stain
figure('Position', [100 100 1500 1000]);
imagesc(zstain); hold on;
% Overlay with annotation
h = imagesc(overlay);
set(h,'AlphaData',annot * 0.5);
set(gca,'XTick',[]); set(gca,'YTick',[])
c = colorbar; clim([zmin,zmax]);
c.Label.String = 'z-score';
set(gca,'FontSize',fsize)

%%% Add scale bar (microns)
sbar(slen, sth, pix, zstain)

%%% Save figure
fout = fullfile(figdir, 'histogram_matched_stain_z_score_overlay_donuts');
saveas(gca,fout,'png');

%% Overlay z-score w/ donuts (cropped)

%%% Take EPVS annotation and create overlay
epvs = lhe(2).epvs;
se1 = strel('disk',0);
se2 = strel('disk',rad_sm);
inner = imdilate(epvs,se1);
outter = imdilate(epvs,se2);
annot = xor(inner, outter);
% Overlay EPVS onto z_stain
overlay = cat(3, ones(size(zstain)), zeros(size(zstain)), zeros(size(zstain)));

%%% Crop image to desired ROI
% Define bounds
x1 = 2500; x2 = 2900;
y1 = 5950; y2 = 6450;
% Crop zstain and overlay
zstainCropped = zstain(y1:y2, x1:x2);
overlayCropped = overlay(y1:y2, x1:x2, :);
annot = annot(y1:y2, x1:x2, :);

%%% Plot z-stain and then overlay
% z-stain
figure('Position', [100 100 1500 1000]);
imagesc(zstainCropped); hold on;
% Overlay with annotation
h = imagesc(overlayCropped);
set(h,'AlphaData',annot * 0.5);
set(gca,'XTick',[]); set(gca,'YTick',[])
c = colorbar; clim([zmin,zmax]);
c.Label.String = 'z-score';
set(gca,'FontSize',fsize)

%%% Add scale bar
sbar(100, 20, pix, zstainCropped)

%%% Save figure
fout = fullfile(figdir, 'histogram_matched_stain_z_score_overlay_donuts_zoom');
saveas(gca,fout,'png');

