%% Create subfigures for Figure 1 (image analysis workflow)
%{
This figure uses CAA 22 frontal:
    - scattering coefficient (mus)
    - white matter mask
    - vessel donuts
    - EPVS donuts

To Do:
    - extract depth 3
    - fig.1a (mus enface)
    - fig.1b (WM mask enface)
    - fig.1c masked w/ yellow box
    - fig.1d zoomed w/ EPVS and lumen (not highlighted)
    - fig.1e/f vessel segmentation + donut
    - fig.1g/h EPVS segmentation + donut
%}

% TODO:

%{
- Crop TIF of the zoomed in region (grayscale, 4 masks)

- Transfer files from
/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/figures/
Fig1_flowchart

- Write code for overlaying
%}

%% Prepare environment
clc; close all;
% Add top-level directory + subdirectories
parentDir = fileparts(pwd);
fsDir = fullfile(parentDir, 'freesurfer');
cstatDir = fullfile(parentDir, 'CircStat2012a');
addpath(fsDir);
addpath(cstatDir);
addpath(parentDir);

% Voxel dimensions (microns) for all runs
res = [20,20,20]; % resolution in microns

%%% Directories on SCC w/ Matlab struct
ddir = '/projectnb/npbssmic/ns/CAA/figures/image_processing_flowchart/';

%% Load the TIFs for each figure
% Grayscale of the first zoomed inset
t = Tiff(fullfile(ddir,'caa22f_depth0_mus_inset1.tif'),'r');
inset = read(t);

%%% Vessel inset
% vessel mus inset
t = Tiff(fullfile(ddir,'caa22f_depth0_mus_ves_inset.tif'),'r');
ves = read(t);
% vessel lumen segmentation mask
t = Tiff(fullfile(ddir,'caa22f_depth0_ves_mask_inset.tif'),'r');
ves_seg = read(t); ves_seg = ves_seg > 0;
% vessel parenchyma donut
t = Tiff(fullfile(ddir,'caa22f_depth0_ves_donut_inset.tif'),'r');
ves_donut = read(t); ves_donut = ves_donut > 0;

%%% EPVS inset
% vessel mus inset
t = Tiff(fullfile(ddir,'caa22f_depth0_mus_epvs_inset.tif'),'r');
epvs = read(t);
% vessel lumen segmentation mask
t = Tiff(fullfile(ddir,'caa22f_depth0_epvs_mask_inset.tif'),'r');
epvs_seg = read(t); epvs_seg = epvs_seg > 0;
% vessel parenchyma donut
t = Tiff(fullfile(ddir,'caa22f_depth0_epvs_donut_inset.tif'),'r');
epvs_donut = read(t); epvs_donut = epvs_donut > 0;

%% First inset figure

% Scale bar settings
len = 200; % microns
% scale bar thickness
th = 10;

% Vessel segmentation
f = figure; imshow(inset);
sbar(len, th, res(1), inset);
fout = fullfile(ddir, 'inset1.png');
exportgraphics(f,fout,'Resolution', 600); pause(0.5);

%% Overlays

% Transparency
alpha = 0.5;

% COLORBLIND FRIENDLY
cseg = [240, 228, 66]; % yellow
cdon = [240, 0, 128];  % Sky blue

% Overlays
ves_seg_overlay = mask_overlay(ves, ves_seg, cseg, alpha);
ves_don_overlay = mask_overlay(ves, ves_donut, cdon, alpha);
epvs_seg_overlay = mask_overlay(epvs, epvs_seg, cseg, alpha);
epvs_don_overlay = mask_overlay(epvs, epvs_donut, cdon, alpha);

%%% Show figures + add scale bar
% Scale bar settings
len = 200; % microns
% scale bar thickness
th = 4;

% Vessel segmentation
f = figure; imshow(ves_seg_overlay);
sbar(len, th, res(1), ves);
fout = fullfile(ddir, 'ves_lumen_overlay.png');
exportgraphics(f,fout,'Resolution', 600); pause(0.5);

% vessel donut
f = figure; imshow(ves_don_overlay);
sbar(len, th, res(1), ves);
fout = fullfile(ddir, 'ves_donut_overlay.png');
exportgraphics(f,fout,'Resolution', 600); pause(0.5);

% EPVS segmentation
f = figure; imshow(epvs_seg_overlay);
sbar(len, th, res(1), ves);
fout = fullfile(ddir, 'epvs_lumen_overlay.png');
exportgraphics(f,fout,'Resolution', 600); pause(0.5);

%% EPVS donut
f = figure; imshow(epvs_don_overlay);
sbar(len, th, res(1), ves);
fout = fullfile(ddir, 'epvs_donut_overlay.png');
exportgraphics(f,fout,'Resolution', 600); pause(0.5);

%% scale bar function
function sbar(sbar_len, sbar_thick, vox, im)
% Add scale bar to bottom right corner
% INPUTS:
%   - sbar_len (uint): scale bar length (microns)
%   - sbar_thick (uint): scale bar line width (builtin linewidth)
%   - vox (vector): isotropic voxel size (microns)
%   - slice (float matrix): single depth of image to display

%%% Define size
% scalebar length in pixels
sbar_px = sbar_len ./ vox;
% get image size
[imHeight, imWidth] = size(im,[1,2]);

%%%% Position: bottom right margin
% small margin (2% of width)
x_end = imWidth - round(imWidth*0.02);
x_start = x_end - sbar_px;          
% a little above bottom (3% of height)
y_pos = imHeight - round(imHeight*0.03);

%%% Draw scale bar (white line)
hold on;
plot([x_start x_end], [y_pos y_pos], 'w', 'LineWidth', sbar_thick);
hold off;
% Disable x,y ticks
set(gca,'XTick',[]); set(gca,'YTick',[])

end


function rgb_overlay = mask_overlay(grayImg, mask, overlay_color, alpha)
    % grayImg: 8-bit grayscale image (MxN)
    % mask: binary mask (MxN)
    % overlay_color: [R G B] in [0 255]
    % alpha: transparency (0 = original, 1 = full color)
    
    grayImg = uint8(grayImg); % Ensure format
    mask = mask > 0; % Logical
    
    % Start with grayscale converted to RGB
    rgb_overlay = double(repmat(grayImg,1,1,3));
    
    for c = 1:3
        rgb_overlay(:,:,c) = rgb_overlay(:,:,c) .* (~mask) + ...
                             (1-alpha) * rgb_overlay(:,:,c) .* mask + ...
                             alpha * overlay_color(c) * mask;
    end
    rgb_overlay = uint8(rgb_overlay);
end



