%% Measure the size-weighted proximity (SWP)
% 
clear; clc; close all;

%% Directories and Data
data_dir = '/projectnb/npbssmic/ns/CAA/';
swp_dir = '/projectnb/npbssmic/ns/CAA/swp/';
% SWP string base
epvs_base = 'swp_voxelwise_radius_500_exp_2_interpolated_heatmap';
ves_base = 'swp_voxelwise_ves_radius_500_exp_2_interpolated_heatmap';
% Directory to store csv
fig_out = '/projectnb/npbssmic/ns/CAA/figures/fig4_mus_swp_epvs/';
% flag for loading data
load_caa_flag = true;
load_swp_flag = true;

%% Array of subject ID, regions, depth (slice)
subjects = struct();
subjects(1).subject_name = 'caa26';
subjects(1).region = 'front';
subjects(2).subject_name = 'caa6';
subjects(2).region = 'front';
subjects(3).subject_name = 'caa17';
subjects(3).region = 'occip';
subjects(4).subject_name = 'caa22';
subjects(4).region = 'front';
% Set slices for each subject
subjects(1).slice = 10;
subjects(2).slice = 10;
subjects(3).slice = 113;
subjects(4).slice = 10;

%% Figure properties
% Voxel size (microns)
vox = 20;
% scalebar length in microns
scaleBarLength = 5000;


%% Import the CAA structs & SWP

if load_caa_flag
    % CAA 26 - control
    fprintf('Loading CAA6\n')
    caa26 = load(fullfile(data_dir,'/caa26/caa26.mat'));
    caa26 = caa26.caa26;
    fprintf('Finished Loading CAA6\n')
    % CAA 6 - mild
    fprintf('Loading CAA6\n')
    caa6 = load(fullfile(data_dir,'/caa6/caa6.mat'));
    caa6 = caa6.caa6;
    fprintf('Finished Loading CAA6\n')
    % CAA 17 - moderate
    fprintf('Loading CAA17\n')
    caa17 = load(fullfile(data_dir,'/caa17/occip/caa17.mat'));
    caa17 = caa17.caa17;
    fprintf('Finished Loading CAA17\n')
    % CAA 22 - severe
    fprintf('Loading CAA22\n')
    caa22 = load(fullfile(data_dir,'/caa22/caa22.mat'));
    caa22 = caa22.caa22;
    fprintf('Finished Loading CAA22\n')
    % Integrate into data structs
    op = struct();
    op.caa26 = caa26;
    op.caa6 = caa6;
    op.caa17 = caa17;
    op.caa22 = caa22;
end
%% Import swp for each subject
if load_swp_flag
    for ii = 1:length(subjects)
        fprintf('Loading SWP for subject %s\n',subjects(ii).subject_name)
        % Create full file name
        subid = subjects(ii).subject_name;
        reg = subjects(ii).region;
        fname_base = strcat(subid,'_',reg,'_');
        
        %%% EPVS
        epvs_fname = strcat(fname_base, epvs_base);
        epvs_fname = fullfile(swp_dir,subid,reg,epvs_fname);    
        % Import SWP
        epvs_swp = load(epvs_fname);
        epvs_swp = single(epvs_swp.interpolated_volume);
        % Add to subjects struct
        subjects(ii).epvs_swp = epvs_swp;

        %%% Vessel
        ves_fname = strcat(fname_base, ves_base);
        ves_fname = fullfile(swp_dir,subid,reg,ves_fname);    
        % Import SWP
        ves_swp = load(ves_fname);
        ves_swp = single(ves_swp.interpolated_volume);
        % Add to subjects struct
        subjects(ii).ves_swp = ves_swp;
    end
end

%% Add white matter mask, mus, vessels, EPVS to subjects struct
for ii = 1:length(subjects)
    % Retrieve subject name and region
    sub = subjects(ii).subject_name;
    reg = subjects(ii).region;
    % Add mask, vessel, epvs, mus
    subjects(ii).mask = op.(sub).(reg).mask;
    subjects(ii).mask_wm = op.(sub).(reg).mask_wm;
    subjects(ii).seg = op.(sub).(reg).seg;
    subjects(ii).epvs = op.(sub).(reg).epvs;
    subjects(ii).mus = op.(sub).(reg).mus;
    subjects(ii).ret = op.(sub).(reg).ret_full;
    subjects(ii).ori = op.(sub).(reg).orient_rgb;
end

%% Create SWP figures
%%% Find limits to normalize across the three subjects
% Initialize vectors for storing min/max
swp_min.ves = 1;
swp_min.epvs = 1;
swp_max.ves = 40;
swp_max.epvs = 130;

% Iterate over subjects
for ii = 4:length(subjects)
    %%% Extract slice from each matrix
    % Extract slice index
    slice_idx = subjects(ii).slice;
    % Extract WM mask, vessels, EPVS
    mask = subjects(ii).mask_wm(:, :, slice_idx);
    ves = subjects(ii).seg(:, :, slice_idx);
    epvs = subjects(ii).epvs(:, :, slice_idx);
    % Remove EPVS from vessel segmentation
    ves(epvs) = 0;

    % Retrieve Vessel SWP
    ves_heatmap = subjects(ii).ves_swp(:, :, slice_idx);
    epvs_heatmap = subjects(ii).epvs_swp(:, :, slice_idx);
    % Retrieve subject and region names
    subid = subjects(ii).subject_name;
    reg = subjects(ii).region;
    
    % Create heatmap for vessel SWP
    swp_heatmap(ves_heatmap, swp_min.ves, swp_max.ves, mask, ves,...
                scaleBarLength, vox,...
                fig_out, subid, reg, slice_idx, ves_base)

    % Create heatmap for EPVS SWP
    swp_heatmap(epvs_heatmap, swp_min.epvs, swp_max.epvs, mask, epvs,...
                scaleBarLength, vox,...
                fig_out, subid, reg, slice_idx, epvs_base)
end

%% Create zoomed SWP figures

% Initialize x,y dimensions for zoomed subset (voxels)
sz = 400;
% Initialize the zoomed box for SWP heatmaps
coords = struct();
coords.caa26 = [600, 760];
coords.caa6 = [1000, 350];
coords.caa17 = [930, 470];
coords.caa22 = [425, 340];
% Zoom inset scale bar (micron)
zoom_scalebar = 1000;

% Iterate over subjects
for ii = 4:length(subjects)
    %%% Extract slice from each matrix
    % Extract slice index
    slice_idx = subjects(ii).slice;
    % Extract WM mask, vessels, EPVS
    mask = subjects(ii).mask_wm(:, :, slice_idx);
    ves = subjects(ii).seg(:, :, slice_idx);
    epvs = subjects(ii).epvs(:, :, slice_idx);
    % Remove EPVS from vessel segmentation
    ves(epvs) = 0;

    %%% Retrieve Vessel SWP
    ves_heatmap = subjects(ii).ves_swp(:, :, slice_idx);
    epvs_heatmap = subjects(ii).epvs_swp(:, :, slice_idx);

    %%% Extract the zoomed subset coordinates
    xy = coords.(subjects(ii).subject_name);
    % Take subset of swp, mask, seg
    ves_heatmap = ves_heatmap(xy(1):xy(1)+sz, xy(2):xy(2)+sz);
    epvs_heatmap = epvs_heatmap(xy(1):xy(1)+sz, xy(2):xy(2)+sz);
    mask = mask(xy(1):xy(1)+sz, xy(2):xy(2)+sz);
    ves = ves(xy(1):xy(1)+sz, xy(2):xy(2)+sz);
    epvs = epvs(xy(1):xy(1)+sz, xy(2):xy(2)+sz);
    
    %%% Retrieve subject and region names
    subid = subjects(ii).subject_name;
    reg = subjects(ii).region;
    % Create heatmap for vessel SWP
    swp_heatmap(ves_heatmap, swp_min.ves, swp_max.ves, mask, ves,...
                zoom_scalebar, vox,...
                fig_out, subid, reg, slice_idx, strcat(ves_base,'_zoom'))
    % Create heatmap for EPVS SWP
    swp_heatmap(epvs_heatmap, swp_min.epvs, swp_max.epvs, mask, epvs,...
                zoom_scalebar, vox,...
                fig_out, subid, reg, slice_idx, strcat(epvs_base,'_zoom'))
end

%% Create grayscale mus subfigures
%{
% import the TIFs that had agarose manually removed
% CAA6f = caa6_front_depth_93_mus.tif
% CAA17o = caa17_occip_depth_536_mus.tif
% CAA22f = caa22_front_depth_347_mus.tif

mus = struct();

% import caa6f
t = Tiff(fullfile(fig_out,'caa6_front_depth_93_mus.tif'),'r');
mus.caa6f = read(t); close(t);

% import caa17o
t = Tiff(fullfile(fig_out,'caa17_occip_depth_536_mus.tif'),'r');
mus.caa17o = read(t); close(t);

% The CAA22f needs to be generated from the structs
mus_caa22f = caa22.front.mus;
mask_caa22f = caa22.front.mask;
mus_caa22f = mus_caa22f .* mask_caa22f;
mus_caa22f = mus_caa22f(:,:,347);
mus.caa22f = mus_caa22f;

%%% Iterate subjects and plot
subs = fields(mus);
for ii = 1:length(subs)
    % figure('position',[100 100 1500 1500]);
    figure('position',[100 100 1500 1500]);
    imagesc(mus.(subs{ii}));
    colormap('gray');
    clim([0, 24])
    h = colorbar;
    h.Ticks = [0 6 12 18 24];
    h.TickLabels = {'0','6','12','18','24'};
    xticks([]); yticks([]);
    set(gca, 'FontSize', 30);
    % Add scale bar
    hold on
    slice = mus.(subs{ii});
    scalebar_fun(scaleBarLength, vox, slice)
    hold off
    % Export image as PNG
    slice_idx = subjects(ii).slice;
    fname = strcat(subjects(ii).subject_name,'_',subjects(ii).region,'_',...
                    'depth_',num2str(slice_idx),'_mus.png');
    fout = fullfile(fig_out,fname);
    exportgraphics(gcf, fout,"Resolution",600)
    pause(1)
    % Export image as PDF
    fname = strcat(subjects(ii).subject_name,'_',subjects(ii).region,'_',...
                    'depth_',num2str(slice_idx),'_mus.pdf');
    fout = fullfile(fig_out,fname);
    exportgraphics(gcf, fout,"Resolution",600)
    pause(1)
    close;
end
%}

%% Orientation subfigures
%{
% Iterate over subjects
for i = 1:length(subjects)
    figure('position',[500 500 1500 1500]);
    %%% Extract slice from each matrix
    % Extract slice index
    slice_idx = subjects(i).slice;
    % Extract swp slice and mask slice
    slice = squeeze(subjects(i).ori(:, :, slice_idx,:));
    mask = subjects(i).mask(:, :, slice_idx);
    mask_wm = subjects(i).mask_wm(:, :, slice_idx);
    mask = mask + mask_wm;
    % Extract vessels and EPVS
    epvs = subjects(i).epvs(:, :, slice_idx);

    % Convert to RGB and set background to black
    mask = logical(mask);                % make sure logical
    mask_rgb = repmat(mask, 1, 1, 3);         % expand mask to RGB
    slice(~mask_rgb) = 0;                       % set background to black
    % Plot the heatmap of orientation
    imagesc(slice); hold on;
    
    % Plot tissue mask boundary
%     hold on;
%     boundaries = bwboundaries(mask);
%     for k = 1:numel(boundaries)
%         boundary = boundaries{k};
%         plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5);
%     end
%     pause(1)

    %%% Add scale bar to bottom right corner
    % scalebar length in pixels
    scaleBar_px = scaleBarLength / vox;
    % get image size
    [imHeight, imWidth] = size(slice,[1,2]);
    % Position: bottom right margin
    % small margin (2% of width)
    x_end = imWidth - round(imWidth*0.02);
    x_start = x_end - scaleBar_px;          
    % a little above bottom (3% of height)
    y_pos = imHeight - round(imHeight*0.03);
    % Draw scale bar (white line)
    plot([x_start x_end], [y_pos y_pos], 'w', 'LineWidth', 5);
    hold off;
    set(gca,'XTick',[]); set(gca,'YTick',[])

    %%% Save output with high quality
    fname = strcat(subjects(i).subject_name,'_',subjects(i).region,'_',...
        'depth_',num2str(slice_idx),'_ori.png');
    fout = fullfile(fig_out,fname);
    exportgraphics(gcf, fout,"Resolution",600)
    pause(1)
end
%}

%% Function scale bar
function scalebar_fun(sbar_len, vox, slice)
%%% Add scale bar to bottom right corner
% INPUTS
%   sbar_len (uint): scalebar length (micron)
%   vox (uint): isotropic voxel size
%   slice (mxn): enface image to determine dimensions

% scalebar length in pixels
sbar_px = sbar_len / vox;
% get image size
[imHeight, imWidth] = size(slice);
% Position: bottom right margin
% small margin (2% of width)
x_end = imWidth - round(imWidth*0.02);
x_start = x_end - sbar_px;          
% a little above bottom (3% of height)
y_pos = imHeight - round(imHeight*0.03);
% Draw scale bar (white line)
plot([x_start x_end], [y_pos y_pos], 'w', 'LineWidth', 5);
end

%% SWP Heatmap Function
function swp_heatmap(heatmap, swp_min, swp_max, mask, seg,...
    scaleBarLength, vox,...
    fig_out, subname, region, slice_idx, base_fname)
% INPUTS:
%   heatmap (float matrix): swp heatmap matrix
%   swp_min (uint): minimum value for heatmap
%   swp_max (uint): maximum value for heatmap
%   mask (logical): white matter tissue mask
%   seg (logical): vessel or EPVS segmentation (depends on iteration)
%   scaleBarLength (uint): scale bar length in microns
%   vox (uint): isotropic voxel size (microns)
%   fig_out (string): output directory
%   subname (string): subject name string
%   region (string): subject brain region string
%   slice_idx (uint): OCT depth for reference in filename
%   base_fname (string): filename base

%%% set colormap
cmap = parula(256);
% Clip the upper and lower to ease plotting
heatmap(heatmap < swp_min) = swp_min; % Clip
heatmap(heatmap > swp_max) = swp_max; % Clip
% Convert heatmap values to colormap indices
idx = round((heatmap - swp_min) / (swp_max - swp_min) * 255) + 1;
% Set any NaN to minimum value
idx(isnan(idx)) = 1;
idx(idx < 1) = 1;
idx(idx > 256) = 256;

% Create RGB image
figure('position',[100, 100, 1100, 1500]);
rgbImage = ind2rgb(idx, cmap);

%%% Overlay black where ~mask_wm OR segmentation==1
blackMask = (~mask) | (seg == 1);
for c = 1:3
    temp = rgbImage(:,:,c);
    temp(blackMask) = 0;
    rgbImage(:,:,c) = temp;
end
imshow(rgbImage,'InitialMagnification','fit');
% Adjust colorbar + tick marks and labels
colormap(cmap);
cb = colorbar;
tick_vals = linspace(0, 1, 10);
cb.Ticks = tick_vals;
cb.TickLabels = string(round(linspace(swp_min, swp_max, 10)));
cb.FontSize = 15;

%%% Purple boundary around tissue
hold on;
boundaries = bwboundaries(mask);
for k = 1:numel(boundaries)
    boundary = boundaries{k};
    plot(boundary(:,2), boundary(:,1),'m','LineWidth', 1.5);
end

%%% Add scale bar to bottom right corner
scalebar_fun(scaleBarLength, vox, heatmap)
hold off;
pause(1)

%%% Save output with high quality
% Export as PNG
fname = strcat(subname,'_',region,'_','depth_',num2str(slice_idx),'_',base_fname,'.png');
fout = fullfile(fig_out,fname);
exportgraphics(gcf, fout,"Resolution",600)
pause(1)
% Export as PDF
fname = strcat(subname,'_',region,'_','depth_',num2str(slice_idx),'_',base_fname,'.pdf');
fout = fullfile(fig_out,fname);
exportgraphics(gcf, fout,"Resolution",600)
pause(1)
close
end