%% Measure the size-weighted proximity (SWP)
% 
clear; clc; close all;

%% Directories and Data
data_dir = '/projectnb/npbssmic/ns/CAA/';
swp_dir = '/projectnb/npbssmic/ns/CAA/swp/';
% SWP string base
str_base = 'radius_200_exp_2_interpolated_heatmap_log10.mat';
% Directory to store csv
fig_out = '/projectnb/npbssmic/ns/CAA/figures/fig4_mus_swp_epvs/';
% flag for loading data
load_struct_flag = false;
load_swp_flag = true;

%% Array of subject ID and regions
subjects = struct();
subjects(1).subject_name = 'caa6';
subjects(1).region = 'front';
subjects(2).subject_name = 'caa17';
subjects(2).region = 'occip';
subjects(3).subject_name = 'caa22';
subjects(3).region = 'front';

%% Figure properties
% Voxel size (microns)
vox = 20;
% scalebar length in microns
scaleBarLength = 5000;
% Set slices for each subject
subjects(1).slice = 93;
subjects(2).slice = 536;
subjects(3).slice = 347;

%% Import the CAA structs & SWP

if load_flag
    % CAA 6
    fprintf('Loading CAA6\n')
    caa6 = load(fullfile(data_dir,'/caa6/caa6.mat'));
    caa6 = caa6.caa6;
    fprintf('Finished Loading CAA6\n')
    % CAA 17
    fprintf('Loading CAA17\n')
    caa17 = load(fullfile(data_dir,'/caa17/occip/caa17.mat'));
    caa17 = caa17.caa17;
    fprintf('Finished Loading CAA17\n')
    % CAA 22
    fprintf('Loading CAA22\n')
    caa22 = load(fullfile(data_dir,'/caa22/caa22.mat'));
    caa22 = caa22.caa22;
    fprintf('Finished Loading CAA22\n')
end
if load_swp_flag
    %%% Import swp for each subject
    for ii = 1:length(subjects)
        fprintf('Loading subject %s\n',subjects(ii).subject_name)
        % Create full file name
        subid = subjects(ii).subject_name;
        reg = subjects(ii).region;
        fname = strcat(subid,'_',reg,'_');
        fname = strcat(fname, str_base);
        fname = fullfile(swp_dir,subid,reg,fname);    
        % Import SWP
        swp = load(fname);
        swp = single(swp.swp);
        % Add to subjects struct
        subjects(ii).swp = swp;
    end
end
% Integrate into data structs
op = struct();
op.caa6 = caa6;
op.caa17 = caa17;
op.caa22 = caa22;

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
swp_min = 1;
swp_max = 5;
% Iterate over subjects under review
for ii = 1:length(subjects)
    % swp of current iteration
    swp = subjects(ii).swp;
    % Find max
    swp_max = max([swp_max,max(swp(:))]);
end
% Round up to nearest integer
swp_max = ceil(swp_max);

% set colormap
cmap = parula(256);

% Iterate over subjects
for i = 1:length(subjects)
    figure('position',[100 100 1500 1500]);
    %%% Extract slice from each matrix
    % Extract slice index
    slice_idx = subjects(i).slice;
    % Extract swp slice and mask slice
    slice = subjects(i).swp(:, :, slice_idx);
    mask = subjects(i).mask_wm(:, :, slice_idx);
    % Extract vessels and EPVS
    vessel = subjects(i).seg(:, :, slice_idx);
    epvs = subjects(i).epvs(:, :, slice_idx);

    % Map heatmap values to parula colormap
    heatmap = slice;  % Your SWP slice
    heatmap(heatmap < swp_min) = swp_min; % Clip
    heatmap(heatmap > swp_max) = swp_max; % Clip
    
    % Convert heatmap values to colormap indices
    idx = round((heatmap - swp_min) / (swp_max - swp_min) * 255) + 1;
    % If there are NaNs (shouldn't be unless you set them)
    idx(isnan(idx)) = 1;
    idx(idx < 1) = 1;
    idx(idx > 256) = 256;
    
    % Create RGB image
    rgbImage = ind2rgb(idx, cmap);
    
    % Overlay black where ~mask_wm OR epvs==1
    blackMask = (~mask) | (epvs == 1);
    for c = 1:3
        temp = rgbImage(:,:,c);
        temp(blackMask) = 0;
        rgbImage(:,:,c) = temp;
    end
    imshow(rgbImage);
    % Adjust colorbar
    colormap(cmap);
    cb = colorbar;
    cb.Ticks = 0 : 1/(swp_max-swp_min) : 1;
    cb.TickLabels = string(1:10);
    cb.FontSize = 15;
    
    % Plot WM boundary
    hold on;
    boundaries = bwboundaries(mask);
    for k = 1:numel(boundaries)
        boundary = boundaries{k};
        plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 1.5);
    end
    pause(1)

    %%% Add scale bar to bottom right corner
    scalebar_fun(scaleBarLength, vox, slice)
    hold off;

    %%% Save output with high quality
    % Export as PNG
    fname = strcat(subjects(i).subject_name,'_',subjects(i).region,'_',...
        'depth_',num2str(slice_idx),'_swp.png');
    fout = fullfile(fig_out,fname);
    exportgraphics(gcf, fout,"Resolution",600)
    pause(1)
    % Export as PDF
    fname = strcat(subjects(i).subject_name,'_',subjects(i).region,'_',...
        'depth_',num2str(slice_idx),'_swp.pdf');
    fout = fullfile(fig_out,fname);
    exportgraphics(gcf, fout,"Resolution",600)
    pause(1)
    close
end

%% Create grayscale mus subfigures
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