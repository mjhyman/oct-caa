%% Create .AVI video of logarithm of size-weighted proximity (SWP)
% This will be used for supplemental figures
clear; clc; close all;

%% Import SWP .mats
data_dir = '/projectnb/npbssmic/ns/CAA/';
ddir = '/projectnb/npbssmic/ns/CAA/swp/';
subjects = struct();
subjects(1).subject_name = 'caa6';
subjects(1).region = 'front';
subjects(2).subject_name = 'caa6';
subjects(2).region = 'occip';
subjects(3).subject_name = 'caa17';
subjects(3).region = 'occip';
subjects(4).subject_name = 'caa22';
subjects(4).region = 'front';
subjects(5).subject_name = 'caa22';
subjects(5).region = 'occip';
subjects(6).subject_name = 'caa25';
subjects(6).region = 'front';
subjects(7).subject_name = 'caa25';
subjects(7).region = 'occip';
subjects(8).subject_name = 'caa26';
subjects(8).region = 'front';
subjects(9).subject_name = 'caa26';
subjects(9).region = 'occip';

% Common string within filename
fin = '_radius_200_exp_2_interpolated_heatmap_log10.mat';
favi = '_radius_200_exp_2_interpolated_heatmap_log10.avi';

% Flag for importing subject-level structs
import_subjects = true;

%% Figure properties
% Voxel size (microns)
vox = 20;
% scalebar length in microns
scaleBarLength = 5000;

%% IMPORT SWP and subject strut

for ii = 1:length(subjects)
    % Create filepath
    sub = subjects(ii).subject_name;
    reg = subjects(ii).region;
    fname = strcat(sub, '_', reg, fin);
    fpath = fullfile(ddir,sub,reg);
    % Import the SWP .MAT
    try
        swp = load(fullfile(fpath, fname));
        swp = swp.swp;
        fprintf('\nImporting %s\n',strcat(sub, '_', reg))
    catch
        fprintf('\nImporting %s\n',strcat(sub, '_', reg))
        continue
    end
    % Add max to struct
    subjects(ii).max = max(swp(:));
    subjects(ii).swp = swp;
    % Measure the median (exclude zero values)
    subjects(ii).score = median(swp(:),'omitnan');
end

%% Import the subject data structs (to extract the tissue masks)

if import_subjects
    % Struct for top-level struct
    subject_files = struct( ...
        'caa6', '/caa6/caa6.mat', ...
        'caa17', '/caa17/occip/caa17.mat', ...
        'caa22', '/caa22/caa22.mat', ...
        'caa25', '/caa25/caa25.mat', ...
        'caa26', '/caa26/caa26.mat');
    
    % Import the structs
    subs = fields(subject_files);
    for ii = 1:length(subs)
        % import struct
        fprintf('\nImporting subject %s\n',subs{ii})
        load(fullfile(data_dir,subject_files.(subs{ii})));
    end

    % Combine subjects into struct
    oct = struct();
    oct.caa6 = caa6;
    oct.caa17 = caa17;
    oct.caa22 = caa22;
    oct.caa25 = caa25;
    oct.caa26 = caa26;
    
    % Clear individual structs
    clear caa6 caa17 caa22 caa25 caa26
end


%% Add the WM mask, vessels, and EPVS to subject struct

% Initialize masks for WM, vessels, and EPVS
for ii = 1:length(subjects)    
    % Create strings to load oct data
    sid = subjects(ii).subject_name;
    reg = subjects(ii).region;
    fprintf('\nCombining WM mask, vessels, EPVS for %s, %s',sid,reg)

    % Load WM mask, vessels, EPVS from oct struct
    subjects(ii).mask_wm = oct.(sid).(reg).mask_wm;
    subjects(ii).epvs = oct.(sid).(reg).epvs;

    % Clear OCT from memory to reduce overhead
    clear oct
end
fprintf('\nDone combining\n')

%% Export video of 3D stack
% Convert to .AVI

% Set global min and max for clim
swp_min = 1;
swp_max = 10;
% set colormap
cmap = parula(256);

%%% Iterate over subjects
for ii = 1:length(subjects)
    % Retrieve volume
    vol = subjects(ii).swp;
    % Create output filepath
    sub = subjects(ii).subject_name;
    reg = subjects(ii).region;
    fname = strcat(sub, '_', reg, favi);
    fpath = fullfile(ddir,sub,reg,fname);
    % Create VideoWriter object
    vidObj = VideoWriter(fpath, 'MPEG-4');
    vidObj.FrameRate = 30;
    open(vidObj);    
    % Create figure for this subject/region
    % [left, bottom, width, height]
    figure('Position', [100, -400, 2400, 1808]);
    % Set background to black
    axes('Color','k');

    % Extract swp and mask
    swp = subjects(ii).swp;
    mask = subjects(ii).mask_wm;
    
    %%% Create the logical mask for black
    % Extract vessels and EPVS
    epvs = subjects(ii).epvs;
    blackMask = (~mask) | (epvs == 1);

    % Map heatmap values to parula colormap
    heatmap = swp;
    heatmap(heatmap < swp_min) = swp_min;
    heatmap(heatmap > swp_max) = swp_max;

    % Convert heatmap values to colormap indices
    idx = round((heatmap - swp_min) / (swp_max - swp_min) * 255) + 1;
    % Set the NaN to 1
    idx(isnan(idx)) = 1;
    idx(idx < 1) = 1;
    idx(idx > 256) = 256;

    % Extract a slice for computing scale bar length
    slice = swp(:,:,1);

    %%% Iterate over each depth and color plot w/ imagesc
    for z = 1:size(vol,3)
        % Create RGB image
        rgbImage = ind2rgb(idx(:,:,z), cmap);
        
        % Outside of tissue boundary, set image black
        for c = 1:3
            temp = rgbImage(:,:,c);
            temp(blackMask(:,:,z)) = 0;
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
        boundaries = bwboundaries(mask(:,:,z));
        for k = 1:numel(boundaries)
            boundary = boundaries{k};
            plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 1.5);
        end

        % Add scale bar to bottom right corner
        scalebar_fun(scaleBarLength, vox, slice)
        hold off;

        % Create Title String
        tstr = sprintf('%s %s - depth %d/%d',sub,reg,z,size(vol,3));
        title(tstr)

        % Save frame to video writer
        frame = getframe(gcf);
        writeVideo(vidObj, frame);
        fprintf('\n\tFinished writing frame %i/%i',z,size(vol,3));
    end
    close(vidObj);
    close(gcf);
    sprintf('\nFinished %s %s\n',sub, reg)
end


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