%% Take log of size-weighted proximity (SWP) + export TIF
% The purpose of this is to ensure all TIFs can be visualized on the same
% colorbar scale
clear; clc; close all;

%% Import SWP .mats
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
fcom = '_radius_200_exp_2_interpolated_heatmap.mat';
fout = '_radius_200_exp_2_interpolated_heatmap_log10.tif';
favi = '_radius_200_exp_2_interpolated_heatmap_log10.avi';
matout = '_radius_200_exp_2_interpolated_heatmap_log10.mat';

%% Take logarithm & measure min/max

for ii = 1:length(subjects)
    % Create filepath
    sub = subjects(ii).subject_name;
    reg = subjects(ii).region;
    fname = strcat(sub, '_', reg, fcom);
    fpath = fullfile(ddir,sub,reg);
    % Import the .MAT
    try
        swp = load(fullfile(fpath, fname));
        fprintf('\nRunning %s\n',strcat(sub, '_', reg))
    catch
        fprintf('\nSkipping %s\n',strcat(sub, '_', reg))
        continue
    end
    % Log transform
    swp = log10(swp.interpolated_volume+10);
    % Remove the log(10) value (1)
    swp(swp==1) = NaN;
    % Add to struct
    subjects(ii).min = min(swp(:));
    subjects(ii).max = max(swp(:));
    subjects(ii).swp = swp;
    % Measure the median (exclude zero values)
    subjects(ii).score = median(swp(:),'omitnan');
    % Create .MAT of the log(SWP)
    fname = strcat(sub, '_', reg, matout);
    matsave = fullfile(fpath,fname);
    save(matsave,'swp','-v7.3');
    % Add to spreadsheet
    
end

%% Export video of 3D stack
% Convert to .AVI

% Set global min and max for clim
global_min = 0;
global_max = max([subjects(:).max]);

%%% Iterate over subjects
for ii = 1:length(subjects)
    % Retrieve volume
    vol = subjects(ii).swp;
    if isempty(vol)
        continue
    end
    % Create output filepath
    sub = subjects(ii).subject_name;
    reg = subjects(ii).region;
    fname = strcat(sub, '_', reg, favi);
    fpath = fullfile(ddir,sub,reg,fname);
    % Create VideoWriter object
    vidObj = VideoWriter(fpath, 'Motion JPEG AVI');
    vidObj.FrameRate = 5;
    open(vidObj);

    
    %%% Create figure for this subject/region
    % [left, bottom, width, height]
    figure('Position', [100, -400, 1200, 904]);
    for z = 1:size(vol,3)
        imagesc(vol(:,:,z));
        axis image off;
        colormap('jet');
        colorbar;
        clim([global_min, global_max]);

        % Create Title String
        tstr = sprintf('%s %s - depth %d/%d',sub,reg,z,size(vol,3));
        title(tstr)

        % Draw the frame and save to video writer
        drawnow;
        frame = getframe(gcf);
        writeVideo(vidObj, frame);
    end
    close(vidObj);
    close(gcf);
    sprintf('\nFinished %s %s\n',sub, reg)
end