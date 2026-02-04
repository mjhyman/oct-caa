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

%% IMPORT SWP and subject strut
for ii = 1:length(subjects)
    % Create filepath
    sub = subjects(ii).subject_name;
    reg = subjects(ii).region;
    fname = strcat(sub, '_', reg, fcom);
    fpath = fullfile(ddir,sub,reg);
    % Import the SWP .MAT
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
    subjects(ii).max = max(swp(:));
    subjects(ii).swp = swp;
    % Measure the median (exclude zero values)
    subjects(ii).score = median(swp(:),'omitnan');
    % Create .MAT of the log(SWP)
    fname = strcat(sub, '_', reg, matout);
    matsave = fullfile(fpath,fname);
    save(matsave,'swp','-v7.3');   
end