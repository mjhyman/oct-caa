%% Analyze the Optical Properties
%{
The purpose of this script is to measure the optical properties in the
parenchymal tissue surrounding the EPVS and the vasculature. This covers
just the following cases:
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)

Outline:
- IMPORT struct containing:
    - tissue masks (exclude background signal)
    - EPVS
    - segmented vasculature (from Etienne)
    - optical properties (mus + retardance)
- Measure mean optical properties in white matter
    - exclude vessels and EPVS
%}

%% Prepare environment
clc; close all;
% Add top-level directory + subdirectories
addpath(genpath(fullfile(pwd, '..')))
% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/projectnb/npbssmic/ns/CAA/';
% flag for reloading the .MAT struct for each subject
flag_load_caa_structs = false;

%% Set the upper limits for mus and retardance
mus_max = 25;
ret_max = 45;

%% Load each subject's .MAT struct and create WM mask
if flag_load_caa_structs
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
    % CAA 25
    fprintf('Loading CAA25\n')
    caa25 = load(fullfile(data_dir,'/caa25/caa25.mat'));
    caa25 = caa25.caa25;
    fprintf('Finished Loading CAA25\n')
    % CAA 26
    fprintf('Loading CAA26\n')
    caa26 = load(fullfile(data_dir,'/caa26/caa26.mat'));
    caa26 = caa26.caa26;
    fprintf('Finished Loading CAA26\n')
end

% Load subjects into struct
subjects = struct();
subjects.caa6 = caa6;
subjects.caa17 = caa17;
subjects.caa22 = caa22;
subjects.caa25 = caa25;
subjects.caa26 = caa26;

%% Analyze mean white matter for each subject

% retrieve subject IDs
subs = fieldnames(subjects);
% Initialize struct
parench_median = struct();

% Loop over each subject
for i = 1:length(subs)
    % subject ID
    sub = subs{i};  
    % regions for this subject
    regions = fieldnames(subjects.(sub));
    % Loop over each region (e.g., 'front', 'occip')
    for j = 1:length(regions)
        % Print to console
        fprintf('Starting subject %s region %s.\n',sub,reg)
        % Region (front or occip)
        reg = regions{j};  
        % Load the local parameters dynamically
        [mus, ret, seg, seg_wm, ~, epvs] = ...
            load_local_params(subjects.(sub), reg);
        
        %%% Apply white matter mask
        % Retrieve white matter mask
        wm = subjects.(sub).(reg).mask_wm;
        % Set tissue outside WM mask to NaN
        mus(~wm) = NaN;
        ret(~wm) = NaN;

        % Exclude vessels from mus and retardance - ensure parenchyma
        % measurements are only of parenchyma
        mus(seg) = NaN;
        ret(seg) = NaN;
        if ~isempty(epvs)
            epvs = logical(epvs);
            mus(epvs) = NaN;
            ret(epvs) = NaN;
        end
        
        %%% Omit values above threshold then apply 1.5*IQR removal
        mus = omit_outlier(mus,mus_max);
        ret = omit_outlier(ret,ret_max);       

        %%% Median mus, retardance, orientation
        mus_median = median(mus, 'omitnan');
        ret_median = median(ret, 'omitnan');

        %%% Save to struct
        fprintf('Finished subject %s region %s.\n',sub,reg)
        % Medians
        parench_median.(sub).(reg).med.mus = mus_median;
        parench_median.(sub).(reg).med.ret = ret_median;
    end
end

% Clear subjects to reduce memory usage
clear subjects mus ret seg epvs

%% Save to folder

% Retrieve current date
d = string(datetime('now','TimeZone','local','Format','d-MMM-y'));
% Write filename
fout = string("median_white_matter_values_" + d + '.mat');
% Create output filepath
fout = fullfile(data_dir,fout);
% Save output
save(fout, 'parench_median');
% Print to console that job is done
fprintf('Finished saving median values to %s.\n', fout);