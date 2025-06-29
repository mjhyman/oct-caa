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
data_dir = ['/autofs/cluster/octdata3/users/mjhyman/' ...
    'oct_caa_analyses/optical_properties'];
% WM mask directory (from Taylor)
wm_dir = ['/autofs/space/turtle_001/users/xz875/projects/' ...
          'Multi-resolution_Unets_Semi_OCT/prediction'];
% Voxel dimensions (microns) for all runs
res = [20,20,20]; % resolution in microns

%%% flag for reloading the .MAT struct for each subject
flag_load_caa_structs = false;



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
parench = struct();

% Loop over each subject
for i = 1:length(subs)
    % subject ID
    sub = subs{i};  
    % regions for this subject
    regions = fieldnames(subjects.(sub));
    % Loop over each region (e.g., 'front', 'occip')
    for j = 1:length(regions)
        % Region (front or occip)
        loc = regions{j};  
        % Load the local parameters dynamically
        [mus, ret, seg, seg_wm, ori, epvs] = ...
            load_local_params(subjects.(sub), loc);
        
        %%% Apply white matter mask
        % Retrieve white matter mask
        wm = subjects.(sub).(loc).mask_wm;
        % Set tissue outside WM mask to NaN
        mus(~wm) = NaN;
        ret(~wm) = NaN;

        % Exclude vessels from mus and retardance - ensure parenchyma
        % measurements are only of parenchyma
        mus(seg) = NaN;
        ret(seg) = NaN;
        ori(seg) = NaN;
        if ~isempty(epvs)
            epvs = logical(epvs);
            mus(epvs) = NaN;
            ret(epvs) = NaN;
            ori(epvs) = NaN;
        end

        %%% Measure average mus, retardance, orientation
        mus = mean(mus,'omitnan');
        ret = mean(ret,'omitnan');
        % Convert orientation to radians
        ori = rmmissing(deg2rad(ori(wm)));
        ori_std = circ_std(ori(:));

        %%% Save to struct
        fprintf('Finished subject %s region %s.\n',sub,loc)
        parench.(sub).(loc).mus = mus;
        parench.(sub).(loc).ret = ret;
        parench.(sub).(loc).ori = ori;

    end
end
