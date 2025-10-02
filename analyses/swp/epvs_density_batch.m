%% Measure the size-weighted proximity (SWP)
% Batch script that runs on individual subjects
clear; clc; close all;

% Input from bash
sid = str2double(getenv('SGE_TASK_ID'));
if isnan(sid)
    sid=4;
    fprintf('The SGE_TASK_ID was not passed to Matlab\n')
    fprintf('Setting SGE_TASK_ID to %s\n',string(sid))
end

% Exponent for denominator in SWP
p = 2;

%% Array of subject ID and regions
subjects = struct();
subjects(1).subject_name = 'caa6';
subjects(1).region = 'front';
subjects(1).p = p;
subjects(2).subject_name = 'caa6';
subjects(2).region = 'occip';
subjects(2).p = p;
subjects(3).subject_name = 'caa17';
subjects(3).region = 'occip';
subjects(3).p = p;
subjects(4).subject_name = 'caa22';
subjects(4).region = 'front';
subjects(4).p = p;
subjects(5).subject_name = 'caa22';
subjects(5).region = 'occip';
subjects(5).p = p;
subjects(6).subject_name = 'caa25';
subjects(6).region = 'front';
subjects(6).p = p;
subjects(7).subject_name = 'caa25';
subjects(7).region = 'occip';
subjects(7).p = p;
subjects(8).subject_name = 'caa26';
subjects(8).region = 'front';
subjects(8).p = p;
subjects(9).subject_name = 'caa26';
subjects(9).region = 'occip';
subjects(9).p = p;

% Set local subject_name and region
subject_name = subjects(sid).subject_name;
region = subjects(sid).region;
p = subjects(sid).p;
fprintf('subject = %s, region = %s, p = %d\n',subject_name, region,p)

%% Settings
data_dir = '/projectnb/npbssmic/ns/CAA/';
radius = 200;
save_base = '/projectnb/npbssmic/ns/CAA/swp/';

%% Map subject names to file paths
subject_files = struct( ...
    'caa6', '/caa6/caa6.mat', ...
    'caa17', '/caa17/occip/caa17.mat', ...
    'caa22', '/caa22/caa22.mat', ...
    'caa25', '/caa25/caa25.mat', ...
    'caa26', '/caa26/caa26.mat');

% Validate subject name
if ~isfield(subject_files, subject_name)
    error('Subject name "%s" not recognized.', subject_name);
end

% Subfolders for saving
subdir1 = subjects(sid).subject_name;
subdir2 = subjects(sid).region;

% Check if heatmap already exists
fname = fullfile(save_base, subdir1, subdir2,...
                strcat(subject_name, '_', region,...
                        '_radius_',num2str(radius),...
                        '_exp_',num2str(p), ...
                        '_interpolated_heatmap.mat'));
if isfile(fname)
    fprintf(['Heatmap already exists for sub=%s, reg=%s p=%d.' ...
        'Skipping.\n'], subject_name, region, p);
    return
end

% Load the subject data
file_path = fullfile(data_dir, subject_files.(subject_name));
if ~isfile(file_path)
    error('File does not exist\n%s',file_path)
else
    fprintf('Loading subject %s\n', subject_name)
end
tmp = load(file_path);
subject_data = tmp.(subject_name);
fprintf('Finished loading %s\n', subject_name)

% Validate region
if ~isfield(subject_data, region)
    error('Region "%s" not found in subject "%s" data.', region, subject_name);
end

% Ensure EPVS field exists
if isfield(subject_data.(region), 'epvs')
    epvs = subject_data.(region).epvs;
    mask = subject_data.(region).mask_wm;
    fprintf('Starting %s %s \n', subject_name, region);
else
    fprintf('No EPVS data for %s %s. Skipping.\n', subject_name, region);
    return
end

%% Calculate EPVS density
[subsampled_volume, interpolated_volume] = ...
    epvs_density_variable_p(epvs, mask, radius, p);

% Save results to .MAT and .TIF
save_epvs_heatmap(save_base, subject_name, region, ...
    subsampled_volume, interpolated_volume, radius, p);

fprintf('Finished processing %s %s \n', subject_name, region);