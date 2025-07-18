%% Measure the epvs density across each volume
% Measure the EPVS density at each parenchymal voxel
% Batch script that runs on individual subjects


%% Get input argument from Bash
sid = str2double(getenv('SGE_TASK_ID'));
% Array of subject ID and regions
subjects = struct();
subjects(1).subject_name = 'caa17';
subjects(1).region = 'occip';
subjects(1).p = 1;
subjects(2).subject_name = 'caa22';
subjects(2).region = 'front';
subjects(2).p = 1;
subjects(3).subject_name = 'caa22';
subjects(3).region = 'occip';
subjects(3).p = 1;
subjects(4).subject_name = 'caa25';
subjects(4).region = 'front';
subjects(4).p = 1;
subjects(5).subject_name = 'caa25';
subjects(5).region = 'occip';
subjects(5).p = 1;
subjects(6).subject_name = 'caa26';
subjects(6).region = 'occip';
subjects(6).p = 1;
% Add same subjectID + region with p=2
for ii=1:6
    % copy the subjectID and region
    subjects(ii+6) = subjects(ii);
    % set p=2
    subjects(ii+6).p = 2;
end

% Set local subject_name and region
subject_name = subjects(sid).subject_name;
region = subjects(sid).region;
p = subjects(sid).p;
fprintf('subject = %s, region = %s, p = %d\n',subject_name, region,p)

%% Settings
data_dir = '/projectnb/npbssmic/s/mhyman/CAA_data/matlab_structs/';
radius = 200;
fbase = '/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/';

%% Map subject names to file paths
subject_files = struct( ...
    'caa17', 'caa17.mat', ...
    'caa22', 'caa22.mat', ...
    'caa25', 'caa25.mat', ...
    'caa26', 'caa26.mat');

% Validate subject name
if ~isfield(subject_files, subject_name)
    error('Subject name "%s" not recognized.', subject_name);
end

% Check if heatmap already exists
fname = fullfile(fbase, strcat(subject_name, '_', region,...
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
fprintf('Loading subject %s\n', subject_name)
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

% Calculate EPVS density
[subsampled_volume, interpolated_volume] = ...
    epvs_density_variable_p(epvs, mask, radius,p);

% Save results to .MAT and .TIF
save_epvs_heatmap(fbase, subject_name, region, ...
    subsampled_volume, interpolated_volume, radius, p);

fprintf('Finished processing %s %s \n', subject_name, region);