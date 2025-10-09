%% Measure the size-weighted proximity (SWP)
% Batch script that runs on individual subjects
% clear; clc; close all;

%% Settings
data_dir = '/projectnb/npbssmic/ns/CAA/';
swp_dir = '/projectnb/npbssmic/ns/CAA/swp/';
% SWP string base
str_base = 'radius_200_exp_2_interpolated_heatmap_log10.mat';
% Directory to store csv
csv_out = '/projectnb/npbssmic/ns/CAA/metrics/';
% flag for importing top-level
import_all_data = false;



%% Array of subject ID and regions
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

%% Import the data structs

if import_all_data
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
        load(fullfile(data_dir,subject_files.(subs{ii})));
    end
end

%% Import the SWP structs

% Percentiles for thresholding
pcts = [5, 10, 25];
% Preallocate for results
results_subject = {};
results_region  = {};
results_means   = [];

% Import swp for each
for ii = 1:length(subjects)
    % Create full file name
    subid = subjects(ii).subject_name;
    reg = subjects(ii).region;
    fname = strcat(subid,'_',reg,'_');
    fname = strcat(fname, str_base);
    fname = fullfile(swp_dir,subid,reg,fname);    
    % Import SWP
    swp = load(fname);
    swp = swp.swp;
    % Keep the values that are NOT NaN
    swp = single(swp(~isnan(swp)));
    
    %%% Compute mean of upper 5%, 10%, 25%
    
    means_upper = zeros(size(pcts));
    % Sort low to high
    vals_sorted = sort(swp, 'ascend');
    N = numel(swp);
    for j = 1:numel(pcts)
        pct = pcts(j);
        n_upper = ceil(N * (pct/100));
        upper_vals = vals_sorted(end - n_upper + 1 : end);
        means_upper(j) = mean(upper_vals);
        % Store results
        results_subject{end+1,1} = subid;
        results_region{end+1,1}  = reg;
        results_means = [results_means; means_upper(:)'];
    end
end

% Build output table (percentiles are columns)
T = table(results_subject, results_region, ...
    results_means(:,1), results_means(:,2), results_means(:,3), ...
    'VariableNames', {'Subject','Region','MeanUpper5','MeanUpper10','MeanUpper25'});
fout = fullfile(csv_out,'mean_upper_percentiles.csv');

% Write to Excel
writetable(T, fout);

