%% Generate Spreadsheet of Measurements at Each Distance
%{
The purpose of this script is to create a spreadsheet of the optical
properties at each EPVS and Vessel at each distance. This covers just the
following cases:
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)
%}

%% Prepare environment
clear; clc; close all;

%%% Add top-level to path
% Get the current folder path
currentFolder = pwd;
% Go two levels up
[parentDir, ~, ~] = fileparts(currentFolder); % Get the parent directory
% Get all subdirectories from two levels up
subDirs = genpath(parentDir);
% Add the subdirectories to the MATLAB path
addpath(subDirs);

%%% Load 40 um donuts
fprintf('\nImporting Data!\n')
% top level directory
data_dir = '/projectnb/npbssmic/ns/CAA/';
% load the parenchymal optical properties donuts matrix
load(fullfile(data_dir, ...
    'parenchyma_optical_properties_40um_thick_09Oct2025.mat'));
% Load median optical properties
load(fullfile(data_dir,"median_white_matter_values_14-Nov-2025.mat"));
% Output filename for spreadsheet
dt = datetime('now','Format','dd-MM-yyyy');

%%% standardize_flag (int):
% 0 : do not change measurement
% 1 : use the absolute change from median
% 2 : use the percentage change from median
standardize_flag = 1;
if standardize_flag == 0
    spreadsheet_name = strcat('caa_all_radii_40um_donut_',...
                            string(dt),'.xlsx');
elseif standardize_flag == 1
    spreadsheet_name = strcat('caa_all_radii_median_subtracted_40um_donut_',...
                            string(dt),'.xlsx');
else
    spreadsheet_name = strcat('caa_all_radii_percentage_diff_40um_donut_',...
                            string(dt),'.xlsx');    
end
% Radius size
rad = 40;

%%% regions to iterate over
regions = {'front','occip'};

%% Set the upper bounds for mus and retardance
mus_upper = 25; % upper bound for mus
ret_upper = 45; % upper bound for retardance

%% Add the CAA disease stage to the parench struct
parench.caa26.stage = 0;
parench.caa6.stage = 1;
parench.caa17.stage = 2;
parench.caa22.stage = 3;
parench.caa25.stage = 3;

%% Count number of control and experimental values

% Counter for # samples in vessels
n_ves = 0;
% Counter for # samples in EPVS
n_epvs = 0;
% Initialize radii strings
rads = fields(parench.caa6.front);

%%% Iterate over subjects in parench
subs = fields(parench);
for ii = 1:length(subs)
    % retrieve subject
    sub = subs{ii};
    %%% iterate over regions
    for j = 1:numel(regions)
        % retrieve region
        reg = regions{j};
        % confirm that this subject contains this region
        if ~isfield(parench.(sub),reg)
            continue
        end
        %%% Iterate over radii
        for r = 1:numel(rads)
            % Initialize radius for iteration
            rad = rads{r};
             % retrieve optical properties (vessels and epvs)
            segmentations = fields(parench.(sub).(reg).(rad).outter);
            %%% Iterate segmentations & add to vectors
            for k=1:length(segmentations)
                % retrieve segmentation
                seg = segmentations{k};
                [mus,~,~,~] = retrieve_op(parench,sub,reg,rad,seg);
                % Count number of measurements in vessel or EPVS
                if strcmp(seg,'ves')
                    n_ves = n_ves + numel(mus);
                elseif strcmp(seg,'epvs')
                    n_epvs = n_epvs + numel(mus);
                end
            end
        end
    end
end

%% Create table containing all measurements across all distances

% Call function to vectorize data
[tbl_mus, tbl_ret, tbl_sori] = ...
    vectorize_struct(parench, regions, standardize_flag, parench_median,...
                    rad, n_ves, n_epvs,...
                    mus_upper, ret_upper);
% Export the tables to CSV
fout = fullfile(data_dir,spreadsheet_name);
fprintf('\nWriting to spreadsheets!\n')
writetable(tbl_mus,fout,'Sheet','scattering');
writetable(tbl_ret,fout,'Sheet','retardance');
writetable(tbl_sori,fout,'Sheet','orientation');


%% LMM for each test
function [tbl_mus, tbl_ret, tbl_sori] =...
    vectorize_struct(parench, regions, standardize_flag, parench_median,...
                    rad, n_ctl, n_exp, mus_upper, ret_upper)
% VECTORIZE_STRUCT: Vectorize the structure into a format for a spreadsheet
%
% INPUTS:
%   parench (struct): parenchyma optical properties struct
%   regions (cell array): cell array of region strings
%   standardize_flag (int):
%       0 : do not change measurement
%       1 : use the absolute change from median
%       2 : use the percentage change from median
%   parench_median (struct): median white matter values
%   rad (string): string indicating EPVS ring radius for structure
%   n_ctl (int): number of vessel measurements (control)
%   n_exp (int): number of EPVS measurements (experimental)
%   mus_upper (float): upper bound for scattering coefficient (1/cm)
%   ret_upper (float): upper bound for retardance (degrees)
% OUTPUTS:
%   tbl_mus (table): scattering coefficient table
%   tbl_ret (table): retardance table
%   tbl_sori (table): orientation table

% Arrays for control (ctl) and experimental (exp)
ctl_mus = zeros(n_ctl,1);
ctl_ret = zeros(n_ctl,1);
ctl_ori = zeros(n_ctl,1);
ctl_reg = cell(n_ctl,1);
ctl_subid = cell(n_ctl,1);
ctl_rad = zeros(n_ctl,1);
ctl_stage = zeros(n_ctl,1);
exp_mus = zeros(n_exp,1);
exp_ret = zeros(n_exp,1);
exp_ori = zeros(n_exp,1);
exp_reg = cell(n_exp,1);
exp_subid = cell(n_exp,1);
exp_rad = zeros(n_exp,1);
exp_stage = zeros(n_exp,1);
% Index to track location in arrays of data (ctl_* and exp_*)
cidx = 1;
eidx = 1;
% Subject ID vectors
ctl_tissue_id = zeros(n_ctl,1);
exp_tissue_id = zeros(n_exp,1);
% Counter for tissue ID
tiss_idx = 1;
% Initialize radii strings
rads = fields(parench.caa6.front);

%%% Iterate over subjects in parench
subs = fields(parench);
for ii = 1:length(subs)
    % update console
    fprintf('Subject %i out of %i\n',ii,length(subs))
    % retrieve subject
    sub = subs{ii};
    % Retrieve subject CAA stage
    stage = parench.(sub).stage;
    %%% iterate over regions
    for j = 1:numel(regions)
        % retrieve region
        reg = regions{j};
        % confirm that this subject contains this region
        if ~isfield(parench.(sub),reg)
            continue
        end
        % retrieve optical properties (vessels and epvs)
        segmentations = fields(parench.(sub).(reg).(rad).outter);
        %%% Iterate over radii
        for r = 1:numel(rads)
            % Initialize radius for iteration
            rad = rads{r};
            %%% Iterate segmentations & add to vectors
            for k=1:length(segmentations)
                % retrieve segmentation
                seg = segmentations{k};
                [mus,ret,ori,n] = retrieve_op(parench,sub,reg,rad,seg);
                %%% Calculate percentage change from median
                if standardize_flag == 1
                    % Retrieve median mus & ret from this subject/region
                    med_mus = parench_median.(sub).(reg).med.mus;
                    med_ret = parench_median.(sub).(reg).med.ret;
                    % subtract median then divide by median to find d%
                    mus = mus - med_mus;
                    ret = ret - med_ret;
                elseif standardize_flag == 2
                    % Retrieve median mus & ret from this subject/region
                    med_mus = parench_median.(sub).(reg).med.mus;
                    med_ret = parench_median.(sub).(reg).med.ret;
                    % subtract median then divide by median to find d%
                    mus = ((mus - med_mus) ./ med_mus) .* 100;
                    ret = ((ret - med_ret) ./ med_ret) .* 100;
                end
                %%% Vessel
                if strcmp(seg,'ves')
                    % copy radius (in microns) to radii vector
                    radius_val = sscanf(rad, 'rad%d');
                    distance = radius_val * 20;
                    ctl_rad(cidx:cidx+n-1,1) = ones(n,1) .* distance;
                    % copy optical properties, subid, region
                    [ctl_mus, ctl_ret, ctl_ori, ctl_tissue_id, ctl_subid,...
                        ctl_reg, ctl_stage, cidx] =...
                        copy_to_vector(mus, ret, ori,...
                                       ctl_mus, ctl_ret, ctl_ori,...
                                       ctl_tissue_id, ctl_subid, ctl_reg,...
                                       ctl_stage,...
                                       cidx, n, tiss_idx, sub, reg, stage);
                %%% EPVS
                elseif strcmp(seg,'epvs')
                    % copy radius (in microns) to radii vector
                    radius_val = sscanf(rad, 'rad%d');
                    distance = radius_val * 20;
                    exp_rad(eidx:eidx+n-1,1) = ones(n,1) .* distance;
                    % copy optical properties, subid, region
                    [exp_mus, exp_ret, exp_ori, exp_tissue_id, exp_subid,...
                        exp_reg, exp_stage, eidx] =...
                        copy_to_vector(mus, ret, ori,...
                                       exp_mus, exp_ret, exp_ori,...
                                       exp_tissue_id, exp_subid, exp_reg,...
                                       exp_stage,...
                                       eidx, n, tiss_idx, sub, reg, stage);
                end
            end
        end
        % Iterate tissue sample counter
        tiss_idx = tiss_idx + 1;
    end
end

%% Create a table for fitting the LME model
% Column 1 = group label
% Column 2 = region
% Column 3 = Tissue ID (TODO)
% Column 4 = Subject ID
% column 5 = distance (micron)
% Column 6 = value of measurement (mus, ret, ori)

%%% Create copies of regions, subject IDs, and distances
% These are for assigning the tables
% mus
mus_exp_reg = exp_reg;
mus_ctl_reg = ctl_reg;
mus_exp_tid = exp_tissue_id;
mus_ctl_tid = ctl_tissue_id;
mus_exp_sid = exp_subid;
mus_ctl_sid = ctl_subid;
mus_exp_rad = exp_rad;
mus_ctl_rad = ctl_rad;
mus_exp_stage = exp_stage;
mus_ctl_stage = ctl_stage;
% retardance
ret_exp_reg = exp_reg;
ret_ctl_reg = ctl_reg;
ret_exp_tid = exp_tissue_id;
ret_ctl_tid = ctl_tissue_id;
ret_exp_sid = exp_subid;
ret_ctl_sid = ctl_subid;
ret_exp_rad = exp_rad;
ret_ctl_rad = ctl_rad;
ret_exp_stage = exp_stage;
ret_ctl_stage = ctl_stage;

%%% Identify indices of values to keep
exp_mus_idx = exp_mus < mus_upper;
exp_ret_idx = exp_ret < ret_upper;
ctl_mus_idx = ctl_mus < mus_upper;
ctl_ret_idx = ctl_ret < ret_upper;

%%% Remove outliers from region, subjectID, radii, mus, retardance
% mus
mus_exp_reg = mus_exp_reg(exp_mus_idx);
mus_ctl_reg = mus_ctl_reg(ctl_mus_idx);
mus_exp_tid = mus_exp_tid(exp_mus_idx);
mus_ctl_tid = mus_ctl_tid(ctl_mus_idx);
mus_exp_sid = mus_exp_sid(exp_mus_idx);
mus_ctl_sid = mus_ctl_sid(ctl_mus_idx);
mus_exp_rad = mus_exp_rad(exp_mus_idx);
mus_ctl_rad = mus_ctl_rad(ctl_mus_idx);
mus_exp_stage = mus_exp_stage(exp_mus_idx);
mus_ctl_stage = mus_ctl_stage(ctl_mus_idx);
% retardance
ret_exp_reg = ret_exp_reg(exp_ret_idx);
ret_ctl_reg = ret_ctl_reg(ctl_ret_idx);
ret_exp_tid = ret_exp_tid(exp_ret_idx);
ret_ctl_tid = ret_ctl_tid(ctl_ret_idx);
ret_exp_sid = ret_exp_sid(exp_ret_idx);
ret_ctl_sid = ret_ctl_sid(ctl_ret_idx);
ret_exp_rad = ret_exp_rad(exp_ret_idx);
ret_ctl_rad = ret_ctl_rad(ctl_ret_idx);
ret_exp_stage = ret_exp_stage(exp_ret_idx);
ret_ctl_stage = ret_ctl_stage(ctl_ret_idx);
% Experimental mus & ret
exp_mus = exp_mus(exp_mus_idx);
exp_ret = exp_ret(exp_ret_idx);
% Control mus & ret
ctl_mus = ctl_mus(ctl_mus_idx);
ctl_ret = ctl_ret(ctl_ret_idx);

%%% Optical property (don't change orientation, it doesn't have outliers)
% Scattering coefficient
mus = [exp_mus; ctl_mus];
% Retardance
ret = [exp_ret; ctl_ret];
% Circular standard deviation of orientation
ori = [exp_ori; ctl_ori];

%%% group labels
% mus
g_exp = repmat({'experimental'},[length(exp_mus),1]);
g_cnrtl = repmat({'control'},[length(ctl_mus),1]);
mus_labels = vertcat(g_exp, g_cnrtl);
% retardance
g_exp = repmat({'experimental'},[length(exp_ret),1]);
g_cnrtl = repmat({'control'},[length(ctl_ret),1]);
ret_labels = vertcat(g_exp, g_cnrtl);
% Orientation
g_exp = repmat({'experimental'},[length(exp_ori),1]);
g_cnrtl = repmat({'control'},[length(ctl_ori),1]);
ori_labels = vertcat(g_exp, g_cnrtl);

%%% Mus Region, Subject ID, Distance, 
% Combine the brain region indices into column vector
mus_region = vertcat(mus_exp_reg, mus_ctl_reg);
% Combine subject ID into column vector
mus_subid = vertcat(mus_exp_sid, mus_ctl_sid);
% Combine the subject IDs into column vector
mus_tids = vertcat(mus_exp_tid, mus_ctl_tid);
% Distance vector
mus_radii = [mus_exp_rad; mus_ctl_rad];
% CAA stage vector
mus_stage = [mus_exp_stage; mus_ctl_stage];

%%% retardance Region, Subject ID, Distance, 
% Combine the brain region indices into column vector
ret_region = vertcat(ret_exp_reg, ret_ctl_reg);
% Combine subject ID into column vector
ret_subid = vertcat(ret_exp_sid, ret_ctl_sid);
% Combine the subject IDs into column vector
ret_tids = vertcat(ret_exp_tid, ret_ctl_tid);
% Distance vector
ret_radii = [ret_exp_rad; ret_ctl_rad];
% CAA stage vector
ret_stage = [ret_exp_stage; ret_ctl_stage];

%%% Orientation Region, Subject ID, Distance, 
% Combine the brain region indices into column vector
ori_region = vertcat(exp_reg, ctl_reg);
% Combine subject ID into column vector
ori_subid = vertcat(exp_subid, ctl_subid);
% Combine the subject IDs into column vector
ori_tids = vertcat(exp_tissue_id, ctl_tissue_id);
% Distance vector
ori_radii = [exp_rad;ctl_rad];
% CAA stage vector
ori_stage = [exp_stage; ctl_stage];

%%% Create Table for each optical property
varnames = {'Groups','subjectID','Stage','Region','tissueID','distance','OpticalProperty'};
tbl_mus = table(mus_labels,mus_subid,mus_stage,mus_region,mus_tids,mus_radii,mus,...
                'VariableNames',varnames);
tbl_ret = table(ret_labels,ret_subid,ret_stage,ret_region,ret_tids,ret_radii,ret,...
                'VariableNames',varnames);
tbl_sori = table(ori_labels,ori_subid,ori_stage,ori_region,ori_tids,ori_radii,ori,...
                'VariableNames',varnames);
% Take real part of complex numbers
tbl_sori.OpticalProperty = real(tbl_sori.OpticalProperty);
end