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
- IMPORT struct containing measurements at various distances
- Create scatterplot (property vs. distance)
- Create box/whisker plot

%}

%% Prepare environment
clear; clc; close all;
% Add top-level directory
current_dir = pwd;
addpath(fullfile(current_dir));

%%% Input directory (Martinos or SCC)
% SCC
data_dir = '/projectnb/npbssmic/ns/CAA/';

%%% Generate output directories
scat_out = fullfile(data_dir,'/scatter_plots');
bw_out = fullfile(data_dir,'/bw_plots');

%%% Load parenchyma struct
% This file contains the measurements at each radius
load(fullfile(data_dir, ...
    "parenchyma_optical_properties_40um_thick_5-Mar-2026.mat"));

%%% Load median white matter value of each subject/region
% Load median white matter values for each subject
load(fullfile(data_dir, 'median_white_matter_values_14-Nov-2025.mat'));

%%% Figure parameters
% Scatter plot dot size
scat_size = 100;
% Isotropic voxel size (microns)
vox = 20;
% Boolean for errorbar plots
err_flag = false;
% Date string for filename
dt = datetime('now','Format','dd-MM-yyyy');

%%% Setting flag for removing outliers
% Flag to remove outliers
rm_outlier = true;
% Create filename 
if rm_outlier
    substr = strcat('_max_removed_OutlierRm_40um_donut_',string(dt));
else
    substr = strcat('_max_removed_40um_donut_03Nov2025',string(dt));
end

%%% Setting flag for subtracting baseline median from vectors
% If this is true, then the following will be achieved:
%   1) remove values above threshold
%   2) remove values outside 1.5*IQR
%   3) subtract the median WM optical property from vector
subtract_median = true;
if subtract_median
    substr = strcat(['_max_removed__outlier_removed__' ...
                    'median_subtracted__40um_donut_'],string(dt));
    standardize_flag = 1;
end

%%% Setting flag for calculating % change from median
% If this is true, then the following will be achieved:
%   1) remove values above threshold
%   2) remove values outside 1.5*IQR
%   4) calculate percentage change from tissue median
%       - subtract the median WM optical property from vector
%       - divide each value by the median WM optical property of tissue
use_percentage = true;
if use_percentage
    substr = strcat(['_max_removed__outlier_removed__median_percentage_diff__' ...
              '40um_donut_'],string(dt));
    standardize_flag = 2;
end


%% Figure axes and limits

%%% Extract number of radii
radii = fieldnames(parench.caa22.occip);
x = zeros(length(radii),1);
% Retrieve inner radii (this will be x-axis)
for ii = 1:length(radii)
    x(ii) = str2num(radii{ii}(4:end));
end
% Convert x-axis to microns
x = x.*vox;
% Scatterplot x-axis ticks
xt = 0:100:500;
% Scatterplot point size
psize = 1000;

%%% Upper threshold for keeping mus and retardance values
mus_max = 25; % 1/cm
ret_max = 45; % degrees

%%% Scatterplot limits for each optical property (no errorbars)
% Set limits for the case of subtracting offset
if use_percentage
    % y-axis limits
    mus_yl = [-5, 10];
    ret_yl = [-10, 10];
    % scatterplot y-axis tick marks
    mus_yt = min(mus_yl) : 2.5 : max(mus_yl);
    ret_yt = min(ret_yl) : 5 : max(ret_yl);
elseif subtract_median
    % y-axis limits
    mus_yl = [-0.1, 0.8];
    ret_yl = [-3, 0.5];
    % scatterplot y-axis tick marks
    mus_yt = min(mus_yl) : 0.2 : max(mus_yl);
    ret_yt = min(ret_yl) : 0.5 : max(ret_yl);
else
    % y-axis limits
    mus_yl = [6, 14];
    ret_yl = [19, 33];
    % scatterplot y-axis tick marks
    mus_yt = min(mus_yl) : 2 : max(mus_yl);
    ret_yt = min(ret_yl) : 2 : max(ret_yl);
end
% Set limits for orientation (unaffected by offset subtraction)
ori_yl = [0, 1];
ori_yt = min(ori_yl) : 0.1 : max(ori_yl);

% place figure properties in a struct for consistency
fprops = struct();
fprops.mus_yl = mus_yl;
fprops.mus_yt = mus_yt;
fprops.ret_yl = ret_yl;
fprops.ret_yt = ret_yt;
fprops.ori_yl = ori_yl;
fprops.ori_yt = ori_yt;
fprops.xt = xt;
fprops.x = x;
fprops.scat_out = scat_out;
fprops.subdir = '/all_subjects/';
fprops.radii = radii;
fprops.parench = parench;
fprops.err_flag = err_flag;
fprops.mus_max = mus_max;
fprops.ret_max = ret_max;
fprops.substr = substr;
fprops.psize = psize;
fprops.rm_outlier = rm_outlier;

%% Scatterplot of optical property vs. distance (disjoint ring)
% y-axis = average optical property  for epvs or vessel
% x-axis = distance of ring from edge of epvs or vessel
% repeat this for all 3 properties (mus, ret, ori)
% This will only use the "outer" ring, which is disjoint from the annotated
% vessel/epvs

% retrieve subject IDs of all subjects
subjects = fieldnames(parench);

%%% call function to create scatter plots
%{
%   standardize_flag (int):
%                  0 : do not change measurement
%                  1 : use the absolute change from parench_median
%                  2 : use the percentage change from parench_median
scatter_main(subjects, fprops, parench_median, standardize_flag)
%}

%% Scatterplot optical property vs. distance (separate subjects)
% y-axis = average optical property  for epvs or vessel
% x-axis = distance of ring from edge of epvs or vessel
% Create for each subject

%%% Rest the scatterplot limits
% Set limits for the case of subtracting offset
if use_percentage
    % y-axis limits
    mus_yl = [-20, 30];
    ret_yl = [-20, 20];
    % scatterplot y-axis tick marks
    mus_yt = min(mus_yl) : 5 : max(mus_yl);
    ret_yt = min(ret_yl) : 5 : max(ret_yl);
elseif subtract_median
    % y-axis limits
    mus_yl = [-0.1, 0.8];
    ret_yl = [-3, 0.5];
    % scatterplot y-axis tick marks
    mus_yt = min(mus_yl) : 0.2 : max(mus_yl);
    ret_yt = min(ret_yl) : 0.5 : max(ret_yl);
else
    % y-axis limits
    mus_yl = [6, 14];
    ret_yl = [19, 33];
    % scatterplot y-axis tick marks
    mus_yt = min(mus_yl) : 2 : max(mus_yl);
    ret_yt = min(ret_yl) : 2 : max(ret_yl);
end
% Add to struct
fprops.mus_yl = mus_yl;
fprops.mus_yt = mus_yt;
fprops.ret_yl = ret_yl;
fprops.ret_yt = ret_yt;

% CAA 6
subjects = {'caa6'};
fprops.subdir = '/caa6/';
scatter_main(subjects,fprops,parench_median,standardize_flag);

% CAA 17
subjects = {'caa17'};
fprops.subdir = '/caa17/';
scatter_main(subjects,fprops,parench_median,standardize_flag);

% CAA 22
subjects = {'caa22'};
fprops.subdir = '/caa22/';
scatter_main(subjects,fprops,parench_median,standardize_flag);

% CAA 25
subjects = {'caa25'};
fprops.subdir = '/caa25/';
scatter_main(subjects,fprops,parench_median,standardize_flag);

% CAA 26
subjects = {'caa26'};
fprops.subdir = '/caa26/';
scatter_main(subjects,fprops,parench_median,standardize_flag);

%% Box/whisker Plots

%{
% Reorder the subjects list in terms of severity
subjects = {'caa6','caa26','caa17','caa25','caa22'};

%%% Iterate distance, subjects, optical property
% iterate distance from edge of epvs/vessel
for ii = 1:length(radii)
    % local radii name
    rname = string(radii(ii));
    % cell arrays to store the values grouped by subject
    ves_mus = {}; ves_ret = {}; ves_ori = {}; ves_name = {};
    epvs_mus = {}; epvs_ret = {}; epvs_ori = {}; epvs_name = {};
    % cell arrays to store names for each observation
    ves_mus_name = {}; ves_ret_name = {}; ves_ori_name = {};
    epvs_mus_name = {}; epvs_ret_name = {}; epvs_ori_name = {};
    % counter for tracking volume index and epvs index
    vol_idx = 1;
    epvs_idx = 1;
    % iterate subjects
    for j = 1:length(subjects)
        % local subject ID
        sub = subjects{j};
        % retrieve regions for this subject
        regions = fieldnames(parench.(sub));
        % iterate regions
        for k = 1:length(regions)
            % region for this iteration
            reg = regions{k};
            % Retrieve vessel measurements (mus, ret, ori)
            ves = parench.(sub).(reg).(rname).outter.ves;
            
            %%% Retrieve mus, ret, ori
            ves_mus{vol_idx} = rmmissing(ves.pmus);
            ves_ret{vol_idx} = rmmissing(ves.pret);
            ves_ori{vol_idx} = real(rmmissing(ves.pori));
            
            %%% Create x-axis names for each datapoint
            % Create x-axis name for boxplot
            xname = [sub,'_',reg];
            % # observations for each optical property
            nmus = length(ves_mus{vol_idx});
            nret = length(ves_ret{vol_idx});
            nori = length(ves_ori{vol_idx});
            % replicate x-name for each observation
            mus_name = repmat({xname},nmus,1);
            ret_name = repmat({xname},nret,1);
            ori_name = repmat({xname},nori,1);
            % Add vector of names to main cell array
            ves_mus_name{vol_idx} = mus_name;
            ves_ret_name{vol_idx} = ret_name;
            ves_ori_name{vol_idx} = ori_name;
            % Iterate counter
            vol_idx = vol_idx + 1;
            % retrieve EPVS measurements (if exist)
            if isfield(parench.(sub).(reg).(rname).outter, 'epvs')
                epvs = parench.(sub).(reg).(rname).outter.epvs;
                % Create x-axis name for boxplot
                epvs_name{vol_idx} = [sub,'_',reg];
                % Take average of mus, ret, ori
                epvs_mus{epvs_idx} = rmmissing(epvs.pmus);
                epvs_ret{epvs_idx} = rmmissing(epvs.pret);
                epvs_ori{epvs_idx} = real(rmmissing(epvs.pori));
               
                %%% Create x-axis names for each datapoint
                % Create x-axis name for boxplot
                xname = [sub,'_',reg];
                % # observations for each optical property
                nmus = length(epvs_mus{epvs_idx});
                nret = length(epvs_ret{epvs_idx});
                nori = length(epvs_ori{epvs_idx});
                % replicate x-name for each observation
                mus_name = repmat({xname},nmus,1);
                ret_name = repmat({xname},nret,1);
                ori_name = repmat({xname},nori,1);
                % Add vector of names to main cell array
                epvs_mus_name{epvs_idx} = mus_name;
                epvs_ret_name{epvs_idx} = ret_name;
                epvs_ori_name{epvs_idx} = ori_name;
                % Iterate counter
                epvs_idx = epvs_idx + 1;
            end
        end
    end
    
    %%% Box/Whisker Plot (mus)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(ves_mus_name{:});
    % Create y-axis of all concatenated observations
    y = [ves_mus{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,50]); set(bx,'LineWidth',2);
    ylabel('\mus (cm^-^1)');
    title({'Vessels - \mus',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('ves_mus_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;

    %%% Box/Whisker Plot (ret)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(ves_ret_name{:});
    % Create y-axis of all concatenated observations
    y = [ves_ret{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,50]); set(bx,'LineWidth',2);
    ylabel('Retardance (degrees)');
    title({'Vessels - Retardance',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('ves_ret_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;

    %%% Box/Whisker Plot (ori)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(ves_ori_name{:});
    % Create y-axis of all concatenated observations
    y = [ves_ori{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,3]); set(bx,'LineWidth',2);
    ylabel('\sigma_{orientation} (radians)');
    title({'Vessels - Orientation',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('ves_ori_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;
    
    %%% EPVS Box/Whisker Plot (mus)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(epvs_mus_name{:});
    % Create y-axis of all concatenated observations
    y = [epvs_mus{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,30]); set(bx,'LineWidth',2);
    ylabel('\mus (cm^-^1)');
    title({'EPVS - \mus',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('epvs_mus_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;

    %%% EPVS Box/Whisker Plot (ret)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(epvs_ret_name{:});
    % Create y-axis of all concatenated observations
    y = [epvs_ret{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,50]); set(bx,'LineWidth',2);
    ylabel('Retardance (degrees)');
    title({'EPVS - Retardance',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('epvs_ret_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;

    %%% EPVS Box/Whisker Plot (ori)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(epvs_ori_name{:});
    % Create y-axis of all concatenated observations
    y = [epvs_ori{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,3]); set(bx,'LineWidth',2);
    ylabel('\sigma_{orientation} (radians)');
    title({'EPVS - Orientation',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('epvs_ori_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;
end

%}

%% Student-t's 95% confidence interval at each distance
% Each data point in the scatter plot represents the mean of each subject
% and region. The N for each point in the scatter plot is small (4-9), so
% the standard error of the mean will be large. Instead the confidence
% interval provides an alternative for this situation

function ci_half = calc_student_t_ci(vector, alpha)
% Calculate the standard error of the mean
n = length(vector);
sem = std(real(vector)) ./ sqrt(n);
% Calculate the critical t-value for 95% confidence interval
df = n-1;
t_crit = tinv(1 - alpha/2,df);
% Multiple t_crit with the standard error of the mean
ci_half = t_crit .* sem;
end

%% Bootstrap 95% confidence interval at each distance
% Each data point in the scatter plot represents the mean of each subject
% and region. The N for each point in the scatter plot is small (4-9), so
% the standard error of the mean will be large. Instead bootstrap the
% vector at each distance to estimate the CI

function [ci_half] = bootstrap_ci(vector, alpha)
% Sample data
x = vector;
% Number of bootstrap samples
B = 10000;
N = length(x);
boot_means = zeros(B,1);

%%% Iterate the bootstrap sampling - take mean of each group of samples
for b = 1:B
    % Sample with replacement
    sample_b = datasample(x, N);
    % Compute mean
    boot_means(b) = mean(sample_b);
end

%%% Calculate 95% CI
% Calculate indices for confidence interval (CI)
lower = alpha/2;
upper = 1 - alpha/2;
% Compute 95% percentile CI (2.5th and 97.5th percentiles)
lower_idx = round(lower * B);
upper_idx = round(upper * B);
% Handle edge cases
if lower_idx < 1, lower_idx = 1; end
if upper_idx > B, upper_idx = B; end
% Sort bootstrap means
boot_means = sort(boot_means);
% Extract lower/upper CI
ci_lower = boot_means(lower_idx);
ci_upper = boot_means(upper_idx);
ci_half = (ci_upper - ci_lower)/2;
end


%% Function to measure the mean optical property across a single subject
function [ves_mus, ves_ret, ves_ori, epvs_mus, epvs_ret, epvs_ori] =...
    meas_sub(parench,sub,reg)

%%% Initialize arrays for storing optical properties measurements
% Retrieve radii within struct
radii = fieldnames(parench.(sub).(reg));
% Initialize array for pairs of optical proprety vs. distance
ves_mus = zeros(length(radii),1);
ves_ret = zeros(length(radii),1);
ves_ori = zeros(length(radii),1);
epvs_mus = zeros(length(radii),1);
epvs_ret = zeros(length(radii),1);
epvs_ori = zeros(length(radii),1);

%%% Iterate distance & optical property
% iterate distance from edge of epvs/vessel
for ii = 1:length(radii)
    % local radii name
    rname = string(radii(ii));
    
    %%% Vessel measurements
    % Retrieve vessel measurements (mus, ret, ori)
    ves = parench.(sub).(reg).(rname).outter.ves;
    % Take average of mus, ret, ori
    ves_mus(ii) = mean(ves.pmus,'omitnan');
    ves_ret(ii)  = mean(ves.pret,'omitnan');
    ves_ori(ii) = mean(ves.pori,'omitnan');

    %%% EPVS measurements (if exist)
    if isfield(parench.(sub).(reg).(rname).outter, 'epvs')
        epvs = parench.(sub).(reg).(rname).outter.epvs;
        % Take average of mus, ret, ori
        epvs_mus(ii) = mean(epvs.pmus,'omitnan');
        epvs_ret(ii) = mean(epvs.pret,'omitnan');
        epvs_ori(ii) = real(mean(epvs.pori,'omitnan'));
    end
end
end