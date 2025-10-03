%% Analyze the EPVS density heatmap
% Import the EPVS density heatmaps
% Import the .MAT of the optical properties
% Measure correlation between EPVS density and optical properties

%% Add top-level directory of code repository to path
clear; clc; close all;
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
ncores = feature('numcores');
maxNumCompThreads(ncores);
% Flag for loading CAA structs (false if already in environment)
flag_load_caa_structs = false;
% Directory containing seg, mus, ret, mask, epvs structs
mat_dir = '/projectnb/npbssmic/s/mhyman/CAA_data/matlab_structs';
% Heatmap directory
heat_dir = '/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps';

%% Load heat maps

%%% subsampled EPVS heat maps
%{
caa17o = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa17/occip/caa17_occip_subsample_heatmap.mat']);
caa17o = caa17o.subsampled_volume;
caa22f = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa22/front/caa22_front_subsample_heatmap.mat']);
caa22f = caa22f.subsampled_volume;
caa22o = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa22/occip/caa22_occip_subsample_heatmap.mat']);
caa22o = caa22o.subsampled_volume;
caa25f = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa25/front/caa25_front_subsample_heatmap.mat']);
caa25f = caa25f.subsampled_volume;
caa25o = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa25/occip/caa25_occip_subsample_heatmap.mat']);
caa25o = caa25o.subsampled_volume;
caa26o = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa26/occip/caa26_occip_subsample_heatmap.mat']);
caa26o = caa26o.subsampled_volume;
% Output filename
stats_out = fullfile(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'epvs_heatmap_stats.mat']);
%}

%%% Fully interpolated EPVS heat maps
caa17o = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa17/occip/caa17_occip_interpolated_heatmap.mat']);
caa17o = caa17o.interpolated_volume;
caa22f = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa22/front/caa22_front_interpolated_heatmap.mat']);
caa22f = caa22f.interpolated_volume;
caa22o = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa22/occip/caa22_occip_interpolated_heatmap.mat']);
caa22o = caa22o.interpolated_volume;
caa25f = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa25/front/caa25_front_interpolated_heatmap.mat']);
caa25f = caa25f.interpolated_volume;
caa25o = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa25/occip/caa25_occip_interpolated_heatmap.mat']);
caa25o = caa25o.interpolated_volume;
caa26o = load(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'caa26/occip/caa26_occip_interpolated_heatmap.mat']);
caa26o = caa26o.interpolated_volume;
% Output filename
stats_out = fullfile(['/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/' ...
    'interpolated_epvs_heatmap_stats.mat']);

%% Load matlab structs
fprintf('Loading CAA17\n')
caa17 = load(fullfile(mat_dir,"caa17.mat"));
fprintf('Finished loading CAA17\n')

fprintf('Loading CAA22\n')
caa22 = load(fullfile(mat_dir,"caa22.mat"));
fprintf('Finished loading CAA22\n')

fprintf('Loading CAA25\n')
caa25 = load(fullfile(mat_dir,"caa25.mat"));
fprintf('Finished loading CAA25\n')

fprintf('Loading CAA26\n')
caa26 = load(fullfile(mat_dir,"caa26.mat"));
fprintf('Finished loading CAA26\n')

% Remove top-level struct
caa17 = caa17.caa17;
caa22 = caa22.caa22;
caa25 = caa25.caa25;
caa26 = caa26.caa26;

%% Create 2D arrays of optical property vs. EPVS density
% 2xN Matrix for each optical property:
%   - top row = EPVS density from heatmap
%   - bottom row = optical property from heatmap
%   - column = pairwise observations from same tissue volume

fprintf('Creating 2D arrays of EPVS density vs. optical prop\n')

% struct for storing pairs
heat_pairs = struct();

% struct for iterating over subjects
subjects = struct();
subjects.caa17 = caa17;
subjects.caa22 = caa22;
subjects.caa25 = caa25;
subjects.caa26 = caa26;

% subject names
subs = fields(subjects);
% add EPVS heatmaps to struct
subjects.caa17.occip.epvs_heat = caa17o;
subjects.caa22.front.epvs_heat = caa22f;
subjects.caa22.occip.epvs_heat = caa22o;
subjects.caa25.front.epvs_heat = caa25f;
subjects.caa25.occip.epvs_heat = caa25o;
subjects.caa26.occip.epvs_heat = caa26o;

for ii = 1:length(fields(subjects))
    sub = subs{ii};
    regions = fields(subjects.(subs{ii}));
    for j = 1:length(regions)
        % retrieve local properties
        reg = regions{j};
        if isfield(subjects.(sub).(reg), 'epvs_heat')
            epvs = subjects.(sub).(reg).epvs_heat;
            mus = subjects.(sub).(reg).mus;
            ret = subjects.(sub).(reg).ret_full;
            % Remove vessels from the epvs heatmap
            ves = subjects.(sub).(reg).seg;
            epvs(ves) = NaN;
            % call function to create pairs
            [mus_pair, ret_pair] = create_pair(epvs,mus,ret);
            % add pairs to heatmap struct
            heat_pairs.(sub).(reg).mus_pair = mus_pair;
            heat_pairs.(sub).(reg).ret_pair = ret_pair;
        else
            continue
        end
    end
end

%% Separate/combine heat_pairs across subject and region

fprintf('Separating heatmap pairs by subject and region\n')

%%% combine across occip + frontal
[mus_combined,ret_combined]=combine_subjects_regions(heat_pairs);

%%% Combine subjects, split regions
region_data = combine_subjects(heat_pairs);
% Frontal pairs
front_mus = region_data.front.mus;
front_ret = region_data.front.ret;
% Occip pairs
occip_mus = region_data.occip.mus;
occip_ret = region_data.occip.ret;

% Average optical property across sliding window

% Window size for EPVS density
window_size = 1e6;

%%% Combined values
fprintf('Averaging frontal+occipital\n')
[mus_x, mus_y, mus_std, mus_sem] = window_avg(mus_combined, window_size);
[ret_x, ret_y, ret_std, ret_sem] = window_avg(ret_combined, window_size);

%%% Frontal
fprintf('Averaging frontal\n')
[front_mus_x, front_mus_y, front_mus_std, front_mus_sem] =...
    window_avg(front_mus, window_size);
[front_ret_x, front_ret_y, front_ret_std, front_ret_sem] =...
    window_avg(front_ret, window_size);

%%% Occipital
fprintf('Averaging occipital\n')
[occip_mus_x, occip_mus_y, occip_mus_std, occip_mus_sem] =...
    window_avg(occip_mus, window_size);
[occip_ret_x, occip_ret_y, occip_ret_std, occip_ret_sem] =...
    window_avg(occip_ret, window_size);

%% Spearman's rho correlation coefficients of windows averages
%{
% Struct for storing rho and p-values
rhop = struct();

%%% combined subjects and regions
% scattering coefficient
[rho, p] = corr([mus_x',mus_y'],'type','Spearman','rows','complete');
rhop.combined.mus.rho = rho;
rhop.combined.mus.p = p;
% retardance
[rho, p] = corr([ret_x',ret_y'],'type','Spearman','rows','complete');
rhop.combined.ret.rho = rho;
rhop.combined.ret.p = p;

%%% Frontal
% Scattering
[rho, p] = corr([front_mus_x',front_mus_y'],'type','Spearman','rows','complete');
rhop.front.mus.rho = rho;
rhop.front.mus.p = p;
% Retardance
[rho, p] = corr([front_ret_x',front_ret_y'],'type','Spearman','rows','complete');
rhop.front.ret.rho = rho;
rhop.front.ret.p = p;

%%% Occipital
% Scattering
[rho, p] = corr([occip_mus_x',occip_mus_y'],'type','Spearman','rows','complete');
rhop.occip.mus.rho = rho;
rhop.occip.mus.p = p;
% Retardance
[rho, p] = corr([occip_ret_x',occip_ret_y'],'type','Spearman','rows','complete');
rhop.occip.ret.rho = rho;
rhop.occip.ret.p = p;
%}

%% Save the matrices

%{
fprintf('Saving the data struct\n')

data = struct();

%%% All pairs (not averaged)
data.entire.mus_combined = mus_combined;
data.entire.ret_combined = ret_combined;
data.entire.front_mus = front_mus;
data.entire.front_ret = front_ret;
data.entire.occip_mus = occip_mus;
data.entire.occip_ret = occip_ret;

%%% Window averaged
% combined mus
data.window.combined_mus_x = mus_x;
data.window.combined_mus_y = mus_y;
data.window.combined_mus_std = mus_std;
data.window.combined_mus_sem = mus_sem;
% combined ret
data.window.combined_ret_x = ret_x;
data.window.combined_ret_y = ret_y;
data.window.combined_ret_std = ret_std;
data.window.combined_ret_sem = ret_sem;
% front mus
data.window.front_mus_x = front_mus_x;
data.window.front_mus_y = front_mus_y;
data.window.front_mus_std = front_mus_std;
data.window.front_mus_sem = front_mus_sem;
% front ret
data.window.front_ret_x = front_ret_x;
data.window.front_ret_y = front_ret_y;
data.window.front_ret_std = front_ret_std;
data.window.front_ret_sem = front_ret_sem;
% occip mus
data.window.occip_mus_x = occip_mus_x;
data.window.occip_mus_y = occip_mus_y;
data.window.occip_mus_std = occip_mus_std;
data.window.occip_mus_sem = occip_mus_sem;
% occip ret
data.window.occip_ret_x = occip_ret_x;
data.window.occip_ret_y = occip_ret_y;
data.window.occip_ret_std = occip_ret_std;
data.window.occip_ret_sem = occip_ret_sem;

%%% Save output
save(stats_out,'data','rhop','-v7.3');
%}

%% Plot optical property vs. EPVS density

%%% Minimum threshold for EPVS density
th = 5*10^5;

%%% combined subjects and regions
% scattering coefficient
xlab = 'EPVS Volume Fraction (percentage occupied by EPVS)';
ylab = '\mu_s (cm^-^1)';
tit = 'Combined -- \mu_s vs. EPVS Density';
scatter_op(mus_combined, xlab, ylab, tit, th)
% retardance
xlab = 'EPVS Volume Fraction (percentage occupied by EPVS)';
ylab = 'retardance (degrees)';
tit = 'Combined -- Retardance vs. EPVS Density';
scatter_op(ret_combined, xlab, ylab, tit, th)

%%% Frontal
% Scattering
xlab = 'EPVS Volume Fraction (percentage occupied by EPVS)';
ylab = '\mu_s (cm^-^1)';
tit = 'Frontal -- \mu_s vs. EPVS Density';
scatter_op(front_mus, xlab, ylab, tit, th)
% Retardance
xlab = 'EPVS Volume Fraction (percentage occupied by EPVS)';
ylab = 'retardance (degrees)';
tit = 'Frontal -- Retardance vs. EPVS Density';
scatter_op(front_ret, xlab, ylab, tit, th)

%%% Occipital
% Scattering
xlab = 'EPVS Volume Fraction (percentage occupied by EPVS)';
ylab = '\mu_s (cm^-^1)';
tit = 'Occipital -- \mu_s vs. EPVS Density';
scatter_op(occip_mus, xlab, ylab, tit, th)
% Retardance
xlab = 'EPVS Volume Fraction (percentage occupied by EPVS)';
ylab = 'retardance (degrees)';
tit = 'Occipital -- Retardance vs. EPVS Density';
scatter_op(occip_ret, xlab, ylab, tit, th)

%%% Plot lower range of EPVS
% 



%%% Plot upper range of epvs
% Find the upper range after 1.5*IQR

%% Remove outliers and replot
outlier_rm = struct();

%% Plot the windowed averages w/ error bars

%%% combined subjects and regions
% scattering coefficient
xlab = 'EPVS Density';
ylab = '\mu_s (cm^-^1)';
tit = 'Combined -- \mu_s vs. EPVS Density';
sliding_errorbar_plot(mus_x, mus_y, mus_sem, xlab, ylab, tit)
% retardance
ylab = 'retardance (degrees)';
tit = 'Combined -- Retardance vs. EPVS Density';
sliding_errorbar_plot(ret_x, ret_y, ret_sem, xlab, ylab, tit)

%%% Frontal
% Scattering
ylab = '\mu_s (cm^-^1)';
tit = 'Frontal -- \mu_s vs. EPVS Density';
sliding_errorbar_plot(front_mus_x, front_mus_y, front_mus_sem, xlab, ylab, tit)
% Retardance
ylab = 'retardance (degrees)';
tit = 'Frontal -- Retardance vs. EPVS Density';
sliding_errorbar_plot(front_ret_x, front_ret_y, front_ret_sem, xlab, ylab, tit)

%%% Occipital
% Scattering
ylab = '\mu_s (cm^-^1)';
tit = 'Occipital -- \mu_s vs. EPVS Density';
sliding_errorbar_plot(occip_mus_x, occip_mus_y, occip_mus_sem, xlab, ylab, tit)
% Retardance
ylab = 'retardance (degrees)';
tit = 'Occipital -- Retardance vs. EPVS Density';
sliding_errorbar_plot(occip_ret_x, occip_ret_y, occip_ret_sem, xlab, ylab, tit)



%% Combine across all sujects/regions

function [mus_combined,ret_combined]=combine_subjects_regions(heat_pairs)
% Combine mus_pair across all subjects and all regions
% Output: all_mus_pairs is a 2xN matrix

% ----------- Pass 1: Count total columns ------------
total_samples = 0;
subjects = fieldnames(heat_pairs);

for i = 1:numel(subjects)
    subj = subjects{i};
    regions = fieldnames(heat_pairs.(subj));
    for j = 1:numel(regions)
        region = regions{j};
        mus = heat_pairs.(subj).(region).mus_pair;
        total_samples = total_samples + size(mus, 1);
    end
end

% ----------- Preallocate ----------------------------
mus_combined = zeros(total_samples,2);
ret_combined = zeros(total_samples,2);
col_idx = 1;

% ----------- Pass 2: Fill data ----------------------
for i = 1:numel(subjects)
    subj = subjects{i};
    regions = fieldnames(heat_pairs.(subj));
    for j = 1:numel(regions)
        region = regions{j};
        mus = heat_pairs.(subj).(region).mus_pair;
        ret = heat_pairs.(subj).(region).ret_pair;

        n = size(mus, 1);
        mus_combined(col_idx:col_idx+n-1,:) = mus;
        ret_combined(col_idx:col_idx+n-1,:) = ret;

        col_idx = col_idx + n;
    end
end
end


%% Combine across all sujects. Separate by regions
function region_data = combine_subjects(heat_pairs)
% Combine mus_pair and ret_pair across all subjects, separated by region
% Output: region_data.(region).mus and region_data.(region).ret are 2xN matrices

region_data = struct();

subjects = fieldnames(heat_pairs);
for i = 1:numel(subjects)
    subj = subjects{i};
    regions = fieldnames(heat_pairs.(subj));
    for j = 1:numel(regions)
        region = regions{j};
        
        % Extract mus_pair and ret_pair
        mus = heat_pairs.(subj).(region).mus_pair;
        ret = heat_pairs.(subj).(region).ret_pair;

        % Initialize or append
        if isfield(region_data, region)
            region_data.(region).mus = [region_data.(region).mus; mus];
            region_data.(region).ret = [region_data.(region).ret; ret];
        else
            region_data.(region).mus = mus;
            region_data.(region).ret = ret;
        end
    end
end
end

%% Scatter plot of combined vectors
function scatter_op(pair, xlab, ylab, tit, th)
% Scatter plot of the optical properts vs. EPVS density
% INPUTS:
%   pair (Nx2 matrix): [EPVS density, optical property]
%   xlab (string): x-axis label
%   ylab (string): y-axis label
%   tit (string): title of figure
%   th (double): minimum threshold for EPVS density

% Remove pairs where EPVS density < threshold (th)
keep_idx = pair(:,1) >= th;
pairf = pair(keep_idx,:);

% scatter plot
figure; scatter(pairf(:,1),pairf(:,2),20,'k','filled');
% Plot features
xlabel(xlab)
ylabel(ylab)
title(tit);
set(gca,'fontsize',20)
end

%% Errorbar plot of windowed averages
function sliding_errorbar_plot(x,y,err, xlab, ylab, tit)
% Scatter plot of the optical properts vs. EPVS density
% INPUTS:
%   pair (Nx2 matrix): [EPVS density, optical property]
%   err (vector): standard error at each window
%   xlab (string): x-axis label
%   ylab (string): y-axis label
%   tit (string): title of figure

% scatter plot
figure
errorbar(x,y,err,'o-')
% Plot features
xlabel(xlab)
ylabel(ylab)
title(tit);
set(gca,'fontsize',20)

end

%% Apply minimum threshold
function [pairf] = threshold_epvs(pair, th)

% Remove pairs where EPVS density < threshold (th)
keep_idx = pair(:,1) >= th;
pairf = pair(keep_idx,:);

end

