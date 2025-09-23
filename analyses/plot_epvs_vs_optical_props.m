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
% Figure outputs
fig_dir = '/projectnb/npbssmic/s/mhyman/CAA_data/figures/';

%% Load heat map data struct
load('/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/epvs_heatmap_stats.mat');

%%% All pairs (not averaged)
combined_mus = data.entire.mus_combined;
combined_ret = data.entire.ret_combined;
front_mus = data.entire.front_mus;
front_ret = data.entire.front_ret;
occip_mus = data.entire.occip_mus;
occip_ret = data.entire.occip_ret;

%% Remove outliers

%%% Combined
% mus
[combined_mus_clean, mus_lb1, mus_lb2, mus_ub1, mus_ub2] =...
    remove_joint_outliers(combined_mus);
% ret
[combined_ret_clean, ret_lb1, ret_lb2, ret_ub1, ret_ub2] =...
    remove_joint_outliers(combined_ret);

%%% Front
% mus
front_mus_clean = cutoff_array(front_mus,mus_lb1,mus_ub1,mus_lb2,mus_ub2);
% ret
front_ret_clean = cutoff_array(front_ret,ret_lb1,ret_ub1,ret_lb2,ret_ub2);

%%% Occip
% mus
occip_mus_clean = cutoff_array(occip_mus,mus_lb1,mus_ub1,mus_lb2,mus_ub2);
% ret
occip_ret_clean = cutoff_array(occip_ret,ret_lb1,ret_ub1,ret_lb2,ret_ub2);

%% Take binned window average
% Size of x-axis window for averaging
window_size = 1e5;

%%% Combined
% mus
[comb_mus_x, comb_mus_y,~,~] = window_avg(combined_mus_clean,window_size);
% ret
[comb_ret_x, comb_ret_y,~,~] = window_avg(combined_ret_clean,window_size);

%%% Front
% mus
[front_mus_x, front_mus_y,~,~] = window_avg(front_mus_clean,window_size);
% ret
[front_ret_x, front_ret_y,~,~] = window_avg(front_ret_clean,window_size);

%%% Occip
% mus
[occip_mus_x, occip_mus_y,~,~] = window_avg(occip_mus_clean,window_size);
% ret
[occip_ret_x, occip_ret_y,~,~] = window_avg(occip_ret_clean,window_size);

%% Spearman's rho correlation coefficients of windows averages

% Struct for storing rho and p-values
rhop = struct();

%%% combined subjects and regions
% scattering coefficient
[rho, p] = corr([comb_mus_x',comb_mus_y'],'type','Spearman','rows','complete');
rhop.combined.mus.rho = rho(1,2);
rhop.combined.mus.p = p(1,2);
% retardance
[rho, p] = corr([comb_ret_x',comb_ret_y'],'type','Spearman','rows','complete');
rhop.combined.ret.rho = rho(1,2);
rhop.combined.ret.p = p(1,2);

%%% Frontal
% Scattering
[rho, p] = corr([front_mus_x',front_mus_y'],'type','Spearman','rows','complete');
rhop.front.mus.rho = rho(1,2);
rhop.front.mus.p = p(1,2);
% Retardance
[rho, p] = corr([front_ret_x',front_ret_y'],'type','Spearman','rows','complete');
rhop.front.ret.rho = rho(1,2);
rhop.front.ret.p = p(1,2);

%%% Occipital
% Scattering
[rho, p] = corr([occip_mus_x',occip_mus_y'],'type','Spearman','rows','complete');
rhop.occip.mus.rho = rho(1,2);
rhop.occip.mus.p = p(1,2);
% Retardance
[rho, p] = corr([occip_ret_x',occip_ret_y'],'type','Spearman','rows','complete');
rhop.occip.ret.rho = rho(1,2);
rhop.occip.ret.p = p(1,2);

%% Plot optical property vs. EPVS density (without window averaging)
%{
%%% Minimum threshold for EPVS density
th = 5*10^5;

%%% combined subjects and regions
% scattering coefficient
xlab = 'EPVS Volume Fraction (percentage occupied by EPVS)';
ylab = '\mu_s (cm^-^1)';
tit = 'Combined -- \mu_s vs. EPVS Density';
scatter_op(combined_mus, xlab, ylab, tit, th)
% retardance
xlab = 'EPVS Volume Fraction (percentage occupied by EPVS)';
ylab = 'retardance (degrees)';
tit = 'Combined -- Retardance vs. EPVS Density';
scatter_op(combined_ret, xlab, ylab, tit, th)

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
%}

%% Plot the windowed averages with error bars
%{
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
%}

%% Plot the windowed averages as scatter plot
% y-axis limits
xlims = [0, 4.2e6];
mus_ylims = [11,13.5];
ret_ylims = [21,32];
fsize = 32;

%%% combined subjects and regions
% scattering coefficient
figure('Position',[200,-400,900,900],'Resize','off');
scatter(comb_mus_x, comb_mus_y,100,'filled')
xlabel('EPVS Density'); ylabel('\mu_s (cm^-^1)');
title({'Combined:','\mu_s vs. EPVS Density'}); set(gca,'fontsize',fsize)
ylim(mus_ylims); xlim(xlims);
fout = fullfile(fig_dir,'combined_mus_vs_epvs.png');
saveas(gcf,fout); pause(0.5)
% retardance
figure('Position',[200,-400,900,900],'Resize','off');
scatter(comb_ret_x, comb_ret_y,100,'filled')
xlabel('EPVS Density'); ylabel('Retardance (\circ)');
title({'Combined:','Ret vs. EPVS Density'}); set(gca,'fontsize',fsize)
ylim(ret_ylims); xlim(xlims);
fout = fullfile(fig_dir,'combined_ret_vs_epvs.png');
saveas(gcf,fout); pause(0.5)

%%% Frontal
% Scattering
figure('Position',[200,-400,900,900],'Resize','off');
scatter(front_mus_x, front_mus_y,100,'filled')
xlabel('EPVS Density'); ylabel('\mu_s (cm^-^1)');
title({'Frontal:','\mu_s vs. EPVS Density'}); set(gca,'fontsize',fsize)
ylim(mus_ylims); xlim(xlims);
fout = fullfile(fig_dir,'front_mus_vs_epvs.png');
saveas(gcf,fout); pause(0.5)
% Retardance
figure('Position',[200,-400,900,900],'Resize','off');
scatter(front_ret_x, front_ret_y,100,'filled')
xlabel('EPVS Density'); ylabel('Retardance (\circ)');
title({'Frontal:','Ret vs. EPVS Density'}); set(gca,'fontsize',fsize)
ylim(ret_ylims); xlim(xlims);
fout = fullfile(fig_dir,'front_ret_vs_epvs.png');
saveas(gcf,fout); pause(0.5)

%%% Occipital
% Scattering
figure('Position',[200,-400,900,900],'Resize','off');
scatter(occip_mus_x, occip_mus_y,100,'filled')
xlabel('EPVS Density'); ylabel('\mu_s (cm^-^1)');
title({'Occipital:','\mu_s vs. EPVS Density'}); set(gca,'fontsize',fsize)
ylim(mus_ylims); xlim(xlims);
fout = fullfile(fig_dir,'occip_mus_vs_epvs.png');
saveas(gcf,fout); pause(0.5)
% Retardance
figure('Position',[200,-400,900,900],'Resize','off');
scatter(occip_ret_x, occip_ret_y,100,'filled')
xlabel('EPVS Density'); ylabel('Retardance (\circ)');
title({'Occipital:','Ret vs. EPVS Density'}); set(gca,'fontsize',fsize)
ylim(ret_ylims); xlim(xlims);
fout = fullfile(fig_dir,'occip_ret_vs_epvs.png');
saveas(gcf,fout); pause(0.5)

%% Apply upper and lower bounds to array
function data_clean = cutoff_array(data,lb1,ub1,lb2,ub2)

% Identify outliers in EPVS and optical property
outlier_col1 = (data(:,1) < lb1) | (data(:,1) > ub1);
outlier_col2 = (data(:,2) < lb2) | (data(:,2) > ub2);

% Combine outlier flags
outliers = outlier_col1 | outlier_col2;

% Remove outlier rows
data_clean = data(~outliers, :);

end

%% Remove outliers with IQR method
function [data_clean,lb1,lb2,ub1,ub2] = remove_joint_outliers(data)
% Remove outliers with interquartile range (IQR)
% INPUTS:
%   data (Nx2 array): data w/ outliers
%       Column 1: EPVS density
%       Column 2: Optical property
% OUTPUTS:
%   clean_data (Nx2) array: outliers removed

% Remove rows with NaNs first (optional but recommended)
data = data(~any(isnan(data),2), :);

% --- Outlier detection in Column 1 (EPVS density) ---
Q1 = prctile(data(:,1), 25);
Q3 = prctile(data(:,1), 75);
IQR = Q3 - Q1;
lb1 = Q1 - 1.5*IQR;
ub1 = Q3 + 1.5*IQR;
outlier_col1 = (data(:,1) < lb1) | (data(:,1) > ub1);

% --- Outlier detection in Column 2 (Optical property) ---
Q1 = prctile(data(:,2), 25);
Q3 = prctile(data(:,2), 75);
IQR = Q3 - Q1;
lb2 = Q1 - 1.5*IQR;
ub2 = Q3 + 1.5*IQR;
outlier_col2 = (data(:,2) < lb2) | (data(:,2) > ub2);

% Combine outlier flags
outliers = outlier_col1 | outlier_col2;

% Remove outlier rows
data_clean = data(~outliers, :);

% Display how many rows were removed
fprintf('Removed %d outliers out of %d samples\n', sum(outliers), size(data,1));

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

%% Sliding Window over x-axis
function [x_cent, y_mean, y_stds, y_sems] = slide_avg(data, window_size)
% Window average over x-axis, averaging y-values within each window.
% Disjoint windows.
%
% Inputs:
%   data        - Nx2 matrix: column 1 = x, column 2 = y
%   window_size - width of the window on the x-axis
%
% Outputs:
%   x_cent   - center x value of each window
%   y_mean   - average y value within each window
%   y_stds   - standard deviation of y within each window
%   y_sems   - standard error of y within each window

x = data(:,1);
y = data(:,2);

x_min = min(x);
x_max = max(x);

x_cent = [];
y_mean = [];
y_stds = [];
y_sems = [];

for start_x = x_min : window_size : x_max
    end_x = start_x + window_size;
    center = (start_x + end_x) / 2;

    in_window = x >= start_x & x < end_x;
    y_window = y(in_window);

    if ~isempty(y_window)
        x_cent(end+1) = center;
        y_mean(end+1) = mean(y_window);
        y_stds(end+1)  = std(y_window);
        y_sems(end+1)  = std(y_window) / sqrt(length(y_window));
    end
end
end
