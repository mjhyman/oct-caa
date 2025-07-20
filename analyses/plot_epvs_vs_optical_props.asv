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
mus_combined = data.entire.mus_combined;
ret_combined = data.entire.ret_combined;
front_mus = data.entire.front_mus;
front_ret = data.entire.front_ret;
occip_mus = data.entire.occip_mus;
occip_ret = data.entire.occip_ret;

%%% Window averaged
mus = data.window.mus_combined;
mus_x = mus(1,:);
mus_y = mus(2,:);
ret = data.window.ret_combined;
ret_x = ret(1,:);
ret_y = ret(2,:);
% front mus
front_mus_x = data.window.front_mus_x;
front_mus_y = data.window.front_mus_y;
front_mus_std = data.window.front_mus_std;
front_mus_sem = data.window.front_mus_sem;
% front ret
front_ret_x = data.window.front_ret_x;
front_ret_y = data.window.front_ret_y;
front_ret_std = data.window.front_ret_std;
front_ret_sem = data.window.front_ret_sem;
% occip mus
occip_mus_x = data.window.occip_mus_x;
occip_mus_y = data.window.occip_mus_y;
occip_mus_std = data.window.occip_mus_std;
occip_mus_sem = data.window.occip_mus_sem;
% occip ret
occip_ret_x = data.window.occip_ret_x;
occip_ret_y = data.window.occip_ret_y;
occip_ret_std = data.window.occip_ret_std;
occip_ret_sem = data.window.occip_ret_sem;

%% Spearman's rho correlation coefficients of windows averages

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

%% Plot optical property vs. EPVS density (without window averaging)

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

%%% combined subjects and regions
% scattering coefficient
figure('Position',[100,100,900,900],'Resize','off');
scatter(mus_x, mus_y)
xlabel('EPVS Density');
ylabel('\mu_s (cm^-^1)');
title('Combined: \mu_s vs. EPVS Density');
set(gca,'fontsize',24)
fout = fullfile(fig_dir,'combined_mus_vs_epvs.png');
saveas(gcf,fout);
% retardance
figure('Position',[100,100,900,900],'Resize','off');
scatter(ret_x, ret_y)
xlabel('EPVS Density');
ylabel('retardance (degrees)');
title('Combined: Retardance vs. EPVS Density');
set(gca,'fontsize',24)
fout = fullfile(fig_dir,'combined_ret_vs_epvs.png');
saveas(gcf,fout);

%%% Frontal
% Scattering
figure('Position',[100,100,900,900],'Resize','off');
scatter(front_mus_x, front_mus_y)
xlabel('EPVS Density');
ylabel('\mu_s (cm^-^1)');
title('Frontal: \mu_s vs. EPVS Density');
set(gca,'fontsize',24)
fout = fullfile(fig_dir,'front_mus_vs_epvs.png');
saveas(gcf,fout);
% Retardance
figure('Position',[100,100,900,900],'Resize','off');
scatter(front_ret_x, front_ret_y)
xlabel('EPVS Density');
ylabel('retardance (degrees)');
title('Frontal: Retardance vs. EPVS Density');
set(gca,'fontsize',24)
fout = fullfile(fig_dir,'front_ret_vs_epvs.png');
saveas(gcf,fout);

%%% Occipital
% Scattering
figure('Position',[100,100,900,900],'Resize','off');
scatter(occip_mus_x, occip_mus_y)
xlabel('EPVS Density');
ylabel('\mu_s (cm^-^1)');
title('Occipital: \mu_s vs. EPVS Density');
set(gca,'fontsize',24)
fout = fullfile(fig_dir,'occip_mus_vs_epvs.png');
saveas(gcf,fout);
% Retardance
figure('Position',[100,100,900,900],'Resize','off');
scatter(occip_ret_x, occip_ret_y)
xlabel('EPVS Density');
ylabel('retardance (degrees)');
title('Occipital: Retardance vs. EPVS Density');
set(gca,'fontsize',24)
fout = fullfile(fig_dir,'occip_ret_vs_epvs.png');
saveas(gcf,fout);


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
