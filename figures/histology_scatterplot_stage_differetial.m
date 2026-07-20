%% Histopathology Mean + SD of all stages combined
%{
This script measures the average stain intensity (units = z-score with
respect to the section's median stain intensity) and standard deviation
across all sections for each region.

This is an alternative to the bayesian statisical model, which had
difficulty converging.

TODO:
- Measure subject level
- Combine the severe
- Subtract control from mild and severe (mild-ctl) & (severe-ctl)
- Plot the difference
- Plot the CD68 separately
%}

%% Top-level settings
clearvars -except caa6 caa17 caa22 caa25 caa26 swp_struct
clc; close all;
% Get the current folder
currentFolder = pwd;
% Move one directory up
parentFolder = fileparts(currentFolder);
% Add the parent folder and all its subfolders to the MATLAB search path
addpath(genpath(parentFolder));

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/';
% Directory to save output figures of averages across sections
figdir = ['/projectnb/npbssmic/ns/CAA/figures/' ...
         'histology_scatterplots/histo_matched/stage_diff/'];
% Directory to save figures of averages within section
subdir = fullfile(figdir, '/subject_level');

%%% Load the stain and segmentation structs
lhe = load(fullfile(hdir,'lhe_rings_21-Nov-2025.mat')); lhe = lhe.lhe;
gfap = load(fullfile(hdir,'gfap_rings_21-Nov-2025.mat')); gfap = gfap.gfap;
cd68 = load(fullfile(hdir,'cd68_rings_21-Nov-2025.mat')); cd68 = cd68.cd68;

%%% measurement outter radii (units = voxels)
radii_sm = [40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486];

%% Take average of severe cases
% Combine subjects but keep regions separate
% Only use the samll radius measurement (rad40)
close all;

%%% Extract the severe sections from each stain
lhe_sev = lhe(2:5);
gfap_sev = gfap(2:5);
cd68_sev = cd68(2:5);

% Define the tissue section codes for frontal and occipital
fcode = "_1";
ocode = "_7";

%%% LHE
[lhe_sev_front] = average_within_region(lhe_sev, fcode, "rad40");
[lhe_sev_occip] = average_within_region(lhe_sev, ocode, "rad40");

%%% GFAP
[gfap_sev_front] = average_within_region(gfap_sev, fcode, "rad40");
[gfap_sev_occip] = average_within_region(gfap_sev, ocode, "rad40");

%%% CD68
[cd68_sev_front] = average_within_region(cd68_sev, fcode, "rad40");
[cd68_sev_occip] = average_within_region(cd68_sev, ocode, "rad40");

%% AVERAGE within subjects -- scatterplot of histology vs. distance
% Create figure for each subject
% This is used for control, mild, moderate
close all;

%%% LHE - extract regions
[lhe_sm] = iterate_subject_radii(lhe, "rad40");
lhe_ctl_front = lhe_sm.CAA_26_1;
lhe_ctl_occip = lhe_sm.CAA_26_7;
lhe_mld_front = lhe_sm.CAA_6_1;
lhe_mld_occip = lhe_sm.CAA_6_7;

%%% GFAP - extract regions
[gfap_sm] = iterate_subject_radii(gfap, "rad40");
gfap_ctl_front = gfap_sm.CAA26_1_GFAP_10x;
gfap_ctl_occip = gfap_sm.CAA26_7_GFAP_10x;
gfap_mld_front = gfap_sm.CAA6_1_GFAP_10x;
gfap_mld_occip = gfap_sm.CAA6_7_GFAP_10x;

%%% CD68 - extract regions
[cd68_sm] = iterate_subject_radii(cd68, "rad40");
cd68_ctl_front = cd68_sm.CAA26_1_CD68_10x;
cd68_ctl_occip = cd68_sm.CAA26_7_CD68_10x;
cd68_mld_occip = cd68_sm.CAA6_7_CD68_10x;

%% LHE - Calculate + Plot differences
% Y-axis limits
ylims = [-0.5, 0.5];
ytick = -0.5 : 0.1 : 0.5;

% Mild Frontal
tstr = 'LHE Frontal (Mild - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(lhe_ctl_front, lhe_mld_front);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'lhe_frontal_mild-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Mild Occipital
tstr = 'LHE Occipital (Mild - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(lhe_ctl_occip, lhe_mld_occip);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'lhe_occipital_mild-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Severe Frontal
tstr = 'LHE Frontal (Severe - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(lhe_ctl_front, lhe_sev_front);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'lhe_frontal_severe-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Severe Occipital
tstr = 'LHE Occipital (Severe - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(lhe_ctl_occip, lhe_sev_occip);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'lhe_occipital_severe-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

%% GFAP - Calculate + Plot differences
% Y-axis limits
ylims = [-1.2, 1.2];
ytick = -1.2 : 0.2 : 1.2;

% Mild Frontal
tstr = 'GFAP Frontal (Mild - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(gfap_ctl_front, gfap_mld_front);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'gfap_frontal_mild-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Mild Occipital
tstr = 'GFAP Occipital (Mild - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(gfap_ctl_occip, gfap_mld_occip);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'gfap_occipital_mild-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Severe Frontal
tstr = 'GFAP Frontal (Severe - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(gfap_ctl_front, gfap_sev_front);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'gfap_frontal_severe-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Severe Occipital
tstr = 'GFAP Occipital (Severe - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(gfap_ctl_occip, gfap_sev_occip);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'gfap_occipital_severe-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

%% CD68 - Calculate + Plot differences
% Y-axis limits
ylims = [-1.2, 1.2];
ytick = -1.2 : 0.2 : 1.2;

% Mild Frontal
% tstr = 'CD68 Frontal (Mild - Control)';
% figure('Position',[100,100,1000,1000]);
% [std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(cd68_ctl_front, cd68_mld_front);
% histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
%                     mean_ves, std_ves, ylims, ytick, tstr)
% fname = fullfile(figdir,'cd68_frontal_mild-ctl_ribbon');
% saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
% exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Mild Occipital
tstr = 'CD68 Occipital (Mild - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(cd68_ctl_occip, cd68_mld_occip);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'cd68_occipital_mild-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Severe Frontal
tstr = 'CD68 Frontal (Severe - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(cd68_ctl_front, cd68_sev_front);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'cd68_frontal_severe-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

% Severe Occipital
tstr = 'CD68 Occipital (Severe - Control)';
figure('Position',[100,100,1000,1000]);
[std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(cd68_ctl_occip, cd68_sev_occip);
histo_mean_std_plot(radii_sm, mean_epvs, std_epvs,...
                    mean_ves, std_ves, ylims, ytick, tstr)
fname = fullfile(figdir,'cd68_occipital_severe-ctl_ribbon');
saveas(gcf,fname,'png'); saveas(gcf,fname,'fig');
exportgraphics(gcf,strcat(fname,'.svg'),'Resolution',600);

%% Measure the mean + StdDev within each subject
function [stats] = iterate_subject_radii(stain, rad_size)
% Measure mean + Std within tissue section
% INPUTS
%   stain (struct): stain struct top level
%   nrad (int): number of radii
%   rad_size (string): sub-field string ("rad40" or "rad100")
% OUTPUTS:
%   mean_epvs (vector): average EPVS
%   std_epvs (vector): std EPVS
%   mean_ves (vector): average Ves
%   std_ves (vector): std Ves

%%% Initialize variables
% Measure number of radii measurements for radius type (40 or 100 um)
nrad = length(fields(stain(1).(rad_size)));
% Initialize vectors
n_epvs = zeros(nrad,1);
n_ves = zeros(nrad,1);
mean_epvs = zeros(nrad,1);
std_epvs = zeros(nrad,1);
mean_ves = zeros(nrad,1);
std_ves = zeros(nrad,1);
% Retrieve names of the fields of small radii
rad_name = fields(stain(1).(rad_size));
% subject names of each section in stain
subs = {stain(:).baseName};
% struct for storing each subject's mean + StdDev
stats = struct();

% Iterate subjects
for s = 1:length(subs)
    % Iterate over the radii
    for ii = 1:nrad
        % EPVS
        epvs = stain(s).(rad_size).(rad_name{ii}).exp;
        n_epvs(ii) = numel(epvs);
        % Vessels
        ves = stain(s).(rad_size).(rad_name{ii}).ctl;
        n_ves(ii) = numel(ves);
        % Take mean & std dev within subject
        mean_epvs(ii) = mean(epvs,'omitnan');
        std_epvs(ii) = std(epvs,'omitnan');
        mean_ves(ii) = mean(ves,'omitnan');
        std_ves(ii) = std(ves,'omitnan');
    end
    % Add each subject vector to struct
    stats.(subs{s}).n_epvs = n_epvs;
    stats.(subs{s}).n_ves = n_ves;
    stats.(subs{s}).mean_epvs = mean_epvs;
    stats.(subs{s}).std_epvs = std_epvs;
    stats.(subs{s}).mean_ves = mean_ves;
    stats.(subs{s}).std_ves = std_ves;
end

end

%% Measure the mean + StdDev within each subject
function [stats] = average_within_region(stain, region, rad_size)
% Measure mean + Std within tissue section. Average within region.
% INPUTS
%   stain (struct): stain struct top level
%   region (string): brain region string
%                       front = "_1"
%                       occip = "_7"
%   rad_size (string): radius size "rad40"
% OUTPUTS:
%   stats (struct): structure containing statistics
%       mean_epvs (vector): average EPVS
%       std_epvs (vector): std EPVS
%       mean_ves (vector): average Ves
%       std_ves (vector): std Ves

%%% Find the tissue sections with the respective brain region
% Extract local variable with basename of each section
names = {stain.baseName};
% Find the index of each tissue section matching region string
if endsWith(names{1},"_10x")
    tf = contains(names,strcat(region,"_"));
else
    tf = endsWith(names,region);
end
% Take the stain indices matching the region
stain_reg = stain(tf);
% subject names of each section in stain
subs = {stain_reg(:).baseName};
% Count total number of subjects within region
nsubs = numel(subs);

%%% Initialize variables
% Measure number of radii measurements for radius type (40 or 100 um)
nrad = length(fields(stain_reg(1).(rad_size)));
% Initialize vectors
n_epvs = zeros(nrad,nsubs);
n_ves = zeros(nrad,nsubs);
mean_epvs = zeros(nrad,nsubs);
std_epvs = zeros(nrad,nsubs);
mean_ves = zeros(nrad,nsubs);
std_ves = zeros(nrad,nsubs);
% Retrieve names of the fields of small radii
rad_name = fields(stain_reg(1).(rad_size));
% struct for storing each subject's mean + StdDev
stats = struct();

% Iterate subjects
for s = 1:length(subs)
    % Iterate over the radii
    for ii = 1:nrad
        % EPVS
        epvs = stain_reg(s).(rad_size).(rad_name{ii}).exp;
        n_epvs(ii,s) = numel(epvs);
        % Vessels
        ves = stain_reg(s).(rad_size).(rad_name{ii}).ctl;
        n_ves(ii,s) = numel(ves);
        % Take mean & std dev within subject
        mean_epvs(ii,s) = mean(epvs,'omitnan');
        std_epvs(ii,s) = std(epvs,'omitnan');
        mean_ves(ii,s) = mean(ves,'omitnan');
        std_ves(ii,s) = std(ves,'omitnan');
    end
end

% Average across subjects
stats.n_epvs = n_epvs(:,1) + n_epvs(:,2);
stats.n_ves = n_ves(:,1) + n_ves(:,2);
stats.std_epvs = std(mean_epvs,0,2,'omitnan');
stats.mean_epvs = mean(mean_epvs,2,'omitnan');
stats.std_ves = std(mean_ves,0,2,'omitnan');
stats.mean_ves = mean(mean_ves,2,'omitnan');
end

%% Function for taking the difference between severities
function [std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(ctl, exp)

% EPVS
mean_epvs = exp.mean_epvs - ctl.mean_epvs;
std_epvs = sqrt((exp.std_epvs).^2./exp.n_epvs +...
                (ctl.std_epvs).^2./ctl.n_epvs );

% Vessel
mean_ves = exp.mean_ves - ctl.mean_ves;
std_ves = sqrt((exp.std_ves).^2./exp.n_ves +...
               (ctl.std_ves).^2./ctl.n_ves);

end