%% Compare the generalized additive models between severities
% Import GAM model results for epvs-swp and ves-swp
% Compare the slive curves
%   These are where the epvs-swp is held constant at a single value while
%   ves-swp is swept. The slices are discrete values, so are not entire
%   story. Figure for each pair of comparisons (ie mild vs severe)
% Create surface difference matrix
%   Displays where joint ves-swp/epvs-swp diverge by severity

%% Initialization
clc; close all;
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
ncores = feature('numcores');
maxNumCompThreads(ncores);

%%% Directories
% Figure output directory (also where the GAM is stored)
swp_dir = '/projectnb/npbssmic/ns/CAA/figures/swp_gam_gmm/';
% Output directory for GAM comparison
fig_dir = '/projectnb/npbssmic/ns/CAA/figures/swp_gam_gmm/gam_comparison/';

%%% Import GAM
% Filename of GAM struct
fname = 'GAM_struct_subjects.mat';
% Load the GAM struct from the specified file
fprintf('\n----Loading GAM struct----\n')
load(fullfile(swp_dir, fname));

%% Map subjectID to severities
% CAA26 = control
% CAA6 = mild
% CAA17o = moderate
% CAA22 = severe

% Create struct for front and occip
gam_front = struct();
gam_occip = struct();

% Assign front structs
gam_front.ctl = gam.pdif.caa26.front;
gam_front.mld = gam.pdif.caa6.front;
gam_front.sev1 = gam.pdif.caa22.front;
gam_front.sev2 = gam.pdif.caa25.front;

% Assign occip structs
gam_occip.ctl = gam.pdif.caa26.occip;
gam_occip.mld = gam.pdif.caa6.occip;
gam_occip.mod = gam.pdif.caa17.occip;
gam_occip.sev1 = gam.pdif.caa22.occip;
gam_occip.sev2 = gam.pdif.caa25.occip;

% Create severity cell array for front + occip
front_sevs = {'ctl','mld','sev1','sev2'};
occip_sevs = {'ctl','mld','mod','sev1','sev2'};

%% Compare slices
fprintf('\n----Running GAM slice comparison----\n')
% Create output directory
dirout = fullfile(fig_dir,'slices');
% Frontal
fprintf('\tGAM slice comparison for frontal\n')
swp_compare_gam_slices(gam_front, front_sevs, fullfile(dirout,'front'))
% Occip
fprintf('\tGAM slice comparison for occip\n')
swp_compare_gam_slices(gam_occip, occip_sevs, fullfile(dirout,'occip'))

%% Compare joint distributions b/w groups
fprintf('\n----Running GAM curve comparison----\n')
% Create output directory
dirout = fullfile(fig_dir,'curves');
% Frontal
fprintf('\tGAM curve comparison for frontal\n')
swp_compare_gam_curves(gam_front, front_sevs, fullfile(dirout,'front'))
% Occip
fprintf('\tGAM curve comparison for occip\n')
swp_compare_gam_curves(gam_occip, occip_sevs, fullfile(dirout,'occip'))