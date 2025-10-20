%% Measure the size-weighted proximity (SWP)
% This script will test iterating the search radius

% clear; clc; close all;

%%% Parameters for SWP calculations
% Exponent for denominator in SWP
p = 2;
% Initial search radius
radius = 200;
% increment for increasing search radius
radius_inc = 4;
data_dir = '/projectnb/npbssmic/ns/CAA/';
save_base = '/projectnb/npbssmic/ns/CAA/swp/';

%%% Load CAA6f (small EPVS)
fprintf('Loading data\n')
load('/projectnb/npbssmic/ns/CAA/caa6/caa6.mat');
fprintf('Finished loading data\n')

%%% Assign local variables
epvs = caa6.front.epvs;
mask = caa6.front.mask_wm;
subject_name = 'caa6_var_rad';
region = 'front';

%% Calculate EPVS density

fprintf('Starting to compute SWP\n')
[subsampled_volume, interpolated_volume] = ...
    swp_var_radius(epvs, mask, radius, radius_inc, p);

% Save results to .MAT and .TIF
fprintf('Saving SWP Results\n')
save_epvs_heatmap(save_base, subject_name, region, ...
    subsampled_volume, interpolated_volume, radius, p);

fprintf('Finished processing %s %s \n', subject_name, region);