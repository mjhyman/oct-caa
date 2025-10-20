%% Measure the size-weighted proximity (SWP)
% Batch script that runs on individual subjects
% clear; clc; close all;

% Exponent for denominator in SWP
p = 2;
data_dir = '/projectnb/npbssmic/ns/CAA/';
radius = 200;
save_base = '/projectnb/npbssmic/ns/CAA/swp/';

%%% Load CAA6f (small EPVS)
fprintf('Loading data\n')
load('/projectnb/npbssmic/ns/CAA/caa6/caa6.mat');
fprintf('Finished loading data\n')

%%% Assign local variables
epvs = caa6.front.epvs;
mask = caa6.front.mask_wm;
subject_name = 'caa6_test_static';
region = 'front';

%% Calculate EPVS density

try
    fprintf('Starting to compute SWP\n')
    [subsampled_volume, interpolated_volume] = ...
        epvs_density_variable_p_mod(epvs, mask, radius, p);
catch ME
    fprintf('A fatal error occurred:\n%s\n', ME.message);
   for k = 1:length(ME.stack)
       fprintf('In %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
   end
end

% Save results to .MAT and .TIF
fprintf('Saving SWP Results\n')
save_epvs_heatmap(save_base, subject_name, region, ...
    subsampled_volume, interpolated_volume, radius, p);

fprintf('Finished processing %s %s \n', subject_name, region);

%% Compare the two volumes
ddir = '/projectnb/npbssmic/ns/CAA/swp/';
swp_orig = load(fullfile(ddir,['caa6/front/' ...
    'caa6_front_radius_200_exp_2_interpolated_heatmap.mat']));
swp_test = load(fullfile(ddir,['caa6_test/front/' ...
    'caa6_test_front_radius_200_exp_2_interpolated_heatmap.mat']));

% Compare the two .MAT
swp_orig = swp_orig.interpolated_volume;
swp_test = swp_test.interpolated_volume;
A = swp_orig;
B = swp_test;

% Set tolerance
tol = 2;
isEqual = isequal(size(A), size(B)) && all(abs(A(:) - B(:)) < tol);
tf = {'false','true'};
fprintf('The latest revisions produced a similar volume: %s\n',tf{isEqual+1})