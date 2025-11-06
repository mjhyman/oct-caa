%% Import the .MAT of the optical properties
clear; clc; close all;
% Set maximum number of threads equal to number of threads for script
ncores = feature('numcores');
maxNumCompThreads(ncores);
% Directory containing seg, mus, ret, mask, epvs structs
mat_dir = '/projectnb/npbssmic/ns/CAA/';

%% Load matlab structs

fprintf('Loading CAA6\n')
caa6 = load(fullfile(mat_dir,"/caa6/caa6.mat"));
fprintf('Finished loading CAA17\n')

fprintf('Loading CAA17\n')
caa17 = load(fullfile(mat_dir,"/caa17/occip/caa17.mat"));
fprintf('Finished loading CAA17\n')

fprintf('Loading CAA22\n')
caa22 = load(fullfile(mat_dir,"/caa22/caa22.mat"));
fprintf('Finished loading CAA22\n')

fprintf('Loading CAA25\n')
caa25 = load(fullfile(mat_dir,"/caa25/caa25.mat"));
fprintf('Finished loading CAA25\n')

fprintf('Loading CAA26\n')
caa26 = load(fullfile(mat_dir,"/caa26/caa26_archive.mat"));
fprintf('Finished loading CAA26\n')

% Remove top-level struct
caa6 = caa6.caa6;
caa17 = caa17.caa17;
caa22 = caa22.caa22;
caa25 = caa25.caa25;
caa26 = caa26.caa26;

% Print out which subjects have been loaded
