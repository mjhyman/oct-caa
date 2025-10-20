%% This script converts the orientation from degrees to an RGB colorwheel
% The orientation in RGB is used for the figures

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
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
ncores = feature('numcores');
maxNumCompThreads(ncores);
% Flag for loading CAA structs (false if already in environment)
flag_load_caa_structs = true;
% Directory containing seg, mus, ret, mask, epvs structs
mat_dir = '/projectnb/npbssmic/ns/CAA/';

%% Load matlab structs

if flag_load_caa_structs
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
    caa26 = load(fullfile(mat_dir,"/caa26/caa26.mat"));
    fprintf('Finished loading CAA26\n')
    
    % Remove top-level struct
    caa6 = caa6.caa6;
    caa17 = caa17.caa17;
    caa22 = caa22.caa22;
    caa25 = caa25.caa25;
    caa26 = caa26.caa26;
end

%% Visual QA

%%% CAA 6
% Frontal
rgb = caa6.front.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 6 Front')
% Occipital
rgb = caa6.occip.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 6 Occip')

%%% CAA 17
% Occipital
rgb = caa17.occip.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 17 Occip')

%%% CAA 22
% Frontal
rgb = caa22.front.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 22 Front')
% Occipital
rgb = caa22.occip.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 22 Occip')

%%% CAA 25
% Frontal
rgb = caa25.front.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 25 Front')
% Occipital
rgb = caa25.occip.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 25 Occip')

%%% CAA 26
% Frontal
rgb = caa26.front.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 26 Front')
% Occipital
rgb = caa26.occip.orient_rgb;
figure; imagesc(squeeze(rgb(:,:,10,:)));
title('CAA 26 Occip')
