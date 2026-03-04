%% Figure 1: compare mus, ret, ori
% All subfigures are made except control retardance
% Control = CAA 26f, top depth
clear; clc; close all;

%% Import CAA 26
ddir = '/projectnb/npbssmic/ns/CAA/caa26/';
% Load the tissue data from the specified directory
caa26 = load(fullfile(ddir, 'caa26.mat'));
caa26 = caa26.caa26;
ret = caa26.front.ret_full;

%% Import Tissue Mask for depth 1 & apply to retardance
% Define path to TIF
mask = fullfile(ddir,"front/","caa26_front_ret_mask_depth1.tif");
% Import TIF
mask = imread(mask);
mask = logical(mask);
% transpose mask to match dimensions
figure; imagesc(mask);

%%% Apply tissue mask
ret1 = ret(:,:,1);
ret1(~mask) = 0;
figure; imagesc(ret1);

%% Export image as TIF
clc;
% Ensure ret1 is uint32
ret1_u32 = uint32(ret1);

tout = fullfile(ddir,'front','caa26_front_ret_depth1.tif');

t = Tiff(tout,'w');
tagStruct.ImageLength = size(ret1_u32,1);
tagStruct.ImageWidth = size(ret1_u32,2);
tagStruct.Photometric = Tiff.Photometric.MinIsBlack; % grayscale
tagStruct.BitsPerSample = 32;
tagStruct.SamplesPerPixel = 1;
tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagStruct.Compression = Tiff.Compression.None;
t.setTag(tagStruct);
t.write(ret1_u32);
t.close();
