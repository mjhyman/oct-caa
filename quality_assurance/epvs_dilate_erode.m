%% Dilate and erode EPVS
% This is to fill in the missing voxels and gaps in the annotations
clear; clc; close all;

%%% Add freesurfer dir to path
% Get parent directory (oct-caa)
currentDir = pwd;
parentDir = fileparts(currentDir);
% Compose path to 'freesurfer'
freesurferDir = fullfile(parentDir, 'freesurfer');
% Add to MATLAB path
addpath(freesurferDir);

%% CAA 17 Occipital
% Import volume
epvs_file = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
        'optical_properties/caa17/occip/epvs.nii'];
epvs = MRIread(epvs_file,0,0);
vol = epvs.vol;

% Dilate
se = strel('sphere',3);
vol = imdilate(vol,se);

% Erode
se = strel('sphere',2);
vol = imerode(vol,se);

% Write back to NIFTI
epvs.vol = vol;
epvs_file = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
        'optical_properties/caa17/occip/epvs_dilate_erode3_2.nii'];
MRIwrite(epvs,epvs_file,'uchar',0);