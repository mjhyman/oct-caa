%% Convert EPVS to uchar
clear; clc; close all;

%%% Add freesurfer dir to path
% Get parent directory (oct-caa)
currentDir = pwd;
parentDir = fileparts(currentDir);
% Compose path to 'freesurfer'
freesurferDir = fullfile(parentDir, 'freesurfer');
% Add to MATLAB path
addpath(freesurferDir);

%% CAA 6 Frontal
% Import volume
epvs_file = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
        'optical_properties/caa6/front/epvs.nii'];
epvs = MRIread(epvs_file,0,0);

% Write back to NIFTI as uchar
epvs_file = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
        'optical_properties/caa6/front/epvs_uchar.nii.gz'];
MRIwrite(epvs,epvs_file,'uchar',0);

%% CAA 6 Occipital
% Import volume
epvs_file = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
        'optical_properties/caa6/occip/epvs.nii'];
epvs = MRIread(epvs_file,0,0);

% Write back to NIFTI as uchar
epvs_file = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
        'optical_properties/caa6/occip/epvs_uchar.nii.gz'];
MRIwrite(epvs,epvs_file,'uchar',0);