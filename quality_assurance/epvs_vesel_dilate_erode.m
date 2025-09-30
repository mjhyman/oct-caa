%% Dilate and erode EPVS
% This is to fill in the missing voxels and gaps in the annotations for all
% volumes.
clear; clc; close all;

%%% Add freesurfer dir to path
% Get parent directory (oct-caa)
currentDir = pwd;
parentDir = fileparts(currentDir);
% Compose path to 'freesurfer'
freesurferDir = fullfile(parentDir, 'freesurfer');
% Add to MATLAB path
addpath(freesurferDir);

%% Cell array of names and directories
topdir = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
            'optical_properties'];
subdirs = {'/caa6/front/','/caa6/occip/',...
        '/caa17/occip/',...
        '/caa22/front/','/caa22/occip/',...
        '/caa25/front/','/caa25/occip/',...
        '/caa26/front/','/caa26/occip/'};

fnames = {'epvs.nii','epvs_uchar.nii.gz',...
    'epvs_dilate_erode5_4.nii',...
    'epvs_manual_latest.nii','epvs.nii',...
    'epvs_updated_MH.mgz','epvs.nii',...
    'epvs.nii','epvs.nii'};

ves_fin = 'seg.nii';

% Dilation/erosion and filenames
ves_r = 4;
ves_fout = 'seg_dilate_erode_4_4.nii.gz';
epvs_r = 4;
epvs_fout = 'epvs_dilate_erode_4_4.nii.gz';

%% Vessel Iterate volumes
for ii = 1:length(subdirs)
    % Retrieve directory and filename
    ffull = fullfile(topdir,subdirs{ii},ves_fin);
    % Import Volume
    ves = MRIread(ffull,0,0);
    vol = logical(ves.vol);
    % Dilate and erode
    vol = dilerode(vol,ves_r);
    % Write back to NIFTI
    ves.vol = logical(vol);
    fout = fullfile(topdir,subdirs{ii},ves_fout);
    MRIwrite(ves,fout,'uchar',0);
end

%% EPVS Iterate volumes
for ii = 1:length(subdirs)
    % Retrieve directory and filename
    ffull = fullfile(topdir,subdirs{ii},fnames{ii});
    % Import Volume
    epvs = MRIread(ffull,0,0);
    vol = logical(epvs.vol);
    % Dilate and erode
    vol = dilerode(vol,epvs_r);
    % Write back to NIFTI
    epvs.vol = logical(vol);
    fout = fullfile(topdir,subdirs{ii},epvs_fout);
    MRIwrite(epvs,fout,'uchar',0);
end

%% Function to dilate and erode
function vol = dilerode(vol,r)
% Create structuring element
se = strel('sphere',r);
% Dilate
vol = imdilate(vol,se);
% Erode
vol = imerode(vol,se);
end