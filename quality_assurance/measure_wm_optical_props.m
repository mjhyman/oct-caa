%% Measure white matter parenchymal Optical Properties
%{
The purpose of this script is to measure the optical properties in the
parenchymal white matter tissue. The average of these values will be
compared against the values in the NIFTI volumes. This is a sanity check.

The full volumes are here:
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)

Outline:
- IMPORT struct containing:
    - tissue masks (exclude background signal)
    - EPVS
    - segmented vasculature (from Etienne)
    - optical properties (mus + retardance)
- MASK
    - Apply WM tissue mask to vessels and optical properties
- Measure average WM value of retardance, scattering, orientation
%}

%% Prepare environment
clear; clc; close all;
% Add top-level directory + subdirectories
addpath(genpath(fullfile(pwd, '..')))
% Directory for loading seg, mus, ret, mask, epvs
data_dir = ['/autofs/cluster/octdata3/users/mjhyman/' ...
    'oct_caa_analyses/optical_properties'];
% WM mask directory (from Taylor)
wm_dir = ['/autofs/space/turtle_001/users/xz875/projects/' ...
          'Multi-resolution_Unets_Semi_OCT/prediction'];
% Voxel dimensions (microns) for all runs
res = [20,20,20]; % resolution in microns

%%% Flag for loading .MAT struct and creating WM masks
% flag for reloading the .MAT struct for each subject
flag_load_caa_structs = true;
% flag for importing updated wm masks
flag_import_wm_mask = true;

%% Load each subject's .MAT struct and create WM mask
%%% Load the .MAT structs
if flag_load_caa_structs
    % CAA 6
    fprintf('Loading CAA6\n')
    caa6 = load(fullfile(data_dir,'/caa6/caa6.mat'));
    caa6 = caa6.caa6;
    fprintf('Finished Loading CAA6\n')
    % CAA 17
    fprintf('Loading CAA17\n')
    caa17 = load(fullfile(data_dir,'/caa17/occip/caa17.mat'));
    caa17 = caa17.caa17;
    fprintf('Finished Loading CAA17\n')
    % CAA 22
    fprintf('Loading CAA22\n')
    caa22 = load(fullfile(data_dir,'/caa22/caa22.mat'));
    caa22 = caa22.caa22;
    fprintf('Finished Loading CAA22\n')
    % CAA 25
    fprintf('Loading CAA25\n')
    caa25 = load(fullfile(data_dir,'/caa25/caa25.mat'));
    caa25 = caa25.caa25;
    fprintf('Finished Loading CAA25\n')
    % CAA 26
    fprintf('Loading CAA26\n')
    caa26 = load(fullfile(data_dir,'/caa26/caa26.mat'));
    caa26 = caa26.caa26;
    fprintf('Finished Loading CAA26\n')
end

if flag_import_wm_mask
    %%% CAA 6 Frontal
    fprintf('Importing caa6 frontal\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa6/front/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa6.front.seg;
    % Apply WM mask to vasculature
    caa6.front.seg_wm = seg .* mask_wm;
    caa6.front.mask_wm = mask_wm;
    
    %%% CAA 6 Occip
    fprintf('Importing caa6 occip\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa6/occip/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa6.occip.seg;
    % Create WM mask and apply to vasculature
    caa6.occip.seg_wm = seg .* mask_wm;
    caa6.occip.mask_wm = mask_wm;
    
    %%% CAA 17 Occip
    fprintf('Importing caa17 occip\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa17/occip/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa17.occip.seg;
    % Apply WM mask tp vasculature
    caa17.occip.seg_wm = seg .* mask_wm;
    caa17.occip.mask_wm = mask_wm;
    
    %%% CAA 22 Frontal
    fprintf('Importing caa22 frontal\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa22/front/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa22.front.seg;
    % Create WM mask and apply to vasculature
    caa22.front.seg_wm = seg .* mask_wm;
    caa22.front.mask_wm = mask_wm;
    
    %%% CAA 22 Occip
    fprintf('Importing caa6 occip\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa22/occip/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa22.occip.seg;
    % Create WM mask and apply to vasculature
    caa22.occip.seg_wm = seg .* mask_wm;
    caa22.occip.mask_wm = mask_wm;
    
    %%% CAA 25 Frontal
    fprintf('Importing caa25 frontal\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa25/front/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa25.front.seg;
    % Create WM mask and apply to vasculature
    caa25.front.seg_wm = seg .* mask_wm;
    caa25.front.mask_wm = mask_wm;
    
    %%% CAA 25 Occip
    fprintf('Importing caa25 occip\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa25/occip/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa25.occip.seg;
    % Create WM mask and apply to vasculature
    caa25.occip.seg_wm = seg .* mask_wm;
    caa25.occip.mask_wm = mask_wm;
    
    %%% CAA 26 Frontal
    fprintf('Importing caa26 frontal\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa26/front/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa26.front.seg;
    % Create WM mask and apply to vasculature
    caa26.front.seg_wm = seg .* mask_wm;
    caa26.front.mask_wm = mask_wm;
    
    %%% CAA 26 Occip
    fprintf('Importing caa26 occip\n')
    % import wm mask
    mask_wm = fullfile(data_dir,'caa26/occip/wm_mask_revised.nii');
    mask_wm = MRIread(mask_wm,0,0);
    % Keep just WM (wm = 1)
    mask_wm = mask_wm.vol;
    mask_wm = logical(mask_wm==1);
    % Flip about the horizontal to align with scattering
    mask_wm = flip(mask_wm,1);
    % import local variables
    seg = caa26.occip.seg;
    % Create WM mask and apply to vasculature
    caa26.occip.seg_wm = seg .* mask_wm;
    caa26.occip.mask_wm = mask_wm;
end

%% Overlay white matter mask + Tissue
% Ensure the coordinates are aligned

%%% CAA 6 - Front
mus = caa6.front.mus;
wm = caa6.front.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')

%%% CAA 6 - Occip
mus = caa6.occip.mus;
wm = caa6.occip.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')

%%% CAA 17 - Occip
mus = caa17.occip.mus;
wm = caa17.occip.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')

%%% CAA 22 - Front
mus = caa22.front.mus;
wm = caa22.front.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')

%%% CAA 22 - Occip
mus = caa22.occip.mus;
wm = caa22.occip.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')

%%% CAA 25 - Front
mus = caa25.front.mus;
wm = caa25.front.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')

%%% CAA 25 - Occip
mus = caa25.occip.mus;
wm = caa25.occip.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')

%%% CAA 26 - Front
mus = caa26.front.mus;
wm = caa26.front.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')

%%% CAA 26 - Occip
mus = caa26.occip.mus;
wm = caa26.occip.mask_wm;
figure; imagesc(wm(:,:,10)); title('White Matter Mask')
figure; imagesc(mus(:,:,10)); title('Scattering')
