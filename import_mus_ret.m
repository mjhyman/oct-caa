%% Import the scattering and retardance volumes for the CAA dataset.
%{
Outline of script:
- import scattering and retardance for each subject
- combine into single data struct
- save struct to: /autofs/cluster/octdata2/users/mjhyman/oct-caa/

Missing Data:
  - CAA 17 frontal

To Do for Each volume:
- Import the tissue mask
- Import segmentation
    - even the ones that were for subsets of the stack
- Scale scattering, retardance, segmentation to tissue mask
    - Create struct: mus, ret, seg, mask
    - Save struct for each subject to individual .MAT structs
%}
%% Prepare environment
% clear; clc; close all;
% Add top-level directory
d = pwd;
addpath(fullfile(pwd));
addpath('/autofs/cluster/octdata2/users/mjhyman/oct-caa/freesurfer');
% Output directory for optical properties
fname = 'optical_properties.mat';
fout_base = ['/autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses/' ...
    'optical_properties'];

%% Parameters for scattering and retardance
%%% Initialize directories
% Top-level directories
base_oct = '/autofs/cluster/octdata2/users/Chao/caa/';
base_omega = '/autofs/space/omega_001/users/caa/';
% Struct to store directory paths
fpaths = struct();
% Filename extension
fext = '.mat';

%%% MRIread parameters
% header only flag
h_flag = 0;
% permute x,y dimensions flag
p_flag = 0;

%%% Voxel dimensions (microns) for all runs
vox_z = 20; 

%%% properties to save MRI
save_nifti_flag = 1; % flag for saving NIFTI of each run
pflag = 0; % permute flag
mus_res = [0.02,0.02,0.02]; % resolution in mm
ret_res = [0.01,0.01,0.1]; % resolution in mm
dtype = 'float'; % float is the equivalent of single

%% CAA 6 front:
%{
%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa6/' ...
    'frontal/derivatives/predictions/caa6-frontal_unet-tissuemask.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa6/' ...
    'frontal/derivatives/predictions/' ...
    'sub-caa6_sample-frontal_chunk-01_label-vessels_mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% Import Scattering Coefficient Maps
fprintf('Starting CAA6 Frontal mus\n')
% Slice thicness (microns) for all runs
slice_th = 100;
% Filepaths
run1_mus = fullfile(base_oct,'caa_6/frontal/process_run1/mus');
run2_mus = fullfile(base_oct,'caa_6/frontal/process_run2/mus');
run3_mus = fullfile(base_oct,'caa_6/frontal/process_run3/mus');
% scattering filename substring
mus_substr = '_mean_xy_20um_iso';
% Import runs
run1 = combine_slices(run1_mus,mus_substr,fext,vox_z,slice_th);
run2 = combine_slices(run2_mus,mus_substr,fext,vox_z,slice_th);
run3 = combine_slices(run3_mus,mus_substr,fext,vox_z,slice_th);
fprintf('Finished CAA6 Frontal mus\n')

%%% mus - pad successive runs to match size b/w runs & create full volume
% Pad run2 w/ zeros
[run2] = pad_runs(run1,run2);
% Pad run3 w/ zeros
[run3] = pad_runs(run1,run3);

%%% Permute mus mask to align with mask (gyri on top)
run1 = permute(run1,[2,1,3]);
run2 = permute(run2,[2,1,3]);
run3 = permute(run3,[2,1,3]);

%%% Align run3 of mus
% Set number of voxels to offset for rows and columns
rdif = 25;
cdif = 32;
% Remove 25 voxels from top rows
run3_mod = run3(rdif+1:end,:,:);
% Add 25 voxels (zeros) to bottom rows
zero_rows = zeros(rdif,size(run3,2),size(run3,3));
run3_mod = cat(1,run3_mod,zero_rows);
% Remove 25 voxels from right columns
run3_mod = run3_mod(:,1:end-cdif,:);
% Add 25 voxels to left columns
zero_col = zeros(size(run3,1),cdif,size(run3,3));
run3_mod = cat(2,zero_col,run3_mod);
% Combine runs 1,2,3
mus = cat(3,run1, run2, run3_mod);

%%% Import retardance
fprintf('Starting CAA6 Frontal retardance\n')
ret = ['/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/StackNii/'...
    'Stacked_Retardance.nii'];
ret = MRIread(ret,h_flag,p_flag);
ret = single(ret.vol);
fprintf('Imported CAA6 Frontal retardance\n')

%%% Scale retardance
% # rows and columns for rescaling
nrow = size(mask,1);
ncol = size(mask,2);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Interpolate retardance in z-axis to match dimensions of mus
% measure dimensions of full volume
full_dims = size(mask);
% Expand retardance to match dimensions fo full volume
[ret_full] = expand_ret(ret,full_dims);

%%% Add aligned mus, ret, seg, mask to struct
% mus individual runs
caa6.front.mus_run1 = run1;
caa6.front.mus_run2 = run2;
caa6.front.mus_run3 = run3;
% mus combined stack
caa6.front.mus = mus;
% Add retardance
caa6.front.ret = ret;
% Add full retardance
caa6.front.ret_full = ret_full;
% Add segmentation
caa6.front.seg = seg;
% Add mask
caa6.front.mask = mask;

%%% Save the aligned files to NIFTI
if save_nifti_flag
    % Save mus runs
    run1_fout = fullfile(fout_base,'/caa6/front/mus_run1.nii');
    run2_fout = fullfile(fout_base,'/caa6/front/mus_run2.nii');
    run3_fout = fullfile(fout_base,'/caa6/front/mus_run3.nii');
    save_mri(run1,run1_fout,mus_res,dtype,pflag);
    save_mri(run2,run2_fout,mus_res,dtype,pflag);
    save_mri(run3,run3_fout,mus_res,dtype,pflag);
    
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa6/front/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned retardance
    ret_path = fullfile(fout_base,'/caa6/front/ret.nii');
    save_mri(ret,ret_path,ret_res,dtype,pflag);

    % Save aligned retardance (full volume)
    ret_path = fullfile(fout_base,'/caa6/front/ret_full.nii');
    save_mri(ret_full,ret_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa6/front/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa6/front/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);
end

clear run1 run2 run3 run3_mod

%% CAA 6 Occip:

%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa6/occipital'...
    '/derivatives/predictions/caa6-occipital_unet-tissuemask.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa6/occipital/'...
    'derivatives/predictions/' ...
    'sub-caa6_sample-occipital_chunk-01_label-vessels_mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% Import Scattering
fprintf('Starting CAA6 Frontal Scattering\n')
mus = fullfile(base_oct,['caa_6/occipital/process_20220209_run2/mus/' ...
    'mus_mean_20um-iso.nii']);
mus = MRIread(mus,h_flag,p_flag);
mus = single(mus.vol);
fprintf('Finished CAA6 Frontal Scattering\n')

%%% Import Retardance
fprintf('Starting CAA6 Frontal Retardance\n')
ret = fullfile(base_oct,['caa_6/occipital/process_20220209_run2/' ...
    'StackNii/Stacked_Retardance.nii']);
ret = MRIread(ret,h_flag,p_flag);
ret = single(ret.vol);
fprintf('Finished CAA6 Frontal Retardance\n')

%%% Rescale retardance with mask
% # rows and columns for rescaling
nrow = size(mask,1);
ncol = size(mask,2);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Interpolate retardance in z-axis to match dimensions of mus
% measure dimensions of full volume
full_dims = size(mask);
% Expand retardance to match dimensions fo full volume
[ret_full] = expand_ret(ret,full_dims);

%%% Add to data struct and save
caa6.occip.mus = mus;
caa6.occip.ret = ret;
caa6.occip.ret_full = ret_full;
caa6.occip.seg = seg;
caa6.occip.mask = mask;
save(fullfile(fout_base,'/caa6/caa6.mat'),'caa6','-v7.3');

%%% Save the aligned files to NIFTI
if save_nifti_flag   
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa6/occip/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned retardance
    ret_path = fullfile(fout_base,'/caa6/occip/ret.nii');
    save_mri(ret,ret_path,ret_res,dtype,pflag);

    % Save aligned full retardance
    ret_path = fullfile(fout_base,'/caa6/occip/ret_full.nii');
    save_mri(ret_full,ret_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa6/occip/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa6/occip/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);
end

%% CAA 17 occip

%%% Import EPVS
epvs = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/reg/' ...
    'EPVS/segmentation/segmentation_07072023.nii'];
epvs = MRIread(epvs,h_flag,p_flag);
epvs = imbinarize(epvs.vol);

%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa17/' ...
    'occipital/derivatives/predictions/caa17_occipital_unet-mask.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% Import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa17/' ...
    'occipital/derivatives/predictions/' ...
    'sub-caa17_sample-occipital_chunk-01_label-vessels-mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% Import Scattering
% Import mus nifti
mus = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/' ...
    'occipital/process_run1/mus/mus_mean_20um-iso.nii'];
mus = MRIread(mus,h_flag,p_flag);
mus = single(mus.vol);

%%% Import retardance
fprintf('Starting CAA17 Occip retardance\n')
ret = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/' ...
    'process_run1/StackNii/Stacked_Retardance.nii'];
ret = MRIread(ret,h_flag,p_flag);
ret = single(ret.vol);
fprintf('Finished CAA17 Occip retardance\n')

%%% Align mus, retardance, segmentation to mask
% Align scattering to mask
mus = permute(mus,[2,1,3]);
% # rows and columns for rescaling
nrow = size(mask,1);
ncol = size(mask,2);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Interpolate retardance in z-axis to match dimensions of mus
% measure dimensions of full volume
full_dims = size(mask);
% Expand retardance to match dimensions fo full volume
[ret_full] = expand_ret(ret,full_dims);

%%% Add to data struct and save
caa17.occip.mus = mus;
caa17.occip.ret = ret;
caa17.occip.ret_full = ret_full;
caa17.occip.seg = seg;
caa17.occip.mask = mask;
caa17.occip.epvs = epvs;
save(fullfile(fout_base,'/caa17/occip/caa17.mat'),'caa17','-v7.3');

%%% Save the aligned files to NIFTI
if save_nifti_flag   
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa17/occip/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned retardance
    ret_path = fullfile(fout_base,'/caa17/occip/ret.nii');
    save_mri(ret,ret_path,ret_res,dtype,pflag);

    % Save aligned full retardance
    ret_path = fullfile(fout_base,'/caa17/occip/ret_full.nii');
    save_mri(ret_full,ret_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa17/occip/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa17/occip/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);

    % Save EPVS
    epvs_path = fullfile(fout_base,'/caa17/occip/epvs.nii');
    save_mri(epvs,epvs_path,mus_res,dtype,pflag);
end

%% CAA 22 front

%%% import EPVS
epvs = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/reg/fron/' ...
    'EPVS/EPVS_segmentation_03262024.mgz'];
epvs = MRIread(epvs,h_flag,p_flag);
epvs = imbinarize(epvs.vol);

%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa22/' ...
    'frontal/derivatives/predictions/caa22-frontal_unet-tissuemask.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% Import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa22/' ...
    'frontal/derivatives/predictions/' ...
    'sub-caa22_sample-frontal_chunk-01_label-vessels-mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% Combine mus 
fprintf('Importing CAA22 Frontal mus\n')
mus = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/' ...
    '20211018/Z_Stitched/mus_mean_20um-iso.nii'];
mus = MRIread(mus,h_flag,p_flag);
mus = single(mus.vol);
fprintf('Finished CAA22 Frontal mus\n')

%%% Import retardance
fprintf('Starting CAA22 Frontal ret\n')
ret = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/' ...
        '20211018/StackNii/Stacked_Retardance.nii'];
ret = MRIread(ret,h_flag,p_flag);
ret = single(ret.vol);
fprintf('Finished CAA22 Frontal ret\n')

%%% Align mus, retardance, segmentation to mask
% # rows and columns for rescaling
nrow = size(mask,1);
ncol = size(mask,2);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Interpolate retardance in z-axis to match dimensions of mus
% measure dimensions of full volume
full_dims = size(mask);
% Expand retardance to match dimensions fo full volume
[ret_full] = expand_ret(ret,full_dims);

%%% Add to data struct and save
caa22.front.mus = mus;
caa22.front.ret = ret;
caa22.front.ret_full = ret_full;
caa22.front.seg = seg;
caa22.front.mask = mask;
caa22.front.epvs = epvs;

%%% Save the aligned files to NIFTI
if save_nifti_flag   
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa22/front/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned retardance
    ret_path = fullfile(fout_base,'/caa22/front/ret.nii');
    save_mri(ret,ret_path,ret_res,dtype,pflag);

    % Save aligned full retardance
    ret_path = fullfile(fout_base,'/caa22/front/ret_full.nii');
    save_mri(ret_full,ret_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa22/front/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa22/front/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);

    % Save EPVS
    epvs_path = fullfile(fout_base,'/caa22/front/epvs.nii');
    save_mri(epvs,epvs_path,mus_res,dtype,pflag);
end

%% CAA 22 occip:

%%% import EPVS
epvs = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/reg/occi/' ...
    'EPVS_segmentation/EPVS_segmentation_02132024_registered.mgz'];
epvs = MRIread(epvs,h_flag,p_flag);
epvs = imbinarize(epvs.vol);

%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa22/' ...
    'occipital/derivatives/predictions/' ...
    'caa22-occipital_unet-tissuemask.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% Import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa22/occipital/'...
    'derivatives/predictions/' ...
    'sub-caa22_sample-occipital_chunk-01_label-vessels-mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% Import scattering
fprintf('Starting CAA22 Occip mus\n')
mus = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/reg/occi/' ...
    'mus_stitch_mean_20um-iso.orient.man.rb_CC.nii'];
mus = MRIread(mus,h_flag,p_flag);
mus = single(mus.vol);
fprintf('Finished CAA22 Occip mus\n')

%%% Import Retardance
fprintf('Starting CAA22 Occip Retardance\n')
ret = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/' ...
    'StackNii/Stacked_Retardance.nii'];
ret = MRIread(ret,h_flag,p_flag);
ret = single(ret.vol);
fprintf('Finished CAA22 Occip Retardance\n')
% Rescale retardanceto mask
nrow = size(mask,1);
ncol = size(mask,2);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Offset retardance
% Add columns to left and rows to bottom
% Set number of voxels to offset for rows and columns
rdif = 46;
cdif = 11;
% Remove redif voxels from top rows
ret = ret(rdif+1:end,:,:);
% Add redif voxels (zeros) to bottom rows
zero_rows = zeros(rdif,size(ret,2),size(ret,3));
ret = cat(1,ret,zero_rows);
% Remove cdif voxels from right columns
ret = ret(:,1:end-cdif,:);
% Add cdif voxels to left columns
zero_col = zeros(size(ret,1),cdif,size(ret,3));
ret = cat(2,zero_col,ret);
% Display retardance
figure; imagesc(ret(:,:,1)); title('retardance - 1');

%%% Interpolate retardance in z-axis to match dimensions of mus
% measure dimensions of full volume
full_dims = size(mask);
% Expand retardance to match dimensions fo full volume
[ret_full] = expand_ret(ret,full_dims);

%%% Add to data struct and save
caa22.occip.mus = mus;
caa22.occip.ret = ret;
caa22.occip.ret_full = ret_full;
caa22.occip.seg = seg;
caa22.occip.mask = mask;
caa22.occip.epvs = epvs;

%%% Save the aligned files to NIFTI
if save_nifti_flag   
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa22/occip/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned retardance
    ret_path = fullfile(fout_base,'/caa22/occip/ret.nii');
    save_mri(ret,ret_path,ret_res,dtype,pflag);

    % Save aligned full retardance
    ret_path = fullfile(fout_base,'/caa22/occip/ret_full.nii');
    save_mri(ret_full,ret_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa22/occip/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa22/occip/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);

    % Save epvs
    epvs_path = fullfile(fout_base,'/caa22/occip/epvs.nii');
    save_mri(epvs,epvs_path,mus_res,dtype,pflag);
end

%%% Save struct
save(fullfile(fout_base,'/caa22/caa22.mat'),'caa22','-v7.3');

%% CAA 25 front:

%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa25/' ...
    'frontal/derivatives/predictions/caa25-frontal_unet-tissuemask.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% Import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa25/' ...
    'frontal/derivatives/predictions/' ...
    'sub-caa25_sample-frontal_chunk-01_label-vessels-mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% Import scattering
fprintf('Starting CAA25 Front mus\n')
mus = ['/autofs/space/omega_001/users/caa/CAA25_Frontal/' ...
    'Process_caa25_frontal/mus/mus_mean_20um-iso.nii'];
mus = MRIread(mus,h_flag,p_flag);
mus = single(mus.vol);
fprintf('Finished CAA25 Front mus\n')

%%% Add retardance
fprintf('Starting CAA25 Front Retardance\n')
ret = ['/autofs/space/omega_001/users/caa/CAA25_Frontal/' ...
    'Process_caa25_frontal/StackNii/Stacked_Retardance.nii'];
ret = MRIread(ret,h_flag,p_flag);
ret = single(ret.vol);

%%% Align mus, retardance, segmentation to mask
% Align scattering to mask
mus = permute(mus,[2,1,3]);
% # rows and columns for rescaling
nrow = size(mask,1);
ncol = size(mask,2);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Interpolate retardance in z-axis to match dimensions of mus
% measure dimensions of full volume
full_dims = size(mask);
% Expand retardance to match dimensions fo full volume
[ret_full] = expand_ret(ret,full_dims);

%%% Add to data struct and save
caa25.front.mus = mus;
caa25.front.ret = ret;
caa25.front.ret_full = ret_full;
caa25.front.seg = seg;
caa25.front.mask = mask;
fprintf('Finished CAA25 Front Retardance\n')

%%% Save the aligned files to NIFTI
if save_nifti_flag   
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa25/front/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned retardance
    ret_path = fullfile(fout_base,'/caa25/front/ret.nii');
    save_mri(ret,ret_path,ret_res,dtype,pflag);

    % Save aligned full retardance
    ret_path = fullfile(fout_base,'/caa25/front/ret_full.nii');
    save_mri(ret_full,ret_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa25/front/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa25/front/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);
end

%% CAA 25 occip

%%% Import EPVS
epvs = ['/autofs/space/omega_001/users/caa/CAA25_Occipital/reg/occi/' ...
    'EPVS/EPVS_mus_segmentation_02052024.mgz'];
epvs = MRIread(epvs,h_flag,p_flag);
epvs = imbinarize(epvs.vol);

%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/' ...
    'caa25/occipital/derivatives/predictions/' ...
    'caa25-occipital_unet-tissuemask.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% Import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa25/occipital/' ...
    'derivatives/predictions/' ...
    'sub-caa25_sample-occipital_chunk-01_label-vessels-mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% import scattering
% Filenames of scattering
mus_substr = '_mean_xy_20um_iso';
% Slice thicness (microns) for all runs
slice_th = 100;
% Filepaths
run1 = fullfile(base_omega,'/CAA25_Occipital/process_caa25_occipital_run1/mus/');
run2 = fullfile(base_omega,'/CAA25_Occipital/process_caa25_occipital_run3/mus/');
run3 = fullfile(base_omega,'/CAA25_Occipital/process_caa25_occipital_run4/mus/');
fprintf('Starting CAA25 Occip mus\n')
run1 = combine_slices(run1,mus_substr,fext,vox_z,slice_th);
run2 = combine_slices(run2,mus_substr,fext,vox_z,slice_th);
run3 = combine_slices(run3,mus_substr,fext,vox_z,slice_th);

%%% Shift run2 and run3 mus down by 30
% Add rows to bottom
rdif = 30;
% Remove rdif voxels from bottom rows
run2_d = run2(rdif+1:end,:,:);
run3_d = run3(rdif+1:end,:,:);
% Add redif voxels (zeros) to bottom rows
run2_zrows = zeros(rdif,size(run2,2),size(run2,3));
run3_zrows = zeros(rdif,size(run3,2),size(run3,3));
run2_d = cat(1,run2_d,run2_zrows);
run3_d = cat(1,run3_d,run3_zrows);
% Display shifted mus
figure; imagesc(run1(:,:,end)); title('mus - run1 - end'); clim([0,50]);
figure; imagesc(run2(:,:,1)); title('mus - run2 - 1'); clim([0,50]);
figure; imagesc(run2_mus_stack(:,:,1)); title('mus - run2 - 1'); clim([0,50]);
figure; imagesc(run2_d(:,:,end)); title('mus - run2 - end'); clim([0,50]);
figure; imagesc(run3_d(:,:,1)); title('mus - run3 - 1'); clim([0,50]);
figure; imagesc(ret(:,:,1)); title('ret - 1'); clim([0,100]);
% concatenate runs
mus = cat(3,run1,run2_d,run3_d);
% transpose scattering to match mask
mus = permute(mus,[2,1,3]);

%%% Retardance
% File paths / names
run1 = fullfile(base_omega,['/CAA25_Occipital/process_caa25_occipital_run1/' ...
    'StackNii/Stacked_Retardance.nii']);
run2 = fullfile(base_omega,['/CAA25_Occipital/process_caa25_occipital_run3/' ...
    'StackNii/Stacked_Retardance.nii']);
run3 = fullfile(base_omega,['/CAA25_Occipital/process_caa25_occipital_run4/' ...
    'StackNii/Stacked_Retardance.nii']);
% Import retardance runs
fprintf('Starting CAA25 Occip ret\n')
run1 = MRIread(run1,h_flag,p_flag); run1 = single(run1.vol);
run2 = MRIread(run2,h_flag,p_flag); run2 = single(run2.vol);
run3 = MRIread(run3,h_flag,p_flag); run3 = single(run3.vol);
% Concatenate retardance runs
ret = cat(3,run1,run2,run3);
fprintf('Finished CAA25 Occip ret\n')

%%% Rescale retardance
% # rows and columns for rescaling
nrow = size(mask,1);
ncol = size(mask,2);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Shift retardance by 30 voxels
% The retardance depths 1-9 are aligned
% The retardance afterwards are shifted 30 voxels up
ret1 = ret(:,:,1:10);
ret2 = ret(:,:,11:end);
% variable of rows to add to bottom
rdif = 30;
% Remove rdif voxels from bottom rows
ret2 = ret2(:,rdif+1:end,:);
% Add redif voxels (zeros) to bottom rows
zcol = zeros(size(ret2,1),rdif,size(ret2,3));
ret2_d = cat(2,ret2,zcol);
% concatenate runs
ret = cat(3,ret1,ret2_d);

%%% Interpolate retardance in z-axis to match dimensions of mus
% measure dimensions of full volume
full_dims = size(mask);
% Expand retardance to match dimensions fo full volume
[ret_full] = expand_ret(ret,full_dims);

% Save to NIFTI
if save_nifti_flag   
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa25/occip/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned retardance
    ret_path = fullfile(fout_base,'/caa25/occip/ret.nii');
    save_mri(ret,ret_path,ret_res,dtype,pflag);

    % Save aligned full retardance
    ret_path = fullfile(fout_base,'/caa25/occip/ret_full.nii');
    save_mri(ret_full,ret_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa25/occip/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa25/occip/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);

    % Save epvs
    epvs_path = fullfile(fout_base,'/caa25/occip/epvs.nii');
    save_mri(epvs,epvs_path,mus_res,dtype,pflag);
end

%%% Add to data struct and save
caa25.occip.mus = mus;
caa25.occip.ret = ret;
caa25.occip.ret_full = ret_full;
caa25.occip.seg = seg;
caa25.occip.mask = mask;
caa25.occip.epvs = epvs;
save(fullfile(fout_base,'/caa25/caa25.mat'),'caa25','-v7.3');

%% CAA 26 front:
% Caroline manually aligned the retardance and then saved it as:
% /autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses/
% optical_properties/caa26/front/ret_reg2mus.nii
% This file contains the full retardance stack (interpolated in z-axis to
% match dimensions of mask, mus).

%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa26/frontal/' ...
    'derivatives/predictions/caa26-frontal_unet-tissuemask.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% Import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa26/frontal/' ...
    'derivatives/predictions/' ...
    'sub-caa26_sample-frontal_chunk-01_label-vessels-mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% Scattering filename and proerties
% Filenames of scattering
mus_substr = '_mean_xy_20um_iso';
% Slice thicness (microns) for all runs
slice_th = 100;
% Filepaths
run1 = fullfile(base_omega,'/CAA26_Frontal/Process_caa26_frontal_run1/mus/');
run2 = fullfile(base_omega,'/CAA26_Frontal/Process_caa26_frontal_run2/mus/');

%%% import scattering
fprintf('Starting CAA26 Front mus\n')
run1 = combine_slices(run1,mus_substr,fext,vox_z,slice_th);
run2 = combine_slices(run2,mus_substr,fext,vox_z,slice_th);
fprintf('Finished CAA26 Front mus\n')

%%% Align Scattering runs 1 and 2      
% Run 1 and 2 are offset in the rows (first dimension). This section of
% code manually aligns the two runs. Run 1 is positioned 138 rows lower
% than run 2.
% Also, run 1 has 1591 columns whereas run 2 has 1452. This code will pad
% run2 with columns of zeros.
% - truncate 138 rows from the top of run 1 
% - add 138 rows to the bottom of run 1
% - add 139 columns to run 2

%%% Align scattering rows b/w runs 1,2
% Manually crop the rows from the top of run 1
row_diff = 138;
run1_minus = run1(row_diff:end-1,:,:);
% Add row_diff to bottom of run 1 to align runs 1&2
run1_pad = zeros(row_diff,size(run1,2),size(run1,3));
run1_pad = cat(1,run1_minus,run1_pad);

%%% Align scattering columns b/w runs 1,2
% Measure column difference
col_diff = size(run1,2) - size(run2,2);
% Create columns of zeros
pad = zeros(size(run2,1),col_diff,size(run2,3));
% Add columns of zeros to the left side of run 2
run2_pad = cat(2,pad,run2);

%%% Combine scattering run1 and run2
% Concatenate 
mus = cat(3, run1_pad, run2_pad);
% Plot overlay of run1 and run2 interface
figure; imshowpair(imadjust(run1_pad(:,:,end)),imadjust(run2_pad(:,:,1)));
title('Run 1 and Run 2 aligned')
% Transpose mus x,y
mus = permute(mus,[2,1,3]);

%%% Shift scattering (mus) by 147 to right (add columns) to match mask
% Manuall set offset
offset = 135;
% Remove zeros to pad scattering
mus_rm = mus(:,1:end-offset,:);
% Create zeros to pad scattering
pad = zeros(size(mus,1),offset,size(mus,3));
% concatenate
mus = cat(2,pad,mus_rm);

%%% Retardance
%}
fprintf('Starting CAA26 Front Retardance\n')
ret = ['/autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses/'...
    'optical_properties/caa26/front/ret_reg2mus.nii'];
ret = MRIread(ret,h_flag,p_flag);
ret = single(ret.vol);
fprintf('Finished CAA26 Front Retardance\n')

%%% Add to struct
caa26.front.mus = mus;
caa26.front.ret_full = ret;
caa26.front.seg = seg;
caa26.front.mask = mask;

%%% Save the aligned files to NIFTI
if save_nifti_flag   
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa26/front/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa26/front/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa26/front/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);
end

%% CAA 26 occip:

%%% Import EPVS
epvs = ['/autofs/space/omega_001/users/caa/CAA26_Occipital/' ...
    'Process_caa26_occipital/reg/EPVS/EPVS_mus_segmentation.mgz'];
epvs = MRIread(epvs,h_flag,p_flag);
epvs = imbinarize(epvs.vol);

%%% import tissue mask
mask = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa26/' ...
    'occipital/sub-caa26_sample-occipital_chunk-01_OCT_TISSUEMASK.nii'];
mask = MRIread(mask,h_flag,p_flag);
mask = imbinarize(mask.vol);

%%% Import vasculature
seg = ['/autofs/cluster/octdata2/users/epc28/data/CAA/caa26/' ...
    'occipital/derivatives/predictions/' ...
    'sub-caa26_sample-occipital_chunk-01_label-vessels-mask.nii'];
seg = MRIread(seg,h_flag,p_flag);
seg = imbinarize(seg.vol);

%%% Import scattering
mus = ['/autofs/space/omega_001/users/caa/CAA26_Occipital/' ...
    'Process_caa26_occipital/mus/mus_mean_20um-iso.nii'];
mus = MRIread(mus,h_flag,p_flag);
mus = single(mus.vol);

%%% Retardance
fprintf('Starting CAA26 Occip Retardance\n')
ret = ['/autofs/space/omega_001/users/caa/CAA26_Occipital/' ...
    'Process_caa26_occipital/StackNii/Stacked_Retardance.nii'];
ret = MRIread(ret,h_flag,p_flag);
ret = single(ret.vol);
fprintf('Finished CAA26 Occip Retardance\n')

%%% Align mus, retardance, segmentation to mask
mus = permute(mus,[2,1,3]);
% # rows and columns for rescaling
nrow = size(mask,1);
ncol = size(mask,2);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Interpolate retardance in z-axis to match dimensions of mus
% measure dimensions of full volume
full_dims = size(mask);
% Expand retardance to match dimensions fo full volume
[ret_full] = expand_ret(ret,full_dims);

%%% Add to struct and save CAA 26
caa26.occip.mus = mus;
caa26.occip.ret = ret;
caa26.occip.ret_full = ret_full;
caa26.occip.seg = seg;
caa26.occip.mask = mask;
caa26.occip.epvs = epvs;
fprintf('Saving CAA26.mat\n')
save(fullfile(fout_base,'/caa26/caa26.mat'),'caa26','-v7.3');
fprintf('Finished saving CAA26.mat\n')

%%% Save the aligned files to NIFTI
if save_nifti_flag 
    fprintf('Saving CAA26-occip niftis\n')
    % Save combined mus runs
    mus_path = fullfile(fout_base,'/caa26/occip/mus.nii');
    save_mri(mus,mus_path,mus_res,dtype,pflag);

    % Save aligned retardance
    ret_path = fullfile(fout_base,'/caa26/occip/ret.nii');
    save_mri(ret,ret_path,ret_res,dtype,pflag);

    % Save aligned full retardance
    ret_path = fullfile(fout_base,'/caa26/occip/ret_full.nii');
    save_mri(ret_full,ret_path,mus_res,dtype,pflag);

    % Save aligned segmentation
    seg_path = fullfile(fout_base,'/caa26/occip/seg.nii');
    save_mri(seg,seg_path,mus_res,dtype,pflag);

    % Save mask
    mask_path = fullfile(fout_base,'/caa26/occip/mask.nii');
    save_mri(mask,mask_path,mus_res,dtype,pflag);

    % Save epvs
    epvs_path = fullfile(fout_base,'/caa26/occip/epvs.nii');
    save_mri(epvs,epvs_path,mus_res,dtype,pflag);
    fprintf('Finished saving CAA26-occip niftis\n')
end

%% Pad runs to match dimensions
function [mv_pad] = pad_runs(fixed,moving)
% PAD_RUNS: add zeros to moving to match dimension of fixed
% INPUTS
%   fixed (matrix): larger image
%   moving (matrix): smaller image

% Pad X-axis
x_pad = abs(size(fixed,2) - size(moving,2));
x_pad = zeros(size(moving,1),x_pad,size(moving,3));
mv_pad = [x_pad, moving];

% Pad Y-axis
y_pad = abs(size(fixed,1) - size(moving,1));
y_pad = zeros(y_pad,size(mv_pad,2),size(moving,3));
mv_pad = [mv_pad; y_pad];

end

%% Import all slices of optical property for each subject
function [stack] = combine_slices(fpath,fpart,fext,zdim,slice_thick)
% COMBINE_SLICES concatenate all slices from a run into a matrix
% INPUTS
%   fpath (string): directory to either scattering or retardance subfolder
%                   containing a 2D array for each slice.
%   fpart (string): part of filename to search for within fpath
%   fext (string): filename extension
%   zdim (array): voxel dimensions (microns)
%   slice_thick (double): slice thickness (microns)
% OUTPUTS
%   stack (matrix): concatenated slices

%%% Extract filenames of slices
% List out directory contents
fext = strcat('*',fext);
files = dir(fullfile(fpath,fext));
% Extract filenames from struct
fnames = {files.name};
fnames = fnames';
% Find files with matching substring (fpart)
tf = contains(fnames,fpart);
% Extract files matching substring
slices = fnames(tf,1);
% Count number of slices for this subject & run
n = size(slices,1);

%%% Calculate # depths (nz) from each slice
% There is an excessive amount of overlap between slices. We only need the
% depth that covers the thickness of the slice. This subsection calculates
% the number of voxels to include from each slice, based upon the voxel
% dimensions and the thickness of the slice.
nz = uint16(slice_thick ./ zdim);

%%% Initialize stack
% Load single slices to find dimensions to initialize matrix
slice = load(fullfile(fpath,slices{1,1}));
slice = slice.Id;
[y,x,~] = size(slice);
% Initialize matrix for storing stack
stack = zeros(y,x,n*nz,'single');

%%% iterate over filenames of each slice and add to stack
z0 = 1;
zf = nz;
for ii = 1:n
    % Retrieve slice
    slice = load(fullfile(fpath,slices{ii,1}));
    slice = slice.Id;
    % Add slice to stack
    stack(:,:,z0:zf) = slice(:,:,1:nz);
    % Iterate nz offset
    z0 = zf + 1;
    zf = zf + nz;
end

end

%% Create stack of retardance maps
function  [ret_full] = expand_ret(ret,full_dims)
% EXPAND_RET copy each retardance mask to recreate a full volume the same
% dimensions as the scattering.
% INPUTS
%   - ret (matrix): retardance matrix
%   - full_dims (matrix): [y,x,z] dimensions of full volume
% OUTPUTS
%   - ret_full (matrix): full retardance matrix

% calculate number of mus depths per retardance slice
z_mus = full_dims(3);
z_ret = size(ret,3);

% Calculate number of depths for each physical slice
n_depth = floor(z_mus ./ z_ret);

% Initialize matrix same size as mask
ret_full = zeros(full_dims);

% Iterate over each retardance map for each section
for ii = 1:size(ret,3)
    % Set the first and last z-axis indices
    z0 = (ii-1)*n_depth + 1;
    zf = ii*n_depth;
    % Expand retardance in z-axis
    ret_full(:,:,z0:zf) = repmat(ret(:,:,ii),[1,1,n_depth]);
end

end