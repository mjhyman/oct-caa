%% Analyze the Optical Properties
%{
The purpose of this script is to measure the optical properties in the
parenchymal tissue surrounding the EPVS and the vasculature. This covers
just the following cases:
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
    - Apply tissue mask to vessels and optical properties
- Create WM mask
- Measure optical properties around EPVS and normal vessels
%}

%% Prepare environment
clc; close all;
% Add top-level directory
current_dir = pwd;
addpath(fullfile(current_dir));
% Directory for loading seg, mus, ret, mask, epvs
data_dir = ['/autofs/cluster/octdata2/users/mjhyman/' ...
    'oct_caa_analyses/optical_properties'];
% Initialize structuring element for dilation
r = 3; % units are voxels from origin to edge
se = strel('disk',r);
% Structure for storing parenchymal optical properties
parench = struct();

%%% MRIread parameters
% Voxel z-dimension (microns)
zvox = 20;
% header only flag
h_flag = 0;
% permute x,y dimensions flag
pflag = 0;
% Voxel dimensions (microns) for all runs
res = [20,20,20]; % resolution in microns
dtype = 'float'; % float is the equivalent of single

%%% Flag for loading .MAT struct and creating WM masks
% flag for reloading the .MAT struct for each subject
flag_load_caa_structs = true;
% flag for creating wm masks
flag_make_wm_mask = true;
% flag for saving wm mask nifti
flag_save_wm_mask = true;

%% Load each subject's .MAT struct and create WM mask

%%% Load the .MAT structs
if flag_load_caa_structs
    % CAA 6
    caa6 = load(fullfile(data_dir,'/caa6/caa6.mat'));
    caa6 = caa6.caa6;
    % CAA 17
    caa17 = load(fullfile(data_dir,'/caa17/occip/caa17.mat'));
    caa17 = caa17.caa17;
    % CAA 22
    caa22 = load(fullfile(data_dir,'/caa22/caa22.mat'));
    caa22 = caa22.caa22;
    % CAA 25
    caa25 = load(fullfile(data_dir,'/caa25/caa25.mat'));
    caa25 = caa25.caa25;
    % CAA 26
    % caa26 = load(fullfile(data_dir,'/caa26/caa26.mat'));
    % caa26 = caa26.caa26;
end

if flag_make_wm_mask
    %%% CAA 6 Frontal
    % import local variables
    ret = caa6.front.ret_full;
    seg = caa6.front.seg;
    mask = caa6.front.mask;
    % retardance white matter threshold
    th = 19;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa6/front/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa6.front.seg_wm = seg_wm;
    caa6.front.mask_wm = mask_wm;
    
    %%% CAA 6 Occip
    % import local variables
    ret = caa6.occip.ret_full;
    seg = caa6.occip.seg;
    mask = caa6.occip.mask;
    % retardance white matter threshold
    th = 22;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa6/occip/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa6.occip.seg_wm = seg_wm;
    caa6.occip.mask_wm = mask_wm;
    
    %%% CAA 17 Occip
    % import local variables
    ret = caa17.occip.ret_full;
    seg = caa17.occip.seg;
    mask = caa17.occip.mask;
    % retardance white matter threshold
    th = 19;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa17/occip/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa17.occip.seg_wm = seg_wm;
    caa17.occip.mask_wm = mask_wm;
    
    %%% CAA 22 Frontal
    % import local variables
    ret = caa22.front.ret_full;
    seg = caa22.front.seg;
    mask = caa22.front.mask;
    % retardance white matter threshold
    th = 16.9;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa22/front/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa22.front.seg_wm = seg_wm;
    caa22.front.mask_wm = mask_wm;
    
    %%% CAA 22 Occip
    % import local variables
    ret = caa22.occip.ret_full;
    seg = caa22.occip.seg;
    mask = caa22.occip.mask;
    % retardance white matter threshold
    th = 22.52;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa22/occip/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa22.occip.seg_wm = seg_wm;
    caa22.occip.mask_wm = mask_wm;
    
    %%% CAA 25 Frontal
    % import local variables
    ret = caa25.front.ret_full;
    seg = caa25.front.seg;
    mask = caa25.front.mask;
    % retardance white matter threshold
    th = 20.5;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa25/front/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa25.front.seg_wm = seg_wm;
    caa25.front.mask_wm = mask_wm;
    
    %%% CAA 25 Occip
    % import local variables
    ret = caa25.occip.ret_full;
    seg = caa25.occip.seg;
    mask = caa25.occip.mask;
    % retardance white matter threshold
    th = 22.5;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa25/occip/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa25.occip.seg_wm = seg_wm;
    caa25.occip.mask_wm = mask_wm;
    
    %%% CAA 26 Frontal
    % import local variables
    ret = caa26.front.ret_full;
    seg = caa26.front.seg;
    mask = caa26.front.mask;
    % retardance white matter threshold
    th = 26;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa26/front/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa26.front.seg_wm = seg_wm;
    caa26.front.mask_wm = mask_wm;
    
    %%% CAA 26 Occip
    % import local variables
    ret = caa26.occip.ret_full;
    seg = caa26.occip.seg;
    mask = caa26.occip.mask;
    % retardance white matter threshold
    th = 26.6;
    % Create WM mask and apply to vasculature
    fout = fullfile(data_dir,'/caa26/occip/wm_mask.nii');
    [mask_wm,seg_wm] = create_wm_mask(mask,ret,seg,th,fout,flag_save_wm_mask);
    caa26.occip.seg_wm = seg_wm;
    caa26.occip.mask_wm = mask_wm;
end

%% CAA 6 Frontal
fprintf('Starting CAA6 frontal\n')
% patient ID
sub = 'caa6';
loc = 'front';

%%% Load local parameters
mus = caa6.front.mus;
ret = caa6.front.ret_full;
seg = caa6.front.seg; % all segmented vessels
seg_wm = caa6.front.seg_wm; % white matter vessels
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;

%%% Measure optical properties of tissue surrounding vessels
[mus_mat, ret_mat] = parenchyma_optical_props(seg_wm,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = mus_mat;
parench.(sub).(loc).ves.pret = ret_mat;
fprintf('Finished CAA6 frontal\n')

%% CAA 6 Occipital
fprintf('Starting CAA6 occipital\n')
% patient ID
sub = 'caa6';
loc = 'occip';

%%% Load local parameters
mus = caa6.occip.mus;
ret = caa6.occip.ret_full;
seg = caa6.occip.seg; % all segmented vessels
seg_wm = caa6.occip.seg_wm;  % white matter vessels
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;

%%% Measure optical properties of tissue surrounding vessels
[mus_mat, ret_mat] = parenchyma_optical_props(seg_wm,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = rmmissing(mus_mat);
parench.(sub).(loc).ves.pret = rmmissing(ret_mat);
fprintf('Finished CAA6 occipital\n')

%% CAA 17 Occipital
fprintf('Starting CAA17 occipital\n')
% patient ID
sub = 'caa17';
loc = 'occip';

%%% Load local parameters
mus = caa17.occip.mus;
ret = caa17.occip.ret_full;
seg = caa17.occip.seg;
seg_wm = caa17.occip.seg_wm;
epvs = caa17.occip.epvs;
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;
% Exclude EPVS from mus and retardance
mus(epvs) = NaN;
ret(epvs) = NaN;

%%% exclude vessels overlapping w/ EPVS
[seg_wm_no_epvs] = exclude_epvs_ves(seg_wm,epvs);

%%% Measure optical properties of tissue surrounding vessels
[mus_mat,ret_mat] = parenchyma_optical_props(seg_wm_no_epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = mus_mat;
parench.(sub).(loc).ves.pret = ret_mat;

%%% Measure EPVS parenchyma
[mus_mat,ret_mat] = parenchyma_optical_props(epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = mus_mat;
parench.(sub).(loc).epvs.pret = ret_mat;
fprintf('Finished CAA17 occipital\n')

%% CAA 22 Frontal
fprintf('Starting CAA22 frontal\n')
% patient ID
sub = 'caa22';
loc = 'front';
% load mus, ret, seg, mask
mus = caa22.front.mus;
ret = caa22.front.ret_full;
seg = caa22.front.seg;
seg_wm = caa22.front.seg_wm;
epvs = caa22.front.epvs;
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;
% Exclude EPVS from mus and retardance
mus(epvs) = NaN;
ret(epvs) = NaN;

%%% exclude vessels overlapping w/ EPVS
[seg_wm_no_epvs] = exclude_epvs_ves(seg_wm,epvs);

%%% Measure optical properties of tissue surrounding vessels
[mus_mat,ret_mat] = parenchyma_optical_props(seg_wm_no_epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = mus_mat;
parench.(sub).(loc).ves.pret = ret_mat;

%%% Measure EPVS parenchyma
[mus_mat,ret_mat] = parenchyma_optical_props(epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = mus_mat;
parench.(sub).(loc).epvs.pret = ret_mat;
fprintf('Finished CAA22 frontal\n')

%% CAA 22 occipital
fprintf('Starting CAA22 occip\n')
% patient ID
sub = 'caa22';
loc = 'occip';
% load mus, ret, seg, mask
mus = caa22.occip.mus;
ret = caa22.occip.ret_full;
seg = caa22.occip.seg;
seg_wm = caa22.occip.seg_wm;
epvs = caa22.occip.epvs;
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;
% Exclude EPVS from mus and retardance
mus(epvs) = NaN;
ret(epvs) = NaN;

%%% exclude vessels overlapping w/ EPVS
[seg_wm_no_epvs] = exclude_epvs_ves(seg_wm,epvs);

%%% Measure optical properties of tissue surrounding vessels
[mus_mat,ret_mat] = parenchyma_optical_props(seg_wm_no_epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = mus_mat;
parench.(sub).(loc).ves.pret = ret_mat;

%%% Measure EPVS parenchyma
[mus_mat,ret_mat] = parenchyma_optical_props(epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = mus_mat;
parench.(sub).(loc).epvs.pret = ret_mat;
fprintf('Finished CAA22 occip\n')

%% CAA 25 front
fprintf('Starting CAA25 front\n')
% patient ID
sub = 'caa25';
loc = 'front';
% load mus, ret, seg, mask
mus = caa25.front.mus;
ret = caa25.front.ret_full;
seg = caa25.front.seg;
seg_wm = caa25.front.seg_wm;
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;

%%% Measure optical properties of tissue surrounding vessels
[mus_mat, ret_mat] = parenchyma_optical_props(seg_wm,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = mus_mat;
parench.(sub).(loc).ves.pret = ret_mat;
fprintf('Finished CAA25 front\n')

%% CAA 25 Occipital
fprintf('Starting CAA25 occip\n')
% patient ID
sub = 'caa25';
loc = 'occip';
% load mus, ret, seg, mask
mus = caa25.occip.mus;
ret = caa25.occip.ret_full;
seg = caa25.occip.seg;
seg_wm = caa25.occip.seg_wm;
epvs = caa25.occip.epvs;
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;
% Exclude EPVS from mus and retardance
mus(epvs) = NaN;
ret(epvs) = NaN;

%%% exclude vessels overlapping w/ EPVS
[seg_wm_no_epvs] = exclude_epvs_ves(seg_wm,epvs);

%%% Measure optical properties of tissue surrounding vessels
[mus_mat,ret_mat] = parenchyma_optical_props(seg_wm_no_epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = mus_mat;
parench.(sub).(loc).ves.pret = ret_mat;

%%% Measure EPVS parenchyma
[mus_mat,ret_mat] = parenchyma_optical_props(epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = mus_mat;
parench.(sub).(loc).epvs.pret = ret_mat;
fprintf('Finished CAA25 occip\n')

%% CAA 26 front
fprintf('Starting CAA26 front\n')
% patient ID
sub = 'caa26';
loc = 'front';
% load mus, ret, seg, mask
mus = caa26.front.mus;
ret = caa26.front.ret_full;
seg = caa26.front.seg;
seg_wm = caa26.front.seg_wm;
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;

%%% Measure optical properties of tissue surrounding vessels
[mus_mat, ret_mat] = parenchyma_optical_props(seg_wm,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = rmmissing(mus_mat);
parench.(sub).(loc).ves.pret = rmmissing(ret_mat);
fprintf('Finished CAA26 front\n')

%% CAA 26 Occipital
fprintf('Starting CAA26 occip\n')
% patient ID
sub = 'caa26';
loc = 'occip';
% load mus, ret, seg, mask
mus = caa26.occip.mus;
ret = caa26.occip.ret_full;
seg = caa26.occip.seg;
seg_wm = caa26.occip.seg_wm;
epvs = caa26.occip.epvs;
% Exclude vessels from mus and retardance
mus(seg) = NaN;
ret(seg) = NaN;
% Exclude EPVS from mus and retardance
mus(epvs) = NaN;
ret(epvs) = NaN;

%%% exclude vessels overlapping w/ EPVS
[seg_wm_no_epvs] = exclude_epvs_ves(seg_wm,epvs);

%%% Measure optical properties of tissue surrounding vessels
[mus_mat,ret_mat] = parenchyma_optical_props(seg_wm_no_epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).ves.pmus = mus_mat;
parench.(sub).(loc).ves.pret = ret_mat;

%%% Measure EPVS parenchyma
[mus_mat,ret_mat] = parenchyma_optical_props(epvs,mus,ret,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = mus_mat;
parench.(sub).(loc).epvs.pret = ret_mat;
fprintf('Finished CAA26 occip\n')

%% Save optical properties of tissue surrounding EPVS
% Backup struct
fout = fullfile(data_dir,'parenchyma_optical_properties.mat');
save(fout,"parench",'-v7.3');