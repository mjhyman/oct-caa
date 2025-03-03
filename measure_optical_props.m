%% Analyze the Optical Properties
%{
The purpose of this script is to analyze the optical properties in the
parenchymal tissue surrounding the EPVS and the vasculature. This covers
just the following cases:
- CAA 17 occipital
- CAA 22 Frontal
- CAA 22 occipital
- CAA 25 Occipital
- CAA 26 Occipital

Outline:
- IMPORT:
    - tissue masks (exclude background signal)
    - EPVS
    - segmented vasculature (from Etienne)
    - optical properties (mus + retardance)
- MASK
    - Apply tissue mask to vessels and optical properties
- Measure optical properties around EPVS and normal vessels

TODO:
- create WM masks from retardance
- mask the vessels

%}

%% Prepare environment
clc; close all;
% Add top-level directory
current_dir = pwd;
addpath(fullfile(current_dir));
% Output directory for optical properties
analysis_dir = '/autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses';
% optical properties struct
opdir = fullfile(analysis_dir,'optical_properties.mat');
% Top-level directory
mask_dir = '/autofs/cluster/octdata2/users/epc28/for/mack/2024-10-10';
% Output of "parench" struct
op_out = fullfile(analysis_dir,'parenchyma.mat');
% Initialize structuring element for dilation
r = 25; % units are voxels from origin to edge
se = strel('disk',r);

%%% MRIread parameters
% header only flag
h_flag = 0;
% permute x,y dimensions flag
p_flag = 0;
% Voxel z-dimension (microns)
zvox = 20;

%%% Structs for storing masks, seg, epvs, optical props
masks = struct();
segs = struct();
epvs = struct();
op = struct();
% Optical properties in tissue surrounding EPVS and vessels
if isfile(op_out)
    load(op_out)
else
    parench = struct();
end

%% CAA 17 Occipital
% patient ID
sub = 'caa17';
loc = 'occip';
% slice thickness (microns)
thick = 50;

%%% Import
% Import tissue masK
tmp = MRIread(fullfile(mask_dir,...
    'caa17-occipital_unet-tissuemask.nii'),h_flag,p_flag);
mask = imbinarize(tmp.vol);
% Import vessel segmentation
seg = MRIread(['/autofs/cluster/octdata2/users/epc28/data/CAA/caa17/' ...
    'occipital/segmentations/caa17_occipital_THRESH-0.5.nii'],h_flag,p_flag);
seg = imbinarize(seg.vol);
ves.caa17.occip = seg;
% Import retardance
oprop = load(opdir,'caa17');
ret = oprop.caa17.occip.ret;
% Import scattering AIP
mus = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/' ...
    'process_run1/StackNii/Stacked_mus.nii'];
mus = MRIread(mus,h_flag,p_flag);
mus = single(mus.vol);
% Import EPVS
f = ['/autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses/EPVS/' ...
    'caa17_occipital_EPVS_segmentation_07072023.nii'];
epvs_ = MRIread(f,h_flag,p_flag);
epvs_ = imbinarize(epvs_.vol);

%%% Scaling
% scale scattering to mask
nrow = size(mask,1);
ncol = size(mask,2);
mus = imresize3(mus,[nrow,ncol,size(mus,3)]);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);
% Save to struct
masks.(sub).(loc) = mask;
segs.(sub).(loc) = seg;
op.(sub).(loc).mus = mus;
op.(sub).(loc).ret = ret;
epvs.(sub).(loc) = epvs_;

%%% Ensure all x,y dimensions match
y = size(mask,1); x = size(mask,2);
assert(size(seg,1) == y); assert(size(seg,2) == x);
assert(size(epvs_,1) == y); assert(size(epvs_,2) == x);
assert(size(ret,1) == y); assert(size(ret,2) == x);
assert(size(mus,1) == y); assert(size(mus,2) == x);

%%% Measure optical properties of tissue surrounding EPVS
% Retrieve all disjoint EPVS groups
cc = bwconncomp(epvs_,26);
% Cell array of lists of pixel indices for each disjoint EPVS group
idxlst = cc.PixelIdxList;
% Initialize arrays to store parenchymal mus and retardance
pmus = zeros(size(idxlst,2),1);
pret = zeros(size(idxlst,2),1);
% Iterate over each EPVS
for ii = 1:length(idxlst)
    % Retrieve EPVS voxel indices
    idx = idxlst{ii};
    % Create zeros matrix the size of EPVS matrix
    emat = false(size(epvs_));
    % Set voxels within EPVS group equal to 1
    emat(idx) = 1;
    % Measure the optical properties
    [pmus(ii),pret(ii)] = meas_optical_props(emat,mus,ret,thick,zvox,se);
    % Print progress
    fprintf('Finished %i of %i\n',ii,length(idxlst))
end
% Add to struct
parench.(sub).(loc).epvs.pmus = rmmissing(pmus);
parench.(sub).(loc).epvs.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');
pause(0.01)

%%% Exclude vessels w/ EPVS
% Apply tissue mask to vessels
seg = seg .* mask;
% Pixel indices of segmentation
cc = bwconncomp(seg,26);
idxlst = cc.PixelIdxList;
% Create matrix of retained vessels
seg_keep = false(size(seg));
% Indices to keep
pmus = parench.(sub).(loc).epvs.pmus;
kp = 1:length(rmmissing(pmus));
% Assign vessels as true
idxlst = idxlst(kp);
seg_keep(cell2mat(idxlst(:))) = true;
% Check if remaining vessels have EPVS
[seg_keep] = exclude_epvs_ves(seg_keep,epvs_);
cc = bwconncomp(seg_keep,26);
seg_keep_lst = cc.PixelIdxList;

%%% Measure optical properties of tissue surrounding non-epvs vessels
[pmus, pret] =...
    iterate_segments(seg_keep_lst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).ves.pmus = rmmissing(pmus);
parench.(sub).(loc).ves.pret = rmmissing(pret);

%% CAA 22 Frontal
% patient ID
sub = 'caa22';
loc = 'front';

%%% Import
% Import mask
tmp = MRIread(fullfile(mask_dir,...
    'caa22-frontal_unet-tissuemask.nii'),h_flag,p_flag);
mask = imbinarize(tmp.vol);
% Import vessel segmentation
seg = MRIread(['/autofs/cluster/octdata2/users/epc28/data/CAA/caa22/'...
    'frontal/segmentations/caa22_frontal_THRESH-0.5.nii'],h_flag,p_flag);
seg = imbinarize(seg.vol);
ves.caa22.front = seg;
% Import retardance
oprop = load(opdir,'caa22');
ret = oprop.caa22.front.ret;
% Import scattering
mus = MRIread(['/autofs/cluster/octdata2/users/Chao/caa/caa_22/' ...
    'processed/20211018/StackNii/Stacked_mus_AIP.nii'],h_flag,p_flag);
mus = mus.vol;
% Import EPVS
epvs_=logical(tiffreadVolume(['/autofs/cluster/octdata2/users/mjhyman/'...
    'oct_caa_analyses/EPVS/caa22_frontal_EPVS_segmentation_03262024.tif']));
epvs_ = permute(epvs_,[2,1,3]);

%%% Scale mus and ret to match masks
% Concatenate scattering matrices
mus = imresize3(mus,[nrow,ncol,size(mus,3)]);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);
% Save to struct
masks.(sub).(loc) = mask;
segs.(sub).(loc) = seg;
op.(sub).(loc).mus = mus;
op.(sub).(loc).ret = ret;
epvs.(sub).(loc) = epvs_;

%%% Ensure all x,y dimensions match
y = size(mask,1); x = size(mask,2);
assert(size(seg,1) == y); assert(size(seg,2) == x);
assert(size(epvs_,1) == y); assert(size(epvs_,2) == x);
assert(size(ret,1) == y); assert(size(ret,2) == x);
assert(size(mus,1) == y); assert(size(mus,2) == x);

%%% Measure optical properties of tissue surrounding EPVS
% Retrieve all disjoint EPVS groups
cc = bwconncomp(epvs_,26);
% Cell array of lists of pixel indices for each disjoint EPVS group
idxlst = cc.PixelIdxList;
[pmus, pret] =...
    iterate_segments(idxlst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = rmmissing(pmus);
parench.(sub).(loc).epvs.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');

%%% Exclude vessels w/ EPVS
% Apply tissue mask to vessels
seg = seg .* mask;
% Pixel indices of segmentation
cc = bwconncomp(seg,26);
idxlst = cc.PixelIdxList;
% Keep first 200 elements
seg_keep_lst = idxlst(1:200);

%%% Measure optical properties of tissue surrounding non-epvs vessels
[pmus, pret] =...
    iterate_segments(seg_keep_lst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).ves.pmus = rmmissing(pmus);
parench.(sub).(loc).ves.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');

%% CAA 22 occipital (TODO: align mus, ret. measure props)
% matching: seg, mask, epvs
% patient ID
sub = 'caa22';
loc = 'occip';

%%% Import
% Import mask
tmp = MRIread(fullfile(mask_dir,...
    'caa22-occipital_unet-tissuemask.nii'),h_flag,p_flag);
mask = imbinarize(tmp.vol);
% Import vessel segmentation
seg = MRIread(['//autofs/cluster/octdata2/users/epc28/data/CAA/caa22/' ...
    'occipital/segmentations/caa22-occipital_prediction-r2-BINARY.nii'], ...
    h_flag,p_flag);
seg = imbinarize(seg.vol);
ves.caa22.occip = seg;
% Import retardance
oprop = load(opdir,'caa22');
ret = oprop.caa22.occip.ret;
% Import scattering
mus = MRIread(['/autofs/cluster/octdata2/users/Chao/caa/caa_22/' ...
    'processed/20211007/StackNii/Stacked_mus_AIP.nii'],h_flag,p_flag);
mus = mus.vol;
% Import EPVS
epvs_ = logical(tiffreadVolume(['/autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses/EPVS/' ...
    'caa22_occipital_EPVS_segmentation_02132024.tif']));
epvs_ = permute(epvs_,[2 1 3]);

%%% Scale mus and ret to match masks
% scale scattering to mask
nrow = size(mask,1);
ncol = size(mask,2);
% Scale scattering to mask
mus = imresize3(mus,[nrow,ncol,size(mus,3)]);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Ensure all x,y dimensions match
y = size(mask,1); x = size(mask,2);
assert(size(seg,1) == y); assert(size(seg,2) == x);
assert(size(epvs_,1) == y); assert(size(epvs_,2) == x);
assert(size(ret,1) == y); assert(size(ret,2) == x);
assert(size(mus,1) == y); assert(size(mus,2) == x);
% Save to struct
masks.(sub).(loc) = mask;
segs.(sub).(loc) = seg;
op.(sub).(loc).mus = mus;
op.(sub).(loc).ret = ret;
epvs.(sub).(loc) = epvs_;

%%% Measure optical properties of tissue surrounding EPVS
% Retrieve all disjoint EPVS groups
cc = bwconncomp(epvs_,26);
% Cell array of lists of pixel indices for each disjoint EPVS group
idxlst = cc.PixelIdxList;
% Select first 200 groups
nmax = min(200,length(idxlst));
idxlst = idxlst(1:nmax);
% Iterate over each group and measure optical properties
[pmus, pret] =...
    iterate_segments(idxlst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = rmmissing(pmus);
parench.(sub).(loc).epvs.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');

%%% Exclude vessels w/ EPVS
% Apply tissue mask to vessels
seg = seg .* mask;
% Pixel indices of segmentation
cc = bwconncomp(seg,26);
idxlst = cc.PixelIdxList;
% NOTE: the following code was commented to expedite plotting
% % Create matrix of retained vessels
% seg_keep = false(size(seg));
% seg_keep(cell2mat(idxlst(:))) = true;
% % Exclude vessels with EPVS
% [seg_keep] = exclude_epvs_ves(seg_keep,epvs_);
% cc = bwconncomp(seg_keep,26);
% seg_keep_lst = cc.PixelIdxList;
% % Take just 200 groups of vessels
pmus = parench.(sub).(loc).epvs.pmus;
kp = 1:length(rmmissing(pmus));
% Assign vessels as true
idxlst = idxlst(kp);
seg_keep_lst = idxlst;

%%% Measure optical properties of tissue surrounding non-epvs vessels
[pmus, pret] =...
    iterate_segments(seg_keep_lst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).ves.pmus = rmmissing(pmus);
parench.(sub).(loc).ves.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');

%% CAA 25 Occipital
% TODO: determine why most of mus and ret return NaN

% patient ID
sub = 'caa25';
loc = 'occip';

%%% Import
% Import mask
tmp = MRIread(fullfile(mask_dir,...
    'caa25-occipital_unet-tissuemask.nii'),h_flag,p_flag);
mask = imbinarize(tmp.vol);
% Import vessel segmentation
seg = MRIread(['/autofs/cluster/octdata2/users/epc28/data/CAA/caa25/' ...
    'occipital/segmentations/caa25_occipital_THRESH-0.5.nii'],h_flag,p_flag);
seg = imbinarize(seg.vol);
ves.caa25.occip = seg;
% Import mus + retardance
mus = MRIread(['/autofs/space/omega_001/users/caa/CAA25_Occipital/' ...
    'process_caa25_occipital_run1/StackNii/Stacked_AIP.nii'],h_flag,p_flag);
mus = mus.vol;
ret = MRIread(['/autofs/space/omega_001/users/caa/CAA25_Occipital/' ...
    'process_caa25_occipital_run1/StackNii/Stacked_Retardance.nii'],...
    h_flag,p_flag);
ret = ret.vol;
% Import EPVS
epvs_ = logical(tiffreadVolume(['/autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses/EPVS/' ...
    'caa25_occipital_EPVS_mus_segmentation_02052024.tif']));
epvs_ = permute(epvs_,[2 1 3]);

%%% TODO Scale mus and ret to match masks
% scale scattering to mask
nrow = size(mask,1);
ncol = size(mask,2);
% Scale scattering to mask
mus = imresize3(mus,[nrow,ncol,size(mus,3)]);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Ensure all x,y dimensions match
y = size(mask,1); x = size(mask,2);
assert(size(seg,1) == y); assert(size(seg,2) == x);
assert(size(epvs_,1) == y); assert(size(epvs_,2) == x);
assert(size(ret,1) == y); assert(size(ret,2) == x);
assert(size(mus,1) == y); assert(size(mus,2) == x);
% Save to struct
masks.(sub).(loc) = mask;
segs.(sub).(loc) = seg;
op.(sub).(loc).mus = mus;
op.(sub).(loc).ret = ret;
epvs.(sub).(loc) = epvs_;

%%% Measure optical properties of tissue surrounding EPVS
% Retrieve all disjoint EPVS groups
cc = bwconncomp(epvs_,26);
% Cell array of lists of pixel indices for each disjoint EPVS group
idxlst = cc.PixelIdxList;
% Select first 200 groups
nmax = min(200,length(idxlst));
idxlst = idxlst(1:nmax);
% Iterate over each group and measure optical properties
[pmus, pret] =...
    iterate_segments(idxlst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = rmmissing(pmus);
parench.(sub).(loc).epvs.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');

%%% Exclude vessels w/ EPVS
% Apply tissue mask to vessels
seg = seg .* mask;
% Pixel indices of segmentation
cc = bwconncomp(seg,26);
idxlst = cc.PixelIdxList;
% NOTE: the following code was commented to expedite plotting
% % Create matrix of retained vessels
% seg_keep = false(size(seg));
% seg_keep(cell2mat(idxlst(:))) = true;
% % Exclude vessels with EPVS
% [seg_keep] = exclude_epvs_ves(seg_keep,epvs_);
% cc = bwconncomp(seg_keep,26);
% seg_keep_lst = cc.PixelIdxList;
% % Take just 200 groups of vessels
pmus = parench.(sub).(loc).epvs.pmus;
kp = 1:length(rmmissing(pmus));
% Assign vessels as true
idxlst = idxlst(kp);
seg_keep_lst = idxlst;

%%% Measure optical properties of tissue surrounding non-epvs vessels
[pmus, pret] =...
    iterate_segments(seg_keep_lst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).ves.pmus = rmmissing(pmus);
parench.(sub).(loc).ves.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');

%% CAA 26 Occipital (no mask available)
% patient ID
sub = 'caa26';
loc = 'occip';

%%% Import
% Import vessel segmentation
seg = MRIread(['/autofs/cluster/octdata2/users/epc28/data/CAA/caa26/' ...
    'occipital/segmentations/caa26_occipital_THRESH-0.5.nii'],h_flag,p_flag);
seg = imbinarize(seg.vol);
ves.caa26.occip = seg;
% Import mus + retardance
mus = MRIread(['/autofs/space/omega_001/users/caa/CAA26_Occipital/' ...
  'Process_caa26_occipital/StackNii/Stacked_mus_AIP.nii'],h_flag,p_flag);
mus = mus.vol;
ret = MRIread(['/autofs/space/omega_001/users/caa/CAA26_Occipital/' ...
  'Process_caa26_occipital/StackNii/Stacked_Retardance.nii'],h_flag,p_flag);
ret = ret.vol;
% Import EPVS
epvs_ = logical(tiffreadVolume(['/autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses/EPVS/' ...
    'caa26_occipital_EPVS_segmentation.tif']));
epvs_ = permute(epvs_,[2 1 3]);

%%% TODO Scale mus and ret to match masks
% scale scattering to mask
nrow = size(seg,1);
ncol = size(seg,2);
% Scale scattering to mask
mus = imresize3(mus,[nrow,ncol,size(mus,3)]);
% Scale retardance to mask
ret = imresize3(ret,[nrow,ncol,size(ret,3)]);

%%% Ensure all x,y dimensions match
y = size(seg,1); x = size(seg,2);
assert(size(seg,1) == y); assert(size(seg,2) == x);
assert(size(epvs_,1) == y); assert(size(epvs_,2) == x);
assert(size(ret,1) == y); assert(size(ret,2) == x);
assert(size(mus,1) == y); assert(size(mus,2) == x);
% Save to struct
segs.(sub).(loc) = seg;
op.(sub).(loc).mus = mus;
op.(sub).(loc).ret = ret;
epvs.(sub).(loc) = epvs_;

%% Measure optical properties of tissue surrounding EPVS
% Retrieve all disjoint EPVS groups
cc = bwconncomp(epvs_,26);
% Cell array of lists of pixel indices for each disjoint EPVS group
idxlst = cc.PixelIdxList;
% Select first 200 groups
nmax = min(200,length(idxlst));
idxlst = idxlst(1:nmax);
% Iterate over each group and measure optical properties
[pmus, pret] =...
    iterate_segments(idxlst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).epvs.pmus = rmmissing(pmus);
parench.(sub).(loc).epvs.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');

%%% Exclude vessels w/ EPVS
% Pixel indices of segmentation
cc = bwconncomp(seg,26);
idxlst = cc.PixelIdxList;
% NOTE: the following code was commented to expedite plotting
% % Create matrix of retained vessels
% seg_keep = false(size(seg));
% seg_keep(cell2mat(idxlst(:))) = true;
% % Exclude vessels with EPVS
% [seg_keep] = exclude_epvs_ves(seg_keep,epvs_);
% cc = bwconncomp(seg_keep,26);
% seg_keep_lst = cc.PixelIdxList;
% % Take just 200 groups of vessels
pmus = parench.(sub).(loc).epvs.pmus;
kp = 1:length(rmmissing(pmus));
% Assign vessels as true
idxlst = idxlst(kp);
seg_keep_lst = idxlst;

%%% Measure optical properties of tissue surrounding non-epvs vessels
[pmus, pret] =...
    iterate_segments(seg_keep_lst,size(seg),mus,ret,thick,zvox,se);
% Add to struct
parench.(sub).(loc).ves.pmus = rmmissing(pmus);
parench.(sub).(loc).ves.pret = rmmissing(pret);
% Backup struct
save(op_out,"parench",'-v7.3');


%% Statistics
% Subjects: caa17_occip, caa_22_front, caa22_occip, caa25_occip,
% caa26_occip
% Notes: excluding caa25 because most values are returned as NaN. Need to
% debug this.

%%% Retrieve data
% Retrieve vessels (scattering)
ves_pmus_17 = parench.caa17.occip.ves.pmus;
ves_pmus_22f = parench.caa22.front.ves.pmus;
ves_pmus_22o = parench.caa22.occip.ves.pmus;
% ves_pmus_25 = parench.caa25.occip.ves.pmus;
ves_pmus_26 = parench.caa26.occip.ves.pmus;
% ves_mus = vertcat(ves_pmus_17(:),ves_pmus_22f(:),ves_pmus_22o(:),...
%     ves_pmus_25(:),ves_pmus_26(:));
ves_mus = vertcat(ves_pmus_17(:),ves_pmus_22f(:),ves_pmus_22o(:),...
                    ves_pmus_26(:));
% Retrieve vessels (retardance)
ves_pret_17 = parench.caa17.occip.ves.pret;
ves_pret_22f = parench.caa22.front.ves.pret;
ves_pret_22o = parench.caa22.occip.ves.pret;
% ves_pret_25 = parench.caa25.occip.ves.pret;
ves_pret_26 = parench.caa26.occip.ves.pret;
% ves_ret = vertcat(ves_pret_17(:),ves_pret_22f(:),ves_pret_22o(:),...
%     ves_pret_25(:),ves_pret_26(:));
ves_ret = vertcat(ves_pret_17(:),ves_pret_22f(:),ves_pret_22o(:),...
    ves_pret_26(:));
% Retrieve EPVS (scattering)
epvs_pmus_17 = parench.caa17.occip.epvs.pmus;
epvs_pmus_22f = parench.caa22.front.epvs.pmus;
epvs_pmus_22o = parench.caa22.occip.epvs.pmus;
% epvs_pmus_25 = parench.caa25.occip.epvs.pmus;
epvs_pmus_26 = parench.caa26.occip.epvs.pmus;
% epvs_mus = vertcat(epvs_pmus_17(:),epvs_pmus_22f(:),epvs_pmus_22o(:),...
%     epvs_pmus_25(:),epvs_pmus_26(:));
epvs_mus = vertcat(epvs_pmus_17(:),epvs_pmus_22f(:),epvs_pmus_22o(:),...
    epvs_pmus_26(:));
% Retrieve EPVS (retardance)
epvs_pret_17 = parench.caa17.occip.epvs.pret;
epvs_pret_22f = parench.caa22.front.epvs.pret;
epvs_pret_22o = parench.caa22.occip.epvs.pret;
% epvs_pret_25 = parench.caa25.occip.epvs.pret;
epvs_pret_26 = parench.caa26.occip.epvs.pret;
% epvs_ret = vertcat(epvs_pret_17(:),epvs_pret_22f(:),epvs_pret_22o(:),...
%     epvs_pret_25(:),epvs_pret_26(:));
epvs_ret = vertcat(epvs_pret_17(:),epvs_pret_22f(:),epvs_pret_22o(:),...
    epvs_pret_26(:));

%%% Bar Chart - Scattering
% Subject mean EPVS scattering coefficient
mean_epvs_mus = [mean(epvs_pmus_17),mean(epvs_pmus_22f),...
                mean(epvs_pmus_22o),mean(epvs_pmus_26)];
% Subject mean vessel scattering coefficient
mean_ves_mus = [mean(ves_pmus_17),mean(ves_pmus_22f),...
                mean(ves_pmus_22o),mean(ves_pmus_26)];
% Barchart
figure;
x = categorical({'CAA 17','CAA 22 Front','CAA 22 Occip','CAA 26 Occip'});
bar(x,[mean_ves_mus',mean_epvs_mus'])
title('Scattering Coefficient')
xtickangle(45)
set(gca,'FontSize',25)

%%% Bar Chart - Retardance
% Subject mean EPVS scattering coefficient
mean_epvs_ret = [mean(epvs_pret_17),mean(epvs_pret_22f),...
                mean(epvs_pret_22o),mean(epvs_pret_26)];
% Subject mean vessel scattering coefficient
mean_ves_ret = [mean(ves_pret_17),mean(ves_pret_22f),...
                mean(ves_pret_22o),mean(ves_pret_26)];
% Barchart
figure;
x = categorical({'CAA 17','CAA 22 Front','CAA 22 Occip','CAA 26 Occip'});
bar(x,[mean_ves_ret',mean_epvs_ret'])
title('Retardance')
xtickangle(45)
set(gca,'FontSize',25)

%%% t-test scattering coefficient (EPVS vs. non-EPVS)
[h_mus,p_mus] = ttest2(ves_mus,epvs_mus,'Alpha',0.05);

%%% t-test retardance (EPVS vs. non-EPVS)
[h_ret,p_ret] = ttest2(ves_ret,epvs_ret,'Alpha',0.05);

%% 

function [pmus, pret] = iterate_segments(idxlst,seg_size,mus,ret,...
                                        thick,zvox,se)
% ITERATE_SEGMENTS measure optical properties for each segment in list
%  Iterate cell array of segment or EPVS voxel indices for each grouping of
%  voxels (segment or EPVS).
%
% INPUTS:
%   idxlst (cell array): cell array of indices of the voxels of the vessel
%           or the EPVS
%   seg_size (uint array): [y,x,z] size of segmentation array
%   mus (single matrix): scattering coefficient matrix
%   ret (single matrix): retardance matrix
%   thick (uint): thickness of slice (microns)
%   zvox (uint): voxel dimensions (microns)
%   se (matrix): structuring element for dilation
% OUTPUTS:
%   pmus (single matrix): scattering coefficients in surrounding
%       parenchymal tissue
%   pret (single matrix): retardance in surrounding
%       parenchymal tissue

% Initialize arrays to store parenchymal mus and retardance
pmus = zeros(size(idxlst,2),1);
pret = zeros(size(idxlst,2),1);
% Iterate over all EPVS
for ii = 1:length(idxlst)
    % Retrieve EPVS voxel indices
    idx = idxlst{ii};
    % Create zeros matrix the size of EPVS or vessel matrix
    emat = false(seg_size);
    % Set voxels within EPVS group equal to 1
    emat(idx) = 1;
    % Measure the optical properties
    [pmus(ii),pret(ii)] = meas_optical_props(emat,mus,ret,thick,zvox,se);
    fprintf('Finished %i of %i\n',ii,length(idxlst))
end

end