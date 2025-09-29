%% Update .MAT structs with updated EPVS and masks
% The EPVS and white matter masks are being updated, and this script will
% update the structs.

% Directory to optical properties
dpath = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
            'optical_properties'];

%{
TODO:
- dilate / erode vessel segmentation
- import/update EPVS for all
- convert all EPVS to logical
- update CAA 26 frontal mask
- apply masks to EPVS as well

After:
Transfer .MAT to SCC
Run donuts on SCC
Run EPVS heatmaps on SCC

%}

%% Load all structures

%%% Load subject structs
% CAA 6
fprintf('Loading CAA6\n')
caa6 = load(fullfile(dpath,'/caa6/caa6.mat'));
caa6 = caa6.caa6;
fprintf('Finished Loading CAA6\n')
% CAA 17
fprintf('Loading CAA17\n')
caa17 = load(fullfile(dpath,'/caa17/occip/caa17.mat'));
caa17 = caa17.caa17;
fprintf('Finished Loading CAA17\n')
% CAA 22
fprintf('Loading CAA22\n')
caa22 = load(fullfile(dpath,'/caa22/caa22.mat'));
caa22 = caa22.caa22;
fprintf('Finished Loading CAA22\n')
% CAA 25
fprintf('Loading CAA25\n')
caa25 = load(fullfile(dpath,'/caa25/caa25.mat'));
caa25 = caa25.caa25;
fprintf('Finished Loading CAA25\n')
% CAA 26
fprintf('Loading CAA26\n')
caa26 = load(fullfile(dpath,'/caa26/caa26.mat'));
caa26 = caa26.caa26;
fprintf('Finished Loading CAA26\n')

%% TODO: Dilate/Erode vessel segmentation


%% TODO: Update EPVS for all subjects
% CAA 22 front
tmp = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
    'optical_properties/caa22/front/epvs_manual_latest.nii'];
epvs = MRIread(tmp,0,0);
epvs = epvs.vol;
caa22.front.epvs = epvs;
fprintf('Updated caa22 front\n')

%% Load updated white matter structs
%%% CAA 6 Frontal
fprintf('Importing caa6 frontal\n')
% import wm mask
mask_wm = fullfile(dpath,'caa6/front/wm_mask_revised.nii');
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
mask_wm = fullfile(dpath,'caa6/occip/wm_mask_revised.nii');
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
mask_wm = fullfile(dpath,'caa17/occip/wm_mask_revised.nii');
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
mask_wm = fullfile(dpath,'caa22/front/wm_mask_revised.nii');
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
mask_wm = fullfile(dpath,'caa22/occip/wm_mask_revised.nii');
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
mask_wm = fullfile(dpath,'caa25/front/wm_mask_revised.nii');
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
mask_wm = fullfile(dpath,'caa25/occip/wm_mask_revised.nii');
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
mask_wm = fullfile(dpath,'caa26/front/wm_mask_revised.nii');
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
mask_wm = fullfile(dpath,'caa26/occip/wm_mask_revised.nii');
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

%% Save .MAT

% CAA 6
fprintf('Saving CAA6/n')
dout = fullfile(dpath,'/caa6/caa6.mat');
save(dout,'caa6','-v7.3');

% CAA 17
fprintf('Saving CAA17/n')
dout = fullfile(dpath,'/caa17/occip/caa17.mat');
save(dout,'caa17','-v7.3');

% CAA 22
fprintf('Saving CAA22/n')
dout = fullfile(dpath,'/caa22/caa22.mat');
save(dout,'caa22','-v7.3');

% CAA 25
fprintf('Saving CAA25/n')
dout = fullfile(dpath,'/caa25/caa25.mat');
save(dout,'caa25','-v7.3');

% CAA 26
fprintf('Saving CAA26/n')
dout = fullfile(dpath,'/caa26/caa26.mat');
save(dout,'caa26','-v7.3');

fprintf('Finished updating volumes')