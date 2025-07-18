%% Co-register Histology b/w Gallyas and OCT
%{
The gallyas stain was performed on the same sections as the OCT sections.
The Gallyas stains were co-registered to OCT with the landmark-based
registration tool. This script will:

- measure the correlation between the density of gallyas surrounding the
EPVS and the respective retardance in the OCT images

%}

%% Top-level settings
clc; close all;
% directory to Gallyas stains
gallyas_dir=['/autofs/cluster/octdata3/users/mjhyman' ...
    '/oct_caa_analyses/histology/Gallyas/'];
% Subfolders for each subject
subdirs = {'CAA6_frontal_61_Gallyas','CAA6_occipital_57_Gallyas',...
    'CAA25_frontal_65_Gallyas_5x',...
    'CAA25_occipital_81_Gallyas_5x','CAA26_frontal_48_Gallyas_5x',...
    'CAA26_occipital_120_Gallyas_5x'};
% Initialize struct array for storing histology
hist = struct();

%% Import

%%% Set filename suffixes
% WM mask
suf_mask = '_5x_ds_mask.tif';
% stain
suf_stain = '_5x_ds_gray_reg.tif';
% EPVS annotations
suf_epvs = '_5x_ds_EPVS_reg.tif';
% vessel annotations
suf_ves = '_5x_ds_VES_reg.tif';
% OCT (retardance)
suf_oct = '_RET.tif';
% string to remove to find basename
suf_rm = '_5x_ds_mask.tif';

% Iterate over subject/region
for ii = 1:length(subdirs)
    % Set the folder for subject/region
    folder = fullfile(gallyas_dir,subdirs{ii});
    
    %%% Extract base filename
    % List all EPVS mask files
    mask_fnames = strcat('*',suf_mask);
    mask_files = dir(fullfile(folder, mask_fnames));
    % Extract base filename
    epvs_filename = mask_files.name;
    base_name = erase(epvs_filename, suf_rm);
    
    %%% Import each image
    % Construct full file paths
    mask_path = fullfile(folder, [base_name suf_mask]);
    stain_path = fullfile(folder, [base_name suf_stain]);
    epvs_path = fullfile(folder, [base_name suf_epvs]);
    ves_path = fullfile(folder, [base_name suf_ves]);
    ret_path = fullfile(folder, [base_name suf_oct]);
    
    % Load the images (assume both files exist)
    mask = logical(imread(mask_path));
    stain = imread(stain_path);
    epvs = logical(imread(epvs_path));
    ves = logical(imread(ves_path));
    ret = imread(ret_path);

    %%% Check if the images have 3 channels, if so take first
    mask = single_ch(mask);
    epvs = single_ch(epvs);
    ves = single_ch(ves);
    ret = single_ch(ret);

    % Store in struct
    hist(ii).baseName = string(base_name);
    hist(ii).mask = mask;
    hist(ii).image = stain;
    hist(ii).epvs = epvs;
    hist(ii).ves = ves;
    hist(ii).ret = ret;
end

%%% Measure the histology

% Define radii of measurements
radii = [0,2,4,6,8];
% Define threshold for increasing donut
th = 5;
% Measure Gallyas
hist = histology_measure_epvs(hist,radii,th,[],[]);


%% Function to keep first channel from logical
% some of the imported images have three channels, but they are all
% logical. this is due to an issue from Fiji. Just keep one of these
% channels

function im = single_ch(im)

if ~ismatrix(im)
    im = im(:,:,1);
end

end