%% Import downsampled stains
%{
These are the stains that our collaborator performed. The grayscale images
can be directly used to compute the stain density surrounding the EPVS
because there is no counterstain. The deconvolved images had a
counterstain and therefore had to be deconvolved to extract the target
stain.

Purpose of this script:
- import stains
    - convert to uint8 grayscale (if not already)
- import EPVS annotations
    - convert to logical
- import non-EPVS Vessel annotations
    - convert to logical
- import tissue border
    - convert to logical
- systematic control measurements:
    - evenly-spaced control measurements (exclude EPVS)
- create struct of results

Grayscale Images:
- LHE
- Gallyas

Deconvolved:
- CD68
- Fibrin
- GFAP

%}

%% Top-level settings
clear; clc; close all;

%%% Directories
% Input directory
hdir='/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/histology';
% Directory to save output
stat_sheet = fullfile(hdir, 'histo_stats_04Sep2025.xlsx');
% Directory to save output figures
figdir = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
            'figures/Histology/ret_gallyas'];

%%% measurement radius (units = pixels)
mrad = 25;

%%% Flags
% histogram matching
hflag = true;
% plot figure flag
pflag = false;

%% LHE
% Stain directory
stain_dir = fullfile(hdir,'LHE/');
% stain suffix
stain_suffix = '_M_BackgroundMask.tif';
% EPVS suffix
epvs_suffix = '_EPVS_Mask.tif';
% Vessel filename suffix
ves_suffix = '_VESSEL_Mask.tif';
% Mask suffix
mask_suffix = '_mask.tif';

%%% Import LHE stain
[lhe] = import_stain(stain_dir,stain_suffix,epvs_suffix, ...
                    ves_suffix,mask_suffix);

%% LHE EPVS + Vessel measurements
fout = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
        'figures/Histology/ret_LHE_myelin'];
lhe = measure_epvs_and_vessel_constant(lhe,mrad,fout,hflag,pflag);

%% LHE Statistics
% One-sided Wilcoxon signed-rank test
% Hypothesize myelin rarefaction around EPVS
%   --> (med(EPVS) - med(control)) < 0
%   --> left tailed
tail = 'left';
% Wilcoxon signed rank
[lhe_stats] = wilcoxon_epvs_ves(lhe,tail);
% Write to specific sheet named "LHE"
writetable(lhe_stats, stat_sheet, 'WriteRowNames', true, 'Sheet', 'LHE');

%% CD68
% Stain directory
stain_dir = fullfile(hdir,'CD68/');
% stain suffix
stain_suffix = '_shrunk_deconv.tif';
% EPVS suffix
epvs_suffix = '_shrunk_EPVS.tif';
% Vessel filename suffix
ves_suffix = '_shrunk_VESSEL.tif';
% Mask suffix
mask_suffix = '_shrunk_mask.tif';

%%% Import LHE stain
[cd68] = import_stain(stain_dir,stain_suffix,epvs_suffix, ...
                    ves_suffix,mask_suffix);

%% CD68 EPVS + Vessel measurements
fout = [];
cd68 = measure_epvs_and_vessel_constant(cd68,mrad,fout,hflag,pflag);

%% CD68 Statistics
% One-sided Wilcoxon signed-rank test
% Hypothesize increased scattering (CD68) around EPVS
%   --> (med(EPVS) - med(control)) > 0
%   --> right tailed
tail = 'right';
% Wilcoxon signed rank
[cd68_stats] = wilcoxon_epvs_ves(cd68,tail);
% Write to specific sheet named "LHE"
writetable(cd68_stats, stat_sheet, 'WriteRowNames', true, 'Sheet', 'CD68');

%% GFAP
% Stain directory
stain_dir = fullfile(hdir,'GFAP/');
% stain suffix
stain_suffix = '_shrunk_deconv.tif';
% EPVS suffix
epvs_suffix = '_shrunk_EPVS.tif';
% Vessel filename suffix
ves_suffix = '_shrunk_VESSEL.tif';
% Mask suffix
mask_suffix = '_shrunk_mask.tif';

%%% Import LHE stain
[gfap] = import_stain(stain_dir,stain_suffix,epvs_suffix, ...
                    ves_suffix,mask_suffix);

%% GFAP EPVS + Vessel measurements
fout = [];
gfap = measure_epvs_and_vessel_constant(gfap,mrad,fout,hflag,pflag);

%% GFAP Statistics
% One-sided Wilcoxon signed-rank test
% Hypothesize increased scattering (GFAP) around EPVS
%   --> (med(EPVS) - med(control)) > 0
%   --> right tailed
tail = 'right';
% Wilcoxon signed rank
[gfap_stats] = wilcoxon_epvs_ves(gfap,tail);
% Write to specific sheet named "LHE"
writetable(gfap_stats, stat_sheet, 'WriteRowNames', true, 'Sheet', 'gfap');

%% Gallyas

%%% Top-Level
% directory to Gallyas stains
gallyas_dir=['/autofs/cluster/octdata3/users/mjhyman' ...
    '/oct_caa_analyses/histology/Gallyas/'];
% Subfolders for each subject
subdirs = {'CAA6_frontal_61_Gallyas','CAA6_occipital_57_Gallyas',...
    'CAA25_frontal_65_Gallyas_5x',...
    'CAA25_occipital_81_Gallyas_5x','CAA26_frontal_48_Gallyas_5x',...
    'CAA26_occipital_120_Gallyas_5x'};
% Initialize struct array for storing histology
gal = struct();

%%% Import
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

    %%% Invert the stain
    % Darker = more stain, but appears as lower number in Matlab
    % We want the opposite trend where higher number denotes more stain.
    stain = 255 - stain;

    %%% Check if the images have 3 channels, if so take first
    mask = single_ch(mask);
    epvs = single_ch(epvs);
    ves = single_ch(ves);
    ret = single_ch(ret);

    % Store in struct
    gal(ii).baseName = string(base_name);
    gal(ii).mask = mask;
    gal(ii).image = stain;
    gal(ii).epvs = epvs;
    gal(ii).ves = ves;
    gal(ii).ret = ret;
end

%% Statistics: histology vs. retardance
% Define constant for dilating donut from the inner radius
r = 2;
% Define vector of radii (pixels)
radii = 0:r:20;
% Measure at each distance and measure correlation
gal_stat = corr_histo_oct(gal,radii,r);
% Scatter plot + stats for each subject
scatter_histo_oct(gal_stat,figdir);

%% Combine all subjects at each distance
% Limits for scatter plots
xl = [-0.2, 0.65];
yl = [30, 32.2];
gal_stats_table = corr_region_distance(gal_stat,xl,yl,figdir);

%% Combine all subjects within each stain (separate EVPS & Vessel)

[epvs, ves] = stain_epvs_ves(lhe, 'LHE');
fprintf('LHE: Median EPVS = %f\n',median(epvs,'omitnan'));
fprintf('LHE: Median Ves = %f\n',median(ves,'omitnan'));
fprintf('LHE: Mean EPVS = %f\n',mean(epvs,'omitnan'));
fprintf('LHE: Mean Ves = %f\n',mean(ves,'omitnan'));
[epvs, ves] = stain_epvs_ves(cd68, 'CD68');
fprintf('CD68: Median EPVS = %f\n',median(epvs,'omitnan'));
fprintf('CD68: Median Ves = %f\n',median(ves,'omitnan'));
fprintf('CD68: Mean EPVS = %f\n',mean(epvs,'omitnan'));
fprintf('CD68: Mean Ves = %f\n',mean(ves,'omitnan'));
[epvs, ves] = stain_epvs_ves(gfap, 'GFAP');
fprintf('GFAP: Median EPVS = %f\n',median(epvs,'omitnan'));
fprintf('GFAP: Median Ves = %f\n',median(ves,'omitnan'));
fprintf('GFAP: Mean EPVS = %f\n',mean(epvs,'omitnan'));
fprintf('GFAP: Mean Ves = %f\n',mean(ves,'omitnan'));

%% Create figure overlaying z-score with donuts

% Create first layer of z-score stain
zstain = lhe(2).z_stain;
mask = lhe(2).mask;
zstain(~mask) = -4;
figure;
imagesc(zstain); hold on;

% Overlay with annotation
epvs = lhe(2).epvs;
se1 = strel('disk',0);
se2 = strel('disk',mrad);
inner = imdilate(epvs,se1);
outter = imdilate(epvs,se2);
annot = xor(inner, outter);
overlay = cat(3, ones(size(zstain)), zeros(size(zstain)), zeros(size(zstain)));
h = imagesc(overlay);
% set(h,'AlphaData',0);
set(h,'AlphaData',annot * 0.5);
hold off
set(gca,'XTick',[]); set(gca,'YTick',[])
colorbar
set(gca,'FontSize',20)

%% Function to keep first channel from logical
% some of the imported images have three channels, but they are all
% logical. this is due to an issue from Fiji. Just keep one of these
% channels
function im = single_ch(im)
if ~ismatrix(im)
    im = im(:,:,1);
end
end

%% Function to create overlapping histograms of individual donuts
% Purpose: compare distributions b/w EPVS and Vessel to determine whether
% there is a difference in distributions. 

function [epvs, ves] = stain_epvs_ves(stain, tstring)
% Overlay histograms of the distributions of the EPVS & Vessel donuts
% INPUTS:
%   stain (struct): contains measurements
%   tstring (string): title of histogram

epvs = [];
ves = [];
for ii = 1:length(stain)
    epvs = [epvs, stain(ii).exp];
    ves = [ves, stain(ii).ctl];
end
pause(0.1)

% Plot the first histogram (red)
figure;
histogram(epvs, 'Normalization', 'pdf', 'FaceColor', [1 0 0],...
                'FaceAlpha', 0.5);
hold on
% Plot the second histogram (blue)
histogram(ves, 'Normalization', 'pdf', 'FaceColor', [0 0 1],...
                'FaceAlpha', 0.5);
hold off
legend('EPVS','Vessel')
xlabel('Z-Score Normalized Value')
ylabel('Probability Density')
title(tstring)
set(gca, 'FontSize', 20)

end