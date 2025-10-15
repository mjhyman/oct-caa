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
% Add top-level directory + subdirectories
parentDir = fileparts(pwd);
fsDir = fullfile(parentDir, 'freesurfer');
cstatDir = fullfile(parentDir, 'CircStat2012a');
addpath(fsDir);
addpath(cstatDir);
addpath(parentDir);

% Voxel dimensions (microns) for all runs
res = [20,20,20]; % resolution in microns

%%% Flag for loading .MAT struct and creating WM masks
% flag for reloading the .MAT struct for each subject
flag_load_caa_structs = true;

%%% Directories on Martinos Center w/ Matlab struct
% data_dir = ['/autofs/cluster/octdata3/users/mjhyman/' ...
%     'oct_caa_analyses/optical_properties'];

%%% Directories on SCC w/ Matlab struct
data_dir = '/projectnb/npbssmic/ns/CAA/';

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

%%% Combine all subjects into single struct for ease
subjects = struct();
subjects.caa6 = caa6;
subjects.caa17 = caa17;
subjects.caa22 = caa22;
subjects.caa25 = caa25;
subjects.caa26 = caa26;

%% Measure from edge of ves -> edge of enclosed cylinder
% The purpose of this is to identify the optimal radius that distinguishes
% a statistical difference between pathological vs non-pathological

%%% Iterate over radius distance for measuring parenchyma
% Thickness of ring in microns
th = 40;
% Maximum distance of ring
dmax = 500;
dmax = floor(dmax/th) .* th;
% radii in units of microns
% radii = [100, 200, 260, 300, 400, 500];
radii = 40 : 40 : dmax;
% radii in units of voxels
radii = radii ./ res(1);
% Radii to include for generating "donut" masks
radii_include = radii;
% Thickness of ring in voxels
th = th ./ res(1);
% Minimum number of voxels per group
n_min = 50;
% Call parench over all radii
parench = parse_caa_measure_parenchyma(subjects,radii,radii_include,...
                                       th,data_dir,res,n_min);
% Backup struct
fout = fullfile(data_dir,'parenchyma_optical_properties_40um_thick_09Oct2025.mat');
fprintf('\nFinished measuring parenchyma\n')
fprintf('\nStarting to Save .MAT to\n%s',fout)
save(fout,"parench",'-v7.3');
fprintf('\nFinished saving .MAT')

%% Function to load data
function [parench] = parse_caa_measure_parenchyma( ...
                        subjects,radii,radii_include,th,dpath,res,n_min)
% Function to process CAA data for multiple subjects and regions.
% INPUTS:
%   subjects (struct): struct containing substruct for each subject
%   radii (array): radii of dilation in terms of voxels
%   radii_include (int vector): these are the radii (in voxels) that will
%                               be used to generate the "donut" masks
%                               (and saved as .nii files)
%   th (double): thickness of segmentation ring
%   dpath (string): file path for saving the donut masks
%   res (array): resolution of voxels (microns)
%   n_min (int): minimum number of elements in volume to retain
% OUTPUTS:
%   parench (struct): the structure containing the optical properties for
%                     each subject/region. For example, it will appear as:
%                   parench
%                       -> subID
%                           -> region (occip/front)
%                               -> radius1
%                                   -> ves/epvs
%                                       -> pmus (scattering)
%                                       -> pret (retardance)
%                                       -> pori (orientation)
%                               -> radius2
%                               -> ...
%                               -> radiusN

% retrieve subject IDs
subs = fieldnames(subjects);
parench = struct();

% Loop over each subject
for i = 1:length(subs)
    % subject ID
    sub = subs{i};  
    fprintf('STARTING sub %s\n',sub)
    % regions for this subject
    regions = fieldnames(subjects.(sub));
    % Loop over each region (e.g., 'front', 'occip')
    for j = 1:length(regions)
        % Region (front or occip)
        loc = regions{j};  
        fprintf('Starting %s %s \n', sub, loc);
        % Load the local parameters dynamically
        [mus, ret, seg, seg_wm, ori, epvs] = ...
            load_local_params(subjects.(sub), loc);
        
        %%% Discard segments and EPVS with fewer than N voxels (these are
        % likely false positives)
        seg = keep_large(seg, n_min);
        epvs = keep_large(epvs, n_min);

        % Exclude ves + EPVS from mus and retardance to ensure
        % parenchyma measurements are only of parenchyma
        mus(seg) = NaN;
        ret(seg) = NaN;
        ori(seg) = NaN;
        if ~isempty(epvs)
            mus(epvs) = NaN;
            ret(epvs) = NaN;
            ori(epvs) = NaN;
        end

        %%% Exclude vessels overlapping with EPVS if epvs exists
        if ~isempty(epvs)
            [seg_wm_no_epvs] = exclude_epvs_ves(seg_wm, epvs);
        else
            seg_wm_no_epvs = seg_wm;
        end
        
        %%% Loop over the radii
        for k = 1:length(radii)
            fprintf('  --starting radius %d of %d\n', k, length(radii));
            %%% Create the structuring elements for two dilations
            % se1 is the inner dilation
            se1 = strel('disk',radii(k));
            % se2 = "th" voxels greater than > inner dilation
            se2 = strel('disk',radii(k)+th);
            % Measure optical properties of tissue surrounding vessels
            [mus_out,ret_out,ori_out,mus_in,ret_in,ori_in,pout,pin] = ...
                parenchyma_optical_props(seg_wm_no_epvs,mus,ret,ori,se1,se2);
            
            %%% Add to the struct
            % radius string
            rad_str = ['rad',num2str(radii(k))];
            % add the inner cylinder to struct
            parench.(sub).(loc).(rad_str).inner.ves.pmus = mus_in;
            parench.(sub).(loc).(rad_str).inner.ves.pret = ret_in;
            parench.(sub).(loc).(rad_str).inner.ves.pori = ori_in;
            % add the outer cylinder to struct
            parench.(sub).(loc).(rad_str).outter.ves.pmus = mus_out;
            parench.(sub).(loc).(rad_str).outter.ves.pret = ret_out;
            parench.(sub).(loc).(rad_str).outter.ves.pori = ori_out;
            % Save the parenchyma ring if its smallest, 250, or largest
            if ismember(radii(k),radii_include)
                %%% Inner
                % set the filename
                fname = strcat(rad_str, '_inner.nii');
                % Set the save path
                fname = fullfile(dpath,sub,loc,'/ves_donuts/',fname);
                % Save as nifti
                fprintf('Saving inner donut to NII: %s, %s, %s\n',...
                    sub,loc,rad_str)
                save_mri(pin,fname,res./1000,'uchar',0);
                %%% Outer
                % set the filename
                fname = strcat(rad_str, '_outer.nii');
                % Set the save path
                fname = fullfile(dpath,sub,loc,'/ves_donuts/',fname);
                % Save as nifti
                fprintf('Saving outer donut to NII: %s, %s, %s\n',...
                    sub,loc,rad_str)
                save_mri(pout,fname,res./1000,'uchar',0);
            end
            
            %%% If EPVS exists, process EPVS parenchyma as well
            if ~isempty(epvs)
                %%% Remove vessel segmentation (ensure independence)
                mus_tmp = mus;
                ret_tmp = ret;
                ori_tmp = ori;
                mus_tmp(pout) = NaN;
                ret_tmp(pout) = NaN;
                ori_tmp(pout) = NaN;

                % Measure optical properties
                [mus_out,ret_out,ori_out,mus_in,ret_in,ori_in,pout,pin] = ...
                    parenchyma_optical_props(epvs,mus_tmp,ret_tmp,ori_tmp,...
                                             se1,se2);
                % add the inner cylinder to struct
                parench.(sub).(loc).(rad_str).inner.epvs.pmus = mus_in;
                parench.(sub).(loc).(rad_str).inner.epvs.pret = ret_in;
                parench.(sub).(loc).(rad_str).inner.epvs.pori = ori_in;
                % add the outer cylinder to struct
                parench.(sub).(loc).(rad_str).outter.epvs.pmus = mus_out;
                parench.(sub).(loc).(rad_str).outter.epvs.pret = ret_out;
                parench.(sub).(loc).(rad_str).outter.epvs.pori = ori_out;
                % Save the parenchyma ring if its smallest, 250, or largest
                if ismember(radii(k),radii_include)
                    %%% Inner
                    % set the filename
                    fname = strcat(rad_str, '_inner.nii');
                    % Set the save path
                    fname = fullfile(dpath,sub,loc,'/epvs_donuts/',fname);
                    % Save as nifti
                    fprintf('Saving inner donut to NII: %s, %s, %s\n',...
                        sub,loc,rad_str)
                    save_mri(pin,fname,res./1000,'uchar',0);
                    %%% Outer
                    % set the filename
                    fname = strcat(rad_str, '_outer.nii');
                    % Set the save path
                    fname = fullfile(dpath,sub,loc,'/epvs_donuts/',fname);
                    % Save as nifti
                    fprintf('Saving outer donut to NII: %s, %s, %s\n',...
                        sub,loc,rad_str)
                    save_mri(pout,fname,res./1000,'uchar',0);
                end
            end
            fprintf('  --finished radius %d of %d\n', k, length(radii));
        end
        % Finished processing for this subject and region
        fprintf('Finished %s %s\n', sub, loc);
    end
end
end

function seg_parsed = keep_large(seg, nmin)
% keepLargeComponents - Keeps connected components with at least nmin voxels (vectorized)
%
% Syntax: seg_parsed = keepLargeComponents(seg, nmin)
%
% Inputs:
%    seg - 3D binary matrix (segmentation)
%    nmin - Minimum number of voxels for a component to be retained
%
% Outputs:
%    seg_parsed - 3D binary matrix with small components removed

% Identify connected components in 3D
CC = bwconncomp(seg, 26); % 26-connectivity for 3D

% Measure component sizes
componentSizes = cellfun(@numel, CC.PixelIdxList);

% Find components that meet size threshold
largeComponentsIdx = componentSizes >= nmin;

% Concatenate indices of large components
voxelsToKeep = vertcat(CC.PixelIdxList{largeComponentsIdx});

% Create the output segmentation
seg_parsed = false(size(seg));
seg_parsed(voxelsToKeep) = true;
end
