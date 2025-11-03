%% Generate Spreadsheet of Measurements at Each Distance
%{
The purpose of this script is to create a spreadsheet of the optical
properties at each EPVS and Vessel at each distance. This covers just the
following cases:
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)
%}

%% Prepare environment
% clear; clc; close all;
% Add top-level directory
current_dir = pwd;
addpath(fullfile(current_dir));
% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/projectnb/npbssmic/ns/CAA/';
%%% 40 um donut
% load the parenchymal optical properties
load(fullfile(data_dir, ...
    'parenchyma_optical_properties_40um_thick_09Oct2025.mat'));
% Output filename for spreadsheet
spreadsheet_name = 'caa_all_radii_40um_donut_03Nov2025.xlsx';

%% Count number of control and experimental values

% Counter for # samples in vessels
n_ves = 0;
% Counter for # samples in EPVS
n_epvs = 0;
% Initialize radii strings
rads = fields(parench.caa6.front);

%%% Iterate over subjects in parench
subs = fields(parench);
for ii = 1:length(subs)
    % retrieve subject
    sub = subs{ii};
    % retrieve regions for subject
    regs = fields(parench.(sub));
    %%% iterate over regions
    for j = 1:length(regs)
        % retrieve region
        reg = regs{j};
        %%% Iterate over radii
        for r = 1:numel(rads)
            % Initialize radius for iteration
            rad = rads{r};
             % retrieve optical properties (vessels and epvs)
            segmentations = fields(parench.(sub).(reg).(rad).outter);
            %%% Iterate segmentations & add to vectors
            for k=1:length(segmentations)
                % retrieve segmentation
                seg = segmentations{k};
                [mus,~,~,~] = retrieve_op(parench,sub,reg,rad,seg);
                % Count number of measurements in vessel or EPVS
                if strcmp(seg,'ves')
                    n_ves = n_ves + numel(mus);
                elseif strcmp(seg,'epvs')
                    n_epvs = n_epvs + numel(mus);
                end
            end
        end
    end
end

%% Create table containing all measurements across all distances

% Call function to vectorize data
[tbl_mus, tbl_ret, tbl_sori] = vectorize_struct(parench, rad, n_ves, n_epvs);
% Export the tables to CSV
fout = fullfile(data_dir,spreadsheet_name);
writetable(tbl_mus,fout,'Sheet','scattering');
writetable(tbl_ret,fout,'Sheet','retardance');
writetable(tbl_sori,fout,'Sheet','orientation');


%% LMM for each test
function [tbl_mus, tbl_ret, tbl_sori] =...
    vectorize_struct(parench, rad, n_ctl, n_exp)
% LMM_TEST Create linear mixed effects model to compare optical props
% This test will measure across all subjects and compare:
%   - vessels w/ EPVS vs. vessels w/o EPVS 
% INPUTS:
%   check_lin (logical): flag for checking linearity assumptions
%   test_idx (int): index for the respective statistical test
%   flag_subs (int): flag for how to increment the subject ID assignment
%       0 -> iterate subject ID for each tissue volume and EPVS/vessel.
%            The vessels and EPVS within the same tissue volume will have
%            different subject IDs.
%       1 -> iterate subject ID for each tissue volume. Both the vessel and
%            EPVS within the same volume will have the same subID
%   parench (struct): parenchyma optical properties struct
%   rad (string): string indicating EPVS ring radius for structure
%   n_ves (int): number of vessel measurements
%   n_epvs (int): number of EPVS measurements
% OUTPUTS:
%   stats (struct): contains the FixedEffects output tables for each
%                   optical property
%   data_exp (struct): experimental data arrays for each optical property
%   data_cntrl (struct): control data arrays for each optical property
%   p (struct): p-values for: mus, retardance, mean orientation, std dev of
%                       orienation

% Arrays for control (ctl) and experimental (exp)
ctl_mus = zeros(n_ctl,1);
ctl_ret = zeros(n_ctl,1);
ctl_ori = zeros(n_ctl,1);
ctl_reg = cell(n_ctl,1);
ctl_rad = zeros(n_ctl,1);
exp_mus = zeros(n_exp,1);
exp_ret = zeros(n_exp,1);
exp_ori = zeros(n_exp,1);
exp_reg = cell(n_exp,1);
exp_rad = zeros(n_exp,1);
% Index to track location in arrays of data (ctl_* and exp_*)
cidx = 1;
eidx = 1;
% Subject ID vectors
ctl_subid = zeros(n_ctl,1);
exp_subid = zeros(n_exp,1);
% Counter for tissue ID
tiss_idx = 1;
% Initialize radii strings
rads = fields(parench.caa6.front);

%%% Iterate over subjects in parench
subs = fields(parench);
for ii = 1:length(subs)
    % update console
    fprintf('Subject %i out of %i\n',ii,length(subs))
    % retrieve subject
    sub = subs{ii};
    % retrieve regions for subject
    regs = fields(parench.(sub));
    %%% iterate over regions
    for j = 1:length(regs)
        % retrieve region
        reg = regs{j};
        % retrieve optical properties (vessels and epvs)
        segmentations = fields(parench.(sub).(reg).(rad).outter);
        %%% Iterate over radii
        for r = 1:numel(rads)
            % Initialize radius for iteration
            rad = rads{r};
            %%% Iterate segmentations & add to vectors
            for k=1:length(segmentations)
                % retrieve segmentation
                seg = segmentations{k};
                [mus,ret,ori,n] = retrieve_op(parench,sub,reg,rad,seg);
                %%% Vessel
                if strcmp(seg,'ves')
                    % copy radius (in microns) to radii vector
                    radius_val = sscanf(rad, 'rad%d');
                    distance = radius_val * 20;
                    ctl_rad(cidx:cidx+n-1,1) = ones(n,1) .* distance;
                    % copy optical properties, subid, region
                    [ctl_mus,ctl_ret,ctl_ori,ctl_subid,ctl_reg,cidx] =...
                        copy_to_vector(mus, ret, ori,...
                                        ctl_mus, ctl_ret,...
                                        ctl_ori, ctl_subid,...
                                        ctl_reg, cidx, n, tiss_idx, reg);
                %%% EPVS
                elseif strcmp(seg,'epvs')
                    % copy radius (in microns) to radii vector
                    radius_val = sscanf(rad, 'rad%d');
                    distance = radius_val * 20;
                    exp_rad(eidx:eidx+n-1,1) = ones(n,1) .* distance;
                    % copy optical properties, subid, region
                    [exp_mus,exp_ret,exp_ori,exp_subid,exp_reg,eidx] =...
                        copy_to_vector(mus, ret, ori,...
                                        exp_mus, exp_ret,...
                                        exp_ori, exp_subid,...
                                        exp_reg, eidx, n, tiss_idx, reg);
                end
            end
        end
        % Iterate tissue sample counter
        tiss_idx = tiss_idx + 1;
    end
end

%% Create a table for fitting the LME model
% Column 1 = group label
% Column 2 = region
% Column 3 = subject ID
% column 4 = distance (micron)
% Column 5 = value of measurement (mus, ret, ori)

%%% group label column
g_exp = repmat({'experimental'},[length(exp_mus),1]);
g_cnrtl = repmat({'control'},[length(ctl_mus),1]);
group = vertcat(g_exp, g_cnrtl);

%%% Optical property
% Scattering coefficient
mus = [exp_mus; ctl_mus];
% Retardance
ret = [exp_ret; ctl_ret];
% Circular standard deviation of orientation
ori = [exp_ori; ctl_ori];

%%% Region, Subject ID, Distance, 
% Combine the brain region indices into column vector
region = vertcat(exp_reg,ctl_reg);
% Combine the subject IDs into column vector
subids = vertcat(exp_subid, ctl_subid);
% Distance vector
radii = [exp_rad;ctl_rad];

%%% Create Table for each optical property
tbl_mus = table(group,region,subids,radii,mus,...
      'VariableNames',{'Groups','Region','subID','distance','OpticalProperty'});
tbl_ret = table(group,region,subids,radii,ret,...
      'VariableNames',{'Groups','Region','subID','distance','OpticalProperty'});
tbl_sori = table(group,region,subids,radii,ori,...
      'VariableNames',{'Groups','Region','subID','distance','OpticalProperty'});
% Take real part of complex numbers
tbl_sori.OpticalProperty = real(tbl_sori.OpticalProperty);
end













