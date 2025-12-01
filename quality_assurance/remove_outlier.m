%% Remove outliers
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
    - optical property donut measurements at each radii
    - inner and outter
    - each radius
- REMOVE outliers
    - mus range = [0,25]
    - retardance range = [0,45]
- Save as .MAT with latest date (03Nov2025) to match the spreadsheet
    "caa_all_radii_40um_donut_03Nov2025.xlsx"
%}

%% Prepare environment
clear; clc; close all;
% Add top-level directory + subdirectories
addpath(genpath(fullfile(pwd, '..')))
% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/projectnb/npbssmic/ns/CAA/';
% Voxel dimensions (microns) for all runs
res = [20,20,20]; % resolution in microns

% Import the spreadsheet
fname = 'parenchyma_optical_properties_40um_thick_09Oct2025.mat';
parench = load(fullfile(data_dir, fname));
parench = parench.parench;

% Output .MAT name
fout = 'parenchyma_optical_properties_40um_thick_03Nov2025.mat';

% Outlier threshold
max_pmus = 25;
max_pret = 45;

%% Remove outliers from mus and retardance

% Iterate subjects
subs = fieldnames(parench);
for ii = 1:numel(subs)
    % Iterate regions (front and occip)
    regions = fieldnames(parench.(subs{ii}));
    for j = 1:numel(regions)
        % Iterate radii
        radii = fieldnames(parench.(subs{ii}).(regions{j}));
        for k = 1:numel(radii)
            % Assign local variables
            sub = subs{ii};
            reg = regions{j};
            rad = radii{k};
            % Vessels
            parench = rm_struct(parench, sub, reg, rad, 'ves',...
                                max_pmus, max_pret);
            % EPVS
            parench = rm_struct(parench, sub, reg, rad, 'epvs',...
                                max_pmus, max_pret);
        end
    end
end

%% Save as .MAT
fout = fullfile(data_dir, fout);
save(fout,"parench",'-mat','-v7.3');

%% Function to identify outliers, remove, add back to struct

function parench = rm_struct(parench, sub, reg, rad, loc, max_pmus, max_pret)
% REMOVE_FROM_STRUCT
% parench (struct): top-level structure
% sub (string): subject
% reg (string): region
% rad (string): radius
% loc (string): ves or epvs
% max_pmus (int): maximum scattering coefficient
% max_pret (int): maximum retardance

% Select mus and retardance for this one
mus = parench.(sub).(reg).(rad).outter.(loc).pmus;
ret = parench.(sub).(reg).(rad).outter.(loc).pret;

% Remove outliers based on defined thresholds
mus_keep = mus < max_pmus;
ret_keep = ret < max_pret;

% Store cleaned data back into the structure
parench.(sub).(reg).(rad).outter.(loc).pmus = mus(mus_keep);
parench.(sub).(reg).(rad).outter.(loc).pret = ret(ret_keep);

% Print the number of removed outliers for each location
fprintf('Subject: %s, Region: %s, Radius: %s, Location: %s - Removed %d outliers from mus and %d from retardance.\n', ...
    sub, reg, rad, loc, sum(~mus_keep), sum(~ret_keep));

end