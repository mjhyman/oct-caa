%% Import downsampled stains
%{
These are the stains that our collaborator performed. The grayscale images
can be directly used to compute the stain density surrounding the EPVS
because there is no counterstain. The deconvolved images had a
counterstain and therefore had to be deconvolved to extract the target
stain.

The purpose of this script is to create the images for the histology
workflow

%}

%% Top-level settings
clearvars -except caa6 caa17 caa22 caa25 caa26 swp_struct
clc; close all;
% Get the current folder
currentFolder = pwd;
% Move one directory up
parentFolder = fileparts(currentFolder);
% Add the parent folder and all its subfolders to the MATLAB search path
addpath(genpath(parentFolder));

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/figures/histology_scatterplots/';

%%% Load the stain and segmentation structs
lhe = load(fullfile(hdir,'lhe_rings_15-Oct-2025.mat')); lhe = lhe.lhe;
gfap = load(fullfile(hdir,'gfap_rings_15-Oct-2025.mat')); gfap = gfap.gfap;
cd68 = load(fullfile(hdir,'cd68_rings_15-Oct-2025.mat')); cd68 = cd68.cd68;

%%% measurement outter radii (units = voxels)
radii_sm = [40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486];
radii_lg = [100, 201, 301, 402, 503];

%% Scatterplot of histology vs. distance
close all;
% LHE
[sm, lg] = mean_std_stain(lhe, radii_sm, radii_lg);
stain_errorbar(sm, lg, radii_sm, radii_lg, 'LHE', figdir)

% Import GFAP
[sm, lg] = mean_std_stain(gfap, radii_sm, radii_lg);
stain_errorbar(sm, lg, radii_sm, radii_lg, 'GFAP', figdir)

% Import CD68
[sm, lg] = mean_std_stain(cd68, radii_sm, radii_lg);
stain_errorbar(sm, lg, radii_sm, radii_lg, 'CD68', figdir)

%% Function to scatterplot
function stain_errorbar(sm, lg, radii_sm, radii_lg, tstr, figdir)
% Offset for EPVS and vessels
offset = 10;
% Errorbar line width
lwidth = 2.5;

%%% SMALL RADII    
% Prepare figure for scatterplot
figure; hold on;
xlabel('Distance (\mum)');
ylabel('Stain Density (z-score)');
title(strcat(tstr, ' Small Radii (40 \mum)'));
% EPVS
errorbar(radii_sm, sm.meanEpvs, sm.stdEpvs, 'o','LineWidth',lwidth,...
        'color','red','DisplayName', 'EPVS');
hold on;
% Vessels
errorbar(radii_sm + offset, sm.meanVes, sm.stdVes, 'o','LineWidth',lwidth,...
        'color','blue','DisplayName', 'vessel');
hold off;
legend show;
% Increase size of dots and bars
set(gca, 'FontSize', 40);
set(findobj(gca, 'Type', 'scatter'), 'SizeData', 100);
xlim([0,550]);
% Save output
fout = fullfile(figdir,strcat(tstr,'_small_radii_scatter'));
pause(1); saveas(gca,fout,'png')

%%% LARGE RADII    
% Prepare figure for scatterplot
figure; hold on;
xlabel('Distance (\mum)');
ylabel('Stain Density (z-score)');
title(strcat(tstr, ' Large Radii (100 \mum)'));
% EPVS
errorbar(radii_lg, lg.meanEpvs, lg.stdEpvs, 'o','LineWidth',lwidth,...
        'color','red','DisplayName', 'EPVS');
hold on;
% Vessels
errorbar(radii_lg + offset, lg.meanVes, lg.stdVes,'LineWidth',lwidth,...
        'color','blue','DisplayName', 'vessel');
hold off;
legend show;
% Increase size of dots and bars
set(gca, 'FontSize', 40);
set(findobj(gca, 'Type', 'scatter'), 'SizeData', 100);
xlim([0,550]);
% Save output
fout = fullfile(figdir,strcat(tstr,'_large_scatter'));
pause(1); saveas(gca,fout,'png')
end

%% Function to iterate stain and take mean + Std Dev across subjects
function [sm, lg] = mean_std_stain(stain, radii_sm, radii_lg)
% MEAN_STD_STAIN measure mean Std. across stains
% INPUTS:
%   stain (struct): stain struct
%   radii_sm (vector): vector of small radii
%   radii_lg (vector): vector of large radii
% OUTPUTS:
%   sm (struct): small radii struct
%   lg (struct): large radii struct


% Extract # radii for each radii
nsmall = length(radii_sm);
nlarge = length(radii_lg);

% Measure small radii
[sm.meanEpvs, sm.stdEpvs, sm.meanVes, sm.stdVes] =...
    iterate_stains(stain, nsmall, 'rad40');

% Measure large radii
[lg.meanEpvs, lg.stdEpvs, lg.meanVes, lg.stdVes] =...
    iterate_stains(stain, nlarge, 'rad100');




function [meanEpvs, stdEpvs, meanVes, stdVes] =...
         iterate_stains(stain, nrad, rad_size)
% ITERATE_STAINS measure mean + Std across stains
% INPUTS
%   stain (struct): stain struct
%   nrad (int): number of radii
%   rad_size (string): sub-field string
% OUTPUTS:
%   meanEpvs (vector): average EPVS
%   stdEpvs (vector): std EPVS
%   meanVes (vector): average Ves
%   stdVes (vector): std Ves

% Initialize vectors
meanEpvs = zeros(nrad,1);
stdEpvs = zeros(nrad,1);
meanVes = zeros(nrad,1);
stdVes = zeros(nrad,1);
% Retrieve names of the fields of small radii
rad_name = fields(stain(1).(rad_size));
% Iterate over the small radii
for ii = 1:nrad
    % EPVS
    epvs = arrayfun(@(s) s.(rad_size).(rad_name{ii}).exp_mean, stain);
    % Vessels
    ves = arrayfun(@(s) s.(rad_size).(rad_name{ii}).ctl_mean, stain);
    % Take mean & std dev across subjects
    meanEpvs(ii) = mean(epvs);
    stdEpvs(ii) = std(epvs);
    meanVes(ii) = mean(ves);
    stdVes(ii) = std(ves);
end

end
end