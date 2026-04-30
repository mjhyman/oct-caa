%% Measure size-weighted proximity (SWP) of stains
% Purpose: plot stain vs. log(SWP) for each section and average
% Overview:
%{
- import each stain's SWP struct
- plot stain vs. SWP (EPVS and vessel) for each stain
- average for each subject
%}

%% Top-level settings
clc; clear; close all;

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/swp';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/histology/gam';

%%% Add top-level directory to path
addpath(genpath('/projectnb/npbssmic/s/mhyman/oct-caa'));

%%% Filenames
lhe_fname = 'histo_swp_hmatched_zscore__LHE_voxelwise_radius4mm_ds4_19-Mar-2026.mat';
gfap_fname = 'histo_swp_hmatched_zscore__GFAP_voxelwise_radius4mm_ds4_20-Mar-2026.mat';
cd68_fname = 'histo_swp_hmatched_zscore__CD68_voxelwise_radius4mm_ds4_20-Mar-2026.mat';

%%% Output figure settings
% x-axis binning number of points
nbins = 100;
% Date for file output
t = datetime('today');
% Font size for figures
fsize = 30;

%% Import stains

% % Import LHE
lhe = load(fullfile(hdir, lhe_fname));
lhe = lhe.lhe;
% Import GFAP
gfap = load(fullfile(hdir, gfap_fname));
gfap = gfap.gfap;
% Import CD68
cd68 = load(fullfile(hdir, cd68_fname));
cd68 = cd68.cd68;

% Add all stains to the struct "stain"
stains = struct();
stains.cd68 = cd68;
stains.gfap = gfap;
stains.lhe = lhe;

%% Create Matrix of ves-swp, epv-swp, stain for each section
% The next code block will combine all sections within the same severity
% This is the order in the matrix:
%   EPVS-swp, ves-swp, stain

%%% Iterate stains
stain_names = fields(stains);
for ii = 1:numel(stain_names)
    % Get the current stain name
    stain = stains.(stain_names{ii});
    % print stain name to console
    fprintf('\nRunning stain %s',stain_names{ii})
    %%% Iterate sections within stain
    for j = 1:numel(stain)
        % print section # to console
        fprintf('\n\tSection %i / %i\n',j,numel(stain))
        % Extract vessel SWP + EPVS SWP
        ves_swp = stain(j).ves_swp;
        epvs_swp = stain(j).epvs_swp;        
        
        %%% Retrieve tissue mask, ves mask, epvs mask
        % Retrieve mask, vessels, EPVS
        mask = stain(j).mask;
        ves = stain(j).ves;
        epvs = stain(j).epvs;
        % Set the EPVS and vessels to 0 within mask
        mask(ves) = 0;
        mask(epvs) = 0;
        
        %%% Apply mask to SWP and stain section
        ves_swp = ves_swp(mask);
        epvs_swp = epvs_swp(mask);
        section = stain(j).z_stain(mask);
        
        %%% Create matrix to store SWP and save
        epv_ves_swp_stain = [epvs_swp, ves_swp, section];
        stains.(stain_names{ii})(j).epv_ves_swp_stain = epv_ves_swp_stain;
    end
end

%% Combine severe (caa22 and caa25)
% Iterate the histo_swp struct + concatenate same severity
% Control: caa26
% Mild: caa6
% Moderate: caa17
% Severe: caa22, caa25

% Initialize struct for storing epv-swp, ves-swp, stain
histo_swp = struct();

% Iterate stains
stain_names = fields(stains);
for ii = 1:numel(stain_names)
    % Get the current stain name
    stain = stains.(stain_names{ii});
    % print stain name to console
    fprintf('\nRunning stain %s',stain_names{ii})

    %%% Find sections for each subject/region
    % Retrieve all sections
    section_names = {stain.baseName};
    % string patterns
    caa22f_txt = {'CAA22_1','CAA_22_1'};
    caa22o_txt = {'CAA22_7','CAA_22_7'};
    caa25f_txt = {'CAA25_1','CAA_25_1'};
    caa25o_txt = {'CAA25_7','CAA_25_7'};
    % Find matching cell
    idx_22f = find(strcmpi(section_names, caa22f_txt));
    idx_25f = find(strcmpi(section_names, caa25f_txt));
    % Find matching cell
    idx_22o = find(strcmpi(section_names, caa22o_txt));
    idx_25o = find(strcmpi(section_names, caa25o_txt));
    
    %%% Combine Frontal
    % Retrieve 22_front and 25_front
    tmp_22 = stain(idx_22f).epv_ves_swp_stain;
    tmp_25 = stain(idx_25f).epv_ves_swp_stain;
    % Concatenate and add to struct
    tmp = [tmp_22; tmp_25];
    n = length(stain);
    stain(n+1).baseName = 'severe_front';
    stain(n+1).epv_ves_swp_stain = tmp;

    %%% Combine Occipital
    % Retrieve 22_occip and 25_occip
    tmp_22 = stain(idx_22o).epv_ves_swp_stain;
    tmp_25 = stain(idx_25o).epv_ves_swp_stain;
    % Concatenate and add to struct
    tmp = [tmp_22; tmp_25];
    n = length(stain);
    stain(n+1).baseName = 'severe_occip';
    stain(n+1).epv_ves_swp_stain = tmp;

    %%% Add stain back to stains parent struct
    stains.(stain_names{ii}) = stain;
end

%% Create table for each subject/region

% Initialize Table
T = struct();

fprintf('Creating table\n')
vnames = {'epv_swp','ves_swp','stain'};
stain_names = fields(stains);
for ii = 1:numel(stain_names)
    % Get the current stain name
    stain = stains.(stain_names{ii});
    fprintf('\tStarting Stain %s\n',stain_names{ii})
    % Name of each section
    section_names = {stain.baseName};
    % Iterate sections and add to table
    for j = 1:numel(section_names)
        T.(stain_names{ii}).(section_names{j}) =...
            array2table(stain(j).epv_ves_swp_stain,'VariableNames',vnames);
    end
end

%% Create GAM for each stain and section

%%% Sampling and Learning
% Set maximum number of samples
n_max = 1e6;
% Set learning parameters (avoid over fitting)
ntrees = 500;
max_splits = 4;
learn_rate = 0.05;

%%% Confidence intervals for slices
nbootstrap = 100;
ci_alpha = 0.05;

%%% Run GAM for each stain and section
% Initialize GAM
histo_gam = struct();

% Iterate Stain
for ii = 1:numel(stain_names)
    % Name of each section
    stain = stains.(stain_names{ii});
    section_names = {stain.baseName};
    % Iterate sections and add to table
    for j = 1:numel(section_names)
        % Set the subdirectory
        subdir = fullfile(figdir,stain_names{ii});
        % Define title string (which is also filename)
        tstr = strcat(stain_names{ii},'_',section_names{j});
        % Run GAM
        histo_gam.(stain).(section) = fit_swp_histo_gam(...
                    T.(stain_names{ii}).(section_names{j}),...
                    'NumTrees',ntrees,...
                    'MaxSplits',max_splits,'LearnRate',learn_rate,...
                    'MaxSample',n_max,'NBootstrap',nbootstrap,...
                    'CIAlpha',ci_alpha,'TitleStr',tstr,...
                    'dirout',subdir,...
                    'stain_name',stain_names{ii},...
                    'section_name',section_names{j});
    end
end
% Save GAM
dt = datetime("now",'TimeZone','local','Format','dd-MMM-yyyy');
fout = strcat('GAM_histo_',string(dt),'.mat');
full_fout = fullfile(figdir,fout);
fprintf('Saving histo_gam data struct to %s\n',full_fout)
save(full_fout,'histo_gam','-v7.3');