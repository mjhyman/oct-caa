%% Analyze the EPVS SWP heatmap
% Import the EPVS SWP heatmaps
% Import the .MAT of the optical properties
% Measure correlation between SVP vs. optical properties

%% Add top-level directory of code repository to path
clearvars -except caa6 caa17 caa22 caa25 caa26 swp_struct
clc; close all;
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
ncores = feature('numcores');
maxNumCompThreads(ncores);
% Directory containing seg, mus, ret, mask, epvs structs
mat_dir = '/projectnb/npbssmic/ns/CAA/';
% Heatmap directory
swp_dir = '/projectnb/npbssmic/ns/CAA/swp';
% Scatterplot directory
plt_dir = '/projectnb/npbssmic/ns/CAA/swp/plots';

%%% Flags for importing data
% Flag for SWP structs
flag_load_swp_structs = false;
% Flag for loading CAA structs (false if already in environment)
flag_load_caa_structs = false;

%% Load EPVS SWP

if flag_load_swp_structs
    % Struct for storing SWP
    swp_struct = struct();
    
    % CAA 6 Frontal
    fprintf('Importing SWP for CAA6 Frontal\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa6/front/' ...
                'caa6_front_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa6f = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa6/front/' ...
                'caa6_front_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa6f = tmp.swp;
    
    % CAA 6 Occipital
    fprintf('Importing SWP for CAA6 Occip\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa6/occip/' ...
                'caa6_occip_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa6o = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa6/occip/' ...
                'caa6_occip_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa6o = tmp.swp;
    
    % CAA 17 Occipital
    fprintf('Importing SWP for CAA17 Occip\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa17/occip/' ...
                'caa17_occip_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa17o = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa17/occip/' ...
                'caa17_occip_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa17o = tmp.swp;
    
    % CAA 22 front
    fprintf('Importing SWP for CAA22 Front\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa22/front/'...
                'caa22_front_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa22f = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa22/front/' ...
                'caa22_front_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa22f = tmp.swp;
    
    % CAA 22 Occip
    fprintf('Importing SWP for CAA22 Occip\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa22/occip/' ...
                'caa22_occip_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa22o = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa22/occip/' ...
                'caa22_occip_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa22o = tmp.swp;
    
    % CAA 25 Front
    fprintf('Importing SWP for CAA25 Front\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa25/front/' ...
                'caa25_front_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa25f = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa25/front/' ...
                'caa25_front_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa25f = tmp.swp;
    
    % CAA 25 Occip
    fprintf('Importing SWP for CAA25 Occip\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa25/occip/' ...
                'caa25_occip_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa25o = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa25/occip/' ...
                'caa25_occip_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa25o = tmp.swp;
    
    % CAA 26 Front
    fprintf('Importing SWP for CAA26 Front\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa26/front/' ...
                'caa26_front_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa26f = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa26/front/' ...
                'caa26_front_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa26f = tmp.swp;
    
    % CAA 26 Occip
    fprintf('Importing SWP for CAA26 Occip\n')
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa26/occip/' ...
                'caa26_occip_radius_200_exp_2_interpolated_heatmap.mat']);
    swp_struct.raw.caa26o = tmp.interpolated_volume;
    tmp = load(['/projectnb/npbssmic/ns/CAA/swp/caa26/occip/' ...
                'caa26_occip_radius_200_exp_2_interpolated_heatmap_log10.mat']);
    swp_struct.log10.caa26o = tmp.swp;
end

%% Load matlab structs

if flag_load_caa_structs
    fprintf('Loading CAA6\n')
    caa6 = load(fullfile(mat_dir,"/caa6/caa6.mat"));
    fprintf('Finished loading CAA17\n')
    
    fprintf('Loading CAA17\n')
    caa17 = load(fullfile(mat_dir,"/caa17/occip/caa17.mat"));
    fprintf('Finished loading CAA17\n')
    
    fprintf('Loading CAA22\n')
    caa22 = load(fullfile(mat_dir,"/caa22/caa22.mat"));
    fprintf('Finished loading CAA22\n')
    
    fprintf('Loading CAA25\n')
    caa25 = load(fullfile(mat_dir,"/caa25/caa25.mat"));
    fprintf('Finished loading CAA25\n')
    
    fprintf('Loading CAA26\n')
    caa26 = load(fullfile(mat_dir,"/caa26/caa26.mat"));
    fprintf('Finished loading CAA26\n')
    
    % Remove top-level struct
    caa6 = caa6.caa6;
    caa17 = caa17.caa17;
    caa22 = caa22.caa22;
    caa25 = caa25.caa25;
    caa26 = caa26.caa26;
end

%% Flag for subtracting subject median offset
rm_offset = true;
if rm_offset
    % Load the parenchyma median offsets
    load(fullfile(mat_dir, "median_white_matter_values_14-Nov-2025.mat"));
end

%% Create 2D arrays of optical property vs. EPVS density
% 2xN Matrix for each optical property:
%   - top row = EPVS density from heatmap
%   - bottom row = optical property from heatmap
%   - column = pairwise observations from same tissue volume

fprintf('Creating 2D arrays of EPVS density vs. optical prop\n')

% struct for storing pairs
heat_pairs = struct();

% struct for iterating over subjects
subjects = struct();
subjects.caa6 = caa6;
subjects.caa17 = caa17;
subjects.caa22 = caa22;
subjects.caa25 = caa25;
subjects.caa26 = caa26;

%%% add EPVS swp to struct
subs = fields(subjects);
subjects.caa6.front.swp.raw = swp_struct.raw.caa6f;
subjects.caa6.occip.swp.raw = swp_struct.raw.caa6o;
subjects.caa17.occip.swp.raw = swp_struct.raw.caa17o;
subjects.caa22.front.swp.raw = swp_struct.raw.caa22f;
subjects.caa22.occip.swp.raw = swp_struct.raw.caa22o;
subjects.caa25.front.swp.raw = swp_struct.raw.caa25f;
subjects.caa25.occip.swp.raw = swp_struct.raw.caa25o;
subjects.caa26.front.swp.raw = swp_struct.raw.caa26f;
subjects.caa26.occip.swp.raw = swp_struct.raw.caa26o;

%%% add EPVS log(swp) to struct
subjects.caa6.front.swp.log10 = swp_struct.log10.caa6f;
subjects.caa6.occip.swp.log10 = swp_struct.log10.caa6o;
subjects.caa17.occip.swp.log10 = swp_struct.log10.caa17o;
subjects.caa22.front.swp.log10 = swp_struct.log10.caa22f;
subjects.caa22.occip.swp.log10 = swp_struct.log10.caa22o;
subjects.caa25.front.swp.log10 = swp_struct.log10.caa25f;
subjects.caa25.occip.swp.log10 = swp_struct.log10.caa25o;
subjects.caa26.front.swp.log10 = swp_struct.log10.caa26f;
subjects.caa26.occip.swp.log10 = swp_struct.log10.caa26o;

for ii = 1:length(fields(subjects))
    sub = subs{ii};
    regions = fields(subjects.(subs{ii}));
    for j = 1:length(regions)
        % retrieve local properties
        reg = regions{j};
        if isfield(subjects.(sub).(reg), 'swp')
            fprintf('Starting %s %s\n', sub, reg)
            
            %%% Import SWP, mus, ret
            % Import raw SWP and log(swp)
            swp = subjects.(sub).(reg).swp.raw;
            swplg = subjects.(sub).(reg).swp.log10;
            % Import mus and ret
            mus = subjects.(sub).(reg).mus;
            ret = subjects.(sub).(reg).ret_full;

            %%% Subtract median mus,ret offset
            if rm_offset
                mus = mus - parench_median.(sub).(reg).med.mus;
                ret = ret - parench_median.(sub).(reg).med.ret;
            end
            
            %%% Remove vessels and EPVS from the epvs heatmap
            ves = subjects.(sub).(reg).seg;
            epvs = subjects.(sub).(reg).epvs;
            swp(ves) = NaN;
            swplg(ves) = NaN;
            swp(epvs) = NaN;
            swplg(epvs) = NaN;

            %%% Create pairs for raw SWP
            [mus_pair, ret_pair] = create_pair(swp,mus,ret);
            % add pairs to heatmap struct
            heat_pairs.(sub).(reg).raw.mus_pair = mus_pair;
            heat_pairs.(sub).(reg).raw.ret_pair = ret_pair;

            %%% Create pairs for log(SWP)
            [mus_pair, ret_pair] = create_pair(swplg,mus,ret);
            % add pairs to heatmap struct
            heat_pairs.(sub).(reg).log10.mus_pair = mus_pair;
            heat_pairs.(sub).(reg).log10.ret_pair = ret_pair;
        else
            continue
        end
    end
end

%% Bin each subject

%%% plotting properties
% struct for storing binned vectors
binned = struct();
% flag for whether to plot the binned data
nbin = 190;
xlab = 'EPVS log(SWP)';

%%% Separate the frontal and occipital
% Iterate subjects
for ii = 1:length(subs)
    % Retrieve regions for this subject
    regs = fields(heat_pairs.(subs{ii}));
    % iterate over regions
    for j = 1:length(regs)
        % Print to console
        fprintf('\nBinning subject %s region %s',subs{ii},regs{j})
        % Retrieve mus and retardance LOG10 pairs
        mus_pair = heat_pairs.(subs{ii}).(regs{j}).log10.mus_pair;
        ret_pair = heat_pairs.(subs{ii}).(regs{j}).log10.ret_pair;
        % scattering coefficient
        ylab = '\mu_s (cm^-^1)';
        tit = string(subs{ii}) + ' ' + string(regs{j});
        xy = bin_swp(mus_pair, nbin);
        binned.(subs{ii}).(regs{j}).mus = xy;
        % retardance
        ylab = 'retardance (degrees)';
        tit = string(subs{ii}) + ' ' + string(regs{j});
        xy = bin_swp(ret_pair, nbin);
        binned.(subs{ii}).(regs{j}).ret = xy;
    end

    %%% Subject-level Combine the frontal and occipital
    % Print status to console
    fprintf('\nCombining Front + Occip for %s',subs{ii})
    % Initialize struct to ensure only combining individual subject
    heat_pair_sub = struct();
    % Retrieve the pairs for each subject
    heat_pair_sub.(subs{ii}) = heat_pairs.(subs{ii});
    
    % Combine frontal and occipital for this subject
    [~,~,log_mus,log_ret] = combine_subjects_regions(heat_pair_sub);
    
    % Bin the data
    log_mus = bin_swp(log_mus, nbin);
    log_ret = bin_swp(log_ret, nbin);

    % Add to struct
    binned.(subs{ii}).comb.mus = log_mus;
    binned.(subs{ii}).comb.ret = log_ret;
end

%% SUBJECT: plot "binned" subset

% Figure porperties
xlab = 'EPVS log(SWP)'; % x-axis label
mus_ylab = '\mu_s (cm^-^1)';
ret_ylab = 'retardance (degrees)';
xlims = [1,10];

% Set limits if removing offset
if rm_offset
    mus_ylim = [-5,5];
    ret_ylim = [-10,5];
    substr = '_median_removed_';
else
    mus_ylim = [5,15];
    ret_ylim = [16,35];
    substr = '';
end

% Iterate subjects
for ii = 1:length(subs)
    % Retrieve regions for this subject
    regs = fields(heat_pairs.(subs{ii}));
    % iterate over regions
    for j = 1:length(regs)
        % Retrieve the binned data
        mus_pair = binned.(subs{ii}).(regs{j}).mus;
        ret_pair = binned.(subs{ii}).(regs{j}).ret;
        % Set title string for both plots
        tit = string(subs{ii}) + ' ' + string(regs{j});
        % scattering coefficient
        fout = strcat(subs{ii},'_',regs{j},substr,'_mus_vs_log_swp');
        swp_scatterplot(mus_pair, xlab, mus_ylab, tit, xlims, mus_ylim,...
                        plt_dir,fout)
        % Retardance
        fout = strcat(subs{ii},'_',regs{j},substr,'_ret_vs_log_swp');
        swp_scatterplot(ret_pair, xlab, ret_ylab, tit, xlims, ret_ylim,...
                        plt_dir,fout)
    end
    
    %%% Plot the combined front + occipital
    % Retrieve the binned data
    mus_pair = binned.(subs{ii}).comb.mus;
    ret_pair = binned.(subs{ii}).comb.ret;
    % Set title string for both plots
    tit = string(subs{ii}) + ' combined';
    % scattering coefficient
    fout = strcat(subs{ii},substr,'_comb_mus_vs_log_swp');
    swp_scatterplot(mus_pair, xlab, mus_ylab, tit, xlims, mus_ylim,...
                    plt_dir,fout)
    % Retardance
    fout = strcat(subs{ii},substr,'_comb_ret_vs_log_swp');
    swp_scatterplot(ret_pair, xlab, ret_ylab, tit, xlims, ret_ylim,...
                    plt_dir,fout)
end

%%% SUBJECT: overlay all subjects
%%% Define Labels (legend) and Colors
labels = struct();
colors = struct();
% Frontal (no CAA 17)
labels.front = {'CAA26 (Control)','CAA6 (Control)',...
               'CAA 25 (Severe)','CAA 22 (Severe)'};
colors.front = {
    [0.3010 0.7450 0.9330];   % Cyan
    [0.0000 0.4470 0.7000];   % Blue
    [0.8500 0.3250 0.0980];   % Orange
    [0.6350 0.0780 0.1840];   % Red
    };   
% Occipital (with CAA 17)
labels.occip = {'CAA26 (Control)','CAA6 (Control)','CAA 17 (Moderate)',...
               'CAA 25 (Severe)','CAA 22 (Severe)'};
colors.occip = {
    [0.3010 0.7450 0.9330];   % Cyan
    [0.0000 0.4470 0.7000];   % Blue
    [0.9290 0.6940 0.1250];   % Yellow
    [0.8500 0.3250 0.0980];   % Orange
    [0.6350 0.0780 0.1840];   % Red
    };   

%%% scattering coefficient
% frontal
xy1 = binned.caa26.front.mus;
xy2 = binned.caa6.front.mus;
xy4 = binned.caa25.front.mus;
xy5 = binned.caa22.front.mus;
xy_cell = {xy1, xy2, xy4, xy5};
tit = 'Frontal: mus vs. log(SWP)';
fname = strcat('overlay_frontal_mus_vs_log_swp',substr);
swp_scatterplot_overlay(xy_cell, labels.front, colors.front, xlab,...
                        mus_ylab, tit, xlims, mus_ylim, plt_dir, fname)
% occipital
xy1 = binned.caa26.occip.mus;
xy2 = binned.caa6.occip.mus;
xy3 = binned.caa17.occip.mus;
xy4 = binned.caa25.occip.mus;
xy5 = binned.caa22.occip.mus;
xy_cell = {xy1, xy2, xy3, xy4, xy5};
tit = 'Occipital: mus vs. log(SWP)';
fname = strcat('overlay_occipital_mus_vs_log_swp',substr);
swp_scatterplot_overlay(xy_cell, labels.occip, colors.occip, xlab,...
                        mus_ylab, tit, xlims, mus_ylim, plt_dir, fname)
% Combined
xy1 = binned.caa26.comb.mus;
xy2 = binned.caa6.comb.mus;
xy3 = binned.caa17.comb.mus;
xy4 = binned.caa25.comb.mus;
xy5 = binned.caa22.comb.mus;
xy_cell = {xy1, xy2, xy3, xy4, xy5};
tit = 'Combined: mus vs. log(SWP)';
fname = strcat('overlay_combined_mus_vs_log_swp',substr);
swp_scatterplot_overlay(xy_cell, labels.occip, colors.occip, xlab,...
                        mus_ylab, tit, xlims, mus_ylim, plt_dir, fname)

%%% Retardance
% frontal
xy1 = binned.caa26.front.ret;
xy2 = binned.caa6.front.ret;
xy4 = binned.caa25.front.ret;
xy5 = binned.caa22.front.ret;
xy_cell = {xy1, xy2, xy4, xy5};
tit = 'Frontal: ret vs. log(SWP)';
fname = strcat('overlay_frontal_ret_vs_log_swp',substr);
swp_scatterplot_overlay(xy_cell, labels.front, colors.front, xlab,...
                        ret_ylab, tit, xlims, ret_ylim, plt_dir, fname)
% occipital
xy1 = binned.caa26.occip.ret;
xy2 = binned.caa6.occip.ret;
xy3 = binned.caa17.occip.ret;
xy4 = binned.caa25.occip.ret;
xy5 = binned.caa22.occip.ret;
xy_cell = {xy1, xy2, xy3, xy4, xy5};
tit = 'Occipital: ret vs. log(SWP)';
fname = strcat('overlay_occipital_ret_vs_log_swp',substr);
swp_scatterplot_overlay(xy_cell, labels.occip, colors.occip, xlab,...
                        ret_ylab, tit, xlims, ret_ylim, plt_dir, fname)
% Combined
xy1 = binned.caa26.comb.ret;
xy2 = binned.caa6.comb.ret;
xy3 = binned.caa17.comb.ret;
xy4 = binned.caa25.comb.ret;
xy5 = binned.caa22.comb.ret;
xy_cell = {xy1, xy2, xy3, xy4, xy5};
tit = 'Combined: ret vs. log(SWP)';
fname = strcat('overlay_combined_ret_vs_log_swp',substr);
swp_scatterplot_overlay(xy_cell, labels.occip, colors.occip, xlab,...
                        ret_ylab, tit, xlims, ret_ylim, plt_dir, fname)

%% Separate/combine heat_pairs across subject and region
fprintf('Separating heatmap pairs by subject and region\n')
% combine across occip + frontal
[~, ~, comb_log_mus, comb_log_ret] = combine_subjects_regions(heat_pairs);

% Combine subjects, split regions
region_data = combine_subjects(heat_pairs);
% Frontal pairs
front_log_mus = region_data.front.log10.mus;
front_log_ret = region_data.front.log10.ret;
% Occip pairs
occip_log_mus = region_data.occip.log10.mus;
occip_log_ret = region_data.occip.log10.ret;

%% ALL SUBJECTS - Bin and Plot optical property vs. log(SWP)
% x-axis label
xlab = 'EPVS log(SWP)';
% Number of bins in x-axis
nbin = 190;

%%% combined subjects and regions
fprintf('Binning for front + occip\n')
% scattering coefficient
xy = bin_swp(comb_log_mus, nbin);
binned.comb_mus = xy;
% retardance
xy = bin_swp(comb_log_ret, nbin);
binned.comb_ret = xy;

%%% Frontal
fprintf('Binning for front\n')
% Scattering
xy = bin_swp(front_log_mus, nbin);
binned.front_mus = xy;
% Retardance
xy = bin_swp(front_log_ret, nbin);
binned.front_ret = xy;

%%% Occipital
fprintf('Binning for occip\n')
% Scattering
xy = bin_swp(occip_log_mus, nbin);
binned.occip_mus = xy;
% Retardance
xy = bin_swp(occip_log_ret, nbin);
binned.occip_ret = xy;

%% Calculate the min/max values across all binned matrices
groups = {'comb_mus','comb_ret','front_mus','front_ret',...
          'occip_mus','occip_ret'};
nfield = length(groups);
xmin = 1;
xmax = 5;
% Set starting limits if removing offsets
if rm_offset
    mus_min = -1;
    mus_max = 1;
    ret_min = -1;
    ret_max = 1;
else
    mus_min = 13;
    mus_max = 13;
    ret_min = 30;
    ret_max = 30;
end

for ii = 1:nfield
    xy = binned.(groups{ii});
    if contains(groups{ii},'mus')
        mus_min = min([mus_min,min(xy(:,2))]);
        mus_max = max([mus_max,max(xy(:,2))]);
    else
        ret_min = min([ret_min,min(xy(:,2))]);
        ret_max = max([ret_max,max(xy(:,2))]);
    end
    % Set x-axis limits across both
    xmax = max([xmax, max(xy(:,1))]);
end

% Consolidate limits
xlims = [xmin,ceil(xmax)];
mus_lims = [ceil(mus_min), ceil(mus_max)];
ret_lims = [ceil(ret_min), ceil(ret_max)];

%%% Plot combined frontal + occipital
% scattering coefficient
ylab = '\mu_s (cm^-^1)';
tit = {'Frontal + Occipital', '\mu_s vs. EPVS log(SWP)'};
xy = binned.comb_mus;
fname = strcat(substr,'comb_mus_swp.png');
swp_scatterplot(xy, xlab, ylab, tit, xlims, mus_lims, plt_dir, fname);
% retardance
ylab = 'retardance (degrees)';
tit = {'Frontal + Occipital', 'Retardance vs. EPVS log(SWP)'};
xy = binned.comb_ret;
fname = strcat(substr,'comb_ret_swp.png');
swp_scatterplot(xy, xlab, ylab, tit, xlims, ret_lims, plt_dir, fname);

%%% Plot frontal
% scattering coefficient
ylab = '\mu_s (cm^-^1)';
tit = {'Frontal', '\mu_s vs. EPVS log(SWP)'};
xy = binned.front_mus;
fname = strcat(substr,'front_mus_swp.png');
swp_scatterplot(xy, xlab, ylab, tit, xlims, mus_lims, plt_dir, fname);
% retardance
ylab = 'retardance (degrees)';
tit = {'Frontal', 'Retardance vs. EPVS log(SWP)'};
xy = binned.front_ret;
fname = strcat(substr,'front_ret_swp.png');
swp_scatterplot(xy, xlab, ylab, tit, xlims, ret_lims, plt_dir, fname);

%%% Plot occipital
% scattering coefficient
ylab = '\mu_s (cm^-^1)';
tit = {'Occipital', '\mu_s vs. EPVS log(SWP)'};
xy = binned.occip_mus;
fname = strcat(substr,'occip_mus_swp.png');
swp_scatterplot(xy, xlab, ylab, tit, xlims, mus_lims, plt_dir, fname);
% retardance
ylab = 'retardance (degrees)';
tit = {'Occipital', 'Retardance vs. EPVS log(SWP)'};
xy = binned.occip_ret;
fname = strcat(substr,'occip_ret_swp.png');
swp_scatterplot(xy, xlab, ylab, tit, xlims, ret_lims, plt_dir, fname);

%% Plot median optical property vs. median SWP

%%% Create arrays for storing pairs
% frontal + occipital
logswp_mus_med = zeros(9:2);
logswp_ret_med = zeros(9:2);
% frontal
logswp_mus_med_front = zeros(5:2);
logswp_ret_med_front = zeros(5:2);
% occip
logswp_mus_med_occip = zeros(4:2);
logswp_ret_med_occip = zeros(4:2);

%%% Counters for iterations
tot_cnt = 1;
front_cnt = 1;
occip_cnt = 1;

% Iterate subject names
subs = fields(subjects);
for ii = 1:length(subs)
    sub = subs{ii};
    regions = fields(subjects.(sub));
    % Iterate regions
    for j=1:length(regions)
        reg = regions{j};
        fprintf('Median for %s %s\n',sub,reg)
        % Check if the SWP field exists
        if isfield(subjects.(sub).(reg), 'swp')
            %%% Median optical properties from parench_median
            mus = parench_median.(sub).(reg).med.mus;
            ret = parench_median.(sub).(reg).med.ret;
            % Calculate median log(swp)
            seg = subjects.(sub).(reg).seg;
            epvs = subjects.(sub).(reg).epvs;
            wm = subjects.(sub).(reg).mask_wm;
            swp_log = subjects.(sub).(reg).swp.log10;
            swp_log(seg) = NaN;
            swp_log(epvs) = NaN;
            swp_log = omit_outlier(swp_log(wm),10);
            swp_log = median(swp_log,'omitnan');
    
            %%% Add to pairs
            logswp_mus_med(tot_cnt,:) = [swp_log; mus];
            logswp_ret_med(tot_cnt,:) = [swp_log; ret];
            % Iterate counter
            tot_cnt = tot_cnt + 1;

            %%% Frontal
            if strcmp(reg,'front')
                logswp_mus_med_front(front_cnt,:) = [swp_log; mus];
                logswp_ret_med_front(front_cnt,:) = [swp_log; ret];
                front_cnt = front_cnt + 1;
            %%% Occipital
            else
                logswp_mus_med_occip(occip_cnt,:) = [swp_log; mus];
                logswp_ret_med_occip(occip_cnt,:) = [swp_log; ret];
                occip_cnt = occip_cnt + 1;
            end
        end
    end
end

%%% Combined front + occip
% mus
figure('Position', [956  -274   909   844]);
scatter(logswp_mus_med(:,1),logswp_mus_med(:,2),200,'b','filled');
title('Median \mu_s vs. Median log(swp) - Combined')
xlabel('SWP'); ylabel('\mu_s');set(gca,'fontsize',25)
fout = fullfile(plt_dir, 'median_mus_vs_median_log_swp_combined.png');
saveas(gcf,fout,'png'); pause(1); close;
% ret
figure('Position', [956  -274   909   844]);
scatter(logswp_ret_med(:,1),logswp_ret_med(:,2),200,'b','filled');
title('Median Ret vs. Median log(swp) - Combined')
xlabel('SWP'); ylabel('Ret');set(gca,'fontsize',25)
fout = fullfile(plt_dir, 'median_ret_vs_median_log_swp_combined.png');
saveas(gcf,fout,'png'); pause(1); close;

%%% Frontal
% mus
figure('Position', [956  -274   909   844]);
scatter(logswp_mus_med_front(:,1),logswp_mus_med_front(:,2),200,'b','filled');
title('Median \mu_s vs. Median log(swp) - Front')
xlabel('SWP'); ylabel('\mu_s');set(gca,'fontsize',25)
fout = fullfile(plt_dir, 'median_mus_vs_median_log_swp_front.png');
saveas(gcf,fout,'png'); pause(1); close;
% ret
figure('Position', [956  -274   909   844]);
scatter(logswp_ret_med_front(:,1),logswp_ret_med_front(:,2),200,'b','filled');
title('Median Ret vs. Median log(swp) - Front')
xlabel('SWP'); ylabel('Ret');set(gca,'fontsize',25)
fout = fullfile(plt_dir, 'median_ret_vs_median_log_swp_front.png');
saveas(gcf,fout,'png'); pause(1); close;

%%% Occipital
% mus
figure('Position', [956  -274   909   844]);
scatter(logswp_mus_med_occip(:,1),logswp_mus_med_occip(:,2),200,'b','filled');
title('Median \mu_s vs. Median log(swp) - Occipital')
xlabel('SWP'); ylabel('\mu_s');set(gca,'fontsize',25)
fout = fullfile(plt_dir, 'median_mus_vs_median_log_swp_occip.png');
saveas(gcf,fout,'png'); pause(1); close;
% ret
figure('Position', [956  -274   909   844]);
scatter(logswp_ret_med_occip(:,1),logswp_ret_med_occip(:,2),200,'b','filled');
title('Median Ret vs. Median log(swp) - Occipital')
xlabel('SWP'); ylabel('Ret');set(gca,'fontsize',25)
fout = fullfile(plt_dir, 'median_ret_vs_median_log_swp_occip.png');
saveas(gcf,fout,'png'); pause(1); close;

%% Export OP/SWP pairs to spreadsheets

%%% OV vs. SWP
% Retrieve current date for filename
dt = string(datetime('now','Format','d-MMM-y'));
% construct filename
field_names = fields(binned);
swp_filename = strcat('op_vs_swp_',substr,dt','.xlsx');
swp_filename = fullfile(swp_dir, swp_filename);
% Iterate over stains
for idx = 1:length(fields(binned))
    % retrieve field name
    structName = field_names{idx};
    
    % set the variable name for the spreadsheet
    if contains(structName, 'mus')
        var_name = 'mus';
    else
        var_name = 'ret';
    end

    % Retrieve matrix
    xy = binned.(structName);
    
    % Create Table for all measurements from stain
    T = array2table(xy, 'VariableNames',...
                   {'log(SWP)',var_name});
    writetable(T, swp_filename, 'Sheet', structName);
end

%%% Median OP vs. Median log(SWP)
dt = string(datetime('now','Format','d-MMM-y'));
median_swp_filename = strcat('med_op_vs_med_swp_',substr,dt','.xlsx');
median_swp_filename = fullfile(swp_dir, median_swp_filename);
% Write the mus sheet
xy = logswp_mus_med;
% Create Table for all measurements from stain
T = array2table(xy, 'VariableNames',...
               {'log(SWP)','mus'});
writetable(T, median_swp_filename, 'Sheet', 'mus');
% Write the retardance sheet
xy = logswp_ret_med;
% Create Table for all measurements from stain
T = array2table(xy, 'VariableNames',...
               {'log(SWP)','ret'});
writetable(T, median_swp_filename, 'Sheet', 'ret');

%% Combine across all sujects/regions

function [raw_mus,raw_ret,log_mus,log_ret]=...
    combine_subjects_regions(heat_pairs)
% Combine mus_pair across all subjects and all regions
% Output: all_mus_pairs is a 2xN matrix

% ----------- Pass 1: Count total columns ------------
raw_samples = 0;
log_samples = 0;
subjects = fieldnames(heat_pairs);

for i = 1:numel(subjects)
    subj = subjects{i};
    regions = fieldnames(heat_pairs.(subj));
    for j = 1:numel(regions)
        region = regions{j};
        % Raw
        mus = heat_pairs.(subj).(region).raw.mus_pair;
        raw_samples = raw_samples + size(mus, 1);
        % log
        mus = heat_pairs.(subj).(region).log10.mus_pair;
        log_samples = log_samples + size(mus, 1);
    end
end

% ----------- Preallocate ----------------------------
% Raw SWP
raw_mus = zeros(raw_samples,2);
raw_ret = zeros(raw_samples,2);
% log(swp)
log_mus = zeros(log_samples,2);
log_ret = zeros(log_samples,2);

% ----------- Pass 2: Fill data ----------------------
% heat_pairs.(sub).(reg).log10.mus_pair = mus_pair;
raw_idx = 1;
log_idx = 1;
for i = 1:numel(subjects)
    subj = subjects{i};
    regions = fieldnames(heat_pairs.(subj));
    for j = 1:numel(regions)
        region = regions{j};
        %%% Raw pairs
        mus = heat_pairs.(subj).(region).raw.mus_pair;
        ret = heat_pairs.(subj).(region).raw.ret_pair;
        n = size(mus, 1);
        raw_mus(raw_idx:raw_idx+n-1,:) = mus;
        raw_ret(raw_idx:raw_idx+n-1,:) = ret;
        % Iterate column index
        raw_idx = raw_idx + n;
        
        %%% log pairs
        mus = heat_pairs.(subj).(region).log10.mus_pair;
        ret = heat_pairs.(subj).(region).log10.ret_pair;
        n = size(mus, 1);
        log_mus(log_idx:log_idx+n-1,:) = mus;
        log_ret(log_idx:log_idx+n-1,:) = ret;
        % Iterate column index
        log_idx = log_idx + n;
    end
end
end


%% Combine across all sujects. Separate by regions
function region_data = combine_subjects(heat_pairs)
% Combine mus_pair and ret_pair across all subjects, separated by region
% Output: region_data.(region).mus and region_data.(region).ret are 2xN matrices

region_data = struct();
dtype = {'raw','log10'};

subjects = fieldnames(heat_pairs);
for i = 1:numel(subjects)
    subj = subjects{i};
    regions = fieldnames(heat_pairs.(subj));
    for j = 1:numel(regions)
        for k = 1:length(dtype)
            region = regions{j};
            % Extract mus_pair and ret_pair
            mus = heat_pairs.(subj).(region).(dtype{k}).mus_pair;
            ret = heat_pairs.(subj).(region).(dtype{k}).ret_pair;
            % Initialize or append
            if isfield(region_data, (region))
                if isfield(region_data.(region),dtype{k})
                    region_data.(region).(dtype{k}).mus =...
                        [region_data.(region).(dtype{k}).mus; mus];
                    region_data.(region).(dtype{k}).ret =...
                        [region_data.(region).(dtype{k}).ret; ret];
                else
                    region_data.(region).(dtype{k}).mus = mus;
                    region_data.(region).(dtype{k}).ret = ret;
                end
            else
                region_data.(region).(dtype{k}).mus = mus;
                region_data.(region).(dtype{k}).ret = ret;
            end
    
    %         %%% Iterate over log10
    %         % Extract mus_pair and ret_pair
    %         mus = heat_pairs.(subj).(region).log10.mus_pair;
    %         ret = heat_pairs.(subj).(region).log10.ret_pair;
    %         % Initialize or append
    %         if isfield(region_data, (region))
    %             if isfield(region_data.(region),'log10')
    %                 region_data.(region).log10.mus =...
    %                     [region_data.(region).log10.mus; mus];
    %                 region_data.(region).log10.ret =...
    %                     [region_data.(region).log10.ret; ret];
    %             else
    %                 region_data.(region).log10.mus = mus;
    %                 region_data.(region).log10.ret = ret;
    %             end
    %         end
        end
    end
end
end

%% Scatter plot of combined vectors
function scatter_op(pair, xlab, ylab, tit, th)
% Scatter plot of the optical properts vs. EPVS density
% INPUTS:
%   pair (Nx2 matrix): [EPVS density, optical property]
%   xlab (string): x-axis label
%   ylab (string): y-axis label
%   tit (string): title of figure
%   th (double): minimum threshold for EPVS density

% Remove pairs where EPVS density < threshold (th)
keep_idx = pair(:,1) >= th;
pairf = pair(keep_idx,:);

% scatter plot
figure; scatter(pairf(:,1),pairf(:,2),20,'k','filled');
% Plot features
xlabel(xlab)
ylabel(ylab)
title(tit);
set(gca,'fontsize',20)
end

%% Errorbar plot of windowed averages
function sliding_errorbar_plot(x,y,err, xlab, ylab, tit)
% Scatter plot of the optical properts vs. EPVS density
% INPUTS:
%   pair (Nx2 matrix): [EPVS density, optical property]
%   err (vector): standard error at each window
%   xlab (string): x-axis label
%   ylab (string): y-axis label
%   tit (string): title of figure

% scatter plot
figure
errorbar(x,y,err,'o-')
% Plot features
xlabel(xlab)
ylabel(ylab)
title(tit);
set(gca,'fontsize',20)

end

%% Apply minimum threshold
function [pairf] = threshold_epvs(pair, th)

% Remove pairs where EPVS density < threshold (th)
keep_idx = pair(:,1) >= th;
pairf = pair(keep_idx,:);

end