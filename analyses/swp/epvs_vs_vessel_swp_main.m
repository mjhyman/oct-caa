%% Analyze the EPVS SWP and Respective Vessel SWP heatmap
% Import the EPVS SWP heatmaps
% Import the .MAT of the optical properties
% Bin the data points along x-axis
% Column 1 = x-axis = EPVS SWP
% Column 2 = y-axis = Vessel SWP

%% Initialization
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

%%% Directories
% Directory containing seg, mus, ret, mask, epvs structs
mat_dir = '/projectnb/npbssmic/ns/CAA/';
% Heatmap directory
swp_dir = '/projectnb/npbssmic/ns/CAA/swp';
% Scatterplot directory
plt_dir = '/projectnb/npbssmic/ns/CAA/figures/epvs_swp_vs_ves_swp/';

%%% Common filename for importing
epvs_fbase = '_swp_voxelwise_radius_500_exp_2_interpolated_heatmap.mat';
ves_fbase = '_swp_voxelwise_ves_radius_500_exp_2_interpolated_heatmap.mat';

%%% Flags for importing data
% Flag for SWP structs
flag_load_swp_structs = true;
% Flag for importing .MAT structs
flag_load_caa_structs = true;

%%% Number of bins for the x-axis along SWP
nbin = 100;

%%% Plot Labels
xlab = 'EPVS SWP';
ylab = 'Vessel SWP';

%%% Initialize struct for subject/region
subjects = struct();
subjects.caa6.front = [];
subjects.caa6.occip = [];
subjects.caa17.occip = [];
subjects.caa22.front = [];
subjects.caa22.occip = [];
subjects.caa25.front = [];
subjects.caa25.occip = [];
subjects.caa26.front = [];
subjects.caa26.occip = [];

%% Load SWP, Vessel SWP

if flag_load_swp_structs
    % Struct for storing SWP
    swp_struct = struct();
    % Iterate subjects
    subs = fields(subjects);
    for ii = 1:numel(fields(subjects))
        sub = subs{ii};
        regs = fields(subjects.(sub));
        for j = 1:numel(regs)
            % retrieve subject and region
            reg = regs{j};
            fprintf('Importing SWP for %s %s\n',sub,reg)
            % EPVS
            tmp = load(fullfile(swp_dir, sub, reg, strcat(sub, '_', reg, epvs_fbase)));
            swp_struct.epvs.(sub).(reg) = tmp.interpolated_volume;
            % Vessels
            tmp = load(fullfile(swp_dir, sub, reg, strcat(sub, '_', reg, ves_fbase)));
            swp_struct.ves.(sub).(reg) = tmp.interpolated_volume;
        end
    end
end

%% Load CAA structs (OCT data)
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
    caa.caa6 = caa6.caa6;
    caa.caa17 = caa17.caa17;
    caa.caa22 = caa22.caa22;
    caa.caa25 = caa25.caa25;
    caa.caa26 = caa26.caa26;
end

%% Create 2D arrays of EPVS SWP vs. Vessel SWP
% 2xN Matrix for each optical property:
%   - column 1 = EPVS SWP from heatmap
%   - column 2 = Vessel SWP from heatmap
%   - row = pairwise observations from same tissue volume

fprintf('Creating 2D arrays of EPVS SWP vs. vessel SWP\n')

% struct for storing pairs
heat_pairs = struct();

% Iterate subjects
for ii = 1:numel(fields(subjects))
    sub = subs{ii};
    regs = fields(subjects.(sub));
    for j = 1:numel(regs)
        reg = regs{j};
        %%% Retrieve SWP for EPVS and Vessel
        swp_epvs = swp_struct.epvs.(sub).(reg);
        swp_ves = swp_struct.ves.(sub).(reg);

        %%% Remove vessels and EPVS from the EPVS and Vessel SWP 
        % Retrieve white matter mask
        mask_wm = caa.(sub).(reg).mask_wm;
        % Remove vesselsa and EPVS from white matter mask
        ves = caa.(sub).(reg).seg;
        epvs = caa.(sub).(reg).epvs;
        mask_wm(ves) = 0;
        mask_wm(epvs) = 0;

        %%% Apply WM mask to SWP for both EPVS and vessel
        swp_ves = swp_ves(mask_wm);
        swp_epvs = swp_epvs(mask_wm);

        %%% Create pairs for raw SWP
        swp_pair = [swp_epvs(:), swp_ves(:)];

        % add pairs to heatmap struct
        heat_pairs.(sub).(reg) = swp_pair;
    end
end
fprintf('Finished making 2D arrays of EPVS density vs. optical prop\n')

%%% Combine subjects, split regions
region_data = combine_subjects(heat_pairs);

%%% plotting properties
% struct for storing binned vectors
binned = struct();

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
        pair = heat_pairs.(subs{ii}).(regs{j});
        % scattering coefficient
        [xy,se] = bin_swp(pair, nbin);
        binned.(subs{ii}).(regs{j}).pair = xy;
        binned.(subs{ii}).(regs{j}).se = se;
    end

    %%% Subject-level Combine the frontal and occipital
    % Print status to console
    fprintf('\nCombining Front + Occip for %s',subs{ii})
    % Initialize struct to ensure only combining individual subject
    heat_pair_sub = struct();
    % Retrieve the pairs for each subject
    heat_pair_sub.(subs{ii}) = heat_pairs.(subs{ii});
    
    % Combine frontal and occipital for this subject
    [comb] = combine_subjects_regions(heat_pair_sub);
    
    % Bin the combined regions
    [comb,comb_se] = bin_swp(comb, nbin);

    % Add to struct
    binned.(subs{ii}).comb      = comb;
    binned.(subs{ii}).comb_se   = comb_se;
end

% Combine CAA22 or CAA25 for severe
heat_pairs.sev.front = [heat_pairs.caa22.front; heat_pairs.caa25.front];
heat_pairs.sev.occip = [heat_pairs.caa22.occip; heat_pairs.caa25.occip];

%%% Separate/combine heat_pairs across subject and region
fprintf('Separating heatmap pairs by subject and region\n')
% combine occip + frontal across all subjects
[comb] = combine_subjects_regions(heat_pairs);

%%% combined subjects and regions
fprintf('Binning for front + occip\n')
% scattering coefficient
[xy,se] = bin_swp(comb, nbin);
binned.comb = xy;
binned.comb_se = se;
fprintf('Finished binning for front + occip\n')

%% 2D Joint Density Heatmap
% TODO: finish this

% Settings
groups = {'caa26','caa6','caa17','sev'};
severities = {'control', 'mild', 'moderate', 'severe'};
n_bins     = 100;

%%% Find histogram edges for epvs-swp and ves-swp (front, occip separate)
% lims.ves.front = [0,5];
% lims.ves.occip = [0,5];
% lims.epv.front = [0,5];
% lims.epv.occip = [0,5];
% % First find min/max for each region and ves-swp and epvs-swp
% for ii = 1:length(groups)
%     % Find number of regions
%     regs = fields(heat_pairs.(groups{ii}));
%     for j = 1:length(regs)
%         % Retrieve [epv, ves] pairs and limits
%         xy = heat_pairs.(groups{ii}).(regs{j});
%         epv_lims = lims.epv.(regs{j});
%         ves_lims = lims.ves.(regs{j});
%         % Compare limits to vector
%         lims.epv.(regs{j}) = [0,ceil(max(epv_lims(2),max(xy(:,1))))];
%         lims.ves.(regs{j}) = [0,ceil(max(ves_lims(2),max(xy(:,2))))];
%         pause(1)
%     end
% end
% % Calculate edges
% type = {'ves','epv'}; regs = {'front','occip'};
% % iterate ves/epvs
% for ii = 1:2
%     % iterate front/occip
%     for j = 1:2
%         edges.(type{ii}).(regs{j}) = linspace(min(lims.(type{ii}).(regs{j})),...
%                                               max(lims.(type{ii}).(regs{j})),...
%                                               n_bins+1);
%     end
% end

%%% Edges based on visual inspection
% control = [0,50] for both ves,epvs for front/occip
% mild = [0,50] for both ves,epvs for front/occip
% moderate
%   epvs = [0,150] + ves = [0,100]
% severe
%   front: epvs = [0,500] + ves = [0,100]
%   occip: epvs = [0,100] + ves = [0,100]
low_lim = 25;
edges.epv.caa26.front = linspace(0,low_lim,n_bins+1);
edges.epv.caa26.occip = linspace(0,low_lim,n_bins+1);
edges.ves.caa26.front = linspace(0,low_lim,n_bins+1);
edges.ves.caa26.occip = linspace(0,low_lim,n_bins+1);
edges.epv.caa6 = edges.epv.caa26;
edges.ves.caa6 = edges.ves.caa26;
edges.epv.caa17.occip = linspace(0,100,n_bins+1);
edges.ves.caa17.occip = linspace(0,50,n_bins+1);
edges.epv.sev.front = linspace(0,400,n_bins+1);
edges.epv.sev.occip = linspace(0,50,n_bins+1);
edges.ves.sev.front = linspace(0,50,n_bins+1);
edges.ves.sev.occip = linspace(0,50,n_bins+1);

%%% Iterate severities
for ii = 1:length(groups)
    % Find number of regions
    regs = fields(heat_pairs.(groups{ii}));
    for j = 1:length(regs)
        % Retrieve [epv, ves] pairs and limits
        epv_ves = heat_pairs.(groups{ii}).(regs{j});
        epv_edge = edges.epv.(groups{ii}).(regs{j});
        ves_edge = edges.ves.(groups{ii}).(regs{j});
        
        % Create 2D histogram (X = epv-swp, Y = ves-swp)
        histo = histcounts2(epv_ves(:,1), epv_ves(:,2),epv_edge,ves_edge,...
                            'Normalization','percentage');
        
        % Generate figure of histogram
        figure;
        imagesc(epv_edge, ves_edge, histo'); set(gca,'YDir','normal');
        colorbar; clim([0,0.15])
        xlabel('EPVS SWP'); ylabel('Vessel SWP');
        title(sprintf('%s %s',severities{ii}, regs{j}));
    end
end
% sgtitle('Joint distribution of SWP and EPVS SWP by severity');

%% ALL SUBJECTS - Bin Frontal and Occipital separately
%%% Frontal
fprintf('Binning for front\n')
% Scattering
[xy,se] = bin_swp(region_data.front, nbin);
binned.front = xy;
binned.front_se = se;

%%% Occipital
fprintf('Binning for occip\n')
% Scattering
[xy,se] = bin_swp(region_data.occip, nbin);
binned.occip = xy;
binned.occip_se = se;

%%% Calculate the min/max values across all binned matrices
% x-axis = EPVS SWP
% y-axis = Vessel SWP
% Initialize
groups = {'comb','front','occip'};
nfield = length(groups);
xmin = 0;
xmax = 5;
% Set starting limits if removing offsets
ymin = 0;
ymax = 2;
for ii = 1:nfield
    xy = binned.(groups{ii});
    ymin = min([ymin,min(xy(:,2))]);
    ymax = max([ymax,max(xy(:,2))]);
    % Set x-axis limits across both
    xmax = max([xmax, max(xy(:,1))]);
end
% Consolidate limits
xlims = [xmin,ceil(xmax)];
ylims = [ceil(ymin), ceil(ymax)];

%% Plot combined, frontal, occipital
%%% Initialize the labels
ylab = 'Vessel SWP';
xlab = 'EPVS SWP';

%%% Combined
tit = {'Frontal + Occipital', 'Vessel SWP vs. EPVS SWP'};
xy = binned.comb;
se = binned.comb_se;
fname = strcat('comb_ves_swp_vs_epvs_swp');
swp_scatterplot_subject(xy, xlab, ylab, tit, plt_dir, fname);
fname = strcat('comb_ves_swp_vs_epvs_swp_errorbar');
swp_scatterplot_errorbars(xy, se, xlab, ylab, tit, xlims, ylims, plt_dir, fname);

%%% Plot frontal
% scattering coefficient
tit = {'Frontal', 'Vessel SWP vs. EPVS SWP'};
xy = binned.front;
se = binned.front_se;
fname = strcat('front_ves_swp_vs_epvs_swp');
swp_scatterplot_subject(xy, xlab, ylab, tit, plt_dir, fname);
fname = strcat('front_ves_swp_vs_epvs_swp_errorbar');
swp_scatterplot_errorbars(xy, se, xlab, ylab, tit, xlims, ylims, plt_dir, fname);

%%% Plot occipital
% scattering coefficient
tit = {'Occipital', 'Vessel SWP vs. EPVS SWP'};
xy = binned.occip;
se = binned.occip_se;
fname = strcat('occip_ves_swp_vs_epvs_swp');
swp_scatterplot_subject(xy, xlab, ylab, tit, plt_dir, fname);
fname = strcat('occip_ves_swp_vs_epvs_swp_errorbar');
swp_scatterplot_errorbars(xy, se, xlab, ylab, tit, xlims, ylims, plt_dir, fname);


%% Export EPVS-SWP/Vessel-SWP pairs to spreadsheets
%{
%%% OV vs. SWP
% Retrieve current date for filename
dt = string(datetime('now','Format','d-MMM-y'));
% construct filename
field_names = {'front_mus','front_ret','occip_mus','occip_ret'};
swp_filename = strcat('op_vs_swp_',substr,dt','.xlsx');
swp_filename = fullfile(swp_dir, swp_filename);
% Iterate over stains
for idx = 1:length(field_names)
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

    %%% Measure Spearman's rho correlation
    [rho,p] = corr(xy(:,1), xy(:,2), 'Type','Spearman');
    % convert to table
    rho_t = array2table([rho,p],'VariableNames',{'rho','p-value'});
    % write table
    sheet_name = strcat(structName,'_spearman');
    writetable(rho_t, swp_filename, 'Sheet', sheet_name);
    
    %%% Create Table for all measurements from stain
    T = array2table(xy, 'VariableNames',...
                   {'log(SWP)',var_name});
    writetable(T, swp_filename, 'Sheet', structName);
end
%}

%% STAGING: plot each stage / region

% Map subject to severity
binned.ctl = binned.caa26;
binned.mld = binned.caa6;
binned.mod = binned.caa17;

%%% Severe: combine CAA22 and CAA25 to 
% Frontal
x1 = binned.caa22.front.pair(:,1);
x2 = binned.caa25.front.pair(:,1);
x = mean([x1,x2],2);
y1 = binned.caa22.front.pair(:,2);
y2 = binned.caa25.front.pair(:,2);
y = mean([y1,y2],2);
y1_se = binned.caa22.front.se;
y2_se = binned.caa25.front.se;
y_se = mean([y1_se,y2_se],2);
% Frontal - place into struct
binned.sev.front.pair = [x,y];
binned.sev.front.se = y_se;
% Occipital
x1 = binned.caa22.occip.pair(:,1);
x2 = binned.caa25.occip.pair(:,1);
x = mean([x1,x2],2);
y1 = binned.caa22.occip.pair(:,2);
y2 = binned.caa25.occip.pair(:,2);
y = mean([y1,y2],2);
y1_se = binned.caa22.occip.se;
y2_se = binned.caa25.occip.se;
y_se = mean([y1_se,y2_se],2);
% Frontal - place into struct
binned.sev.occip.pair = [x,y];
binned.sev.occip.se = y_se;

%%% Set the min/max values across all binned matrices
% x-axis = EPVS SWP
% y-axis = Vessel SWP
xlims.front = [0,450];
xlims.occip = [0,300];
ylims.front = [0,25];
ylims.occip = [0,42];

%%% Plot each severity
sevs = {'ctl','mld','mod','sev'};
% Iterate severities
for ii = 1:length(sevs)
    % Check which regions are present
    if isfield(binned.(sevs{ii}),'front')
        regs = {'front','occip'};
    else
        regs = {'occip'};
    end
    % Iterate regions
    for j = 1:length(regs)
        % Retrieve xy pair
        xy = binned.(sevs{ii}).(regs{j}).pair;
        % set limits
        xl = xlims.(regs{j});
        yl = ylims.(regs{j});
        % Plot
        tstr = strcat(sevs{ii},' ',regs{j});
        fname = strcat(sevs{ii},'_',regs{j},'__ves_swp_vs_epvs_swp');
        swp_scatterplot(xy,xlab,ylab,tstr,xl,yl,plt_dir,fname);
    end
end


%% SUBJECT: plot "binned" subset separately

% Iterate subjects
for ii = 1:length(subs)
    % Retrieve regions for this subject
    regs = fields(binned.(subs{ii}));
    % If occip is part, then add to regs
    if ismember('front',regs)
        regs = {'front','occip'};
    else
        regs = {'occip'};
    end
    % iterate over regions
    for j = 1:length(regs)
        % Retrieve the binned data
        pair = binned.(subs{ii}).(regs{j}).pair;
        % Set title string for both plots
        tit = string(subs{ii}) + ' ' + string(regs{j});
        % scattering coefficient
        fout = strcat(subs{ii},'_',regs{j},'_ves_swp_vs_epvs_swp');
        swp_scatterplot_subject(pair, xlab, ylab, tit,plt_dir,fout)
    end
    
    %%% Plot the combined front + occipital
    % Retrieve the binned data
    pair = binned.(subs{ii}).comb;
    % Set title string for both plots
    tit = string(subs{ii}) + ' combined';
    % scattering coefficient
    fout = strcat(subs{ii},'_comb_ves_swp_vs_epvs_swp');
    swp_scatterplot_subject(pair, xlab, ylab, tit,plt_dir,fout)
end

%% SUBJECT: overlay all subjects
% Define Labels (legend) and Colors
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
xy1 = binned.caa26.front.pair;
xy2 = binned.caa6.front.pair;
xy4 = binned.caa25.front.pair;
xy5 = binned.caa22.front.pair;
xy_cell = {xy1, xy2, xy4, xy5};
tit = 'Frontal: Vessel SWP vs. EPVS SWP';
fname = strcat('overlay_frontal_ves_swp_vs_epvs_swp');
swp_scatterplot_overlay(xy_cell, labels.front, colors.front, xlab,...
                        ylab, tit, xlims, ylims, plt_dir, fname)
% occipital
xy1 = binned.caa26.occip.pair;
xy2 = binned.caa6.occip.pair;
xy3 = binned.caa17.occip.pair;
xy4 = binned.caa25.occip.pair;
xy5 = binned.caa22.occip.pair;
xy_cell = {xy1, xy2, xy3, xy4, xy5};
tit = 'Occipital: Vessel SWP vs. EPVS SWP';
fname = strcat('overlay_occipital_ves_swp_vs_epvs_swp');
swp_scatterplot_overlay(xy_cell, labels.occip, colors.occip, xlab,...
                        ylab, tit, xlims, ylims, plt_dir, fname)
% Combined
xy1 = binned.caa26.comb;
xy2 = binned.caa6.comb;
xy3 = binned.caa17.comb;
xy4 = binned.caa25.comb;
xy5 = binned.caa22.comb;
xy_cell = {xy1, xy2, xy3, xy4, xy5};
tit = 'Combined: Vessel SWP vs. EPVS SWP';
fname = strcat('overlay_combined_ves_swp_vs_epvs_swp');
swp_scatterplot_overlay(xy_cell, labels.occip, colors.occip, xlab,...
                        ylab, tit, xlims, ylims, plt_dir, fname)

%% Combine across all sujects/regions
function [comb] = combine_subjects_regions(heat_pairs)
% Combine mus_pair across all subjects and all regions
% Output: all_mus_pairs is a 2xN matrix

% ----------- Pass 1: Count total columns ------------
raw_samples = 0;
subjects = fieldnames(heat_pairs);

for i = 1:numel(subjects)
    subj = subjects{i};
    regions = fieldnames(heat_pairs.(subj));
    for j = 1:numel(regions)
        region = regions{j};
        pair = heat_pairs.(subj).(region);
        raw_samples = raw_samples + size(pair, 1);
    end
end

% ----------- Preallocate ----------------------------
% Raw SWP
comb = zeros(raw_samples,2);

% ----------- Pass 2: Fill data ----------------------
% heat_pairs.(sub).(reg).log10.mus_pair = mus_pair;
idx = 1;
for i = 1:numel(subjects)
    subj = subjects{i};
    regions = fieldnames(heat_pairs.(subj));
    for j = 1:numel(regions)
        region = regions{j};
        %%% Raw pairs
        pair = heat_pairs.(subj).(region);
        n = size(pair, 1);
        comb(idx:idx+n-1,:) = pair;
        % Iterate column index
        idx = idx + n;
    end
end
end


%% Combine across all subjects. Separate by regions
function region_data = combine_subjects(heat_pairs)
% Combine heatpars across all subjects, separated by region
% Output: region_data.(region) = 2xN matrix

    region_data = struct();
    subjects = fieldnames(heat_pairs);
    
    % First pass: Determine size for preallocation
    total_per_region = containers.Map('KeyType', 'char', 'ValueType', 'int32'); 
    for i = 1:numel(subjects)
        subj = subjects{i};
        regions = fieldnames(heat_pairs.(subj));
        for j = 1:numel(regions)
            region = regions{j};
            pair = heat_pairs.(subj).(region);
            if isKey(total_per_region, region)
                total_per_region(region) = total_per_region(region) + size(pair, 1);
            else
                total_per_region(region) = size(pair, 1);
            end
        end
    end
    
    % Initialize region_data with preallocated matrices
    for region = keys(total_per_region)
        region_name = region{1};
        n_pairs = total_per_region(region_name);
        % Assuming pairs have the same number of columns, use a placeholder pair to get column size
        pair_col_size = size(heat_pairs.(subjects{1}).(region_name), 2); 
        region_data.(region_name) = zeros(n_pairs, pair_col_size); 
    end

    % Second pass: Fill the preallocated matrices
    for i = 1:numel(subjects)
        subj = subjects{i};
        regions = fieldnames(heat_pairs.(subj));
        for j = 1:numel(regions)
            region = regions{j};
            pair = heat_pairs.(subj).(region);
            % Find the current starting index for filling
            if isfield(region_data, region)
                start_idx = size(region_data.(region), 1) - total_per_region(region) + 1;
            else
                start_idx = 1;
            end
            end_idx = start_idx + size(pair, 1) - 1;
            region_data.(region)(start_idx:end_idx, :) = pair;
        end
    end
end
