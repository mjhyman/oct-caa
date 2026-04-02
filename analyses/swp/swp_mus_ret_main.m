%% Analyze the EPVS SWP heatmap
% Import the EPVS and Vessel SWP heatmaps (without log transform)
% Import the .MAT of the optical properties
% Bin the data points along x-axis and y-axis
%   x-axis = EPVS SWP
%   y-axis = Optical Property
%{
TODO:
- 
%}


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
% Figure output directory
fig_dir = '/projectnb/npbssmic/ns/CAA/figures/swp_gam_gmm/';

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

%% Load EPVS SWP & Vessel SWP
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
    caa.caa6 = caa6.caa6;
    caa.caa17 = caa17.caa17;
    caa.caa22 = caa22.caa22;
    caa.caa25 = caa25.caa25;
    caa.caa26 = caa26.caa26;
end

%% Import medians for calculating percentage difference
% This is faster than calculating medians each iteration
load(fullfile(mat_dir, "median_white_matter_values_14-Nov-2025.mat"));

%% Create 4xN matrix (EPVS-SWP, vessel-SWP, mus, ret)
% Each parnechyma voxel has a value for EPVS-swp, vessel-SWP, mus, ret.
% They will be organized according to the following Nx4 Matrix:
%   - column 1 = mus
%   - column 2 = retardance
%   - column 3 = Vessel SWP
%   - column 4 = EPVS SWP

% struct for storing pairs
swp_op = struct();

% Nest for implementation
if ~exist('caa','var')
    caa.caa6 = caa6.caa6;
    caa.caa17 = caa17.caa17;
    caa.caa22 = caa22.caa22;
    caa.caa25 = caa25.caa25;
    caa.caa26 = caa26.caa26;
end

% Iterate subjects
fprintf('Creating 4xN matrices\n')
subs = fields(subjects);
for ii = 1:numel(fields(subjects))
    sub = subs{ii};
    regs = fields(subjects.(sub));
    for j = 1:numel(regs)
        %%% Retrieve region and optical properties        
        reg = regs{j};
        fprintf('\tSubject %s region %s\n',sub,reg)
        mus = caa.(sub).(reg).mus;
        ret = caa.(sub).(reg).ret_full;

        %%% Percentage difference
        % retrieve median
        mus_med = parench_median.(sub).(reg).med.mus;
        ret_med = parench_median.(sub).(reg).med.ret;
        % Take percentage difference from median
        mus_pdif = (mus - mus_med) ./ mus_med;
        ret_pdif = (ret - ret_med) ./ ret_med;

        %%% Retrieve SWP for EPVS and Vessel
        swp_epvs = swp_struct.epvs.(sub).(reg);
        swp_ves = swp_struct.ves.(sub).(reg);

        %%% Exclude vessels and EPVS from white matter mask (mask_wm)
        % Retrieve white matter mask
        mask_wm = caa.(sub).(reg).mask_wm;
        % Retrieve vessels and EPVS
        ves = caa.(sub).(reg).seg;
        epvs = caa.(sub).(reg).epvs;
        % Remove vessels and EPVS from white matter logical mask
        mask_wm(ves) = 0;
        mask_wm(epvs) = 0;
        
        %%% Apply WM mask
        swp_epvs = swp_epvs(mask_wm);
        swp_ves = swp_ves(mask_wm);
        mus = mus(mask_wm);
        ret = ret(mask_wm);
        mus_pdif = mus_pdif(mask_wm);
        ret_pdif = ret_pdif(mask_wm);

        %%% Create matrix
        swp = [mus(:), ret(:), swp_ves(:), swp_epvs(:)];
        swp_pdif = [mus_pdif(:), ret_pdif(:), swp_ves(:), swp_epvs(:)];

        % add pairs to heatmap struct
        swp_op.raw.(sub).(reg) = swp;
        swp_op.pdif.(sub).(reg) = swp_pdif;
    end
end
fprintf('Finished making 2D arrays of EPVS density vs. optical prop\n')

%% Consolidate regions across all subjects
fprintf('Combining regions\n')
swp_op_region_raw = combine_subjects(swp_op.raw);
swp_op_region_pdif = combine_subjects(swp_op.pdif);
fprintf('Finished combining regions\n')

%%% Convert to tables
vnames = {'scattering','retardance','SWP','EPVS_SWP'};
% Raw
T_front_raw = array2table(swp_op_region_raw.front,'VariableNames',vnames);
T_occip_raw = array2table(swp_op_region_raw.occip,'VariableNames',vnames);
% Percent Diff
T_front_pdif = array2table(swp_op_region_pdif.front,'VariableNames',vnames);
T_occip_pdif = array2table(swp_op_region_pdif.occip,'VariableNames',vnames);
% Add to struct
T.raw.front = T_front_raw;
T.raw.occip = T_occip_raw;
T.pdif.front = T_front_pdif;
T.pdif.occip = T_occip_pdif;

%% Create table for each subject/region
% Iterate subjects
fprintf('Creating 4xN matrices\n')
subs = fields(subjects);
for ii = 1:numel(fields(subjects))
    sub = subs{ii};
    regs = fields(subjects.(sub));
    for j = 1:numel(regs)
        %%% Retrieve region and optical properties        
        reg = regs{j};
        fprintf('\tSubject %s region %s\n',sub,reg)
        % Raw data
        tmp = array2table(swp_op.raw.(sub).(reg),'VariableNames',vnames);
        T.raw.(sub).(reg) = tmp;
        % PDiff data
        tmp = array2table(swp_op.pdif.(sub).(reg),'VariableNames',vnames);
        T.pdif.(sub).(reg) = tmp;
    end
end

%% Generalized Additive Model (GAM) - By Region (combined subjects)
%{
% Initialize GAM output struct
gam = struct();

%%% Sampling and Learning
% Set maximum number of samples
n_max = 1e6;
% Set learning parameters (avoid over fitting)
ntrees = 100;
max_splits = 4;
learn_rate = 0.05;
% Set learning parameters (sharper contrast near EPVS)
% ntrees = 200;
% max_splits = 10;
% learn_rate = 0.1;
%%% Plotting confidence intervals settings
nbootstrap = 100;
ci_alpha = 0.05;

%%% Iterate over raw/pdif and front/occip
dtype = {'raw','pdif'};
regs = {'front','occip'};
% Iterate raw/pdif
for ii = 1:2
    % Iterate front/occip
    for j = 1:2
        % Select data type and region
        d = dtype{ii};
        reg = regs{j};
        % Create title substring from data type and region
        tstr = string(sprintf('%s %s',string(d),string(reg)));
        fprintf('Running GAM for %s %s\n',d,reg)
        gam.(d).(reg) = fit_swp_psoct_gam(T.(d).(reg),'NumTrees',ntrees,...
                        'MaxSplits',max_splits,'LearnRate',learn_rate,...
                        'MaxSample',n_max,'NBootstrap',nbootstrap,...
                        'CIAlpha',ci_alpha,'TitleStr',tstr,...
                        'dirout',fig_dir);
    end
end
save(fullfile(fig_dir,'GAM_struct.mat'),'gam','-v7.3');
%}

%% Gausian Mixture Model (GMM) Clustering - By Region (combined subjects)
%{
% Initialize struct for GMM
gmm = struct();

%%% Initialization Parameters
% Set the number of quanitle bins for stratifying by the EPVS SWP
n_quantile = 4;
% Maximum number of samples per GMM
n_max = 1e6;

%%% Iterate over raw/pdif and front/occip
dtype = {'raw','pdif'};
regs = {'front','occip'};
% Iterate raw/pdif
for ii = 1:2
    % Iterate front/occip
    for j = 1:2
        % Select data type and region
        d = dtype{ii};
        reg = regs{j};
        % Create title substring from data type and region
        tstr = string(sprintf('%s %s',string(d),string(reg)));
        % set z-score normalization boolean
        if strcmp(d,'raw')
            zscore = true;
        else
            zscore = false;
        end
        % Run GMM clustering
        gmm.(d).(reg) = cluster_psoct_gmm(T_front_raw,...
                                         'NQuantileBins',n_quantile,...
                                         'MaxSample',n_max,'ZScore',zscore,...
                                         'TitleStr',tstr,...
                                         'dirout',fig_dir);
    end
end

% Save GMM
full_fout = fullfile(fig_dir,'GMM_struct.mat');
sprintf('Saving GMM data structs to %s\n',full_fout)
save(full_fout,'gmm','-v7.3');
%}

%% Generalized Additive Model (GAM) - By Subject/Region
% Initialize GAM output struct
gam = struct();

%%% Sampling and Learning
% Set maximum number of samples
n_max = 1e6;
% Set learning parameters (avoid over fitting)
ntrees = 100;
max_splits = 4;
learn_rate = 0.05;
% Set learning parameters (sharper contrast near EPVS)
% ntrees = 200;
% max_splits = 10;
% learn_rate = 0.1;
%%% Plotting confidence intervals settings
nbootstrap = 100;
ci_alpha = 0.05;

% Initialize data type (raw/pdif)
dtype = {'raw','pdif'};
% Iterate subjects
subs = fields(subjects);
for ii = 1:numel(fields(subjects))
    sub = subs{ii};
    regs = fields(subjects.(sub));
    % Update output directory with subject name sub directory
    subdir = fullfile(fig_dir, sub);
    if ~exist(subdir,'dir'); mkdir(subdir); end
    % Iterate regions
    for j = 1:numel(regs)
        reg = regs{j};
        for k = 1:numel(dtype)
            %%% Retrieve optical properties        
            fprintf('\tSubject %s region %s\n',sub,reg)
            d = dtype{k};
            % Create title substring from data type and region
            tstr = string(sprintf('%s %s %s',string(d),string(sub),...
                                  string(reg)));
            % set z-score normalization boolean
            if strcmp(d,'raw')
                zscore = true;
            else
                zscore = false;
            end            
            % Run GMM clustering
            t = T.(d).(sub).(reg);
            fprintf('Running GAM for %s %s %s\n',d,sub,reg)
            gam.(d).(sub).(reg) = fit_swp_psoct_gam(t,...
                            'NumTrees',ntrees,...
                            'MaxSplits',max_splits,'LearnRate',learn_rate,...
                            'MaxSample',n_max,'NBootstrap',nbootstrap,...
                            'CIAlpha',ci_alpha,'TitleStr',tstr,...
                            'dirout',subdir);
        end
    end
end
pause(1); close all;
% Save GAM
full_fout = fullfile(fig_dir,'GAM_struct_subjects.mat');
sprintf('Saving GAM data structs to %s\n',full_fout)
save(full_fout,'gam','-v7.3');

%%% Gausian Mixture Model (GMM) Clustering - By Subject/Region
% Initialize struct for GMM
gmm = struct();

%%% Initialization Parameters
% Set the number of quanitle bins for stratifying by the EPVS SWP
n_quantile = 4;
% Maximum number of samples per GMM
n_max = 1e6;

% Initialize data type (raw/pdif)
dtype = {'raw','pdif'};
% Iterate subjects
subs = fields(subjects);
for ii = 1:numel(fields(subjects))
    sub = subs{ii};
    regs = fields(subjects.(sub));
    % Update output directory with subject name sub directory
    subdir = fullfile(fig_dir, sub);
    if ~exist(subdir,'dir'); mkdir(subdir); end
    % Iterate regions
    for j = 1:numel(regs)
        reg = regs{j};
        for k = 1:numel(dtype)
            %%% Retrieve optical properties        
            fprintf('\tSubject %s region %s\n',sub,reg)
            d = dtype{k};
            % Create title substring from data type and region
            tstr = string(sprintf('%s %s %s',string(d),string(sub),...
                                  string(reg)));
            % set z-score normalization boolean
            if strcmp(d,'raw')
                zscore = true;
            else
                zscore = false;
            end
            % Run GMM clustering
            t = T.(d).(sub).(reg);
            gmm.(d).(sub).(reg) = cluster_psoct_gmm(t,...
                                             'NQuantileBins',n_quantile,...
                                             'MaxSample',n_max,...
                                             'ZScore',zscore,...
                                             'TitleStr',tstr,...
                                             'dirout',subdir);
        end
    end
end

% Save GMM
full_fout = fullfile(fig_dir,'GMM_struct_subjects.mat');
sprintf('Saving GMM data structs to %s\n',full_fout)
save(full_fout,'gmm','-v7.3');

%% Combine across all subjects. Separate by regions
function region_data = combine_subjects(heat_pairs)
% Combine heatpars across all subjects, separated by region
% Output: region_data.(region) = Nx4 matrix

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