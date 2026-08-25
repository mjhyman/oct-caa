%% Generate Generalized Additive Model (GAM) for EPVS/Ves-SWP, mus, ret
% Import the EPVS and Vessel SWP heatmaps (without log transform)
% Import the .MAT of the optical properties
% Bin the data points along x-axis and y-axis
%   x-axis = EPVS SWP
%   y-axis = Optical Property
% TODO: combine CAA6 into control

%% Initialization
clearvars -except caa caa6 caa17 caa22 caa25 caa26 swp_struct
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
    
    %%% Remove top-level struct, consolidate, remove extra struct
    caa.caa6 = caa6.caa6;       clear caa6;
    caa.caa17 = caa17.caa17;    clear caa17;
    caa.caa22 = caa22.caa22;    clear caa22;
    caa.caa25 = caa25.caa25;    clear caa25;
    caa.caa26 = caa26.caa26;    clear caa26;

    %%% Remove orientation & convert retardance to single
    subjects = fieldnames(caa);
    for i = 1:numel(subjects)
        subj = subjects{i};
        regions = fieldnames(caa.(subj));
    
        for j = 1:numel(regions)
            reg = regions{j};
    
            % remove orient_rgb if it exists
            if isfield(caa.(subj).(reg), 'orient_rgb')
                caa.(subj).(reg) = rmfield(caa.(subj).(reg), 'orient_rgb');
            end
            if isfield(caa.(subj).(reg), 'orient')
                caa.(subj).(reg) = rmfield(caa.(subj).(reg), 'orient');
            end
        end
    end
end

%% Load EPVS SWP & Vessel SWP
if flag_load_swp_structs
    % Struct for storing SWP
    % swp_struct = struct();
    
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
    clear tmp;
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
fprintf('Finished making 4D arrays of EPVS density vs. optical prop\n')
% Remove SWP to save data
clear swp_struct swp_ves swp_epvs swp swp_pdif mus_pdif ret_pdif


%% Combine severe (CAA22 + CAA25) then control (CAA26 + CAA6)

%%% Severe
% Create copy of raw and pdif
swp_op_raw = swp_op.raw;
swp_op_pdif = swp_op.pdif;
% Define fields to keep
keepFields = {'caa22', 'caa25'};
% Remove all other fields
swp_op_raw_severe = rmfield(swp_op_raw, setdiff(fieldnames(swp_op_raw), keepFields));
swp_op_pdif_severe = rmfield(swp_op_pdif, setdiff(fieldnames(swp_op_pdif), keepFields));
% Combine caa22 and caa25
swp_op_raw_severe = combine_subjects(swp_op_raw_severe);
swp_op_pdif_severe = combine_subjects(swp_op_pdif_severe);

%%% Control
% Create copy of raw and pdif
swp_op_raw = swp_op.raw;
swp_op_pdif = swp_op.pdif;
% Define fields to keep
keepFields = {'caa26', 'caa6'};
% Remove all other fields
swp_op_raw_control = rmfield(swp_op_raw, setdiff(fieldnames(swp_op_raw), keepFields));
swp_op_pdif_control = rmfield(swp_op_pdif, setdiff(fieldnames(swp_op_pdif), keepFields));
% Combine caa22 and caa25
swp_op_raw_control = combine_subjects(swp_op_raw_control);
swp_op_pdif_control = combine_subjects(swp_op_pdif_control);

%% Create table for each subject/region
fprintf('Creating 4xN Tables\n')
vnames = {'scattering','retardance','ves_swp','epv_swp'};
subs = fields(subjects);
% Iterate subjects
for ii = 1:numel(fields(subjects))
    sub = subs{ii};
    regs = fields(subjects.(sub));
    % Iterate regions for subject
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
% Add the combined control
fprintf('\tAdding the combined control to the table\n')
T.raw.ctl.front = array2table(swp_op_raw_control.front,'VariableNames',vnames);
T.raw.ctl.occip = array2table(swp_op_raw_control.occip,'VariableNames',vnames);
T.pdif.ctl.front = array2table(swp_op_pdif_control.front,'VariableNames',vnames);
T.pdif.ctl.occip = array2table(swp_op_pdif_control.occip,'VariableNames',vnames);
% Add the combined severe
fprintf('\tAdding the combined severe to the table\n')
T.raw.sev.front = array2table(swp_op_raw_severe.front,'VariableNames',vnames);
T.raw.sev.occip = array2table(swp_op_raw_severe.occip,'VariableNames',vnames);
T.pdif.sev.front = array2table(swp_op_pdif_severe.front,'VariableNames',vnames);
T.pdif.sev.occip = array2table(swp_op_pdif_severe.occip,'VariableNames',vnames);

%% Isosurface Comparison — Severe vs Control (shared grid, tuned per group)
%%% Bootstrap settings for exploration vs. figure quality
% The higher bootsrap leads to cleaner plots (less wobble)
%   -- true while iterating
%   -- false for manuscript figures
EXPLORATORY = true;
if EXPLORATORY
    nbootstrap = 100;     % fast, seed-jittery — NOT for final figures
    n_max_cmp  = 2e5;     % subsample for speed while tuning
else
    nbootstrap = 500;     % reproducible percentile-CI band for the paper
    n_max_cmp  = 1e6;     % your full sample
end
% output directory
cmp_dir = fullfile(fig_dir,'severe_vs_control');
if ~exist(cmp_dir,'dir'); mkdir(cmp_dir); end

% Regions
regs = {'front','occip'};
ci_alpha = 0.05;
gam_cmp = struct();
for ii = 1:numel(regs)
    fprintf('\nStart GAM for severe vs. control for %s\n',regs{ii})
    reg = regs{ii};
    tstr = sprintf('pdif_severe_vs_control_%s', reg);
    fprintf('Comparing severe vs control: %s\n', reg);
    gam_cmp.pdif.(reg) = compare_isosurfaces_severe_control(...
        T.pdif.sev.(reg), T.pdif.ctl.(reg), ...
        'TuneParams', true, ...
        'MaxSplitsGrid', [4 10 20], 'NumTreesGrid', [100 200 500], ...
        'MaxSample', n_max_cmp, ...
        'NBootstrap', nbootstrap, 'CIAlpha', ci_alpha, ...
        'GridRange', 'intersection', 'SharedPercentiles', true, ...
        'TitleStr', tstr, 'dirout', cmp_dir);
    pause(1); close all;
    fprintf('\nFinished GAM for severe vs. control for %s\n',regs{ii})
end
dt = datetime("now",'TimeZone','local','Format','dd-MMM-yyyy');
save(fullfile(fig_dir, strcat('GAM_compare_severe_control_',string(dt),'.mat')), ...
     'gam_cmp','-v7.3');


%%% Isosurface per Subject/Region (tuned per dataset)
%%% Bootstrap settings for exploration vs. figure quality
% The higher bootsrap leads to cleaner plots (less wobble)
%   -- true while iterating
%   -- false for manuscript figures
EXPLORATORY = true;
if EXPLORATORY
    nbootstrap = 100;
    n_max_iso  = 2e5;
else
    nbootstrap = 500;
    n_max_iso  = 1e6;
end
% Initialize subjects struct
subjects.sev.front = []; subjects.sev.occip = [];
subjects.ctl.front = []; subjects.ctl.occip = [];
% Output directory
iso_dir_base = fullfile(fig_dir,'by_subject');

% Initialize variables for looping
ci_alpha = 0.05;
dtype = {'pdif'};              % or {'raw','pdif'}
gam_iso = struct();
subs = fields(subjects);
for ii = 1:numel(subs)
    sub = subs{ii};
    regs = fields(subjects.(sub));
    subdir = fullfile(iso_dir_base, sub);
    if ~exist(subdir,'dir'); mkdir(subdir); end
    for j = 1:numel(regs)
        reg = regs{j};
        fprintf('\nStart GAM for %s %s\n',subs{ii},regs{j})
        for k = 1:numel(dtype)
            d = dtype{k};
            if ~isfield(T.(d),sub) || ~isfield(T.(d).(sub),reg); continue; end
            tstr = sprintf('%s_%s_%s', d, sub, reg);
            fprintf('Isosurface for %s\n', tstr);
            gam_iso.(d).(sub).(reg) = fit_isosurface_by_dataset(...
                T.(d).(sub).(reg), ...
                'TuneParams', true, ...
                'MaxSplitsGrid', [4 10 20], 'NumTreesGrid', [100 200 500], ...
                'MaxSample', n_max_iso, ...
                'NBootstrap', nbootstrap, 'CIAlpha', ci_alpha, ...
                'TitleStr', tstr, 'dirout', subdir);
            pause(1); close all;
        end
        fprintf('\nFinished GAM for %s %s\n',subs{ii},regs{j})
    end
end
% Save the gam_iso
dt = datetime("now",'TimeZone','local','Format','dd-MMM-yyyy');
save(fullfile(fig_dir, strcat('GAM_subject_region_',string(dt),'.mat')), ...
     'gam_iso','-v7.3');

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
