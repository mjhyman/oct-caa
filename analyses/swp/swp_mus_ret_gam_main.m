%% Generate Generalized Additive Model (GAM) for EPVS/Ves-SWP, mus, ret
% Import the EPVS and Vessel SWP heatmaps (without log transform)
% Import the .MAT of the optical properties
% Bin the data points along x-axis and y-axis
%   x-axis = EPVS SWP
%   y-axis = Optical Property
%
% Bug fixes (Oct 2026):
%   1. combine_subjects: every subject was written from row 1, overwriting
%      the previous subject and leaving min(n_A, n_B) all-zero rows in the
%      combined control/severe matrices. Now concatenates correctly.
%   2. save() calls: variable name and '-v7.3' were passed inside fullfile,
%      saving the whole workspace to a malformed path. Fixed.
%   3. 'subjects' was overwritten with a cell array in the orientation-
%      removal loop and only defined inside the SWP-loading block, so the
%      script failed when flag_load_swp_structs = false. Now defined once.
%   4. Optional equal weighting of subjects in the combined groups
%      (EQUAL_WEIGHT_SUBJECTS); off by default to match previous behaviour.
%   5. Fixed RNG seed for reproducible subsampling/bootstrap.

%% Initialization
clearvars -except caa swp_struct
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
% Reproducibility
rng(1);

%%% Directories
% Directory containing seg, mus, ret, mask, epvs structs
mat_dir = '/projectnb/npbssmic/ns/CAA/';
% Heatmap directory
swp_dir = '/projectnb/npbssmic/ns/CAA/swp';
% Figure output directory
fig_dir = '/projectnb/npbssmic/ns/CAA/figures/swp_gam_gmm/';
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

%%% Common filename for importing
epvs_fbase = '_swp_voxelwise_radius_500_exp_2_interpolated_heatmap.mat';
ves_fbase  = '_swp_voxelwise_ves_radius_500_exp_2_interpolated_heatmap.mat';

%%% Flags for importing data (set false to reuse what is in the workspace)
flag_load_swp_structs = ~exist('swp_struct','var');
flag_load_caa_structs = ~exist('caa','var');

%%% Number of bins for the x-axis along SWP
nbin = 100; %#ok<NASGU>

%%% Flag for rough GAM or high fidelity
EXPLORATORY = false;

%%% Weight subjects equally in the combined control/severe groups by
% subsampling each subject to the smallest subject's voxel count (per
% region). false = voxel-weighted concatenation (previous behaviour).
EQUAL_WEIGHT_SUBJECTS = false;

%%% Subject/region layout (defined once; used throughout)
subjects = struct();
subjects.caa6.front  = [];
subjects.caa6.occip  = [];
subjects.caa17.occip = [];
subjects.caa22.front = [];
subjects.caa22.occip = [];
subjects.caa25.front = [];
subjects.caa25.occip = [];
subjects.caa26.front = [];
subjects.caa26.occip = [];

%% Load matlab structs
if flag_load_caa_structs
    fprintf('Loading CAA6\n')
    caa6 = load(fullfile(mat_dir,"/caa6/caa6.mat"));
    fprintf('Finished loading CAA6\n')

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
    caa.caa6  = caa6.caa6;      clear caa6;
    caa.caa17 = caa17.caa17;    clear caa17;
    caa.caa22 = caa22.caa22;    clear caa22;
    caa.caa25 = caa25.caa25;    clear caa25;
    caa.caa26 = caa26.caa26;    clear caa26;

    %%% Remove orientation
    caa_subs = fieldnames(caa);
    for i = 1:numel(caa_subs)
        subj = caa_subs{i};
        caa_regs = fieldnames(caa.(subj));
        for j = 1:numel(caa_regs)
            reg = caa_regs{j};
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
    swp_struct = struct();
    subs = fieldnames(subjects);
    for ii = 1:numel(subs)
        sub = subs{ii};
        regs = fieldnames(subjects.(sub));
        for j = 1:numel(regs)
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

%% Create Nx4 matrix (mus, ret, vessel-SWP, EPVS-SWP)
% Each parenchyma voxel has a value for EPVS-SWP, vessel-SWP, mus, ret:
%   - column 1 = mus
%   - column 2 = retardance
%   - column 3 = Vessel SWP
%   - column 4 = EPVS SWP
swp_op = struct();
fprintf('Creating Nx4 matrices\n')
subs = fieldnames(subjects);
for ii = 1:numel(subs)
    sub = subs{ii};
    regs = fieldnames(subjects.(sub));
    for j = 1:numel(regs)
        %%% Retrieve region and optical properties
        reg = regs{j};
        fprintf('\tSubject %s region %s\n',sub,reg)
        mus = caa.(sub).(reg).mus;
        ret = caa.(sub).(reg).ret_full;

        %%% Percentage difference (fractional) from stored WM median
        mus_med = parench_median.(sub).(reg).med.mus;
        ret_med = parench_median.(sub).(reg).med.ret;
        mus_pdif = (mus - mus_med) ./ mus_med;
        ret_pdif = (ret - ret_med) ./ ret_med;

        %%% Retrieve SWP for EPVS and Vessel
        swp_epvs = swp_struct.epvs.(sub).(reg);
        swp_ves  = swp_struct.ves.(sub).(reg);

        %%% Exclude vessels and EPVS from white matter mask (mask_wm)
        mask_wm = caa.(sub).(reg).mask_wm;
        ves  = caa.(sub).(reg).seg;
        epvs = caa.(sub).(reg).epvs;
        mask_wm(ves)  = 0;
        mask_wm(epvs) = 0;

        %%% Apply WM mask
        swp_epvs = swp_epvs(mask_wm);
        swp_ves  = swp_ves(mask_wm);
        mus      = mus(mask_wm);
        ret      = ret(mask_wm);
        mus_pdif = mus_pdif(mask_wm);
        ret_pdif = ret_pdif(mask_wm);

        %%% Create matrix
        swp_op.raw.(sub).(reg)  = [mus(:), ret(:), swp_ves(:), swp_epvs(:)];
        swp_op.pdif.(sub).(reg) = [mus_pdif(:), ret_pdif(:), swp_ves(:), swp_epvs(:)];
    end
end
fprintf('Finished making Nx4 arrays of SWP vs. optical properties\n')
% Remove SWP to save memory (re-run will reload it)
clear swp_struct swp_ves swp_epvs mus_pdif ret_pdif mus ret mask_wm ves epvs

%% Combine severe (CAA22 + CAA25) then control (CAA26 + CAA6)

%%% Severe
keepFields = {'caa22', 'caa25'};
swp_op_raw_severe  = rmfield(swp_op.raw,  setdiff(fieldnames(swp_op.raw),  keepFields));
swp_op_pdif_severe = rmfield(swp_op.pdif, setdiff(fieldnames(swp_op.pdif), keepFields));
swp_op_raw_severe  = combine_subjects(swp_op_raw_severe,  EQUAL_WEIGHT_SUBJECTS);
swp_op_pdif_severe = combine_subjects(swp_op_pdif_severe, EQUAL_WEIGHT_SUBJECTS);

%%% Control
keepFields = {'caa26', 'caa6'};
swp_op_raw_control  = rmfield(swp_op.raw,  setdiff(fieldnames(swp_op.raw),  keepFields));
swp_op_pdif_control = rmfield(swp_op.pdif, setdiff(fieldnames(swp_op.pdif), keepFields));
swp_op_raw_control  = combine_subjects(swp_op_raw_control,  EQUAL_WEIGHT_SUBJECTS);
swp_op_pdif_control = combine_subjects(swp_op_pdif_control, EQUAL_WEIGHT_SUBJECTS);

%%% Sanity check: combined row counts and no padding rows
check_combined(swp_op.raw, {'caa22','caa25'}, swp_op_raw_severe,  'severe',  EQUAL_WEIGHT_SUBJECTS);
check_combined(swp_op.raw, {'caa26','caa6'},  swp_op_raw_control, 'control', EQUAL_WEIGHT_SUBJECTS);

%% Create table for each subject/region
fprintf('Creating Nx4 Tables\n')
vnames = {'scattering','retardance','ves_swp','epv_swp'};
T = struct();
subs = fieldnames(subjects);
for ii = 1:numel(subs)
    sub = subs{ii};
    regs = fieldnames(subjects.(sub));
    for j = 1:numel(regs)
        reg = regs{j};
        fprintf('\tSubject %s region %s\n',sub,reg)
        T.raw.(sub).(reg)  = array2table(swp_op.raw.(sub).(reg),  'VariableNames',vnames);
        T.pdif.(sub).(reg) = array2table(swp_op.pdif.(sub).(reg), 'VariableNames',vnames);
    end
end
% Add the combined control
fprintf('\tAdding the combined control to the table\n')
T.raw.ctl.front  = array2table(swp_op_raw_control.front,  'VariableNames',vnames);
T.raw.ctl.occip  = array2table(swp_op_raw_control.occip,  'VariableNames',vnames);
T.pdif.ctl.front = array2table(swp_op_pdif_control.front, 'VariableNames',vnames);
T.pdif.ctl.occip = array2table(swp_op_pdif_control.occip, 'VariableNames',vnames);
% Add the combined severe
fprintf('\tAdding the combined severe to the table\n')
T.raw.sev.front  = array2table(swp_op_raw_severe.front,  'VariableNames',vnames);
T.raw.sev.occip  = array2table(swp_op_raw_severe.occip,  'VariableNames',vnames);
T.pdif.sev.front = array2table(swp_op_pdif_severe.front, 'VariableNames',vnames);
T.pdif.sev.occip = array2table(swp_op_pdif_severe.occip, 'VariableNames',vnames);

%% Run tag for output filenames
dt = datetime("now",'TimeZone','local','Format','dd-MMM-yyyy');
if EXPLORATORY
    run_tag = sprintf('nboot100_EXPLORATORY_%s', string(dt));
else
    run_tag = sprintf('nboot500_%s', string(dt));
end
if EQUAL_WEIGHT_SUBJECTS
    run_tag = ['eqw_' run_tag];
end

%% Isosurface Comparison — Severe vs Control (shared grid, tuned per group)
%%% Bootstrap settings for exploration vs. figure quality
if EXPLORATORY
    nbootstrap = 100;     % fast, seed-jittery — NOT for final figures
    n_max_cmp  = 2e5;     % subsample for speed while tuning
else
    nbootstrap = 500;     % reproducible percentile-CI band for the paper
    n_max_cmp  = 1e6;     % full sample
end
cmp_dir = fullfile(fig_dir,'severe_vs_control');
if ~exist(cmp_dir,'dir'); mkdir(cmp_dir); end

regs = {'front','occip'};
ci_alpha = 0.05;
gam_cmp = struct();
for ii = 1:numel(regs)
    reg = regs{ii};
    fprintf('\nStart GAM for severe vs. control for %s\n', reg)
    tstr = sprintf('pdif_severe_vs_control_%s', reg);
    gam_cmp.pdif.(reg) = compare_isosurfaces_severe_control(...
        T.pdif.sev.(reg), T.pdif.ctl.(reg), ...
        'TuneParams', true, ...
        'MaxSplitsGrid', [4 10 20], 'NumTreesGrid', [100 200 500], ...
        'MaxSample', n_max_cmp, ...
        'NBootstrap', nbootstrap, 'CIAlpha', ci_alpha, ...
        'GridRange', 'intersection', 'SharedPercentiles', true, ...
        'TitleStr', tstr, 'dirout', cmp_dir);
    pause(1); close all;
    fprintf('\nFinished GAM for severe vs. control for %s\n', reg)
end
fout = sprintf('GAM_compare_severe_control_%s.mat', run_tag);
save(fullfile(fig_dir, fout), 'gam_cmp', '-v7.3');
fprintf('Saved %s\n', fullfile(fig_dir, fout));

%% Isosurface per Subject/Region (tuned per dataset)
if EXPLORATORY
    nbootstrap = 100;
    n_max_iso  = 2e5;
else
    nbootstrap = 500;
    n_max_iso  = 1e6;
end
% Add the combined groups to the fitting list
iso_subjects = subjects;
iso_subjects.sev.front = []; iso_subjects.sev.occip = [];
iso_subjects.ctl.front = []; iso_subjects.ctl.occip = [];
iso_dir_base = fullfile(fig_dir,'by_subject');

ci_alpha = 0.05;
dtype = {'pdif'};              % or {'raw','pdif'}
gam_iso = struct();
subs = fieldnames(iso_subjects);
for ii = 1:numel(subs)
    sub = subs{ii};
    regs = fieldnames(iso_subjects.(sub));
    subdir = fullfile(iso_dir_base, sub);
    if ~exist(subdir,'dir'); mkdir(subdir); end
    for j = 1:numel(regs)
        reg = regs{j};
        fprintf('\nStart GAM for %s %s\n', sub, reg)
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
        fprintf('\nFinished GAM for %s %s\n', sub, reg)
    end
end
fout = sprintf('GAM_subject_region_%s.mat', run_tag);
save(fullfile(fig_dir, fout), 'gam_iso', '-v7.3');
fprintf('Saved %s\n', fullfile(fig_dir, fout));


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function region_data = combine_subjects(heat_pairs, equal_weight)
% Concatenate subjects within each region.
% Output: region_data.(region) = Nx4 matrix
% equal_weight = true subsamples each subject (without replacement) to the
% smallest subject's row count in that region, so every subject contributes
% equally to the combined fit.
    if nargin < 2; equal_weight = false; end
    region_data = struct();
    subj = fieldnames(heat_pairs);

    % Collect per-region blocks
    blocks = struct();
    for i = 1:numel(subj)
        regs = fieldnames(heat_pairs.(subj{i}));
        for j = 1:numel(regs)
            reg = regs{j};
            if ~isfield(blocks, reg); blocks.(reg) = {}; end
            blocks.(reg){end+1} = heat_pairs.(subj{i}).(reg);
        end
    end

    % Concatenate (optionally equal-weighted)
    regs = fieldnames(blocks);
    for j = 1:numel(regs)
        reg = regs{j};
        B = blocks.(reg);
        if equal_weight
            n_min = min(cellfun(@(x) size(x,1), B));
            for b = 1:numel(B)
                idx  = randperm(size(B{b},1), n_min);
                B{b} = B{b}(idx, :);
            end
        end
        region_data.(reg) = vertcat(B{:});
    end
end

function check_combined(per_subject, keep, combined, label, equal_weight)
% Verify the combined matrix has the expected number of rows and no
% all-zero padding rows (the symptom of the old combine_subjects bug).
    regs = fieldnames(combined);
    for j = 1:numel(regs)
        reg = regs{j};
        n_sub = [];
        for k = 1:numel(keep)
            if isfield(per_subject, keep{k}) && isfield(per_subject.(keep{k}), reg)
                n_sub(end+1) = size(per_subject.(keep{k}).(reg), 1); %#ok<AGROW>
            end
        end
        if equal_weight
            n_exp = numel(n_sub) * min(n_sub);
        else
            n_exp = sum(n_sub);
        end
        n_got  = size(combined.(reg), 1);
        n_zero = sum(all(combined.(reg) == 0, 2));
        fprintf('\t%s %s: %d rows (expected %d), all-zero rows: %d\n', ...
                label, reg, n_got, n_exp, n_zero);
        assert(n_got == n_exp, 'Combined %s %s has %d rows, expected %d.', ...
               label, reg, n_got, n_exp);
    end
end