%% Create a GAM per stain x group (control / severe) x region
% Overview:
%{
Mirrors the optical-property vs SWP GAM pipeline:
  - import each stain's SWP struct
  - build the predictor matrix [epv_swp, ves_swp, ves_epv_swp, stain] per
    section, where ves_epv_swp is the EXPLICIT product interaction term
    (epv_swp .* ves_swp) -- the same third predictor the optical GAMs use, so
    predictions / isosurfaces stay internally consistent
  - combine subjects into the corrected 2x2 groups (frontal / occipital):
        control = CAA_26 + CAA_6
        severe  = CAA_22 + CAA_25
    (CAA_17 "moderate" is excluded -- not part of the corrected 2x2 design;
     the old "mild = CAA_6" label was erroneous, CAA_6 is a second control)
  - fit one gradient-boosted GAM per group x region, leaving control_* and
    severe_* fits side by side for the downstream severe-vs-control
    isosurface / bootstrap-difference comparison.
  - SUBJECT-LEVEL: additionally fit one GAM per subject x region (all five
    subjects, CAA_17 included) and evaluate every subject in a region on ONE
    shared prediction grid, so the supplementary subject panels are directly
    comparable. Saved as histo_gam_subj.(stain).(<subject>_<regioncode>) in
    GAM_histo_subjects_<tag>_<date>.mat, consumed by
    plot_histo_swp_gam_stages_subjects.m.
%}

%% Top-level settings
clc; clear; close all;

% EXPLORATORY toggle: 100 bootstraps for quick runs, 500 for manuscript figures
EXPLORATORY = false;

%%% Directories (SCC)
% Input directory
hdir = '/projectnb/npbssmic/ns/CAA/histology/swp';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/histology/gam';

%%% Add top-level directory to path
addpath(genpath('/projectnb/npbssmic/s/mhyman/oct-caa'));

%%% Filenames
lhe_fname  = 'histo_swp_hmatched_zscore__LHE_voxelwise_radius4mm_ds4_19-Mar-2026.mat';
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

% Import LHE
lhe = load(fullfile(hdir, lhe_fname));
lhe = lhe.lhe;
% Import GFAP
gfap = load(fullfile(hdir, gfap_fname));
gfap = gfap.gfap;
% Import CD68
cd68 = load(fullfile(hdir, cd68_fname));
cd68 = cd68.cd68;

% Add all stains to the struct "stains"
stains = struct();
stains.cd68 = cd68;
stains.gfap = gfap;
stains.lhe  = lhe;

%% Build predictor matrix per section: [epv_swp, ves_swp, ves_epv_swp, stain]
% ves_epv_swp is the explicit product interaction term. Computed AFTER masking
% so it equals epv_swp .* ves_swp element-for-element on the tissue voxels.
stain_names = fieldnames(stains);
for ii = 1:numel(stain_names)
    % Get the current stain name
    stain = stains.(stain_names{ii});
    % print stain name to console
    fprintf('\nRunning stain %s', stain_names{ii})
    %%% Iterate sections within stain
    for j = 1:numel(stain)
        % print section # to console
        fprintf('\n\tSection %i / %i\n', j, numel(stain))
        % Extract vessel SWP + EPVS SWP
        ves_swp  = stain(j).ves_swp;
        epvs_swp = stain(j).epvs_swp;

        %%% Retrieve tissue mask, ves mask, epvs mask
        mask = stain(j).mask;
        ves  = stain(j).ves;
        epvs = stain(j).epvs;
        % Set the EPVS and vessels to 0 within mask
        mask(ves)  = 0;
        mask(epvs) = 0;

        %%% Apply mask to SWP and stain section
        ves_swp  = ves_swp(mask);
        epvs_swp = epvs_swp(mask);
        section  = stain(j).z_stain(mask);

        %%% Explicit interaction term (must equal epv_swp .* ves_swp)
        ves_epv_swp = epvs_swp .* ves_swp;

        %%% Store matrix: [epv_swp, ves_swp, ves_epv_swp, stain]
        epv_ves_swp_stain = [epvs_swp, ves_swp, ves_epv_swp, section];
        stains.(stain_names{ii})(j).epv_ves_swp_stain = epv_ves_swp_stain;
    end
end

%% Combine subjects into the corrected 2x2 groups (control / severe) x region
% Patterns match both baseName conventions (CAA22_1 and CAA_22_1). Frontal = _1,
% occipital = _7. control = CAA_26 + CAA_6, severe = CAA_22 + CAA_25.
groups = struct( ...
    'control_front', {{'CAA26_1','CAA_26_1','CAA6_1','CAA_6_1'}}, ...
    'control_occip', {{'CAA26_7','CAA_26_7','CAA6_7','CAA_6_7'}}, ...
    'severe_front',  {{'CAA22_1','CAA_22_1','CAA25_1','CAA_25_1'}}, ...
    'severe_occip',  {{'CAA22_7','CAA_22_7','CAA25_7','CAA_25_7'}});

grouped    = struct();
group_flds = fieldnames(groups);
fprintf('\nCombining subjects into 2x2 groups\n')
for ii = 1:numel(stain_names)
    sname = stain_names{ii};
    stain = stains.(sname);
    for g = 1:numel(group_flds)
        gname = group_flds{g};
        [M, n_sec] = combine_group(stain, groups.(gname));
        if isempty(M)
            fprintf('  %s / %s: no matching sections, skipping\n', sname, gname);
            continue
        end
        grouped.(sname).(gname) = M;
        fprintf('  %s / %s: combined %d section(s), %d rows\n', ...
                sname, gname, n_sec, size(M,1));
    end
end

%% Create a table for each stain x group
vnames = {'epv_swp','ves_swp','ves_epv_swp','stain'};
T = struct();
fprintf('\nCreating tables\n')
for ii = 1:numel(stain_names)
    sname = stain_names{ii};
    if ~isfield(grouped, sname), continue; end
    gflds = fieldnames(grouped.(sname));
    for g = 1:numel(gflds)
        gname = gflds{g};
        T.(sname).(gname) = array2table(grouped.(sname).(gname), ...
                                        'VariableNames', vnames);
    end
end

%% Create GAM for each stain and group

%%% Sampling and Learning
% Set maximum number of samples
n_max = 1e6;
% Set learning parameters (avoid over fitting)
ntrees     = 500;
max_splits = 4;
learn_rate = 0.05;

%%% Confidence intervals for slices (bootstrap count set by EXPLORATORY toggle)
if EXPLORATORY
    nbootstrap = 100;   % quick exploratory runs
else
    nbootstrap = 500;   % manuscript figures
end
ci_alpha = 0.05;

%%% Run GAM for each stain and group
histo_gam   = struct();
group_names = {'control_front','control_occip','severe_front','severe_occip'};

for ii = 1:numel(stain_names)
    sname = stain_names{ii};
    if ~isfield(T, sname), continue; end
    for g = 1:numel(group_names)
        gname = group_names{g};
        % Skip groups with no data for this stain
        if ~isfield(T.(sname), gname)
            fprintf('  Skipping %s / %s (no data)\n', sname, gname);
            continue
        end
        % Set the subdirectory
        subdir = fullfile(figdir, sname);
        % Define title string (which is also filename)
        tstr = strcat(sname, '_', gname);
        % NOTE: fit_swp_histo_gam must fit
        %       stain ~ epv_swp + ves_swp + ves_epv_swp
        % i.e. include the explicit product interaction as a third predictor,
        % matching tune_swp_psoct_gam in the optical pipeline. If the function
        % currently uses only the first two columns, update it to pick up
        % ves_epv_swp as well.
        histo_gam.(sname).(gname) = fit_swp_histo_gam( ...
            T.(sname).(gname), ...
            'NumTrees', ntrees, ...
            'MaxSplits', max_splits, 'LearnRate', learn_rate, ...
            'MaxSample', n_max, 'NBootstrap', nbootstrap, ...
            'CIAlpha', ci_alpha, 'TitleStr', tstr, ...
            'dirout', subdir, ...
            'stain_name', sname, ...
            'section_name', gname);
    end
end

% Save GAM (tag exploratory vs manuscript so runs don't overwrite each other)
dt = datetime("now", 'TimeZone', 'local', 'Format', 'dd-MMM-yyyy');
if EXPLORATORY, run_tag = 'exploratory_'; else, run_tag = 'manuscript_'; end
fout      = strcat('GAM_histo_', run_tag, string(dt), '.mat');
full_fout = fullfile(figdir, fout);
fprintf('Saving histo_gam data struct to %s\n', full_fout)
save(full_fout, 'histo_gam', '-v7.3');

%% Subject-level GAMs (one per subject x region) on a shared prediction grid
% Fits the same model (stain ~ epv_swp + ves_swp + ves_epv_swp) to each
% subject's section separately, including CAA_17 (moderate; occipital only).
% For each stain x region, every subject is then evaluated on ONE grid spanning
% the union of the subjects' [q_lo, q_hi] SWP ranges, with predictions outside
% a subject's own convex hull set to NaN. Shared grid edges = comparable panels.
%
% Output entry (histo_gam_subj.(stain).(<subject>_<regioncode>)):
%   fit         (ngrid x ngrid) predicted z-score, NaN outside support
%   ci_lo/ci_hi (ngrid x ngrid) row-bootstrap CI on the fit (NaN if disabled)
%   epv_vec, ves_vec, EPV, VES, support, n_rows, subject, region, stage,
%   model       the fitted RegressionGAM
%   table       long-format grid table
%
% NOTE on the subject-level CI: rows within one section are voxels of a single
% tissue volume, so a row bootstrap describes fitting uncertainty for that
% volume only. It says nothing about between-subject variability (that is what
% the side-by-side panels are for), and it must not be read as a group CI.
fprintf('\nFitting subject-level GAMs\n')

%%% Subject definitions (patterns match both baseName conventions)
subjects = struct( ...
    'CAA26', struct('stage','control',  'front',{{'CAA26_1','CAA_26_1'}}, 'occip',{{'CAA26_7','CAA_26_7'}}), ...
    'CAA6',  struct('stage','control',  'front',{{'CAA6_1', 'CAA_6_1' }}, 'occip',{{'CAA6_7', 'CAA_6_7' }}), ...
    'CAA17', struct('stage','moderate', 'front',{{'CAA17_1','CAA_17_1'}}, 'occip',{{'CAA17_7','CAA_17_7'}}), ...
    'CAA22', struct('stage','severe',   'front',{{'CAA22_1','CAA_22_1'}}, 'occip',{{'CAA22_7','CAA_22_7'}}), ...
    'CAA25', struct('stage','severe',   'front',{{'CAA25_1','CAA_25_1'}}, 'occip',{{'CAA25_7','CAA_25_7'}}));
subj_ids = fieldnames(subjects);

%%% Subject-level settings
ngrid_subj        = 60;      % grid points per axis (same as the diff grid)
q_lo_subj         = 0.01;    % per-axis range percentiles
q_hi_subj         = 0.99;
do_subj_bootstrap = false;   % true -> row-bootstrap CI on each subject's fit
nboot_subj        = nbootstrap;

subjdir = fullfile(figdir, 'subjects');
if ~exist(subjdir, 'dir'), mkdir(subjdir); end

histo_gam_subj = struct();
region_codes   = {'front','occip'};

for ii = 1:numel(stain_names)
    sname = stain_names{ii};
    stain = stains.(sname);

    for rr = 1:numel(region_codes)
        rcode = region_codes{rr};

        %%% 1) Collect each subject's table for this stain x region
        tbls = struct();
        for s = 1:numel(subj_ids)
            sid = subj_ids{s};
            [M, n_sec] = combine_group(stain, subjects.(sid).(rcode));
            if isempty(M)
                fprintf('  %s / %s %s: no section, skipping\n', sname, sid, rcode);
                continue
            end
            if n_sec > 1
                fprintf('  %s / %s %s: %d sections matched, concatenating\n', ...
                        sname, sid, rcode, n_sec);
            end
            tbls.(sid) = array2table(M, 'VariableNames', vnames);
        end
        present = fieldnames(tbls);
        if isempty(present), continue; end

        %%% 2) Shared grid for this stain x region (union of subject ranges)
        epv_lo = inf; epv_hi = -inf; ves_lo = inf; ves_hi = -inf;
        for s = 1:numel(present)
            tb = tbls.(present{s});
            epv_lo = min(epv_lo, quantile(tb.epv_swp, q_lo_subj));
            epv_hi = max(epv_hi, quantile(tb.epv_swp, q_hi_subj));
            ves_lo = min(ves_lo, quantile(tb.ves_swp, q_lo_subj));
            ves_hi = max(ves_hi, quantile(tb.ves_swp, q_hi_subj));
        end
        epv_vec = double(linspace(epv_lo, epv_hi, ngrid_subj));
        ves_vec = double(linspace(ves_lo, ves_hi, ngrid_subj));
        [EPV, VES] = meshgrid(epv_vec, ves_vec);
        VE = EPV .* VES;                           % interaction = epv_swp .* ves_swp
        Xg = [EPV(:), VES(:), VE(:)];              % column order matches training

        %%% 3) Fit + predict each subject on the shared grid
        for s = 1:numel(present)
            sid = present{s};
            tb  = tbls.(sid);
            fprintf('  %s / %s %s: fitting (%d rows)\n', sname, sid, rcode, height(tb));

            mdl = fit_gam_local(tb, ntrees, max_splits, learn_rate, n_max);

            F = reshape(predict(mdl, Xg), size(EPV));
            SUP = reshape(in_hull_2d(tb.epv_swp, tb.ves_swp, EPV(:), VES(:)), size(EPV));
            F(~SUP) = NaN;

            lo = nan(size(EPV)); hi = lo;
            if do_subj_bootstrap
                B = nan(numel(Xg(:,1)), nboot_subj);
                for b = 1:nboot_subj
                    mb = fit_gam_local(boot_tbl(tb), ntrees, max_splits, learn_rate, n_max);
                    B(:, b) = predict(mb, Xg);
                end
                a  = ci_alpha / 2;
                lo = reshape(quantile(B, a,     2), size(EPV));
                hi = reshape(quantile(B, 1 - a, 2), size(EPV));
                lo(~SUP) = NaN; hi(~SUP) = NaN;
            end

            tab = table(EPV(:), VES(:), VE(:), F(:), lo(:), hi(:), SUP(:), ...
                'VariableNames', {'epv_swp','ves_swp','ves_epv_swp', ...
                                  'fit','ci_lo','ci_hi','in_support'});

            key = sprintf('%s_%s', sid, rcode);        % e.g. CAA22_front
            histo_gam_subj.(sname).(key) = struct( ...
                'subject', sid, 'region', rcode, 'stage', subjects.(sid).stage, ...
                'epv_vec', epv_vec, 'ves_vec', ves_vec, 'EPV', EPV, 'VES', VES, ...
                'fit', F, 'ci_lo', lo, 'ci_hi', hi, 'support', SUP, ...
                'do_bootstrap', do_subj_bootstrap, 'n_rows', height(tb), ...
                'model', mdl, 'table', tab);

            writetable(tab, fullfile(subjdir, sprintf('%s_%s.csv', sname, key)));
        end
    end
end

% Save the subject-level struct
subj_fout = fullfile(figdir, strcat('GAM_histo_subjects_', run_tag, string(dt), '.mat'));
fprintf('Saving histo_gam_subj data struct to %s\n', subj_fout)
save(subj_fout, 'histo_gam_subj', '-v7.3');

%% Difference: severe - control on a common prediction grid (per stain x region)
% Evaluates the two already-fitted surfaces on one grid, restricted to the
% region of (epv_swp, ves_swp) space where BOTH groups have data (intersection
% of their percentile ranges AND of their convex hulls), then differences them
% (severe - control). ves_epv_swp on the grid is set to epv_swp .* ves_swp so
% the interaction column stays internally consistent with the training tables.

%%% Comparison settings
ngrid        = 60;      % grid points per axis
q_lo         = 0.01;    % lower percentile for the per-axis common range
q_hi         = 0.99;    % upper percentile
do_bootstrap = false;   % true -> bootstrap CIs on the difference (refits GAMs;
                        % first confirm fit_gam_local matches fit_swp_histo_gam)

diffdir = fullfile(figdir, 'severe_minus_control');
if ~exist(diffdir, 'dir'), mkdir(diffdir); end

region_suffix  = {'front','occip'};
histo_gam_diff = struct();

for ii = 1:numel(stain_names)
    sname = stain_names{ii};
    if ~isfield(histo_gam, sname), continue; end
    for rr = 1:numel(region_suffix)
        reg   = region_suffix{rr};
        cgrp  = ['control_' reg];
        sgrp  = ['severe_'  reg];

        % Need both a control and a severe fit for this stain/region
        if ~isfield(T.(sname), cgrp) || ~isfield(T.(sname), sgrp) || ...
           ~isfield(histo_gam.(sname), cgrp) || ~isfield(histo_gam.(sname), sgrp)
            fprintf('  %s %s: missing control or severe, skipping diff\n', sname, reg);
            continue
        end

        tbl_ctl = T.(sname).(cgrp);
        tbl_sev = T.(sname).(sgrp);

        % Reuse the stored fitted models (fallback: refit from the table)
        gam_ctl = get_gam_model(histo_gam.(sname).(cgrp), tbl_ctl, ...
                                ntrees, max_splits, learn_rate, n_max);
        gam_sev = get_gam_model(histo_gam.(sname).(sgrp), tbl_sev, ...
                                ntrees, max_splits, learn_rate, n_max);

        % Common-grid difference (severe - control)
        D = gam_common_grid_diff(gam_ctl, gam_sev, tbl_ctl, tbl_sev, ...
                ngrid, q_lo, q_hi, do_bootstrap, nbootstrap, ci_alpha, ...
                ntrees, max_splits, learn_rate, n_max);
        if isempty(fieldnames(D)), continue; end

        histo_gam_diff.(sname).(reg) = D;

        % Save the long table and the surface figure
        outstem = fullfile(diffdir, sprintf('%s_%s_severe-control', sname, reg));
        writetable(D.table, [outstem '.csv']);
        plot_diff_surface(D, sprintf('%s %s: severe - control', upper(sname), reg), outstem);
        fprintf('  %s %s: wrote severe-control surface -> %s\n', sname, reg, outstem);
    end
end

% Save the difference struct
diff_fout = fullfile(diffdir, strcat('GAM_histo_diff_', run_tag, string(dt), '.mat'));
fprintf('Saving histo_gam_diff data struct to %s\n', diff_fout)
save(diff_fout, 'histo_gam_diff', '-v7.3');


%% ------------------------------------------------------------------------
%% Local function: concatenate epv_ves_swp_stain for all sections in a group
function [M, n_sec] = combine_group(stain, patterns)
% Concatenate the [epv_swp, ves_swp, ves_epv_swp, stain] matrices for every
% section whose baseName matches any pattern in `patterns`.
names = {stain.baseName};
idx   = find(contains(lower(names), lower(patterns)));
M = [];
for k = idx
    M = [M; stain(k).epv_ves_swp_stain]; %#ok<AGROW>
end
n_sec = numel(idx);
end


%% Local function: severe - control on a common prediction grid
function D = gam_common_grid_diff(gam_ctl, gam_sev, tbl_ctl, tbl_sev, ...
        ngrid, q_lo, q_hi, do_bootstrap, nboot, ci_alpha, ...
        ntrees, max_splits, learn_rate, max_sample)

D = struct();

epv_c = tbl_ctl.epv_swp; ves_c = tbl_ctl.ves_swp;
epv_s = tbl_sev.epv_swp; ves_s = tbl_sev.ves_swp;

% Per-axis common range = intersection of the two groups' [q_lo, q_hi] ranges
epv_lo = max(quantile(epv_c, q_lo), quantile(epv_s, q_lo));
epv_hi = min(quantile(epv_c, q_hi), quantile(epv_s, q_hi));
ves_lo = max(quantile(ves_c, q_lo), quantile(ves_s, q_lo));
ves_hi = min(quantile(ves_c, q_hi), quantile(ves_s, q_hi));
if epv_hi <= epv_lo || ves_hi <= ves_lo
    warning('gam_common_grid_diff: no overlapping predictor range; skipping.');
    return
end

% Cast to double: quantiles of single-precision SWP columns stay single, but
% convhull (support mask) and downstream math want double.
epv_vec = double(linspace(epv_lo, epv_hi, ngrid));
ves_vec = double(linspace(ves_lo, ves_hi, ngrid));
[EPV, VES] = meshgrid(epv_vec, ves_vec);
VE = EPV .* VES;                       % interaction term (= epv_swp .* ves_swp)

% Predict as a numeric matrix so predictor-name mismatches can't bite. Column
% order must match the training tables: [epv_swp, ves_swp, ves_epv_swp].
Xg = [EPV(:), VES(:), VE(:)];
yc = predict(gam_ctl, Xg);
ys = predict(gam_sev, Xg);
dvec = ys - yc;                        % severe - control

% Common support: inside the convex hull of BOTH groups' (epv, ves) clouds
inC = in_hull_2d(epv_c, ves_c, EPV(:), VES(:));
inS = in_hull_2d(epv_s, ves_s, EPV(:), VES(:));
SUP = reshape(inC & inS, size(EPV));

YC   = reshape(yc,   size(EPV));
YS   = reshape(ys,   size(EPV));
DIFF = reshape(dvec, size(EPV));
DIFF(~SUP) = NaN;                      % report only where both groups have data

% Optional bootstrap CIs on the difference
lo = nan(size(EPV)); hi = lo; sd = lo; excl0 = false(size(EPV));
if do_bootstrap
    B = nan(numel(dvec), nboot);
    for b = 1:nboot
        gcb = fit_gam_local(boot_tbl(tbl_ctl), ntrees, max_splits, learn_rate, max_sample);
        gsb = fit_gam_local(boot_tbl(tbl_sev), ntrees, max_splits, learn_rate, max_sample);
        B(:, b) = predict(gsb, Xg) - predict(gcb, Xg);
    end
    a  = ci_alpha / 2;
    lo = reshape(quantile(B, a,     2), size(EPV));
    hi = reshape(quantile(B, 1 - a, 2), size(EPV));
    sd = reshape(std(B, 0, 2, 'omitnan'), size(EPV));
    lo(~SUP) = NaN; hi(~SUP) = NaN; sd(~SUP) = NaN;
    excl0 = (lo > 0 | hi < 0) & SUP;   % CI excludes zero (in-support only)
end

% Long-format table of the grid
tab = table(EPV(:), VES(:), VE(:), YC(:), YS(:), DIFF(:), SUP(:), ...
    'VariableNames', {'epv_swp','ves_swp','ves_epv_swp', ...
                      'yhat_control','yhat_severe', ...
                      'diff_severe_minus_control','in_support'});
if do_bootstrap
    tab.diff_lo           = lo(:);
    tab.diff_hi           = hi(:);
    tab.diff_sd           = sd(:);
    tab.ci_excludes_zero  = excl0(:);
end

D = struct('epv_vec', epv_vec, 'ves_vec', ves_vec, 'EPV', EPV, 'VES', VES, ...
           'yhat_control', YC, 'yhat_severe', YS, 'diff', DIFF, 'support', SUP, ...
           'diff_lo', lo, 'diff_hi', hi, 'diff_sd', sd, ...
           'ci_excludes_zero', excl0, 'do_bootstrap', do_bootstrap, 'table', tab);
end


%% Local function: get the stored RegressionGAM, or refit from the table
function mdl = get_gam_model(stored, tbl, ntrees, max_splits, learn_rate, max_sample)
mdl = [];
if isa(stored, 'RegressionGAM')
    mdl = stored; return;
elseif isstruct(stored)
    fn = fieldnames(stored);
    for i = 1:numel(fn)
        if isa(stored.(fn{i}), 'RegressionGAM'); mdl = stored.(fn{i}); return; end
        if isstruct(stored.(fn{i}))                     % look one level deeper
            fn2 = fieldnames(stored.(fn{i}));
            for m = 1:numel(fn2)
                if isa(stored.(fn{i}).(fn2{m}), 'RegressionGAM')
                    mdl = stored.(fn{i}).(fn2{m}); return;
                end
            end
        end
    end
end
warning('get_gam_model: no RegressionGAM found in stored result; refitting from table.');
mdl = fit_gam_local(tbl, ntrees, max_splits, learn_rate, max_sample);
end


%% Local function: fit a fitrgam from a predictor table
function mdl = fit_gam_local(tbl, ntrees, max_splits, learn_rate, max_sample)
% NOTE: match these fitrgam name-value pairs to whatever fit_swp_histo_gam uses
% internally, so bootstrap refits are drawn from the same model family as the
% point-estimate fits. Only needed when do_bootstrap = true or when the stored
% object does not expose a RegressionGAM.
n = height(tbl);
if n > max_sample
    tbl = tbl(randsample(n, max_sample), :);
end
X = double([tbl.epv_swp, tbl.ves_swp, tbl.ves_epv_swp]);   % order matches training
y = double(tbl.stain);
mdl = fitrgam(X, y, ...
    'NumTreesPerPredictor',          ntrees, ...
    'MaxNumSplitsPerPredictor',      max_splits, ...
    'InitialLearnRateForPredictors', learn_rate);
end


%% Local function: bootstrap resample a table (rows, with replacement)
function tb = boot_tbl(tb0)
n  = height(tb0);
tb = tb0(randi(n, n, 1), :);
end


%% Local function: is (xq, yq) inside the convex hull of (x, y)?
function tf = in_hull_2d(x, y, xq, yq)
% convhull requires double; SWP columns are stored single, so cast here.
x  = double(x(:));  y  = double(y(:));
xq = double(xq(:)); yq = double(yq(:));
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) > 50000                    % subsample for speed; extremes dominate hull
    s = randsample(numel(x), 50000);
    x = x(s); y = y(s);
end
k  = convhull(x, y);
tf = inpolygon(xq, yq, x(k), y(k));
end


%% Local function: diverging blue-white-red colormap (self-contained)
function cmap = bluewhitered_local(n)
if nargin < 1, n = 256; end
b = [0.23 0.30 0.75]; w = [1 1 1]; r = [0.75 0.15 0.15];
half = floor(n/2);
top  = [linspace(b(1),w(1),half)', linspace(b(2),w(2),half)', linspace(b(3),w(3),half)'];
bot  = [linspace(w(1),r(1),n-half)', linspace(w(2),r(2),n-half)', linspace(w(3),r(3),n-half)'];
cmap = [top; bot];
end


%% Local function: plot and save the severe - control surface
function plot_diff_surface(D, tstr, outstem)
if ~isfield(D, 'diff') || isempty(D.diff), return; end
figure('Position', [100,100,900,750]);
imagesc(D.epv_vec, D.ves_vec, D.diff, 'AlphaData', ~isnan(D.diff));
set(gca, 'YDir', 'normal'); axis tight; hold on;

% Diverging color scale centred on zero
m = max(abs(D.diff(:)), [], 'omitnan');
if isempty(m) || ~isfinite(m) || m == 0, m = 1; end
caxis([-m m]);
colormap(bluewhitered_local());
cb = colorbar; cb.Label.String = 'severe - control (z-score)';

% Outline where the bootstrap CI excludes zero
if isfield(D, 'do_bootstrap') && D.do_bootstrap
    [~, hc] = contour(D.epv_vec, D.ves_vec, double(D.ci_excludes_zero), [0.5 0.5]);
    set(hc, 'LineColor', 'k', 'LineWidth', 1.2);
end

xlabel('epv\_swp'); ylabel('ves\_swp'); title(tstr);
set(gca, 'fontname', 'Arial');
saveas(gcf, outstem, 'png');
saveas(gcf, outstem, 'fig');
exportgraphics(gcf, strcat(outstem, '.svg'), 'Resolution', 600);
end