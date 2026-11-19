%% Between-group (severe vs control) violins of histology stain intensity
%% in EPVS-SWP x vessel-SWP quantile cells, with subject-pair effect sizes
% Histology counterpart of the OCT between-group violin script
% (between_group_violins.m). Same cells, edge modes, pooled effects,
% subject-pair matrix, cell verdicts, group tests and pair-summary
% figures; the measured quantity is the z-scored stain intensity (CD68,
% GFAP, LHE) instead of mus / retardance.
%
% Per section (as in the histology GAM and within-subject violin scripts):
%   tissue mask (stain(j).mask) with vessels and EPVS removed
%   N x 3 matrix [stain_z, ves_swp, epv_swp]
% Sections are concatenated per subject x region (frontal = _1,
% occipital = _7) using the GAM script's baseName patterns. The matrices
% are cached as histo_op (shared with plot_histo_swp_violins.m).
%
% Interpretation notes specific to histology:
%   - Stains are z-scored within each section, so between-group level
%     differences (pooled delta, cell LME beta_group, Welch difference)
%     compare positions relative to each section's own mean, not absolute
%     stain intensity. The contrast models (beta_group:inA) are the
%     comparison the normalisation supports.
%   - Histology SWP uses a 4 mm radius (ds4) vs 500 um for OCT.
%   - No CD68 section exists for CAA6 frontal: control frontal CD68 is
%     CAA26 alone (no within-control pair; LME group terms have 1 df;
%     no Welch test, which needs >= 2 subjects per group).
%   - CD68 is near-binary, so its violins show quantisation steps and its
%     effect sizes mostly reflect the fraction of positive voxels.
%
% Y-axes are shared across regions per stain: frontal and occipital
% violins of one stain use the same y-range (min / max over subjects and
% both regions of each subject's 0.5-99.5th percentiles). YLIM_LOOKUP
% overrides this per stain (both regions) or per stain x region.
%
% Subject-median markers are kept inside their violin (offset <= 75% of the
% local half-width); 0 < n < MIN_N: hollow marker; n = 0: labelled.
%
% Group tests (subjects as the unit of replication):
%   Linear mixed models on voxels (fitlme, REML). df: Satterthwaite; if
%   fixedEffects returns NaN Satterthwaite df (its numerical Hessian of the
%   covariance parameters fails, e.g. at very large N or a near-zero
%   variance estimate), between-within df = n_subjects - 2 are used for
%   every term (recorded in df_method). Each subject contributes at most
%   LME_MAX_N voxels per cell (Inf = all of its voxels in the cell, after
%   EQUAL_WEIGHT_SUBJECTS subsampling). With a subject random intercept
%   the group SE is set by between-subject variance, so the estimate, df
%   and p are essentially insensitive to this cap (unlike the KS p).
%   Cell model     : y ~ group + (1 | subject), per region x stain x cell
%                    beta_group = severe - control mean (z units).
%                    df_collapsed flags a fit whose subject variance went
%                    to ~0: its p is pseudoreplicated.
%   Contrast model : y ~ group * inA + (1 + inA | subject), diagonal
%                    covariance, cells A and B only (inA = 1 for cell A).
%                    beta_inA          = A - B contrast in control
%                    beta_group:inA    = severe - control difference in
%                                        that within-subject contrast
%                    If the random-slope SD collapses to ~0 (flag
%                    slope_sd_singular) the interaction p is
%                    pseudoreplicated and should not be reported.
%   Welch t-test   : ttest2(..., 'Vartype','unequal') on subject-level
%                    summaries (median and mean of all of each subject's
%                    voxels in the cell), per region x stain x cell.
%                    Robustness check on the LME; the median version is
%                    insensitive to the tails of the voxel distributions.
%   A voxel-level ttest2 is not used: it treats voxels as independent
%   replicates of the group difference (pseudoreplication; p ~ 0).
% With 2 subjects per group the group terms have ~2 df (Welch: 1-2 df);
% the exact subject-label permutation floor is p = 1/3 (two-sided).
%
% Outputs in a per-run folder <out_root>/<tag>_<HHmmss>/ (violins/,
% pair_summary/, csv/, run_params.txt). The folder carries the full
% descriptive tag; file names carry a compact tag so every full path stays
% <= MAX_PATH_LEN characters (Illustrator):
%   tag  = zscore_<qLL-HH>_<edgeMode>Edges[_eqW]_<EXPL|MS>_<yyyymmdd>
%   ftag = zscore_<qLL-HH>_<wsE|shE>[_eqW]_<EXPL|MS>_<yyyymmdd>
%   violin_<stain>_<region>_<ftag>.{png,svg,fig}
%   pair_{delta,ksD}_<ftag>.{png,svg,fig}   (rows = stains)
%   {bin_summary,pooled_effects,pairwise,cell_verdict,
%    lme_cell,lme_contrast,welch_cell}_<ftag>.csv
% Set PAIR_SUMMARY_ONLY = true and PAIR_CSV to redraw the pair summaries
% from an existing pairwise_*.csv (or older between_pairwise_*.csv)
% without loading any data.

%% Initialization
clearvars -except histo_op caa swp_op
clc; close all;
addpath(genpath('/projectnb/npbssmic/s/mhyman/oct-caa'));
maxNumCompThreads(feature('numcores'));
rng(1);

%% Parameters
EXPLORATORY           = false;           % EXPL = smaller KDE / effect-size subsamples
FORCE_REBUILD         = false;           % true = rebuild histo_op even if in workspace
Q_EDGES               = [0.25 0.75];     % low: x <= q(1); high: x >= q(2)
EDGE_MODE             = 'withinSubject'; % 'shared' | 'withinSubject'
EQUAL_WEIGHT_SUBJECTS = true;            % subsample subjects to equal n per stain x region
MIN_N                 = 50;              % min voxels to draw a violin / compute effects
FIG_VISIBLE           = 'off';
MAX_PATH_LEN          = 250;             % max full-path length of any output file

%%% Optional per-stain y-limit lookup (overrides the shared automatic limits):
%   YLIM_LOOKUP.gfap       = [-2 4];     % both regions (keeps them matched)
%   YLIM_LOOKUP.gfap.front = [-2 4];     % one region only
YLIM_LOOKUP           = struct();

%%% Pair summary figures
SPARSE_FRAC           = 0.01;
PAIR_SUMMARY_ONLY     = false;
PAIR_CSV              = '';

%%% Subject display order and markers
SUBS_ORDER  = {'CAA26','CAA6','CAA22','CAA25'};
SUB_MARKERS = {'o','s','^','d'};

if ~any(strcmp(EDGE_MODE, {'shared','withinSubject'}))
    error(['EDGE_MODE must be ''shared'' or ''withinSubject''; group-specific ' ...
           'edges make cells non-comparable across groups.']);
end
if EXPLORATORY
    run_tag = 'EXPL'; KDE_MAX_N = 2e4; CLIFF_MAX_N = 5e5; LME_MAX_N = 5e3;
else
    run_tag = 'MS';   KDE_MAX_N = 1e5; CLIFF_MAX_N = 5e6; LME_MAX_N = Inf;
end
% LME_MAX_N: voxels per subject per cell entering each mixed model (Inf = all)
date_str  = char(datetime('now','Format','yyyyMMdd'));
q_tag     = sprintf('q%02d-%02d', round(100*Q_EDGES));
if EQUAL_WEIGHT_SUBJECTS, eqw_tag = '_eqW'; else, eqw_tag = ''; end
VALUE_TYPE = 'zscore';
edge_abbr = struct('shared','shE', 'withinSubject','wsE');
tag  = sprintf('%s_%s_%sEdges%s_%s_%s', VALUE_TYPE, q_tag, EDGE_MODE, eqw_tag, run_tag, date_str);
ftag = sprintf('%s_%s_%s%s_%s_%s', VALUE_TYPE, q_tag, edge_abbr.(EDGE_MODE), eqw_tag, run_tag, date_str);

%%% Inputs (as in the histology GAM script)
hdir = '/projectnb/npbssmic/ns/CAA/histology/swp';
stain_files = struct( ...
    'cd68', 'histo_swp_hmatched_zscore__CD68_voxelwise_radius4mm_ds4_20-Mar-2026.mat', ...
    'gfap', 'histo_swp_hmatched_zscore__GFAP_voxelwise_radius4mm_ds4_20-Mar-2026.mat', ...
    'lhe',  'histo_swp_hmatched_zscore__LHE_voxelwise_radius4mm_ds4_19-Mar-2026.mat');
stain_names = fieldnames(stain_files);

%%% Output folders: one subfolder per run
out_root = '/projectnb/npbssmic/ns/CAA/histology/figures/histo_swp_between_group/';
time_str = char(datetime('now','Format','HHmmss'));

if PAIR_SUMMARY_ONLY
    csv_dir  = fileparts(PAIR_CSV);
    [parent, leaf] = fileparts(csv_dir);
    if strcmp(leaf, 'csv'), run_dir = parent; else, run_dir = csv_dir; end
    sum_dir  = fullfile(run_dir, sprintf('pair_summary_redraw_%s_%s', date_str, time_str));
    mkdir(sum_dir);
else
    run_dir = fullfile(out_root, sprintf('%s_%s', tag, time_str));
    vio_dir = fullfile(run_dir, 'violins');
    sum_dir = fullfile(run_dir, 'pair_summary');
    csv_dir = fullfile(run_dir, 'csv');

    %%% Path-length check before any computation (longest file of each kind)
    lstn    = stain_names{find(cellfun(@numel, stain_names) == max(cellfun(@numel, stain_names)), 1)};
    longest = {fullfile(vio_dir, sprintf('violin_%s_occip_%s.png', lstn, ftag)), ...
               fullfile(sum_dir, sprintf('pair_delta_%s.png', ftag)), ...
               fullfile(csv_dir, sprintf('pooled_effects_%s.csv', ftag))};
    for i = 1:numel(longest), check_path(longest{i}, MAX_PATH_LEN); end

    for d = {run_dir, vio_dir, sum_dir, csv_dir}
        if ~exist(d{1},'dir'), mkdir(d{1}); end
    end
    fid = fopen(fullfile(run_dir, 'run_params.txt'), 'w');
    fprintf(fid, 'script                : %s\n', mfilename('fullpath'));
    fprintf(fid, 'started               : %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    fprintf(fid, 'tag (folder)          : %s\n', tag);
    fprintf(fid, 'ftag (file names)     : %s\n', ftag);
    fprintf(fid, 'EXPLORATORY           : %d (%s)\n', EXPLORATORY, run_tag);
    fprintf(fid, 'stains                : %s\n', strjoin(stain_names', ', '));
    for si = 1:numel(stain_names)
        fprintf(fid, '  %-5s file           : %s\n', stain_names{si}, stain_files.(stain_names{si}));
    end
    fprintf(fid, 'Q_EDGES               : [%g %g]\n', Q_EDGES);
    fprintf(fid, 'EDGE_MODE             : %s\n', EDGE_MODE);
    fprintf(fid, 'EQUAL_WEIGHT_SUBJECTS : %d\n', EQUAL_WEIGHT_SUBJECTS);
    fprintf(fid, 'MIN_N                 : %d\n', MIN_N);
    fprintf(fid, 'KDE_MAX_N             : %g\n', KDE_MAX_N);
    fprintf(fid, 'CLIFF_MAX_N           : %g\n', CLIFF_MAX_N);
    fprintf(fid, 'LME_MAX_N             : %g\n', LME_MAX_N);
    fprintf(fid, 'LME                   : fitlme, REML, Satterthwaite df (between-within fallback)\n');
    fprintf(fid, 'Welch                 : ttest2 unequal var, subject medians and means\n');
    fprintf(fid, 'y-axes                : shared across regions per stain (YLIM_LOOKUP overrides)\n');
    fprintf(fid, 'SPARSE_FRAC           : %g\n', SPARSE_FRAC);
    fprintf(fid, 'MAX_PATH_LEN          : %d\n', MAX_PATH_LEN);
    fprintf(fid, 'rng seed              : 1\n');
    fclose(fid);
end
fprintf('Run folder: %s\n', run_dir);

scfg = struct('out_dir',sum_dir, 'fig_visible',FIG_VISIBLE, ...
              'sparse_frac',SPARSE_FRAC, 'subs_order',{SUBS_ORDER}, ...
              'props',{stain_names'}, 'prop_disp',{upper(stain_names')}, ...
              'max_path',MAX_PATH_LEN);

%% Pair summary only: redraw from an existing CSV (no data loading)
if PAIR_SUMMARY_ONLY
    T_pw = readtable(PAIR_CSV, 'TextType','char');
    [~, nm] = fileparts(PAIR_CSV);
    plot_pair_summary(T_pw, 'delta', regexprep(nm, '^(between_)?pairwise_', 'pair_delta_'), scfg);
    plot_pair_summary(T_pw, 'ks',    regexprep(nm, '^(between_)?pairwise_', 'pair_ksD_'),   scfg);
    fprintf('Pair summaries written to %s\n', sum_dir);
    return
end

%%% Subjects (baseName patterns from the GAM script) and groups
subjects = struct( ...
    'caa26', struct('front',{{'CAA26_1','CAA_26_1'}}, 'occip',{{'CAA26_7','CAA_26_7'}}), ...
    'caa6',  struct('front',{{'CAA6_1', 'CAA_6_1' }}, 'occip',{{'CAA6_7', 'CAA_6_7' }}), ...
    'caa22', struct('front',{{'CAA22_1','CAA_22_1'}}, 'occip',{{'CAA22_7','CAA_22_7'}}), ...
    'caa25', struct('front',{{'CAA25_1','CAA_25_1'}}, 'occip',{{'CAA25_7','CAA_25_7'}}), ...
    'caa17', struct('front',{{'CAA17_1','CAA_17_1'}}, 'occip',{{'CAA17_7','CAA_17_7'}}));
groups.control = {'caa26','caa6'};
groups.severe  = {'caa22','caa25'};
grp_order = {'control','severe'};
regions   = {'front','occip'};
reg_disp  = struct('front','Frontal','occip','Occipital');

%% Build N x 3 matrices [stain_z, ves_swp, epv_swp] (skipped if cached)
if ~exist('histo_op','var') || FORCE_REBUILD
    histo_op = struct();
    subs = fieldnames(subjects);
    for si = 1:numel(stain_names)
        stn = stain_names{si};
        fprintf('Loading %s\n', upper(stn));
        tmp = load(fullfile(hdir, stain_files.(stn)));
        Sx  = tmp.(stn); clear tmp
        names = lower({Sx.baseName});
        for s = 1:numel(subs)
            sub = subs{s};
            for r = 1:numel(regions)
                reg = regions{r};
                idx = find(contains(names, lower(subjects.(sub).(reg))));
                if isempty(idx), continue; end
                M = [];
                for k = idx
                    mask = Sx(k).mask;
                    mask(Sx(k).ves)  = 0;
                    mask(Sx(k).epvs) = 0;
                    M = [M; double(Sx(k).z_stain(mask)), ...
                            double(Sx(k).ves_swp(mask)), ...
                            double(Sx(k).epvs_swp(mask))]; %#ok<AGROW>
                end
                histo_op.(stn).(sub).(reg) = M;
                fprintf('\t%s %s: %d section(s), %d voxels\n', upper(sub), reg, numel(idx), size(M,1));
            end
        end
        clear Sx
    end
else
    fprintf('histo_op found in workspace; skipping import\n');
end

%% Plot configuration (stain-specific fields set inside the stain loop)
cfg = struct();
cfg.prop_col    = 1;            % stain_z
cfg.ves_col     = 2;
cfg.epv_col     = 3;
cfg.yscale      = 1;
cfg.units       = 'z-score';
cfg.beta_fmt    = '%+.2f';      % beta / Welch difference label precision
cfg.vt          = VALUE_TYPE;
cfg.Q           = Q_EDGES;
cfg.edge_mode   = EDGE_MODE;
cfg.eqw         = EQUAL_WEIGHT_SUBJECTS;
cfg.min_n       = MIN_N;
cfg.kde_max_n   = KDE_MAX_N;
cfg.cliff_max_n = CLIFF_MAX_N;
cfg.lme_max_n   = LME_MAX_N;
cfg.contr = { 'Ves-SWP within EPVS-SWP', 'EPVS-SWP low',  2, 1;
              'Ves-SWP within EPVS-SWP', 'EPVS-SWP high', 4, 3;
              'anti-diagonal',           'none',          3, 2;
              'EPVS-SWP within Ves-SWP', 'Ves-SWP low',   3, 1;
              'EPVS-SWP within Ves-SWP', 'Ves-SWP high',  4, 2;
              'diagonal',                'none',          4, 1 };
cfg.xc          = [1 2.2 3.9 5.1];
cfg.dx          = 0.22;
cfg.hw          = 0.19;
cfg.med_off     = 0.07;
cfg.med_inside  = 0.75;
cfg.gcol        = [0.25 0.25 0.25; 0.80 0.15 0.45];
cfg.glab        = {'Control','Severe'};
cfg.vbin        = [1 2 1 2];
cfg.vlab        = {'Ves-SWP low','Ves-SWP high'};
cfg.ebin_lab    = {'low','low','high','high'};
cfg.vbin_lab    = {'low','high','low','high'};
cfg.cell_lab    = {'EPVS-lo/Ves-lo','EPVS-lo/Ves-hi','EPVS-hi/Ves-lo','EPVS-hi/Ves-hi'};
cfg.subs_order  = SUBS_ORDER;
cfg.markers     = SUB_MARKERS;
cfg.fig_visible = FIG_VISIBLE;

all_bin  = [];
all_pool = [];
all_pw   = [];
all_lme  = [];
all_lmec = [];
all_wel  = [];
all_subs = [groups.control, groups.severe];

%% Stain loop
for si = 1:numel(stain_names)
stn = stain_names{si};
if ~isfield(histo_op, stn), continue; end
D = histo_op.(stn);
cfg.props     = {stn};
cfg.prop_disp = {upper(stn)};
cfg.ylab      = {sprintf('%s (z-score)', upper(stn))};
fprintf('=== %s ===\n', upper(stn));

%%% Shared y-limits for this stain, pooled over regions (full data,
% before any equal-weight subsampling)
ylim_stain = [inf -inf];
for r = 1:numel(regions)
    for s = 1:numel(all_subs)
        sb = all_subs{s};
        if ~(isfield(D, sb) && isfield(D.(sb), regions{r})), continue; end
        x = D.(sb).(regions{r})(:, cfg.prop_col);
        q = cfg.yscale * prctile(x(isfinite(x)), [0.5 99.5]);
        ylim_stain = [min(ylim_stain(1), q(1)), max(ylim_stain(2), q(2))];
    end
end
clear x
fprintf('y-limits %s (both regions): [%.4g, %.4g]\n', stn, ylim_stain);

%% Region loop
for r = 1:numel(regions)
    reg = regions{r};

    %%% Subjects in this region (control first, then severe)
    S = struct('name',{}, 'group',{}, 'g',{}, 'M',{}, 'C',{}, 'eE',{}, 'eV',{});
    for g = 1:numel(grp_order)
        gs = groups.(grp_order{g});
        for s = 1:numel(gs)
            if isfield(D, gs{s}) && isfield(D.(gs{s}), reg)
                S(end+1) = struct('name',gs{s}, 'group',grp_order{g}, 'g',g, ...
                                  'M',clean_rows(D.(gs{s}).(reg)), 'C',[], ...
                                  'eE',[NaN NaN], 'eV',[NaN NaN]); %#ok<SAGROW>
            end
        end
    end
    if ~all(ismember(1:2, [S.g]))
        warning('between:groups', '%s: both groups not present; skipped.', reg);
        continue
    end
    fprintf('Region %s: %s\n', reg, strjoin(upper({S.name}), ', '));

    %%% Optional equal weighting: subsample all subjects to the smallest
    if EQUAL_WEIGHT_SUBJECTS
        n_min_eq = min(arrayfun(@(x) size(x.M,1), S));
        for s = 1:numel(S)
            n = size(S(s).M,1);
            S(s).M = S(s).M(sort(randperm(n, n_min_eq)),:);
        end
    end

    %%% Edges and cell masks
    switch EDGE_MODE
        case 'shared'
            allSWP = cell2mat(arrayfun(@(x) x.M(:,[cfg.ves_col cfg.epv_col]), S(:), 'UniformOutput', false));
            eV = swp_edges(allSWP(:,1), Q_EDGES, sprintf('shared %s Ves-SWP',  reg));
            eE = swp_edges(allSWP(:,2), Q_EDGES, sprintf('shared %s EPVS-SWP', reg));
            clear allSWP
            for s = 1:numel(S), S(s).eE = eE; S(s).eV = eV; end
            edge_title = {edge_line(eE, eV)};
        case 'withinSubject'
            eE = NaN(1,2); eV = NaN(1,2);
            tE = cell(1, numel(S)); tV = tE;
            for s = 1:numel(S)
                S(s).eE = swp_edges(S(s).M(:,cfg.epv_col), Q_EDGES, sprintf('%s %s EPVS-SWP', S(s).name, reg));
                S(s).eV = swp_edges(S(s).M(:,cfg.ves_col), Q_EDGES, sprintf('%s %s Ves-SWP',  S(s).name, reg));
                tE{s} = sprintf('%s \\leq %.3g / \\geq %.3g', upper(S(s).name), S(s).eE(1), S(s).eE(2));
                tV{s} = sprintf('%s \\leq %.3g / \\geq %.3g', upper(S(s).name), S(s).eV(1), S(s).eV(2));
            end
            iC = [S.g] == 1;  iS = [S.g] == 2;
            edge_title = { ['EPVS-SWP edges, control:  ' strjoin(tE(iC), ';    ')], ...
                           ['EPVS-SWP edges, severe:  '  strjoin(tE(iS), ';    ')], ...
                           ['Ves-SWP edges, control:  '  strjoin(tV(iC), ';    ')], ...
                           ['Ves-SWP edges, severe:  '   strjoin(tV(iS), ';    ')] };
    end
    for s = 1:numel(S)
        S(s).C = cell_masks(S(s).M(:,cfg.epv_col), S(s).M(:,cfg.ves_col), S(s).eE, S(s).eV);
    end

    %%% Y-limits: shared across regions for this stain, unless looked up
    ylims = ylim_stain;
    if isfield(YLIM_LOOKUP, stn)
        ylk = YLIM_LOOKUP.(stn);
        if isnumeric(ylk)
            ylims = ylk;
        elseif isstruct(ylk) && isfield(ylk, reg)
            ylims = ylk.(reg);
        end
    end

    ctl_names = strjoin(upper({S([S.g] == 1).name}), ', ');
    sev_names = strjoin(upper({S([S.g] == 2).name}), ', ');

    for p = 1:numel(cfg.props)
        col  = cfg.prop_col(p);
        prop = cfg.props{p};
        base = sprintf('violin_%s_%s_%s', prop, reg, ftag);
        fprintf('\t%s\n', prop);

        %%% Subject-level cells: summary rows and median overlays
        ov = struct('name',{}, 'g',{}, 'med',{}, 'n',{}, 'marker',{}, 'off',{});
        for g = 1:2
            idx = find([S.g] == g);
            if numel(idx) == 1, offs = 0; else, offs = linspace(-cfg.med_off, cfg.med_off, numel(idx)); end
            for j = 1:numel(idx)
                s    = idx(j);
                ntot = size(S(s).M, 1);
                med  = NaN(1,4); nk = zeros(1,4);
                for k = 1:4
                    st = cell_stats(cfg.yscale * S(s).M(S(s).C(:,k), col));
                    nk(k) = st.n;
                    if st.n > 0, med(k) = st.med; end      % low-n medians drawn hollow
                    if st.n < cfg.min_n
                        fprintf('\t\t%s %s: n = %d (< MIN_N = %d)\n', ...
                                upper(S(s).name), cfg.cell_lab{k}, st.n, cfg.min_n);
                    end
                    all_bin = app(all_bin, bin_row('subject', upper(S(s).name), S(s).group, ...
                                  reg, prop, k, S(s).eE, S(s).eV, st, ntot, base, cfg));
                end
                mi = find(strcmpi(cfg.subs_order, S(s).name), 1);
                if isempty(mi), mk = 'p'; else, mk = cfg.markers{mi}; end
                ov(end+1) = struct('name',upper(S(s).name), 'g',g, 'med',med, 'n',nk, ...
                                   'marker',mk, 'off',offs(j)); %#ok<SAGROW>
            end
        end

        %%% Pooled group cells
        Yg = cell(2,4); stg = repmat(cell_stats([]), 2, 4); ntot_g = zeros(1,2);
        for g = 1:2
            idx = find([S.g] == g);
            ntot_g(g) = sum(arrayfun(@(x) size(x.M,1), S(idx)));
            for k = 1:4
                Yg{g,k} = cell2mat(arrayfun(@(x) cfg.yscale * x.M(x.C(:,k), col), ...
                                            S(idx)', 'UniformOutput', false));
                stg(g,k) = cell_stats(Yg{g,k});
                all_bin = app(all_bin, bin_row('group', grp_order{g}, grp_order{g}, ...
                              reg, prop, k, eE, eV, stg(g,k), ntot_g(g), base, cfg));
            end
        end

        %%% Pooled effects (A = severe, B = control)
        dpool = NaN(1,4); Dpool = NaN(1,4); ppool = NaN(1,4);
        for k = 1:4
            na_u = 0; nb_u = 0;
            if stg(2,k).n >= cfg.min_n && stg(1,k).n >= cfg.min_n
                [a, b]   = cap_pair(Yg{2,k}, Yg{1,k}, cfg.cliff_max_n);
                dpool(k) = cliffs_delta(a, b);
                Dpool(k) = ks_stat(a, b);
                ppool(k) = ks_pval(Dpool(k), numel(a), numel(b));
                na_u = numel(a); nb_u = numel(b);
                clear a b
            end
            all_pool = app(all_pool, struct( ...
                'region',reg, 'property',prop, 'value_type',cfg.vt, ...
                'edge_mode',cfg.edge_mode, 'equal_weight',cfg.eqw, ...
                'cell',cfg.cell_lab{k}, 'cell_idx',k, ...
                'n_control',stg(1,k).n, 'n_severe',stg(2,k).n, ...
                'n_control_used',nb_u, 'n_severe_used',na_u, ...
                'frac_control',stg(1,k).n/ntot_g(1), 'frac_severe',stg(2,k).n/ntot_g(2), ...
                'median_control',stg(1,k).med, 'median_severe',stg(2,k).med, ...
                'units',cfg.units, 'cliffs_delta',dpool(k), 'ks_D',Dpool(k), ...
                'ks_p',ppool(k), ...
                'figure',base));
        end

        %%% Subject-pair matrix and per-cell verdicts
        pw = pairwise_cell_effects(S, col, prop, reg, cfg);
        all_pw = app(all_pw, pw);
        pwT = struct2table(pw(:));
        for k = 1:4
            vk = cell_verdict(pwT(pwT.cell_idx == k, :));
            if k == 1, V = vk; else, V(k) = vk; end %#ok<SAGROW>
        end

        %%% Mixed model per cell: y ~ group + (1 | subject), voxels as observations
        for k = 1:4
            meta = struct('region',reg, 'property',prop, 'value_type',cfg.vt, ...
                          'edge_mode',cfg.edge_mode, 'equal_weight',cfg.eqw, ...
                          'cell',cfg.cell_lab{k}, 'cell_idx',k, 'units',cfg.units, ...
                          'figure',base);
            lk = merge_structs(meta, lme_group_cell(S, col, k, cfg));
            if ~isempty(lk.note)
                warning('between:lme', '%s %s %s: LME cell model: %s', ...
                        reg, prop, cfg.cell_lab{k}, lk.note);
            end
            if k == 1, L = lk; else, L(k) = lk; end %#ok<SAGROW>
        end
        all_lme = app(all_lme, L);

        %%% Welch t-test on subject-level summaries (median shown in the figure)
        W = [];
        for k = 1:4
            meta = struct('region',reg, 'property',prop, 'value_type',cfg.vt, ...
                          'edge_mode',cfg.edge_mode, 'equal_weight',cfg.eqw, ...
                          'cell',cfg.cell_lab{k}, 'cell_idx',k, 'units',cfg.units, ...
                          'figure',base);
            for sm = {'median','mean'}
                wk = merge_structs(meta, welch_cell(S, col, k, cfg, sm{1}));
                all_wel = app(all_wel, wk);
                if strcmp(sm{1}, 'median'), W = app(W, wk); end
            end
        end

        %%% Mixed model per within-subject contrast: y ~ group * inA + (1 + inA | subject)
        for c = 1:size(cfg.contr,1)
            lc = lme_group_contrast(S, col, cfg.contr{c,3}, cfg.contr{c,4}, cfg);
            if ~isempty(lc.note) && isnan(lc.beta_interaction)
                warning('between:lme', '%s %s %s vs %s: LME contrast model failed: %s', ...
                        reg, prop, cfg.cell_lab{cfg.contr{c,3}}, cfg.cell_lab{cfg.contr{c,4}}, lc.note);
            end
            meta = struct('region',reg, 'property',prop, 'value_type',cfg.vt, ...
                          'edge_mode',cfg.edge_mode, 'equal_weight',cfg.eqw, ...
                          'contrast_type',cfg.contr{c,1}, 'held_fixed',cfg.contr{c,2}, ...
                          'cell_A',cfg.cell_lab{cfg.contr{c,3}}, ...
                          'cell_B',cfg.cell_lab{cfg.contr{c,4}}, 'units',cfg.units);
            all_lmec = app(all_lmec, merge_structs(meta, lc));
        end

        %%% Figure
        P = struct();
        P.ylims  = ylims;
        P.ylabel = cfg.ylab{p};
        P.title  = [ {sprintf('Control (%s) vs Severe (%s)  |  %s  |  %s', ...
                              ctl_names, sev_names, reg_disp.(reg), cfg.prop_disp{p})}, ...
                     reshape(edge_title, 1, []), ...
                     {sprintf('quantiles %.2f / %.2f, %s edges%s;   \\delta_{Cliff} > 0: severe larger', ...
                              cfg.Q(1), cfg.Q(2), EDGE_MODE, strrep(eqw_tag, '_', '\_')), ...
                      'bold: subject pairs concordant and outside within-group range', ...
                      'D, p: KS statistic and asymptotic p on voxels (pseudoreplicated; not a test of the group difference)', ...
                      'LME: y ~ group + (1|subject) on voxels, \beta = severe - control, Satterthwaite or between-within df (subjects as replicates)', ...
                      'Welch: t-test on subject medians, \Delta = severe - control (1-2 df)', ...
                      'z-scored within section: level differences are relative to each section''s mean'} ];
        fh = plot_between_violins(Yg, stg, dpool, Dpool, ppool, V, L, W, ov, P, cfg);
        save_fig(fh, vio_dir, base, MAX_PATH_LEN);
        clear Yg
    end
    clear S
end

end   % stain loop

%% Export CSVs
T_bin  = struct2table(all_bin(:));
T_pool = struct2table(all_pool(:));
T_pw   = struct2table(all_pw(:));

%%% Cell verdicts
keys  = unique(T_pw(:, {'region','property','cell_idx'}), 'rows');
vrows = [];
for i = 1:height(keys)
    Tk = T_pw(strcmp(T_pw.region, keys.region{i}) & strcmp(T_pw.property, keys.property{i}) ...
              & T_pw.cell_idx == keys.cell_idx(i), :);
    v   = cell_verdict(Tk);
    row = struct('region',keys.region{i}, 'property',keys.property{i}, ...
                 'value_type',VALUE_TYPE, 'edge_mode',EDGE_MODE, ...
                 'equal_weight',EQUAL_WEIGHT_SUBJECTS, ...
                 'cell',Tk.cell{1}, 'cell_idx',keys.cell_idx(i));
    fn = fieldnames(v);
    for j = 1:numel(fn), row.(fn{j}) = v.(fn{j}); end
    vrows = app(vrows, row);
end
T_ver = struct2table(vrows(:));

write_csv(T_bin,  fullfile(csv_dir, sprintf('bin_summary_%s.csv',    ftag)), MAX_PATH_LEN);
write_csv(T_pool, fullfile(csv_dir, sprintf('pooled_effects_%s.csv', ftag)), MAX_PATH_LEN);
write_csv(T_pw,   fullfile(csv_dir, sprintf('pairwise_%s.csv',       ftag)), MAX_PATH_LEN);
write_csv(T_ver,  fullfile(csv_dir, sprintf('cell_verdict_%s.csv',   ftag)), MAX_PATH_LEN);
write_csv(struct2table(all_lme(:)),  fullfile(csv_dir, sprintf('lme_cell_%s.csv',     ftag)), MAX_PATH_LEN);
write_csv(struct2table(all_lmec(:)), fullfile(csv_dir, sprintf('lme_contrast_%s.csv', ftag)), MAX_PATH_LEN);
write_csv(struct2table(all_wel(:)),  fullfile(csv_dir, sprintf('welch_cell_%s.csv',   ftag)), MAX_PATH_LEN);

%% Pair summary figures
plot_pair_summary(T_pw, 'delta', sprintf('pair_delta_%s', ftag), scfg);
plot_pair_summary(T_pw, 'ks',    sprintf('pair_ksD_%s',   ftag), scfg);
fprintf('Done. Outputs in %s\n', run_dir);


%% ======================= Local functions ===============================

function rows = pairwise_cell_effects(S, col, prop, reg, cfg)
% Cliff's delta and KS D for every subject pair, per cell.
% Between-group pairs: A = severe, B = control (delta > 0: severe larger).
% Within-group pairs: A/B in S order; compare on |delta|.
rows = [];
P    = nchoosek(1:numel(S), 2);
ntot = arrayfun(@(x) size(x.M,1), S);
for k = 1:4
    for i = 1:size(P,1)
        a = P(i,1); b = P(i,2);
        if ~strcmp(S(a).group, S(b).group)
            ptype = 'between';
            if strcmp(S(a).group, 'control'), [a, b] = deal(b, a); end
        else
            ptype = ['within_' S(a).group];
        end
        ya = cfg.yscale * S(a).M(S(a).C(:,k), col);
        yb = cfg.yscale * S(b).M(S(b).C(:,k), col);
        d = NaN; Dks = NaN; pks = NaN; na_u = 0; nb_u = 0;
        if numel(ya) >= cfg.min_n && numel(yb) >= cfg.min_n
            [ya_s, yb_s] = cap_pair(ya, yb, cfg.cliff_max_n);
            d    = cliffs_delta(ya_s, yb_s);
            Dks  = ks_stat(ya_s, yb_s);
            pks  = ks_pval(Dks, numel(ya_s), numel(yb_s));
            na_u = numel(ya_s); nb_u = numel(yb_s);
        end
        rows = app(rows, struct( ...
            'region',reg, 'property',prop, 'value_type',cfg.vt, ...
            'edge_mode',cfg.edge_mode, 'equal_weight',cfg.eqw, ...
            'cell',cfg.cell_lab{k}, 'cell_idx',k, 'pair_type',ptype, ...
            'subj_A',upper(S(a).name), 'group_A',S(a).group, ...
            'subj_B',upper(S(b).name), 'group_B',S(b).group, ...
            'cliffs_delta',d, 'abs_delta',abs(d), 'ks_D',Dks, 'ks_p',pks, ...
            'n_A',numel(ya), 'n_B',numel(yb), 'n_A_used',na_u, 'n_B_used',nb_u, ...
            'frac_A',numel(ya)/ntot(a), 'frac_B',numel(yb)/ntot(b)));
    end
end
end


function [y, inA, gc, sc, used] = stack_voxels(S, col, ka, kb, cfg)
% Voxels of cell ka (inA = 1) and, if kb is non-empty, cell kb (inA = 0)
% for every subject with >= cfg.min_n voxels in each, at most
% cfg.lme_max_n per subject per cell (Inf = all). Numeric group / subject
% codes; no cell arrays, so ~1e7-1e8 rows stay affordable in memory.
ks   = [ka kb];
used = false(1, numel(S));
nv   = zeros(numel(S), numel(ks));
for s = 1:numel(S)
    nk      = sum(S(s).C(:,ks), 1);
    used(s) = all(nk >= cfg.min_n);
    nv(s,:) = min(nk, cfg.lme_max_n);
end
N   = sum(sum(nv(used,:)));
y   = zeros(N,1); inA = zeros(N,1); gc = zeros(N,1); sc = zeros(N,1);
i0  = 0;
for s = find(used)
    for j = 1:numel(ks)
        v  = cap_one(cfg.yscale * S(s).M(S(s).C(:,ks(j)), col), cfg.lme_max_n);
        ii = i0 + (1:numel(v))';
        y(ii) = v;  inA(ii) = (j == 1);  gc(ii) = S(s).g;  sc(ii) = s;
        i0 = i0 + numel(v);
    end
end
end


function r = lme_group_cell(S, col, k, cfg)
% y ~ group + (1 | subject) for cell k, voxels as observations.
% beta_group = severe - control (units of y), Satterthwaite df.
r = struct('n_subjects',0, 'n_control',0, 'n_severe',0, 'n_voxels',0, ...
           'mean_control',NaN, 'beta_group',NaN, 'se_group',NaN, 't_group',NaN, ...
           'df_group',NaN, 'p_group',NaN, 'ci_low',NaN, 'ci_high',NaN, ...
           'sd_subject',NaN, 'sd_resid',NaN, 'df_method','', 'df_collapsed',false, 'note','');
[y, ~, gc, sc, used] = stack_voxels(S, col, k, [], cfg);
ng = [sum(used & [S.g] == 1), sum(used & [S.g] == 2)];
r.n_subjects = sum(used); r.n_control = ng(1); r.n_severe = ng(2); r.n_voxels = numel(y);
if any(ng < 1) || sum(ng) < 3
    r.note = 'insufficient subjects'; return
end
tbl = table(y, categorical(gc, [1 2], {'control','severe'}), categorical(sc), ...
            'VariableNames', {'y','group','subject'});
clear y gc sc
try
    lme = fitlme(tbl, 'y ~ group + (1|subject)', 'FitMethod','REML');
    [st, r.df_method] = fe_table(lme, r.n_subjects - 2);
    nm = cellstr(st.Name);
    i  = find(strcmp(nm, 'group_severe'), 1);
    i0 = find(strcmp(nm, '(Intercept)'), 1);
    r.mean_control = st.Estimate(i0);
    r.beta_group   = st.Estimate(i);
    r.se_group     = st.SE(i);
    r.t_group      = st.tStat(i);
    r.df_group     = st.DF(i);
    r.p_group      = st.pValue(i);
    r.ci_low       = st.Lower(i);
    r.ci_high      = st.Upper(i);
    [psi, mse]     = covarianceParameters(lme);
    r.sd_subject   = sqrt(psi{1}(1,1));
    r.sd_resid     = sqrt(mse);
    % expect df ~2 with 2 vs 2; a ~0 subject SD means the SE is voxel-level
    r.df_collapsed = (strcmp(r.df_method, 'satterthwaite') && r.df_group > r.n_subjects) || ...
                     r.sd_subject < 1e-3 * r.sd_resid;
    if r.df_collapsed
        r.note = 'subject variance ~0: group p is pseudoreplicated';
    end
catch ME
    r.note = ME.message;
end
end


function r = lme_group_contrast(S, col, ka, kb, cfg)
% y ~ group * inA + (1 + inA | subject), diagonal random-effects
% covariance, using cells ka (inA = 1) and kb (inA = 0) only.
% beta_inA_control  : A - B contrast in control
% beta_interaction  : severe - control difference in that contrast
r = struct('n_subjects',0, 'n_control',0, 'n_severe',0, 'n_voxels',0, ...
           'beta_inA_control',NaN, 'p_inA_control',NaN, ...
           'beta_interaction',NaN, 'se_interaction',NaN, 't_interaction',NaN, ...
           'df_interaction',NaN, 'p_interaction',NaN, 'ci_low',NaN, 'ci_high',NaN, ...
           'sd_subject_intercept',NaN, 'sd_subject_slope',NaN, 'sd_resid',NaN, ...
           'df_method','', 'slope_sd_singular',false, 'note','');
[y, inA, gc, sc, used] = stack_voxels(S, col, ka, kb, cfg);
ng = [sum(used & [S.g] == 1), sum(used & [S.g] == 2)];
r.n_subjects = sum(used); r.n_control = ng(1); r.n_severe = ng(2); r.n_voxels = numel(y);
if any(ng < 1) || sum(ng) < 3
    r.note = 'insufficient subjects'; return
end
tbl = table(y, inA, categorical(gc, [1 2], {'control','severe'}), categorical(sc), ...
            'VariableNames', {'y','inA','group','subject'});
clear y inA gc sc
try
    lme = fitlme(tbl, 'y ~ group*inA + (1 + inA | subject)', 'FitMethod','REML', ...
                 'CovariancePattern','Diagonal');
    [st, r.df_method] = fe_table(lme, r.n_subjects - 2);
    nm = cellstr(st.Name);
    ia = find(strcmp(nm, 'inA'), 1);
    ix = find(contains(nm, 'group_severe') & contains(nm, 'inA'), 1);
    r.beta_inA_control = st.Estimate(ia);
    r.p_inA_control    = st.pValue(ia);
    r.beta_interaction = st.Estimate(ix);
    r.se_interaction   = st.SE(ix);
    r.t_interaction    = st.tStat(ix);
    r.df_interaction   = st.DF(ix);
    r.p_interaction    = st.pValue(ix);
    r.ci_low           = st.Lower(ix);
    r.ci_high          = st.Upper(ix);
    [psi, mse] = covarianceParameters(lme);
    r.sd_subject_intercept = sqrt(psi{1}(1,1));
    r.sd_subject_slope     = sqrt(psi{1}(2,2));
    r.sd_resid             = sqrt(mse);
    r.slope_sd_singular    = r.sd_subject_slope < 1e-3 * r.sd_resid;
    if r.slope_sd_singular
        r.note = 'random-slope SD ~0: interaction p is pseudoreplicated';
    end
catch ME
    r.note = ME.message;
end
end


function r = welch_cell(S, col, k, cfg, summ)
% Welch t-test (unequal variances) on subject-level summaries of cell k.
% summ = 'median' | 'mean', computed on all of the subject's voxels in the
% cell (subjects with < cfg.min_n voxels excluded). diff = severe - control.
r = struct('summary',summ, 'n_control',0, 'n_severe',0, ...
           'subj_control','', 'subj_severe','', 'vals_control','', 'vals_severe','', ...
           'mean_control',NaN, 'mean_severe',NaN, 'diff',NaN, 't',NaN, 'df',NaN, ...
           'p',NaN, 'ci_low',NaN, 'ci_high',NaN, 'note','');
v = NaN(1, numel(S));
for s = 1:numel(S)
    y = cfg.yscale * S(s).M(S(s).C(:,k), col);
    if numel(y) < cfg.min_n, continue; end
    if strcmp(summ, 'median'), v(s) = median(y); else, v(s) = mean(y); end
end
ic = find([S.g] == 1 & ~isnan(v));  is = find([S.g] == 2 & ~isnan(v));
xc = v(ic);  xs = v(is);
r.n_control    = numel(xc);  r.n_severe = numel(xs);
r.subj_control = strjoin(upper({S(ic).name}), ';');
r.subj_severe  = strjoin(upper({S(is).name}), ';');
r.vals_control = strjoin(arrayfun(@(x) sprintf('%.5g', x), xc, 'UniformOutput', false), ';');
r.vals_severe  = strjoin(arrayfun(@(x) sprintf('%.5g', x), xs, 'UniformOutput', false), ';');
r.mean_control = mean(xc);   r.mean_severe = mean(xs);
r.diff         = r.mean_severe - r.mean_control;
if numel(xc) < 2 || numel(xs) < 2
    r.note = 'fewer than 2 subjects in a group'; return
end
[~, p, ci, st] = ttest2(xs, xc, 'Vartype','unequal');
r.t = st.tstat;  r.df = st.df;  r.p = p;  r.ci_low = ci(1);  r.ci_high = ci(2);
end


function [st, method] = fe_table(lme, df_bw)
% Fixed-effect table (Name, Estimate, SE, tStat, DF, pValue, Lower, Upper).
% Satterthwaite df if fixedEffects returns them; otherwise between-within
% df (df_bw = n_subjects - number of between-subject fixed effects) for
% every term, with t, p and 95% CI recomputed from the model SE.
[~, ~, ds] = fixedEffects(lme, 'DFMethod','satterthwaite');
method = 'satterthwaite';
if any(isnan(ds.DF))
    [~, ~, ds] = fixedEffects(lme, 'DFMethod','none');
    method = sprintf('between-within (df = %d)', df_bw);
end
st = struct('Name',{cellstr(ds.Name)}, 'Estimate',double(ds.Estimate), ...
            'SE',double(ds.SE), 'tStat',double(ds.tStat), 'DF',double(ds.DF), ...
            'pValue',double(ds.pValue), 'Lower',double(ds.Lower), 'Upper',double(ds.Upper));
if ~strcmp(method, 'satterthwaite')
    if df_bw < 1
        st.DF(:) = NaN; st.pValue(:) = NaN; st.Lower(:) = NaN; st.Upper(:) = NaN;
        method = 'none (fewer than 3 subjects)';
        return
    end
    tc         = tinv(0.975, df_bw);
    st.DF(:)   = df_bw;
    st.tStat   = st.Estimate ./ st.SE;
    st.pValue  = 2 * tcdf(-abs(st.tStat), df_bw);
    st.Lower   = st.Estimate - tc * st.SE;
    st.Upper   = st.Estimate + tc * st.SE;
end
end


function v = cap_one(v, cap)
v = v(:);
if numel(v) > cap, v = v(randperm(numel(v), cap)); end
end


function a = merge_structs(a, b)
fn = fieldnames(b);
for i = 1:numel(fn), a.(fn{i}) = b.(fn{i}); end
end


function v = cell_verdict(Tk)
% Descriptive separation criteria for one region x property x cell.
isb = strcmp(Tk.pair_type, 'between');
db  = Tk.cliffs_delta(isb);   dw = Tk.abs_delta(~isb);
kb  = Tk.ks_D(isb);           kw = Tk.ks_D(~isb);
v = struct();
v.n_between           = sum(isb);
v.n_between_valid     = sum(~isnan(db));
v.n_within_valid      = sum(~isnan(dw));
v.delta_between_min   = min_finite(db);
v.delta_between_max   = max_finite(db);
v.absdelta_within_max = max_finite(dw);
v.ksD_between_min     = min_finite(kb);
v.ksD_within_max      = max_finite(kw);
if ismember('ks_p', Tk.Properties.VariableNames)
    v.ksp_between_max = max_finite(Tk.ks_p(isb));
else
    v.ksp_between_max = NaN;
end
complete     = v.n_between > 0 && v.n_between_valid == v.n_between && v.n_within_valid > 0;
v.concordant = v.n_between > 0 && v.n_between_valid == v.n_between && ...
               (all(db > 0) || all(db < 0));
v.sep_delta  = complete && v.concordant && min(abs(db)) > v.absdelta_within_max;
v.sep_ks     = complete && all(~isnan(kb)) && v.ksD_between_min > v.ksD_within_max;
end


function fh = plot_between_violins(Yg, stg, dpool, Dpool, ppool, V, L, W, ov, P, cfg)
fh = figure('Visible',cfg.fig_visible, 'Color','w', 'Position',[100 100 920 980], ...
            'CreateFcn','set(gcbo,''Visible'',''on'')');   % .fig reopens visible
ax = axes(fh); hold(ax,'on');
yl = P.ylims; sp = diff(yl);
xg = @(k,g) cfg.xc(k) + (2*g - 3)*cfg.dx;   % g = 1 control (left), 2 severe (right)

hg   = gobjects(1,2);
prof = cell(2,4);                            % violin profiles {y grid, half-width}
for k = 1:4
    for g = 1:2
        x = xg(k,g);
        prof{g,k} = {[], []};
        if stg(g,k).n >= cfg.min_n
            [h, yi, fw] = draw_violin(ax, x, Yg{g,k}, cfg.gcol(g,:), cfg.hw, cfg.kde_max_n, yl);
            prof{g,k} = {yi, fw};
            if ~isgraphics(hg(g)), hg(g) = h; end
            plot(ax, [x x], [stg(g,k).q1 stg(g,k).q3], 'k-', 'LineWidth', 3);
            plot(ax, x, stg(g,k).med, 'o', 'MarkerFaceColor','w', ...
                 'MarkerEdgeColor','k', 'MarkerSize',6);
        else
            text(ax, x, yl(1) + 0.5*sp, sprintf('n < %d', cfg.min_n), ...
                 'HorizontalAlignment','center', 'Rotation',90, 'Color',[.4 .4 .4]);
        end
    end
end

%%% Subject median overlays, kept inside the violin. The offset from the
% violin centre is capped at cfg.med_inside x the violin half-width at the
% median, so markers in a narrow tail move toward the centreline.
% 0 < n < MIN_N: hollow marker (median shown, no effect sizes).
% n = 0: label at the base of the cell.
hs = gobjects(1, numel(ov));
for s = 1:numel(ov)
    g = ov(s).g;
    for k = 1:4
        x0 = xg(k, g);
        if ov(s).n(k) == 0
            text(ax, x0 + ov(s).off, yl(1) + 0.02*sp, sprintf('%s n=0', ov(s).name), ...
                 'Rotation',90, 'FontSize',7, 'Color',cfg.gcol(g,:), ...
                 'HorizontalAlignment','left', 'VerticalAlignment','middle');
            continue
        end
        y   = ov(s).med(k);
        hwk = local_halfwidth(prof{g,k}, y);
        off = sign(ov(s).off) * min(abs(ov(s).off), cfg.med_inside * hwk);
        if ov(s).n(k) < cfg.min_n
            fc = 'w'; ec = cfg.gcol(g,:);
        else
            fc = cfg.gcol(g,:); ec = 'w';
        end
        h = plot(ax, x0 + off, y, ov(s).marker, 'LineStyle','none', 'MarkerSize',6.5, ...
                 'MarkerFaceColor',fc, 'MarkerEdgeColor',ec, 'LineWidth',0.8, ...
                 'DisplayName',ov(s).name);
        if ~isgraphics(hs(s)) && ov(s).n(k) >= cfg.min_n, hs(s) = h; end
    end
    if ~isgraphics(hs(s))                    % legend entry if every cell was low-n / empty
        hs(s) = plot(ax, NaN, NaN, ov(s).marker, 'LineStyle','none', 'MarkerSize',6.5, ...
                     'MarkerFaceColor',cfg.gcol(g,:), 'MarkerEdgeColor','w', ...
                     'DisplayName',ov(s).name);
    end
end
any_low = any(arrayfun(@(o) any(o.n > 0 & o.n < cfg.min_n), ov));

%%% EPVS separator
xline(ax, mean(cfg.xc(2:3)), ':', 'Color',[.5 .5 .5]);

%%% Between-group bracket and effect labels per cell
h0 = 0.04; tk = 0.015;
for k = 1:4
    xa = xg(k,1); xb = xg(k,2);
    yb = yl(2) + h0*sp;
    plot(ax, [xa xa xb xb], [yb - tk*sp, yb, yb, yb - tk*sp], 'k-', 'LineWidth', 0.8);
    if isnan(dpool(k))
        lines = {'\delta = n/a'};
    else
        lines = {sprintf('\\delta = %+.3f', dpool(k)), ...
                 sprintf('D = %.3f,  p %s', Dpool(k), fmt_p(ppool(k)))};
        if ~isnan(L(k).beta_group)
            lb = sprintf(['LME \\beta = ' cfg.beta_fmt ',  p %s'], L(k).beta_group, fmt_p(L(k).p_group));
            if L(k).df_collapsed, lb = [lb ' (df collapsed)']; end
            lines{end+1} = lb;
        elseif ~isempty(L(k).note)
            lines{end+1} = 'LME: fit failed (see CSV note)';
        end
        if ~isnan(W(k).p)
            lines{end+1} = sprintf(['Welch_{med} \\Delta = ' cfg.beta_fmt ',  p %s'], W(k).diff, fmt_p(W(k).p));
        end
        if ~isnan(V(k).delta_between_min)
            l3 = sprintf('pairs [%+.2f, %+.2f]', V(k).delta_between_min, V(k).delta_between_max);
            l4 = sprintf('within |\\delta| \\leq %.2f', V(k).absdelta_within_max);
            if V(k).sep_delta
                l3 = ['{\bf ' l3 '}'];
                l4 = ['{\bf ' l4 '}'];
            end
            lines = [lines, {l3, l4}]; %#ok<AGROW>
        end
    end
    text(ax, (xa + xb)/2, yb + 0.01*sp, lines, 'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', 'FontSize',7.5);
end
ytop_br = yl(2) + 0.33*sp;

%%% EPVS bin labels above the brackets
pairs = {[1 2], [3 4]};
elab  = {'Low EPVS-SWP', 'High EPVS-SWP'};
for g = 1:2
    text(ax, mean(cfg.xc(pairs{g})), ytop_br, elab{g}, ...
         'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontWeight','bold');
end

%%% Legend
okg = isgraphics(hg);
hl  = [hg(okg), hs];
lab = [cellfun(@(x) [x ' (pooled)'], cfg.glab(okg), 'UniformOutput', false), {ov.name}];
if any_low
    hl(end+1)  = plot(ax, NaN, NaN, 'o', 'LineStyle','none', 'MarkerSize',6.5, ...
                      'MarkerFaceColor','w', 'MarkerEdgeColor',[.3 .3 .3]);
    lab{end+1} = sprintf('subject n < %d', cfg.min_n);
end
legend(ax, hl, lab, 'Location','southoutside', 'Orientation','horizontal', ...
       'Box','off', 'AutoUpdate','off');

%%% Axes
ylim(ax, [yl(1), ytop_br + 0.08*sp]);
yt = ax.YTick; ax.YTick = yt(yt <= yl(2));
xlim(ax, [cfg.xc(1) - 0.7, cfg.xc(4) + 0.7]);
tl = cell(1,4);
for k = 1:4
    tl{k} = sprintf('%s\\newline{\\fontsize{8}n_C %s / n_S %s}', cfg.vlab{cfg.vbin(k)}, ...
                    fmt_n(stg(1,k).n), fmt_n(stg(2,k).n));
end
ax.XTick = cfg.xc;
ax.XTickLabel = tl;
ax.TickLabelInterpreter = 'tex';
ylabel(ax, P.ylabel);
title(ax, P.title, 'FontWeight','normal', 'FontSize',9.5);
set(ax, 'Box','off', 'TickDir','out', 'FontSize',11);
end


function hw = local_halfwidth(pr, y)
% Violin half-width at y (0 outside the drawn support or with no violin).
yi = pr{1}; fw = pr{2};
if isempty(yi) || isnan(y), hw = 0; return; end
hw = interp1(yi, fw, y, 'linear', 0);
end


function plot_pair_summary(T, metric, base, scfg)
% Subject-pair effect per cell. Rows = property (scfg.props), columns = region.
% metric = 'delta' (signed Cliff's delta) or 'ks' (KS statistic D).
if isempty(T)
    warning('between:summary', 'No pairwise rows for the pair summary.');
    return
end
props     = scfg.props(ismember(scfg.props, unique(T.property)));
pdisp_all = containers.Map(scfg.props, scfg.prop_disp);
np_rows   = numel(props);
regions   = {'front','occip'};
reg_disp  = {'Frontal','Occipital'};
cell_lab  = {'EPVS-lo/Ves-lo','EPVS-lo/Ves-hi','EPVS-hi/Ves-lo','EPVS-hi/Ves-hi'};

%%% Between-pair identity: marker = severe subject, colour = control subject
isb_all = strcmp(T.pair_type, 'between');
sev = order_names(unique(T.subj_A(isb_all)), scfg.subs_order);
ctl = order_names(unique(T.subj_B(isb_all)), scfg.subs_order);
mk_sev = {'^','d','v','p','h'};
ccol   = [0.55 0.05 0.30; 0.95 0.50 0.70; 0.75 0.25 0.50; 0.40 0.00 0.20];
pairs  = cell(0,4);
for i = 1:numel(sev)
    for j = 1:numel(ctl)
        pairs(end+1,:) = {sev{i}, ctl{j}, mk_sev{min(i,end)}, ccol(min(j,end),:)}; %#ok<AGROW>
    end
end
np = size(pairs,1);
if np == 1, dyp = 0; else, dyp = linspace(-0.16, 0.16, np); end

wcol = struct('within_control',[0.35 0.35 0.35], 'within_severe',[0.80 0.55 0.68]);
wy   = struct('within_control',-0.34, 'within_severe',0.34);

switch metric
    case 'delta'
        xl   = [-1 1];
        xlab = 'Cliff''s \delta   (between pairs: severe vs control; > 0: severe larger)';
        ttl  = 'Subject-pair Cliff''s \delta per cell';
    case 'ks'
        xl   = [0 1];
        xlab = 'KS statistic D';
        ttl  = 'Subject-pair KS D per cell';
    otherwise
        error('metric must be ''delta'' or ''ks''.');
end
ylim_ = [0.4 4.6];

fh = figure('Visible',scfg.fig_visible, 'Color','w', ...
            'Position',[100 100 1000 330*np_rows + 60], ...
            'CreateFcn','set(gcbo,''Visible'',''on'')');
tl = tiledlayout(fh, np_rows, 2, 'TileSpacing','compact', 'Padding','compact');
for p = 1:np_rows
    for r = 1:2
        ax = nexttile(tl); hold(ax, 'on');
        if strcmp(metric, 'delta')
            for xr = [-0.474 -0.33 -0.147 0.147 0.33 0.474]
                plot(ax, [xr xr], ylim_, ':', 'Color',[0.7 0.7 0.7]);
            end
            plot(ax, [0 0], ylim_, 'k-', 'LineWidth', 0.8);
        end
        plot(ax, xl, [2.5 2.5], '-', 'Color',[0.85 0.85 0.85]);   % EPVS lo | hi

        sel = strcmp(T.property, props{p}) & strcmp(T.region, regions{r});
        for k = 1:4
            Tk = T(sel & T.cell_idx == k, :);
            if isempty(Tk), continue; end
            v   = cell_verdict(Tk);
            isb = strcmp(Tk.pair_type, 'between');

            %%% Within-group envelope
            if strcmp(metric, 'delta')
                w = v.absdelta_within_max;
                if ~isnan(w)
                    patch(ax, [-w w w -w], k + [-0.44 -0.44 0.44 0.44], [0.90 0.90 0.90], ...
                          'EdgeColor','none');
                end
            else
                w = v.ksD_within_max;
                if ~isnan(w)
                    patch(ax, [0 w w 0], k + [-0.44 -0.44 0.44 0.44], [0.90 0.90 0.90], ...
                          'EdgeColor','none');
                end
            end

            %%% Within-group pairs
            for j = find(~isb)'
                pt = Tk.pair_type{j};
                if isfield(wy, pt), yy = k + wy.(pt); c = wcol.(pt); else, yy = k; c = [.5 .5 .5]; end
                if strcmp(metric, 'delta')
                    a = Tk.abs_delta(j);
                    if isnan(a), continue; end
                    plot(ax, [-a a], [yy yy], '-', 'Color',c, 'LineWidth',1.5);
                    plot(ax, [-a a], [yy yy], '|', 'Color',c, 'MarkerSize',6, 'LineWidth',1.5);
                else
                    a = Tk.ks_D(j);
                    if isnan(a), continue; end
                    plot(ax, a, yy, 'o', 'MarkerSize',5, 'MarkerEdgeColor',c, 'MarkerFaceColor',c);
                end
            end

            %%% Between-group pairs
            for j = find(isb)'
                ip = find(strcmp(pairs(:,1), Tk.subj_A{j}) & strcmp(pairs(:,2), Tk.subj_B{j}), 1);
                if isempty(ip), continue; end
                if strcmp(metric, 'delta'), xv = Tk.cliffs_delta(j); else, xv = Tk.ks_D(j); end
                if isnan(xv), continue; end
                if min(Tk.frac_A(j), Tk.frac_B(j)) < scfg.sparse_frac, fc = 'w'; else, fc = pairs{ip,4}; end
                plot(ax, xv, k + dyp(ip), pairs{ip,3}, 'MarkerSize',7, ...
                     'MarkerEdgeColor',pairs{ip,4}, 'MarkerFaceColor',fc, 'LineWidth',1);
            end

            %%% Verdict
            if (strcmp(metric,'delta') && v.sep_delta) || (strcmp(metric,'ks') && v.sep_ks)
                plot(ax, xl(2) - 0.04*diff(xl), k, 'p', 'MarkerSize',10, ...
                     'MarkerFaceColor','k', 'MarkerEdgeColor','k');
            end
        end

        set(ax, 'YDir','reverse', 'YLim',ylim_, 'XLim',xl, 'YTick',1:4, ...
                'Box','off', 'TickDir','out', 'FontSize',10, 'TickLabelInterpreter','tex');
        if r == 1, ax.YTickLabel = cell_lab; else, ax.YTickLabel = {}; end
        title(ax, sprintf('%s  |  %s', reg_disp{r}, pdisp_all(props{p})), 'FontWeight','normal');
        if p == np_rows, xlabel(ax, xlab); end
    end
end

%%% Shared legend (dummy handles)
hl = gobjects(0); lab = {};
for ip = 1:np
    hl(end+1) = plot(ax, NaN, NaN, pairs{ip,3}, 'LineStyle','none', 'MarkerSize',7, ...
                     'MarkerEdgeColor',pairs{ip,4}, 'MarkerFaceColor',pairs{ip,4}); %#ok<AGROW>
    lab{end+1} = sprintf('%s vs %s', pairs{ip,1}, pairs{ip,2}); %#ok<AGROW>
end
wf = fieldnames(wcol);
for i = 1:numel(wf)
    if ~any(strcmp(T.pair_type, wf{i})), continue; end
    wp = unique(strcat(T.subj_A(strcmp(T.pair_type, wf{i})), '-', T.subj_B(strcmp(T.pair_type, wf{i}))));
    if strcmp(metric, 'delta')
        hl(end+1) = plot(ax, NaN, NaN, '-', 'Color',wcol.(wf{i}), 'LineWidth',1.5); %#ok<AGROW>
        lab{end+1} = sprintf('within: %s (\\pm|\\delta|)', strjoin(wp, ', ')); %#ok<AGROW>
    else
        hl(end+1) = plot(ax, NaN, NaN, 'o', 'LineStyle','none', 'MarkerSize',5, ...
                         'MarkerEdgeColor',wcol.(wf{i}), 'MarkerFaceColor',wcol.(wf{i})); %#ok<AGROW>
        lab{end+1} = sprintf('within: %s', strjoin(wp, ', ')); %#ok<AGROW>
    end
end
hl(end+1) = patch(ax, NaN, NaN, [0.90 0.90 0.90], 'EdgeColor','none');
lab{end+1} = 'within-group envelope';
hl(end+1) = plot(ax, NaN, NaN, 'o', 'LineStyle','none', 'MarkerSize',7, ...
                 'MarkerEdgeColor',[.3 .3 .3], 'MarkerFaceColor','w');
lab{end+1} = sprintf('cell < %g%% of voxels', 100*scfg.sparse_frac);
hl(end+1) = plot(ax, NaN, NaN, 'p', 'LineStyle','none', 'MarkerSize',10, ...
                 'MarkerEdgeColor','k', 'MarkerFaceColor','k');
if strcmp(metric, 'delta')
    lab{end+1} = 'concordant & min|\delta_{btw}| > max|\delta_{within}|';
else
    lab{end+1} = 'min D_{btw} > max D_{within}';
end
lg = legend(ax, hl, lab, 'Orientation','horizontal', 'Box','off', 'AutoUpdate','off', ...
            'NumColumns',4, 'FontSize',9);
lg.Layout.Tile = 'south';

title(tl, sprintf('%s  (%s edges, %s)', ttl, T.edge_mode{1}, T.value_type{1}), 'FontSize',12);
save_fig(fh, scfg.out_dir, base, scfg.max_path);
end


function [h, yi, f] = draw_violin(ax, xc, y, col, hw, kde_max_n, yl)
% KDE on a subsample; display support clipped to the cell's 0.5-99.5
% percentiles and to the shared y-limits. Returns the graphics handle and
% the drawn profile (y grid, half-width) for placing subject markers.
if numel(y) > kde_max_n
    y = y(randperm(numel(y), kde_max_n));
end
lims = prctile(y, [0.5 99.5]);
lims = [max(lims(1), yl(1)), min(lims(2), yl(2))];
if lims(2) <= lims(1)
    h = plot(ax, xc + [-hw hw], lims(1)*[1 1], 'Color', col, 'LineWidth', 2);
    yi = []; f = [];
    return
end
yi = linspace(lims(1), lims(2), 200);
f  = ksdensity(y, yi);
f  = f ./ max(f) * hw;
h  = patch(ax, [xc - f, fliplr(xc + f)], [yi, fliplr(yi)], col, ...
           'FaceAlpha',0.55, 'EdgeColor',col*0.7, 'LineWidth',0.8);
end


function [a, b] = cap_pair(a, b, cap)
% Random subsample of each sample to at most cap values.
a = a(:); b = b(:);
if numel(a) > cap, a = a(randperm(numel(a), cap)); end
if numel(b) > cap, b = b(randperm(numel(b), cap)); end
end


function d = cliffs_delta(a, b)
% delta = P(a > b) - P(a < b), exact via ranks (ties count 1/2).
a = a(:); b = b(:);
n1 = numel(a); n2 = numel(b);
if n1 == 0 || n2 == 0, d = NaN; return; end
r = tiedrank([a; b]);
d = 2 * (mean(r(1:n1)) - (n1 + 1)/2) / n2 - 1;
end


function D = ks_stat(a, b)
% Two-sample KS statistic, max |F_a - F_b| over all observed values.
a = a(:); b = b(:);
if isempty(a) || isempty(b), D = NaN; return; end
xs    = unique([a; b]);
edges = [xs; Inf];
Fa = cumsum(histcounts(a, edges)) / numel(a);   % F_a(xs_i) = P(a <= xs_i)
Fb = cumsum(histcounts(b, edges)) / numel(b);
D  = max(abs(Fa - Fb));
end


function p = ks_pval(D, na, nb)
% Asymptotic two-sided KS p-value (Kolmogorov distribution with Stephens'
% small-sample correction), the large-sample formula used by kstest2.
if isnan(D) || na == 0 || nb == 0, p = NaN; return; end
ne  = na * nb / (na + nb);
lam = (sqrt(ne) + 0.12 + 0.11/sqrt(ne)) * D;
if lam < 0.2, p = 1; return; end            % Q_KS(0.2) = 1 to double precision
j = (1:101)';
p = 2 * sum((-1).^(j - 1) .* exp(-2 * lam^2 * j.^2));
p = min(max(p, 0), 1);                       % underflows to 0 for large lam
end


function s = fmt_p(p)
% TeX-formatted p for figure labels.
if isnan(p)
    s = '= n/a';
elseif p < realmin
    s = '< 10^{-307}';
elseif p < 1e-3
    e = floor(log10(p));
    s = sprintf('= %.1f\\times10^{%d}', p / 10^e, e);
else
    s = sprintf('= %.3f', p);
end
end


function C = cell_masks(epv, ves, eE, eV)
El = epv <= eE(1);  Eh = epv >= eE(2);
Vl = ves <= eV(1);  Vh = ves >= eV(2);
C  = [El & Vl, El & Vh, Eh & Vl, Eh & Vh];
end


function r = bin_row(level, name, grp, reg, prop, k, eE, eV, st, ntot, base, cfg)
r = struct('level',level, 'name',name, 'group',grp, 'region',reg, ...
           'property',prop, 'value_type',cfg.vt, 'edge_mode',cfg.edge_mode, ...
           'equal_weight',cfg.eqw, 'cell',cfg.cell_lab{k}, ...
           'epvs_bin',cfg.ebin_lab{k}, 'ves_bin',cfg.vbin_lab{k}, ...
           'q_low',cfg.Q(1), 'q_high',cfg.Q(2), ...
           'epvs_edge_low',eE(1), 'epvs_edge_high',eE(2), ...
           'ves_edge_low',eV(1), 'ves_edge_high',eV(2), ...
           'n',st.n, 'frac_of_voxels',st.n/ntot, 'units',cfg.units, ...
           'median',st.med, 'q1',st.q1, 'q3',st.q3, 'figure',base);
end


function X = app(X, r)
% Append struct(s) r to column struct array X (X may start as []).
X = [X; r(:)];
end


function out = order_names(names, order)
names = reshape(names, 1, []);
out   = [order(ismember(order, names)), setdiff(names, order, 'stable')];
end


function m = min_finite(x)
x = x(~isnan(x));
if isempty(x), m = NaN; else, m = min(x); end
end


function m = max_finite(x)
x = x(~isnan(x));
if isempty(x), m = NaN; else, m = max(x); end
end


function s = edge_line(eE, eV)
s = sprintf(['EPVS-SWP: low \\leq %.3g, high \\geq %.3g    ' ...
             'Ves-SWP: low \\leq %.3g, high \\geq %.3g'], eE(1), eE(2), eV(1), eV(2));
end


function e = swp_edges(x, Q, label)
e = quantile(x, Q);
if e(1) >= e(2)
    warning('violin:edges', ...
        '%s: low/high quantile edges coincide (%.4g); bins overlap due to ties.', label, e(1));
end
end


function st = cell_stats(y)
if isempty(y)
    st = struct('n',0, 'med',NaN, 'q1',NaN, 'q3',NaN);
else
    q  = quantile(y, [0.25 0.5 0.75]);
    st = struct('n',numel(y), 'med',q(2), 'q1',q(1), 'q3',q(3));
end
end


function M = clean_rows(M)
M = M(all(isfinite(M), 2), :);
end


function s = fmt_n(n)
if n >= 1e6
    s = sprintf('%.2fM', n/1e6);
elseif n >= 1e4
    s = sprintf('%.1fk', n/1e3);
else
    s = sprintf('%d', n);
end
end


function check_path(f, max_len)
% Error if the full output path exceeds max_len characters.
if numel(f) > max_len
    error('between:path', 'Output path is %d characters (> %d):\n%s', ...
          numel(f), max_len, f);
end
end


function write_csv(T, f, max_len)
check_path(f, max_len);
writetable(T, f);
end


function save_fig(fh, out_dir, base, max_len)
ext = {'.png','.svg','.fig'};
for i = 1:numel(ext), check_path(fullfile(out_dir, [base ext{i}]), max_len); end
exportgraphics(fh, fullfile(out_dir, [base '.png']), 'Resolution', 300);
exportgraphics(fh, fullfile(out_dir, [base '.svg']), 'ContentType', 'vector');
savefig(fh, fullfile(out_dir, [base '.fig']));
close(fh);
end