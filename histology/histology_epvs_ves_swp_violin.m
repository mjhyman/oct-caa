%% Violin plots of histology stain intensity in EPVS-SWP x vessel-SWP quantile bins
% Histology counterpart of plot_mus_ret_swp_violins.m: same quantile
% binning, four-violin layout, Cliff's delta contrasts, CSVs and delta-
% summary panel. The measured quantity is the z-scored stain intensity
% (CD68, GFAP, LHE) instead of mus / retardance.
%
% Per section (as in the histology GAM pipeline):
%   tissue mask with vessels and EPVS removed
%   N x 3 matrix [stain_z, ves_swp, epv_swp]
% Sections are concatenated per subject x region (frontal = _1,
% occipital = _7) using the same baseName patterns as the GAM script.
%
% Cells (per stain, per subject or group, per region):
%   x = 1.0 : low EPVS-SWP  & low vessel-SWP
%   x = 2.0 : low EPVS-SWP  & high vessel-SWP
%   x = 3.6 : high EPVS-SWP & low vessel-SWP
%   x = 4.6 : high EPVS-SWP & high vessel-SWP
%   low : x <= quantile(x, Q_EDGES(1));  high : x >= quantile(x, Q_EDGES(2))
%   Each SWP is binned on its own marginal distribution.
% Subject figures : edges from that subject's own SWP distribution.
% Group figures   : GROUP_EDGE_MODE = 'withinSubject' | 'group' | 'shared'
%                   (see plot_mus_ret_swp_violins.m). CAA17 is subject-level
%                   only (not part of the 2x2 design).
%
% Cliff's delta (descriptive; voxels are pseudoreplicates), all six pairs:
%   delta = P(A > B) - P(A < B), A = right-hand cell, B = left-hand cell
%
% Outputs (date-stamped) in out_dir:
%   violin_<stain>_<unit>_<region>[_tags]_<qLL-HH>_<EXPL|MS>_<yyyymmdd>.{png,svg,fig}
%   histo_violin_bin_summary_<qLL-HH>_<EXPL|MS>_<yyyymmdd>.csv
%   histo_violin_cliffs_delta_<qLL-HH>_<EXPL|MS>_<yyyymmdd>.csv
%   histo_delta_summary_<qLL-HH>_<EXPL|MS>_<yyyymmdd>.{png,svg,fig}
%     rows = stain, columns = frontal / occipital; one marker per subject.
%     Hollow = contrast involving a cell with < SPARSE_FRAC of the voxels.
%
% NOTE: histology SWP uses a 4 mm radius (ds4) versus 500 um for OCT, so the
% "neighbourhood" in the two figures is not the same spatial scale.

%% Initialization
clearvars -except histo_op caa swp_op swp_struct
clc; close all;
addpath(genpath('/projectnb/npbssmic/s/mhyman/oct-caa'));
maxNumCompThreads(feature('numcores'));
rng(1);

%% Parameters
EXPLORATORY           = false;           % EXPL = smaller KDE / Cliff subsamples
FORCE_REBUILD         = false;           % true = rebuild histo_op even if in workspace
Q_EDGES               = [0.25 0.75];     % low: x <= q(1); high: x >= q(2)
GROUP_EDGE_MODE       = 'withinSubject'; % 'withinSubject' | 'group' | 'shared'
EQUAL_WEIGHT_SUBJECTS = true;            % subsample group subjects to equal n per region
MAKE_GROUP_FIGS       = true;            % false = subject-level figures only
MIN_N                 = 50;              % min voxels to draw a violin / compute delta
ANNOTATE_PAIRS        = 'all';           % 'all' = 6 brackets; 'within' = 2
FIG_VISIBLE           = 'off';

%%% Optional per-stain y-limit lookup (overrides the automatic limits), e.g.
%   YLIM_LOOKUP.gfap.front = [-2 4];
YLIM_LOOKUP           = struct();

%%% Delta summary panel
SPARSE_FRAC           = 0.01;            % hollow marker if a cell holds < 1% of voxels
SUMMARY_GROUPS        = {'control','severe'};   % add 'moderate' to show CAA17
DELTA_SUMMARY_ONLY    = false;           % true = redraw summary from DELTA_CSV and stop
DELTA_CSV             = '';              % existing histo_violin_cliffs_delta_*.csv

if ~any(strcmp(GROUP_EDGE_MODE, {'withinSubject','group','shared'}))
    error('GROUP_EDGE_MODE must be ''withinSubject'', ''group'' or ''shared''.');
end
if EXPLORATORY
    run_tag = 'EXPL'; KDE_MAX_N = 2e4; CLIFF_MAX_N = 5e5;
else
    run_tag = 'MS';   KDE_MAX_N = 1e5; CLIFF_MAX_N = 5e6;
end
date_str = char(datetime('now','Format','yyyyMMdd'));
q_tag    = sprintf('q%02d-%02d', round(100*Q_EDGES));

%%% Directories / inputs
hdir    = '/projectnb/npbssmic/ns/CAA/histology/swp';
out_dir = '/projectnb/npbssmic/ns/CAA/histology/figures/histo_swp_violin_plots/';
if ~exist(out_dir,'dir'); mkdir(out_dir); end
stain_files = struct( ...
    'cd68', 'histo_swp_hmatched_zscore__CD68_voxelwise_radius4mm_ds4_20-Mar-2026.mat', ...
    'gfap', 'histo_swp_hmatched_zscore__GFAP_voxelwise_radius4mm_ds4_20-Mar-2026.mat', ...
    'lhe',  'histo_swp_hmatched_zscore__LHE_voxelwise_radius4mm_ds4_19-Mar-2026.mat');
stain_names = fieldnames(stain_files);

%%% Subjects (patterns match both baseName conventions) and pooled groups
subjects = struct( ...
    'caa26', struct('group','control',  'front',{{'CAA26_1','CAA_26_1'}}, 'occip',{{'CAA26_7','CAA_26_7'}}), ...
    'caa6',  struct('group','control',  'front',{{'CAA6_1', 'CAA_6_1' }}, 'occip',{{'CAA6_7', 'CAA_6_7' }}), ...
    'caa17', struct('group','moderate', 'front',{{'CAA17_1','CAA_17_1'}}, 'occip',{{'CAA17_7','CAA_17_7'}}), ...
    'caa22', struct('group','severe',   'front',{{'CAA22_1','CAA_22_1'}}, 'occip',{{'CAA22_7','CAA_22_7'}}), ...
    'caa25', struct('group','severe',   'front',{{'CAA25_1','CAA_25_1'}}, 'occip',{{'CAA25_7','CAA_25_7'}}));
groups.control = {'caa26','caa6'};
groups.severe  = {'caa22','caa25'};
regions  = {'front','occip'};
reg_disp = struct('front','Frontal','occip','Occipital');

scfg = struct('out_dir',out_dir, 'fig_visible',FIG_VISIBLE, ...
              'sparse_frac',SPARSE_FRAC, 'groups',{SUMMARY_GROUPS}, ...
              'props',{stain_names'});

%% Delta summary only: redraw from an existing CSV (no data loading)
if DELTA_SUMMARY_ONLY
    T_d = readtable(DELTA_CSV, 'TextType','char');
    [~, nm] = fileparts(DELTA_CSV);
    plot_delta_summary(T_d, strrep(nm, 'histo_violin_cliffs_delta_', 'histo_delta_summary_'), scfg);
    fprintf('Delta summary written to %s\n', out_dir);
    return
end

%% Build N x 3 matrices [stain_z, ves_swp, epv_swp] per stain x subject x region
flag_build = ~exist('histo_op','var') || FORCE_REBUILD;
if flag_build
    histo_op = struct();
    subs = fieldnames(subjects);
    for si = 1:numel(stain_names)
        stn = stain_names{si};
        fprintf('Loading %s\n', upper(stn));
        tmp = load(fullfile(hdir, stain_files.(stn)));
        S   = tmp.(stn); clear tmp
        names = lower({S.baseName});
        for s = 1:numel(subs)
            sub = subs{s};
            for r = 1:numel(regions)
                reg = regions{r};
                idx = find(contains(names, lower(subjects.(sub).(reg))));
                if isempty(idx), continue; end
                M = [];
                for k = idx
                    mask = S(k).mask;
                    mask(S(k).ves)  = 0;
                    mask(S(k).epvs) = 0;
                    M = [M; double(S(k).z_stain(mask)), ...
                            double(S(k).ves_swp(mask)), ...
                            double(S(k).epvs_swp(mask))]; %#ok<AGROW>
                end
                histo_op.(stn).(sub).(reg) = M;
                fprintf('\t%s %s: %d section(s), %d voxels\n', upper(sub), reg, numel(idx), size(M,1));
            end
        end
        clear S
    end
else
    fprintf('histo_op found in workspace; skipping import\n')
end

%% Plot configuration (stain-specific fields set inside the stain loop)
cfg = struct();
cfg.prop_col    = 1;          % stain_z
cfg.ves_col     = 2;
cfg.epv_col     = 3;
cfg.units       = 'z-score';
cfg.Q           = Q_EDGES;
cfg.min_n       = MIN_N;
cfg.kde_max_n   = KDE_MAX_N;
cfg.cliff_max_n = CLIFF_MAX_N;
cfg.xpos        = [1 2 3.6 4.6];
cfg.hw          = 0.38;
cfg.vbin        = [1 2 1 2];
cfg.vcol        = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];   % blue = low, orange = high
cfg.vlab        = {'Ves-SWP low','Ves-SWP high'};
cfg.ebin_lab    = {'low','low','high','high'};
cfg.vbin_lab    = {'low','high','low','high'};
cfg.cell_lab    = {'EPVS-lo/Ves-lo','EPVS-lo/Ves-hi','EPVS-hi/Ves-lo','EPVS-hi/Ves-hi'};
cfg.markers     = {'d','s','^','v'};
cfg.fig_visible = FIG_VISIBLE;
cfg.out_dir     = out_dir;
cfg.run_tag     = run_tag;
cfg.date_str    = date_str;
cfg.q_tag       = q_tag;
% Six contrasts (A = right-hand, B = left-hand): {type, held fixed, A, B, bracket level}
cfg.contr = { 'Ves-SWP within EPVS-SWP', 'EPVS-SWP low',  2, 1, 1;
              'Ves-SWP within EPVS-SWP', 'EPVS-SWP high', 4, 3, 1;
              'anti-diagonal',           'none',          3, 2, 2;
              'EPVS-SWP within Ves-SWP', 'Ves-SWP low',   3, 1, 3;
              'EPVS-SWP within Ves-SWP', 'Ves-SWP high',  4, 2, 4;
              'diagonal',                'none',          4, 1, 5 };
switch ANNOTATE_PAIRS
    case 'all',    cfg.annot = true(size(cfg.contr,1), 1);
    case 'within', cfg.annot = strcmp(cfg.contr(:,1), 'Ves-SWP within EPVS-SWP');
    otherwise,     error('ANNOTATE_PAIRS must be ''all'' or ''within''.');
end
if EQUAL_WEIGHT_SUBJECTS, eqw_tag = '_eqW'; else, eqw_tag = ''; end

all_sum = {};
all_d   = {};
grp_names = fieldnames(groups);

%% Stain loop
for si = 1:numel(stain_names)
    stn = stain_names{si};
    if ~isfield(histo_op, stn), continue; end
    D = histo_op.(stn);
    cfg.props     = {stn};
    cfg.prop_disp = {upper(stn)};
    cfg.ylab      = {sprintf('%s (z-score)', upper(stn))};

    %%% Per-stain y-limits per region: lookup if given, else pooled 0.5-99.5 pct
    subs = fieldnames(D);
    cfg.ylims = struct();
    for r = 1:numel(regions)
        reg = regions{r};
        if isfield(YLIM_LOOKUP, stn) && isfield(YLIM_LOOKUP.(stn), reg)
            cfg.ylims.(reg) = YLIM_LOOKUP.(stn).(reg);
            continue
        end
        lo = inf; hi = -inf;
        for s = 1:numel(subs)
            if ~isfield(D.(subs{s}), reg), continue; end
            q  = prctile(clean_col(D.(subs{s}).(reg)(:,cfg.prop_col)), [0.5 99.5]);
            lo = min(lo, q(1)); hi = max(hi, q(2));
        end
        cfg.ylims.(reg) = [lo hi];
    end

    %%% Subject-level figures
    for s = 1:numel(subs)
        sub  = subs{s};
        regs = fieldnames(D.(sub));
        for r = 1:numel(regs)
            reg = regs{r};
            fprintf('%s | subject %s %s\n', upper(stn), upper(sub), reg);
            M  = clean_rows(D.(sub).(reg));
            eE = swp_edges(M(:,cfg.epv_col), Q_EDGES, sprintf('%s %s %s EPVS-SWP', stn, sub, reg));
            eV = swp_edges(M(:,cfg.ves_col), Q_EDGES, sprintf('%s %s %s Ves-SWP',  stn, sub, reg));
            C  = cell_masks(M(:,cfg.epv_col), M(:,cfg.ves_col), eE, eV);
            grp = subjects.(sub).group;
            meta = struct('level','subject', 'name',upper(sub), 'group',grp, ...
                          'region',reg, 'reg_disp',reg_disp.(reg), ...
                          'disp',sprintf('%s (%s)', upper(sub), grp), ...
                          'stub',sprintf('%s_%s', sub, reg), ...
                          'edge_mode','subject', 'eqw',false, ...
                          'edge_title',{{edge_line(eE, eV)}}, 'edge_detail','');
            [sr, dr] = process_unit(M, C, struct('name',{},'M',{},'C',{}), eE, eV, meta, cfg);
            all_sum = [all_sum, sr]; %#ok<AGROW>
            all_d   = [all_d,   dr]; %#ok<AGROW>
        end
    end
    clear M C

    %%% Group-level figures (control and severe on separate figures)
    if ~MAKE_GROUP_FIGS, continue; end
    for r = 1:numel(regions)
        reg = regions{r};
        Ms = struct();
        for g = 1:numel(grp_names)
            gs = groups.(grp_names{g});
            for s = 1:numel(gs)
                if isfield(D, gs{s}) && isfield(D.(gs{s}), reg)
                    Ms.(gs{s}) = clean_rows(D.(gs{s}).(reg));
                end
            end
        end
        ms = fieldnames(Ms);
        if isempty(ms), continue; end

        if EQUAL_WEIGHT_SUBJECTS
            nmin = min(cellfun(@(f) size(Ms.(f),1), ms));
            for s = 1:numel(ms)
                n = size(Ms.(ms{s}),1);
                Ms.(ms{s}) = Ms.(ms{s})(sort(randperm(n, nmin)),:);
            end
        end

        if strcmp(GROUP_EDGE_MODE, 'shared')
            allSWP = cellfun(@(f) Ms.(f)(:,[cfg.ves_col cfg.epv_col]), ms, 'UniformOutput', false);
            allSWP = vertcat(allSWP{:});
            eV_sh = swp_edges(allSWP(:,1), Q_EDGES, sprintf('%s shared %s Ves-SWP',  stn, reg));
            eE_sh = swp_edges(allSWP(:,2), Q_EDGES, sprintf('%s shared %s EPVS-SWP', stn, reg));
            clear allSWP
        end

        for g = 1:numel(grp_names)
            gname = grp_names{g};
            gs = groups.(gname);
            gs = gs(isfield(Ms, gs));
            if isempty(gs), continue; end
            fprintf('%s | group %s %s\n', upper(stn), gname, reg);

            ovl = struct('name', cellfun(@upper, gs, 'UniformOutput', false), ...
                         'M',    cellfun(@(f) Ms.(f), gs, 'UniformOutput', false), ...
                         'C',    cell(1, numel(gs)));
            Mg = vertcat(ovl.M);

            switch GROUP_EDGE_MODE
                case 'withinSubject'
                    eE = NaN(1,2); eV = NaN(1,2);
                    tE = cell(1, numel(gs)); tV = tE; det = tE;
                    for s = 1:numel(gs)
                        Mi  = ovl(s).M;
                        eEs = swp_edges(Mi(:,cfg.epv_col), Q_EDGES, sprintf('%s %s %s EPVS-SWP', stn, gs{s}, reg));
                        eVs = swp_edges(Mi(:,cfg.ves_col), Q_EDGES, sprintf('%s %s %s Ves-SWP',  stn, gs{s}, reg));
                        ovl(s).C = cell_masks(Mi(:,cfg.epv_col), Mi(:,cfg.ves_col), eEs, eVs);
                        tE{s}  = sprintf('%s \\leq %.3g / \\geq %.3g', upper(gs{s}), eEs(1), eEs(2));
                        tV{s}  = sprintf('%s \\leq %.3g / \\geq %.3g', upper(gs{s}), eVs(1), eVs(2));
                        det{s} = sprintf('%s: EPVS<=%.4g,>=%.4g; Ves<=%.4g,>=%.4g', ...
                                         upper(gs{s}), eEs(1), eEs(2), eVs(1), eVs(2));
                    end
                    C = vertcat(ovl.C);
                    edge_title  = { ['EPVS-SWP (per subject):  ' strjoin(tE, ';    ')], ...
                                    ['Ves-SWP (per subject):  '  strjoin(tV, ';    ')] };
                    edge_detail = strjoin(det, ' | ');
                case {'group','shared'}
                    if strcmp(GROUP_EDGE_MODE, 'shared')
                        eE = eE_sh; eV = eV_sh;
                    else
                        eE = swp_edges(Mg(:,cfg.epv_col), Q_EDGES, sprintf('%s %s %s EPVS-SWP', stn, gname, reg));
                        eV = swp_edges(Mg(:,cfg.ves_col), Q_EDGES, sprintf('%s %s %s Ves-SWP',  stn, gname, reg));
                    end
                    C = cell_masks(Mg(:,cfg.epv_col), Mg(:,cfg.ves_col), eE, eV);
                    for s = 1:numel(gs)
                        ovl(s).C = cell_masks(ovl(s).M(:,cfg.epv_col), ovl(s).M(:,cfg.ves_col), eE, eV);
                    end
                    edge_title  = {edge_line(eE, eV)};
                    edge_detail = '';
            end

            meta = struct('level','group', 'name',gname, 'group',gname, ...
                          'region',reg, 'reg_disp',reg_disp.(reg), ...
                          'disp',sprintf('%s%s (%s)', upper(gname(1)), gname(2:end), ...
                                         strjoin(upper(gs), ', ')), ...
                          'stub',sprintf('%s_%s_%sEdges%s', gname, reg, GROUP_EDGE_MODE, eqw_tag), ...
                          'edge_mode',GROUP_EDGE_MODE, 'eqw',EQUAL_WEIGHT_SUBJECTS, ...
                          'edge_title',{edge_title}, 'edge_detail',edge_detail);
            [sr, dr] = process_unit(Mg, C, ovl, eE, eV, meta, cfg);
            all_sum = [all_sum, sr]; %#ok<AGROW>
            all_d   = [all_d,   dr]; %#ok<AGROW>
        end
        clear Ms Mg ovl C
    end
end

%% Export CSVs
T_sum = struct2table([all_sum{:}]);
T_d   = struct2table([all_d{:}]);
f_sum = fullfile(out_dir, sprintf('histo_violin_bin_summary_%s_%s_%s.csv',  q_tag, run_tag, date_str));
f_d   = fullfile(out_dir, sprintf('histo_violin_cliffs_delta_%s_%s_%s.csv', q_tag, run_tag, date_str));
writetable(T_sum, f_sum);
writetable(T_d,   f_d);

%% Cliff's delta summary (subject level)
[~, nm] = fileparts(f_d);
plot_delta_summary(T_d, strrep(nm, 'histo_violin_cliffs_delta_', 'histo_delta_summary_'), scfg);
fprintf('Done. Outputs in %s\n', out_dir);


%% ======================= Local functions ===============================

function [sum_rows, d_rows] = process_unit(M, C, ovl, eE, eV, meta, cfg)
% One figure per property (here: one stain) for one subject or pooled group.
ntot = size(M,1);
sum_rows = {};
d_rows   = {};
for p = 1:numel(cfg.props)
    col   = cfg.prop_col(p);
    y     = M(:,col);
    cellY = cell(1,4);
    st    = repmat(cell_stats([]), 1, 4);
    for k = 1:4
        cellY{k} = y(C(:,k));
        st(k)    = cell_stats(cellY{k});
    end

    %%% Subject medians (group figures only)
    ov = struct('name',{}, 'med',{}, 'marker',{});
    for s = 1:numel(ovl)
        med = NaN(1,4);
        for k = 1:4
            v = ovl(s).M(ovl(s).C(:,k), col);
            if numel(v) >= cfg.min_n, med(k) = median(v); end
        end
        ov(s) = struct('name',ovl(s).name, 'med',med, 'marker',cfg.markers{s});
    end

    base = sprintf('violin_%s_%s_%s_%s_%s', cfg.props{p}, meta.stub, ...
                   cfg.q_tag, cfg.run_tag, cfg.date_str);

    %%% Cliff's delta, all six pairs
    nc    = size(cfg.contr,1);
    d_all = NaN(nc,1);
    for c = 1:nc
        ka = cfg.contr{c,3}; kb = cfg.contr{c,4};
        if st(ka).n >= cfg.min_n && st(kb).n >= cfg.min_n
            [d_all(c), na_u, nb_u] = cliffs_delta(cellY{ka}, cellY{kb}, cfg.cliff_max_n);
        else
            na_u = 0; nb_u = 0;
        end
        d_rows{end+1} = struct( ...
            'level',meta.level, 'name',meta.name, 'group',meta.group, ...
            'region',meta.region, 'property',cfg.props{p}, ...
            'edge_mode',meta.edge_mode, 'equal_weight',meta.eqw, ...
            'contrast_type',cfg.contr{c,1}, 'held_fixed',cfg.contr{c,2}, ...
            'cell_A',cfg.cell_lab{ka}, 'cell_B',cfg.cell_lab{kb}, ...
            'cliffs_delta',d_all(c), 'n_A',st(ka).n, 'n_B',st(kb).n, ...
            'n_A_used',na_u, 'n_B_used',nb_u, ...
            'frac_A',st(ka).n/ntot, 'frac_B',st(kb).n/ntot, ...
            'on_figure',cfg.annot(c), 'figure',base); %#ok<AGROW>
    end

    %%% Figure
    P = struct();
    P.ylims  = cfg.ylims.(meta.region);
    P.ylabel = cfg.ylab{p};
    P.title  = [ {sprintf('%s  |  %s  |  %s', meta.disp, meta.reg_disp, cfg.prop_disp{p})}, ...
                 reshape(meta.edge_title, 1, []), ...
                 {sprintf('quantiles %.2f / %.2f, %s edges;   \\delta_{Cliff} > 0: right-hand cell larger', ...
                          cfg.Q(1), cfg.Q(2), meta.edge_mode)} ];
    fh = plot_four_violins(cellY, st, d_all, ov, P, cfg);
    save_fig(fh, cfg.out_dir, base);

    %%% Summary rows
    for k = 1:4
        sum_rows{end+1} = struct( ...
            'level',meta.level, 'name',meta.name, 'group',meta.group, ...
            'region',meta.region, 'property',cfg.props{p}, ...
            'edge_mode',meta.edge_mode, 'equal_weight',meta.eqw, ...
            'epvs_bin',cfg.ebin_lab{k}, 'ves_bin',cfg.vbin_lab{k}, ...
            'q_low',cfg.Q(1), 'q_high',cfg.Q(2), ...
            'epvs_edge_low',eE(1), 'epvs_edge_high',eE(2), ...
            'ves_edge_low',eV(1), 'ves_edge_high',eV(2), ...
            'edge_detail',meta.edge_detail, ...
            'n',st(k).n, 'frac_of_voxels',st(k).n/ntot, ...
            'units',cfg.units, ...
            'median',st(k).med, 'q1',st(k).q1, 'q3',st(k).q3, ...
            'figure',base); %#ok<AGROW>
    end
end
end


function fh = plot_four_violins(cellY, st, d_all, ov, P, cfg)
fh = figure('Visible',cfg.fig_visible, 'Color','w', 'Position',[100 100 760 700], ...
            'CreateFcn','set(gcbo,''Visible'',''on'')');
ax = axes(fh); hold(ax,'on');
xp = cfg.xpos;
yl = P.ylims; sp = diff(yl);

for k = 1:4
    col = cfg.vcol(cfg.vbin(k),:);
    if st(k).n >= cfg.min_n
        draw_violin(ax, xp(k), cellY{k}, col, cfg.hw, cfg.kde_max_n, yl);
        plot(ax, [xp(k) xp(k)], [st(k).q1 st(k).q3], 'k-', 'LineWidth', 3);
        plot(ax, xp(k), st(k).med, 'o', 'MarkerFaceColor','w', ...
             'MarkerEdgeColor','k', 'MarkerSize',6);
    else
        text(ax, xp(k), yl(1) + 0.5*sp, sprintf('n < %d', cfg.min_n), ...
             'HorizontalAlignment','center', 'Color',[.4 .4 .4]);
    end
end

%%% Subject median overlays (group figures)
ns = numel(ov);
if ns > 0
    if ns == 1, dx = 0; else, dx = linspace(-0.12, 0.12, ns); end
    hov = gobjects(1, ns);
    for s = 1:ns
        hov(s) = plot(ax, xp + dx(s), ov(s).med, ov(s).marker, 'LineStyle','none', ...
                      'MarkerSize',7, 'MarkerFaceColor',[.25 .25 .25], ...
                      'MarkerEdgeColor','k', 'DisplayName',ov(s).name);
    end
    legend(ax, hov, 'Location','southoutside', 'Orientation','horizontal', ...
           'Box','off', 'AutoUpdate','off');
end

xline(ax, mean(xp(2:3)), ':', 'Color',[.5 .5 .5]);

%%% Cliff's delta brackets, stacked by level
lev  = cell2mat(cfg.contr(:,5));
used = unique(lev(cfg.annot));
h0 = 0.05; dh = 0.10; tk = 0.015;
for c = find(cfg.annot(:))'
    L  = find(used == lev(c));
    yb = yl(2) + (h0 + (L-1)*dh)*sp;
    xa = xp(cfg.contr{c,4});
    xb = xp(cfg.contr{c,3});
    plot(ax, [xa xa xb xb], [yb - tk*sp, yb, yb, yb - tk*sp], 'k-', 'LineWidth', 0.8);
    if isnan(d_all(c)), s = '\delta = n/a'; else, s = sprintf('\\delta = %+.3f', d_all(c)); end
    text(ax, (xa + xb)/2, yb + 0.01*sp, s, 'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', 'FontSize',8.5);
end
ytop_br = yl(2) + (h0 + (numel(used) - 1)*dh + 0.06)*sp;

pairs = {[1 2], [3 4]};
elab  = {'Low EPVS-SWP', 'High EPVS-SWP'};
for g = 1:2
    text(ax, mean(xp(pairs{g})), ytop_br + 0.03*sp, elab{g}, ...
         'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontWeight','bold');
end

ylim(ax, [yl(1), ytop_br + 0.12*sp]);
yt = ax.YTick; ax.YTick = yt(yt <= yl(2));
xlim(ax, [xp(1) - 0.7, xp(4) + 0.7]);
tl = cell(1,4);
for k = 1:4
    tl{k} = sprintf('%s\\newline{\\fontsize{8}n = %s}', cfg.vlab{cfg.vbin(k)}, fmt_n(st(k).n));
end
ax.XTick = xp;
ax.XTickLabel = tl;
ax.TickLabelInterpreter = 'tex';
ylabel(ax, P.ylabel);
title(ax, P.title, 'FontWeight','normal', 'FontSize',10);
set(ax, 'Box','off', 'TickDir','out', 'FontSize',11);
end


function plot_delta_summary(T, base, scfg)
% Subject-level Cliff's delta; rows = stain, columns = frontal / occipital.
T = T(strcmp(T.level, 'subject') & ismember(T.group, scfg.groups), :);
if isempty(T)
    warning('violin:summary', 'No subject-level rows for the delta summary.');
    return
end
has_frac = all(ismember({'frac_A','frac_B'}, T.Properties.VariableNames));

C = { 'EPVS-lo/Ves-hi', 'EPVS-lo/Ves-lo', 'Ves hi vs lo | EPVS lo';
      'EPVS-hi/Ves-hi', 'EPVS-hi/Ves-lo', 'Ves hi vs lo | EPVS hi';
      'EPVS-hi/Ves-lo', 'EPVS-lo/Ves-lo', 'EPVS hi vs lo | Ves lo';
      'EPVS-hi/Ves-hi', 'EPVS-lo/Ves-hi', 'EPVS hi vs lo | Ves hi';
      'EPVS-hi/Ves-lo', 'EPVS-lo/Ves-hi', 'EPVS hi\cdotVes lo vs EPVS lo\cdotVes hi';
      'EPVS-hi/Ves-hi', 'EPVS-lo/Ves-lo', 'Both hi vs both lo' };
nc = size(C,1);

props    = scfg.props(ismember(scfg.props, unique(T.property)));
np       = numel(props);
regions  = {'front','occip'};
reg_disp = {'Frontal','Occipital'};

subs_order = {'CAA26','CAA6','CAA22','CAA25','CAA17'};
names = unique(T.name, 'stable');
subs  = [subs_order(ismember(subs_order, names)), ...
         reshape(setdiff(names, subs_order, 'stable'), 1, [])];
mk    = {'o','s','^','d','v','p','h','>'};
gcol  = struct('control',[0.25 0.25 0.25], 'severe',[0.80 0.15 0.45], ...
               'moderate',[0.45 0.55 0.85]);
if numel(subs) == 1, dy = 0; else, dy = linspace(-0.28, 0.28, numel(subs)); end
scol = zeros(numel(subs), 3); sgrp = cell(1, numel(subs));
for s = 1:numel(subs)
    sgrp{s} = T.group{find(strcmp(T.name, subs{s}), 1)};
    if isfield(gcol, sgrp{s}), scol(s,:) = gcol.(sgrp{s}); else, scol(s,:) = [.5 .5 .5]; end
end

fh = figure('Visible',scfg.fig_visible, 'Color','w', ...
            'Position',[100 100 1000 260*np + 160], ...
            'CreateFcn','set(gcbo,''Visible'',''on'')');
tl = tiledlayout(fh, np, 2, 'TileSpacing','compact', 'Padding','compact');
ylim_ = [0.4, nc + 0.6];
for p = 1:np
    for r = 1:2
        ax = nexttile(tl); hold(ax, 'on');
        patch(ax, [-0.147 0.147 0.147 -0.147], ylim_([1 1 2 2]), [0.92 0.92 0.92], ...
              'EdgeColor','none');
        for xr = [-0.474 -0.33 0.33 0.474]
            plot(ax, [xr xr], ylim_, ':', 'Color',[0.6 0.6 0.6]);
        end
        plot(ax, [0 0], ylim_, 'k-', 'LineWidth', 0.8);
        for yb = [2.5 4.5]
            plot(ax, [-1 1], [yb yb], '-', 'Color',[0.85 0.85 0.85]);
        end

        sel = strcmp(T.property, props{p}) & strcmp(T.region, regions{r});
        for s = 1:numel(subs)
            rs = sel & strcmp(T.name, subs{s});
            if ~any(rs), continue; end
            for c = 1:nc
                i = find(rs & strcmp(T.cell_A, C{c,1}) & strcmp(T.cell_B, C{c,2}), 1);
                if isempty(i) || isnan(T.cliffs_delta(i)), continue; end
                is_sparse = has_frac && min(T.frac_A(i), T.frac_B(i)) < scfg.sparse_frac;
                if is_sparse, fc = 'w'; else, fc = scol(s,:); end
                plot(ax, T.cliffs_delta(i), c + dy(s), mk{s}, 'MarkerSize',7, ...
                     'MarkerEdgeColor',scol(s,:), 'MarkerFaceColor',fc, 'LineWidth',1);
            end
        end

        set(ax, 'YDir','reverse', 'YLim',ylim_, 'XLim',[-1 1], 'YTick',1:nc, ...
                'Box','off', 'TickDir','out', 'FontSize',10, 'TickLabelInterpreter','tex');
        if r == 1, ax.YTickLabel = C(:,3); else, ax.YTickLabel = {}; end
        title(ax, sprintf('%s  |  %s', reg_disp{r}, upper(props{p})), 'FontWeight','normal');
        if p == np, xlabel(ax, 'Cliff''s \delta   (> 0: first-named cell larger)'); end
    end
end

hl = gobjects(0); lab = {};
for s = 1:numel(subs)
    hl(end+1) = plot(ax, NaN, NaN, mk{s}, 'LineStyle','none', 'MarkerSize',7, ...
                     'MarkerEdgeColor',scol(s,:), 'MarkerFaceColor',scol(s,:)); %#ok<AGROW>
    lab{end+1} = sprintf('%s (%s)', subs{s}, sgrp{s}); %#ok<AGROW>
end
if has_frac
    hl(end+1) = plot(ax, NaN, NaN, 'o', 'LineStyle','none', 'MarkerSize',7, ...
                     'MarkerEdgeColor',[.3 .3 .3], 'MarkerFaceColor','w');
    lab{end+1} = sprintf('cell < %g%% of voxels', 100*scfg.sparse_frac);
end
hl(end+1) = patch(ax, NaN, NaN, [0.92 0.92 0.92], 'EdgeColor','none');
lab{end+1} = '|\delta| < 0.147';
lg = legend(ax, hl, lab, 'Orientation','horizontal', 'Box','off', 'AutoUpdate','off');
lg.Layout.Tile = 'south';

title(tl, 'Histology: subject-level Cliff''s \delta, per-subject quantile edges', 'FontSize',12);
save_fig(fh, scfg.out_dir, base);
end


function draw_violin(ax, xc, y, col, hw, kde_max_n, yl)
if numel(y) > kde_max_n
    y = y(randperm(numel(y), kde_max_n));
end
lims = prctile(y, [0.5 99.5]);
lims = [max(lims(1), yl(1)), min(lims(2), yl(2))];
if lims(2) <= lims(1)
    plot(ax, xc + [-hw hw], lims(1)*[1 1], 'Color', col, 'LineWidth', 2);
    return
end
yi = linspace(lims(1), lims(2), 200);
f  = ksdensity(y, yi);
f  = f ./ max(f) * hw;
patch(ax, [xc - f, fliplr(xc + f)], [yi, fliplr(yi)], col, ...
      'FaceAlpha',0.55, 'EdgeColor',col*0.7, 'LineWidth',0.8);
end


function [d, n1, n2] = cliffs_delta(a, b, cap)
% delta = P(a > b) - P(a < b), exact via ranks (ties count 1/2).
a = a(:); b = b(:);
if numel(a) > cap, a = a(randperm(numel(a), cap)); end
if numel(b) > cap, b = b(randperm(numel(b), cap)); end
n1 = numel(a); n2 = numel(b);
if n1 == 0 || n2 == 0, d = NaN; return; end
r = tiedrank([a; b]);
d = 2 * (mean(r(1:n1)) - (n1 + 1)/2) / n2 - 1;
end


function C = cell_masks(epv, ves, eE, eV)
El = epv <= eE(1);  Eh = epv >= eE(2);
Vl = ves <= eV(1);  Vh = ves >= eV(2);
C  = [El & Vl, El & Vh, Eh & Vl, Eh & Vh];
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


function x = clean_col(x)
x = x(isfinite(x));
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


function save_fig(fh, out_dir, base)
exportgraphics(fh, fullfile(out_dir, [base '.png']), 'Resolution', 300);
exportgraphics(fh, fullfile(out_dir, [base '.svg']), 'ContentType', 'vector');
savefig(fh, fullfile(out_dir, [base '.fig']));
close(fh);
end