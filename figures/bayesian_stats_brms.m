%% Figure 3 c-e: optical properties vs distance from the brms OCT results
%{
Reads the CSVs written by fit_parenchyma_shells_brms.R and draws:
  c  scattering (µs, Δ%) vs distance: EPVS and vessel, control and severe
  d  retardance (Δ%) vs distance: same layout
  e  severe - control, per class, for µs and retardance
Panels a/b (segmentation images) are not drawn.

Exports (each as png / svg / fig, set by cfg.formats):
  <out_dir>/<out_name>_main.*            full figure (c, d, e)
  <out_dir>/rows/<out_name>_<c|d|e>.*    one figure per lettered row
  <out_dir>/panels/<out_name>_<letter>_<region>_<stage or property>.*
                                         one figure per individual plot (12)

Lines = posterior mean; ribbons = 95% credible interval.
Hollow points (optional) = distances supported by < min_objects_flag objects.
%}

clear; clc; close all;

%% Settings
cfg.beta_root    = '/projectnb/npbssmic/ns/CAA/beta_stats';
cfg.base_name    = 'parenchyma_shells_2026-10-07_pct';
cfg.analysis_tag = 'main_vsub500';     % folder under brms/ (e.g. diam_covariate_vsub500)
cfg.scale        = 'pct';
cfg.model_set    = 'separate';         % 'separate' (group/epvs + group/vessel) or 'joint' (group/joint)

cfg.regions      = {'front', 'occip'};
cfg.region_names = {'Frontal', 'Occipital'};
cfg.stages       = {'control', 'severe'};
cfg.contrast     = 'severe - control';

% Axis limits ([] = automatic)
cfg.xlim                = [0 540];
cfg.ylim.c              = [-5 30];
cfg.ylim.d              = [-15 10];
cfg.ylim.e.scattering   = [-25 15];
cfg.ylim.e.retardance   = [-15 10];

cfg.colors.EPVS   = [0.12 0.42 0.75];
cfg.colors.vessel = [0.85 0.33 0.10];
cfg.ribbon_alpha  = 0.2;
cfg.show_low_support = true;
cfg.font      = 'Helvetica';
cfg.font_size = 10;

% Export
cfg.out_dir  = fullfile(cfg.beta_root, cfg.base_name, 'brms', cfg.analysis_tag, 'figures');
cfg.out_name = sprintf('scatterplot_%s_%s', cfg.scale, cfg.model_set);
cfg.formats  = {'png', 'svg', 'fig'};  % any of 'png', 'svg', 'pdf', 'fig'
cfg.dpi      = 300;
cfg.export_main   = true;
cfg.export_rows   = true;              % c, d, e as separate figures
cfg.export_panels = true;              % each of the 12 plots on its own
cfg.size_main  = [36 30];              % figure sizes (cm)
cfg.size_row   = [36 10];
cfg.size_panel = [9 7.5];
cfg.legend_on_every_panel = true;      % legend on each single-panel figure

% Supplemental figure: moderate stage (CAA17; occipital volume only)
cfg.export_supp       = true;
cfg.supp.region       = 'occip';
cfg.supp.stages       = {'control', 'moderate', 'severe'};   % moderate shown beside control/severe
cfg.supp.stage_titles = struct('control', 'Control', 'moderate', 'Moderate (CAA17)', 'severe', 'Severe');
cfg.supp.contrast     = 'moderate - control';
cfg.supp.letters      = {'a', 'b', 'c'};                      % µs row, retardance row, difference row
cfg.supp.out_name     = sprintf('supp_moderate_CAA17_%s_%s', cfg.scale, cfg.model_set);
cfg.supp.size_main    = [27 30];
cfg.supp.size_row     = [27 10];

% Supplemental: subject-level fits (one joint EPVS + vessel model per subject)
cfg.export_subjects    = true;
cfg.subj.stage_order   = {'control', 'moderate', 'severe'};   % column order of subjects
cfg.subj.out_name      = sprintf('supp_subjects_%s', cfg.scale);
cfg.subj.cm_per_plot   = [7 7];          % size of each plot in the multi-panel figures
cfg.subj.export_ev     = true;           % extra figure: EPVS - vessel for every subject
cfg.subj.ylim_ev       = [];             % [] = automatic
cfg.subj.stage_colors  = struct('control',  [0.35 0.35 0.35], ...
                                'moderate', [0.93 0.55 0.00], ...
                                'severe',   [0.70 0.13 0.13]);

%% Paths, labels, data
root  = fullfile(cfg.beta_root, cfg.base_name, 'brms', cfg.analysis_tag, cfg.scale);
props = {'scattering', 'retardance'};
lab.ylab  = struct('scattering', '\mu_s (\Delta%)', 'retardance', 'Retardance (\Delta%)');
lab.pname = struct('scattering', '\mu_s', 'retardance', 'Retardance');

% Load every CSV once; all figures draw from this
data = struct();
for p = 1:numel(props)
    for g = 1:numel(cfg.regions)
        data.fitted.(props{p}).(cfg.regions{g}) = load_results(root, cfg, props{p}, cfg.regions{g}, 'fitted');
        data.diff.(props{p}).(cfg.regions{g})   = load_results(root, cfg, props{p}, cfg.regions{g}, 'diff');
        data.subj.(props{p}).(cfg.regions{g})    = load_subject_results(root, cfg, props{p}, cfg.regions{g}, 'fitted');
        data.subj_ev.(props{p}).(cfg.regions{g}) = load_subject_results(root, cfg, props{p}, cfg.regions{g}, 'ev');
    end
end

rows = {'c', 'd', 'e'};
if ~exist(cfg.out_dir, 'dir'), mkdir(cfg.out_dir); end

%% Main figure
if cfg.export_main
    fig   = new_figure(cfg.size_main);
    outer = tiledlayout(fig, 3, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
    for r = 1:3
        rowT = tiledlayout(outer, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        rowT.Layout.Tile = r;
        draw_row(rowT, rows{r}, data, cfg, lab);
    end
    export_all(fig, fullfile(cfg.out_dir, [cfg.out_name '_main']), cfg);
end

%% One figure per lettered row
if cfg.export_rows
    rdir = fullfile(cfg.out_dir, 'rows');
    if ~exist(rdir, 'dir'), mkdir(rdir); end
    for r = 1:3
        fig  = new_figure(cfg.size_row);
        rowT = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        draw_row(rowT, rows{r}, data, cfg, lab);
        export_all(fig, fullfile(rdir, sprintf('%s_%s', cfg.out_name, rows{r})), cfg);
        close(fig);
    end
end

%% One figure per individual plot
if cfg.export_panels
    pdir = fullfile(cfg.out_dir, 'panels');
    if ~exist(pdir, 'dir'), mkdir(pdir); end
    for r = 1:3
        for g = 1:numel(cfg.regions)
            rgn = cfg.regions{g};
            for j = 1:2
                fig = new_figure(cfg.size_panel);
                ax  = axes(fig);
                [~, name] = draw_axes(ax, rows{r}, g, j, data, cfg, lab, ...
                                      cfg.legend_on_every_panel, true);
                export_all(fig, fullfile(pdir, sprintf('%s_%s_%s_%s', ...
                           cfg.out_name, rows{r}, rgn, name)), cfg);
                close(fig);
            end
        end
    end
end

%% Supplemental figure: moderate stage (CAA17), occipital
if cfg.export_supp
    supp_props = props;
    supp_ylim  = struct('scattering', cfg.ylim.c, 'retardance', cfg.ylim.d);
    supp_dylim = cfg.ylim.e;
    supp_label = @(prop, is_diff) lab.ylab.(prop);
    supp_name  = @(prop) lab.pname.(prop);
    run_supplemental(data, cfg, supp_props, supp_ylim, supp_dylim, supp_label, supp_name);
end

%% Supplemental figures: subject-level fits
if cfg.export_subjects
    S_subj = struct('props', {props}, ...
                    'ylim_curve', struct('scattering', cfg.ylim.c, 'retardance', cfg.ylim.d), ...
                    'label_fn', @(prop) lab.ylab.(prop), ...
                    'name_fn',  @(prop) lab.pname.(prop));
    run_subject_supplemental(data, cfg, S_subj);
end

fprintf('Done. Figures in %s\n', cfg.out_dir);

%% ========================================================================
%% Drawing
%% ========================================================================
function draw_row(rowT, letter, data, cfg, lab)
% One lettered row: two regions, each with two plots
for g = 1:numel(cfg.regions)
    regT = tiledlayout(rowT, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    regT.Layout.Tile = g;
    if letter == 'e'
        ttl = {cfg.region_names{g}, ['(' cfg.contrast ')']};
    else
        ttl = cfg.region_names{g};
    end
    title(regT, ttl, 'FontName', cfg.font, 'FontSize', cfg.font_size + 2);
    for j = 1:2
        ax = nexttile(regT);
        show_legend = (g == 2 && j == 1);
        draw_axes(ax, letter, g, j, data, cfg, lab, show_legend, false);
        if g == 1 && j == 1, panel_letter(ax, letter, cfg); end
    end
end
end

function [h, name] = draw_axes(ax, letter, g, j, data, cfg, lab, show_legend, standalone)
% Draw one plot. letter c/d: j = stage index; letter e: j = property index.
% standalone = true adds the region to the title (single-panel figures).
rgn = cfg.regions{g};
switch letter
    case {'c', 'd'}
        prop = ternary(letter == 'c', 'scattering', 'retardance');
        stg  = cfg.stages{j};
        T    = data.fitted.(prop).(rgn);
        if ~isempty(T), T = T(T.stage == stg, :); end
        yl   = cfg.ylim.(letter);
        ttl  = [upper(stg(1)) stg(2:end)];
        name = stg;
        if standalone, ttl = sprintf('%s: %s', cfg.region_names{g}, ttl); end
    case 'e'
        prop = ternary(j == 1, 'scattering', 'retardance');
        T    = data.diff.(prop).(rgn);
        if ~isempty(T), T = T(T.contrast == cfg.contrast, :); end
        yl   = cfg.ylim.e.(prop);
        ttl  = lab.pname.(prop);
        name = prop;
        if standalone
            ttl = sprintf('%s (%s): %s', cfg.region_names{g}, cfg.contrast, lab.pname.(prop));
        end
end
h = [];
if ~isempty(T), h = plot_classes(ax, T, cfg); end
style_axes(ax, cfg, yl, lab.ylab.(prop));
title(ax, ttl, 'FontWeight', 'normal');
if show_legend, add_legend(ax, h); end
end

function h = plot_classes(ax, T, cfg)
% Ribbon + line for EPVS and vessel; returns line handles [EPVS, vessel]
hold(ax, 'on');
classes = {'EPVS', 'vessel'};
h = gobjects(1, 2);
for c = 1:2
    cls = classes{c};
    S = sortrows(T(T.class == cls, :), 'distance');
    if isempty(S), continue; end
    col = cfg.colors.(cls);
    x = S.distance; m = S.mean; lo = S.lower; hi = S.upper;
    fill(ax, [x; flipud(x)], [lo; flipud(hi)], col, 'FaceAlpha', cfg.ribbon_alpha, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
    h(c) = plot(ax, x, m, '-', 'Color', col, 'LineWidth', 1.5, ...
                'Marker', '.', 'MarkerSize', 12);
    if cfg.show_low_support && ismember('low_support', S.Properties.VariableNames)
        ls = S.low_support;
        if any(ls)
            plot(ax, x(ls), m(ls), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', col, ...
                 'MarkerFaceColor', 'w', 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
end
end

function style_axes(ax, cfg, yl, ylabel_str)
yline(ax, 0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
grid(ax, 'on'); box(ax, 'on');
if ~isempty(cfg.xlim), xlim(ax, cfg.xlim); end
if ~isempty(yl), ylim(ax, yl); end
xticks(ax, 0:100:500);
xlabel(ax, 'Distance (\mum)');
ylabel(ax, ylabel_str);
set(ax, 'FontName', cfg.font, 'FontSize', cfg.font_size, 'Layer', 'top', ...
        'TickDir', 'out', 'LineWidth', 0.75);
end

function add_legend(ax, h)
if isempty(h) || ~all(isgraphics(h)), return; end
legend(ax, h, {'EPVS', 'Vessel'}, 'Location', 'best', 'Box', 'on');
end

function panel_letter(ax, letter, cfg)
text(ax, -0.32, 1.22, letter, 'Units', 'normalized', 'FontName', cfg.font, ...
     'FontSize', cfg.font_size + 6, 'FontWeight', 'bold', 'Clipping', 'off');
end

%% ========================================================================
%% Figures and export
%% ========================================================================
function fig = new_figure(size_cm)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 size_cm]);
end

function export_all(fig, base, cfg)
base = char(base);          % string paths would turn [base '.png'] into an array
drawnow;
for k = 1:numel(cfg.formats)
    switch lower(cfg.formats{k})
        case 'png', exportgraphics(fig, [base '.png'], 'Resolution', cfg.dpi);
        case 'svg', exportgraphics(fig, [base '.svg'], 'ContentType', 'vector');
        case 'pdf', exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
        case 'fig', savefig(fig, [base '.fig']);
        otherwise,  warning('Unknown format: %s', cfg.formats{k});
    end
end
fprintf('Saved %s (%s)\n', base, strjoin(cfg.formats, ', '));
end

%% ========================================================================
%% Data
%% ========================================================================
function T = load_results(root, cfg, prop, rgn, kind)
% kind = 'fitted' -> columns stage, class, distance, mean, lower, upper, low_support
% kind = 'diff'   -> columns class, contrast, distance, mean, lower, upper
switch cfg.model_set
    case 'separate', srcs = {'epvs', 'vessel'};
    case 'joint',    srcs = {'joint'};
    otherwise, error('cfg.model_set must be ''separate'' or ''joint''');
end
switch kind
    case 'fitted'
        suffix = '_all_stages_fitted.csv';
        keep   = {'stage', 'class', 'distance', 'mean', 'lower', 'upper', 'low_support'};
    case 'diff'
        suffix = '_all_stages_model_diff.csv';
        keep   = {'class', 'contrast', 'distance', 'mean', 'lower', 'upper'};
end

parts = {};
for k = 1:numel(srcs)
    stem = sprintf('%s_%s_%s_%s', cfg.scale, prop, rgn, srcs{k});
    f = fullfile(root, prop, rgn, 'group', srcs{k}, 'all_stages', [stem suffix]);
    if ~isfile(f)
        warning('Missing file, panel left empty: %s', f);
        continue
    end
    Ti = readtable(f, 'TextType', 'string');
    for v = setdiff(keep, Ti.Properties.VariableNames)   % e.g. low_support absent
        Ti.(v{1}) = false(height(Ti), 1);
    end
    Ti = Ti(:, keep);
    if ismember('low_support', keep), Ti.low_support = to_logical(Ti.low_support); end
    parts{end+1} = Ti; %#ok<AGROW>
end
if isempty(parts), T = table(); else, T = vertcat(parts{:}); end
for v = intersect({'stage', 'class', 'contrast'}, keep)
    if ~isempty(T), T.(v{1}) = string(T.(v{1})); end
end
end

function b = to_logical(x)
if islogical(x), b = x; return; end
if isnumeric(x), b = x ~= 0; return; end
b = ismember(lower(string(x)), ["true", "1"]);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

%% ========================================================================
%% Supplemental figure (moderate stage)
%% ========================================================================
function run_supplemental(data, cfg, props, ylim_curve, ylim_diff, label_fn, name_fn)
% Rows 1..numel(props): control | moderate | severe curves for one property
% Last row: moderate - control, one plot per property
sp  = cfg.supp;
rgn = sp.region;
has_mod = false;
for p = 1:numel(props)
    T = data.fitted.(props{p}).(rgn);
    has_mod = has_mod || (~isempty(T) && any(T.stage == "moderate"));
end
if ~has_mod
    warning('No moderate-stage results for %s: supplemental figure skipped.', rgn);
    return
end
S = struct('props', {props}, 'ylim_curve', ylim_curve, 'ylim_diff', ylim_diff, ...
           'label_fn', label_fn, 'name_fn', name_fn, ...
           'region_name', cfg.region_names{strcmp(cfg.regions, rgn)});
sdir   = fullfile(cfg.out_dir, 'supplemental');
if ~exist(sdir, 'dir'), mkdir(sdir); end
n_rows = numel(props) + 1;

if cfg.export_main
    fig   = new_figure(sp.size_main);
    outer = tiledlayout(fig, n_rows, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
    title(outer, sprintf('%s: moderate CAA (CAA17)', S.region_name), ...
          'FontName', cfg.font, 'FontSize', cfg.font_size + 3);
    for r = 1:n_rows
        rowT = tiledlayout(outer, 1, supp_n_cols(r, S, sp), 'TileSpacing', 'compact', 'Padding', 'compact');
        rowT.Layout.Tile = r;
        supp_draw_row(rowT, r, data, cfg, S);
    end
    export_all(fig, fullfile(sdir, [sp.out_name '_main']), cfg);
end
if cfg.export_rows
    for r = 1:n_rows
        fig  = new_figure(sp.size_row);
        rowT = tiledlayout(fig, 1, supp_n_cols(r, S, sp), 'TileSpacing', 'compact', 'Padding', 'compact');
        supp_draw_row(rowT, r, data, cfg, S);
        export_all(fig, fullfile(sdir, sprintf('%s_%s', sp.out_name, sp.letters{r})), cfg);
        close(fig);
    end
end
if cfg.export_panels
    for r = 1:n_rows
        for j = 1:supp_n_cols(r, S, sp)
            fig = new_figure(cfg.size_panel);
            ax  = axes(fig);
            name = supp_draw_axes(ax, r, j, data, cfg, S, cfg.legend_on_every_panel, true);
            export_all(fig, fullfile(sdir, sprintf('%s_%s_%s', sp.out_name, sp.letters{r}, name)), cfg);
            close(fig);
        end
    end
end
end

function n = supp_n_cols(r, S, sp)
if r <= numel(S.props), n = numel(sp.stages); else, n = numel(S.props); end
end

function supp_draw_row(rowT, r, data, cfg, S)
sp = cfg.supp;
for j = 1:supp_n_cols(r, S, sp)
    ax = nexttile(rowT);
    supp_draw_axes(ax, r, j, data, cfg, S, (r == 1 && j == 1), false);
    if j == 1, panel_letter(ax, sp.letters{r}, cfg); end
end
end

function name = supp_draw_axes(ax, r, j, data, cfg, S, show_legend, standalone)
sp  = cfg.supp;
rgn = sp.region;
if r <= numel(S.props)
    prop = S.props{r};
    stg  = sp.stages{j};
    T    = data.fitted.(prop).(rgn);
    if ~isempty(T), T = T(T.stage == stg, :); end
    yl   = S.ylim_curve.(prop);
    ylab = S.label_fn(prop, false);
    ttl  = sp.stage_titles.(stg);
    if standalone, ttl = sprintf('%s, %s: %s', S.region_name, S.name_fn(prop), ttl); end
    name = sprintf('%s_%s', prop, stg);
else
    prop = S.props{j};
    T    = data.diff.(prop).(rgn);
    if ~isempty(T), T = T(T.contrast == sp.contrast, :); end
    yl   = S.ylim_diff.(prop);
    ylab = S.label_fn(prop, true);
    ttl  = sprintf('%s (%s)', S.name_fn(prop), sp.contrast);
    if standalone, ttl = sprintf('%s: %s', S.region_name, ttl); end
    name = sprintf('%s_moderate_minus_control', prop);
end
h = [];
if ~isempty(T), h = plot_classes(ax, T, cfg); end
style_axes(ax, cfg, yl, ylab);
title(ax, ttl, 'FontWeight', 'normal');
if show_legend, add_legend(ax, h); end
end

%% ========================================================================
%% Supplemental: subject-level fits
%% ========================================================================
function run_subject_supplemental(data, cfg, S)
% Per property: rows = regions, columns = subjects (ordered by stage).
% Each plot = that subject's joint-model EPVS and vessel curves.
% Optional extra figure: EPVS - vessel for every subject, rows = properties.
sp   = cfg.subj;
sdir = fullfile(cfg.out_dir, 'supplemental', 'subjects');
if ~exist(sdir, 'dir'), mkdir(sdir); end
nR = numel(cfg.regions);

for p = 1:numel(S.props)
    prop = S.props{p};
    subj = subject_order(data.subj.(prop), cfg);
    if isempty(subj)
        warning('%s: no subject-level results found (was run_subject = TRUE?), skipped.', prop);
        continue
    end
    nS   = height(subj);
    name = sprintf('%s_%s', sp.out_name, prop);

    if cfg.export_main
        fig   = new_figure([sp.cm_per_plot(1) * nS + 2, sp.cm_per_plot(2) * nR + 2]);
        outer = tiledlayout(fig, nR, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
        title(outer, sprintf('%s: subject-level fits', S.name_fn(prop)), 'FontName', cfg.font, ...
              'FontSize', cfg.font_size + 3);
        for g = 1:nR
            rowT = tiledlayout(outer, 1, nS, 'TileSpacing', 'compact', 'Padding', 'compact');
            rowT.Layout.Tile = g;
            subject_row(rowT, prop, g, subj, data, cfg, S);
        end
        export_all(fig, fullfile(sdir, [name '_main']), cfg);
    end
    if cfg.export_rows
        for g = 1:nR
            fig  = new_figure([sp.cm_per_plot(1) * nS + 2, sp.cm_per_plot(2) + 2]);
            rowT = tiledlayout(fig, 1, nS, 'TileSpacing', 'compact', 'Padding', 'compact');
            subject_row(rowT, prop, g, subj, data, cfg, S);
            export_all(fig, fullfile(sdir, sprintf('%s_%s', name, cfg.regions{g})), cfg);
            close(fig);
        end
    end
    if cfg.export_panels
        for g = 1:nR
            for k = 1:nS
                T = subject_rows(data.subj.(prop).(cfg.regions{g}), subj.subjectID(k));
                if isempty(T), continue; end           % e.g. CAA17 has no frontal volume
                fig = new_figure(cfg.size_panel);
                ax  = axes(fig);
                subject_axes(ax, T, prop, g, subj(k, :), cfg, S, cfg.legend_on_every_panel, true);
                export_all(fig, fullfile(sdir, sprintf('%s_%s_%s', name, ...
                           cfg.regions{g}, subj.subjectID(k))), cfg);
                close(fig);
            end
        end
    end
end

if sp.export_ev
    nP    = numel(S.props);
    fig   = new_figure([sp.cm_per_plot(1) * nR * 1.4 + 2, sp.cm_per_plot(2) * nP + 2]);
    outer = tiledlayout(fig, nP, nR, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(outer, 'EPVS - vessel by subject (subject-level joint models)', ...
          'FontName', cfg.font, 'FontSize', cfg.font_size + 3);
    for p = 1:nP
        for g = 1:nR
            ax = nexttile(outer);
            ev_axes(ax, data.subj_ev.(S.props{p}).(cfg.regions{g}), S.props{p}, g, cfg, S);
        end
    end
    export_all(fig, fullfile(sdir, [sp.out_name '_EPVS_minus_vessel']), cfg);
end
end

function subj = subject_order(subj_by_region, cfg)
% Unique subjects across regions, ordered by stage then ID
parts = {};
for g = 1:numel(cfg.regions)
    T = subj_by_region.(cfg.regions{g});
    if ~isempty(T), parts{end+1} = unique(T(:, {'subjectID', 'stage'}), 'rows'); end %#ok<AGROW>
end
if isempty(parts), subj = table(); return; end
subj = unique(vertcat(parts{:}), 'rows');
[~, so] = ismember(subj.stage, string(cfg.subj.stage_order));
so(so == 0) = numel(cfg.subj.stage_order) + 1;
subj.order = so;
subj = sortrows(subj, {'order', 'subjectID'});
end

function T = subject_rows(T, id)
if ~isempty(T), T = T(T.subjectID == id, :); end
end

function subject_row(rowT, prop, g, subj, data, cfg, S)
title(rowT, cfg.region_names{g}, 'FontName', cfg.font, 'FontSize', cfg.font_size + 2);
for k = 1:height(subj)
    ax = nexttile(rowT);
    T  = subject_rows(data.subj.(prop).(cfg.regions{g}), subj.subjectID(k));
    if isempty(T)                                  % subject not imaged in this region
        axis(ax, 'off');
        text(ax, 0.5, 0.5, sprintf('%s\nno %s volume', upper(subj.subjectID(k)), ...
             lower(cfg.region_names{g})), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
             'FontName', cfg.font, 'FontSize', cfg.font_size, 'Color', [0.5 0.5 0.5]);
        continue
    end
    subject_axes(ax, T, prop, g, subj(k, :), cfg, S, (k == 1), false);
end
end

function subject_axes(ax, T, prop, g, subj_row, cfg, S, show_legend, standalone)
h = plot_classes(ax, T, cfg);
style_axes(ax, cfg, S.ylim_curve.(prop), S.label_fn(prop));
ttl = sprintf('%s (%s)', upper(subj_row.subjectID), subj_row.stage);
if standalone, ttl = sprintf('%s, %s: %s', S.name_fn(prop), cfg.region_names{g}, ttl); end
title(ax, ttl, 'FontWeight', 'normal');
if show_legend, add_legend(ax, h); end
end

function ev_axes(ax, T, prop, g, cfg, S)
% One line + ribbon per subject; colour = stage, line style separates subjects
hold(ax, 'on');
styles = {'-', '--', ':', '-.'};
h = gobjects(0); lbl = {};
if ~isempty(T)
    ids = unique(T.subjectID, 'stable');
    stage_count = struct();
    for k = 1:numel(ids)
        D   = sortrows(T(T.subjectID == ids(k), :), 'distance');
        stg = char(D.stage(1));
        if isfield(stage_count, stg), stage_count.(stg) = stage_count.(stg) + 1;
        else, stage_count.(stg) = 1; end
        col = [0.5 0.5 0.5];
        if isfield(cfg.subj.stage_colors, stg), col = cfg.subj.stage_colors.(stg); end
        x = D.distance;
        fill(ax, [x; flipud(x)], [D.lower; flipud(D.upper)], col, 'FaceAlpha', 0.12, ...
             'EdgeColor', 'none', 'HandleVisibility', 'off');
        h(end+1) = plot(ax, x, D.mean, styles{mod(stage_count.(stg) - 1, numel(styles)) + 1}, ...
                        'Color', col, 'LineWidth', 1.5); %#ok<AGROW>
        lbl{end+1} = sprintf('%s (%s)', upper(ids(k)), stg); %#ok<AGROW>
    end
end
style_axes(ax, cfg, cfg.subj.ylim_ev, ['EPVS - vessel, ' S.label_fn(prop)]);
title(ax, sprintf('%s, %s', S.name_fn(prop), cfg.region_names{g}), 'FontWeight', 'normal');
if ~isempty(h), legend(ax, h, lbl, 'Location', 'best', 'Box', 'on', 'FontSize', cfg.font_size - 2); end
end

function T = load_subject_results(root, cfg, prop, rgn, kind)
% R layout: <root>/<prop>/<region>/subject/all_subjects/<scale>_<prop>_<region>_all_subjects_*.csv
switch kind
    case 'fitted'
        suffix = '_fitted.csv';
        keep   = {'subjectID', 'stage', 'class', 'distance', 'mean', 'lower', 'upper', 'low_support'};
    case 'ev'
        suffix = '_EPVS_minus_vessel.csv';
        keep   = {'subjectID', 'stage', 'distance', 'mean', 'lower', 'upper'};
end
f = fullfile(root, prop, rgn, 'subject', 'all_subjects', ...
             sprintf('%s_%s_%s_all_subjects%s', cfg.scale, prop, rgn, suffix));
if ~isfile(f)
    warning('Missing subject-level file: %s', f);
    T = table(); return
end
T = readtable(f, 'TextType', 'string');
for v = setdiff(keep, T.Properties.VariableNames)
    T.(v{1}) = false(height(T), 1);
end
T = T(:, keep);
if ismember('low_support', keep), T.low_support = to_logical(T.low_support); end
for v = intersect({'subjectID', 'stage', 'class'}, keep)
    T.(v{1}) = string(T.(v{1}));
end
end