%% Histology figure: stain intensity vs distance from the brms results
%{
Same layout and exports as plot_fig3_shells.m (OCT), for the histology
models fitted by fit_histo_shells_brms.R on histo_shells_<date>.csv.

Rows:
  one row per stain (LHE, CD68, GFAP): EPVS and vessel curves,
      frontal control | frontal severe | occipital control | occipital severe
  last row: severe - control per class, one plot per stain, per region

Exports (cfg.formats):
  <out_dir>/<out_name>_main.*                     full figure
  <out_dir>/rows/<out_name>_<letter>.*            one figure per row
  <out_dir>/panels/<out_name>_<letter>_<region>_<stage or stain>.*
                                                  one figure per plot

Lines = posterior mean; ribbons = 95% credible interval.
Hollow points (optional) = distances supported by < min_objects_flag objects.
%}

clear; clc; close all;

%% Settings
cfg.hist_root    = '/projectnb/npbssmic/ns/CAA/histology';
cfg.input_csv    = 'histo_shells_2026-10-08.csv';   % CSV the R models were fitted on
cfg.analysis_tag = 'main_vsub500';                  % folder under brms/
cfg.scale        = 'z';                             % 'z', 'hmatch' or 'raw' (as fitted)
cfg.model_set    = 'separate';                      % 'separate' or 'joint'
cfg.test_mode    = false;                           % true = plot the brms_TEST/ results

cfg.stains       = {'LHE', 'CD68', 'GFAP'};
cfg.regions      = {'front', 'occip'};
cfg.region_names = {'Frontal', 'Occipital'};
cfg.stages       = {'control', 'severe'};
cfg.contrast     = 'severe - control';
cfg.letters      = {'a', 'b', 'c', 'd'};            % one per stain row, then the difference row

% Axis limits ([] = automatic); per stain for curve rows and difference row
cfg.xlim = [0 540];

% Y-axis limits, set separately for curve plots and difference plots.
% Modes:
%   'fixed' : use the per-stain limits below
%   'data'  : computed from the plotted curves and their 95% CrI
%             (+ padding, rounded outward to round_step), same for all
%             panels of a stain, so main/supplemental/subject plots match
%   'auto'  : MATLAB default, each plot scaled on its own
cfg.ylim.curve_mode = 'fixed';
cfg.ylim.curve = struct( ...
    'LHE',  [-0.30 0.40], ...
    'CD68', [-0.15 0.10], ...
    'GFAP', [-0.40 0.60]);

cfg.ylim.diff_mode = 'data';
cfg.ylim.diff = struct( ...          % used when diff_mode = 'fixed'
    'LHE',  [-0.30 0.30], ...
    'CD68', [-0.15 0.15], ...
    'GFAP', [-0.40 0.40]);

% Options for 'data' mode
cfg.ylim.pad_frac   = 0.05;          % padding as a fraction of the data range
cfg.ylim.round_step = 0.05;          % round limits outward to this step ([] = no rounding)
cfg.ylim.include_subjects = true;    % curve limits also cover the subject-level curves

cfg.colors.EPVS   = [0.12 0.42 0.75];
cfg.colors.vessel = [0.85 0.33 0.10];
cfg.ribbon_alpha  = 0.2;
cfg.show_low_support = true;
cfg.font      = 'Helvetica';
cfg.font_size = 10;

% Export
[~, base_name] = fileparts(cfg.input_csv);
brms_folder = 'brms';
if cfg.test_mode, brms_folder = 'brms_TEST'; end
cfg.results_dir = fullfile(cfg.hist_root, 'beta_stats', base_name, brms_folder, cfg.analysis_tag);
cfg.out_dir  = fullfile(cfg.results_dir, 'figures');
cfg.out_name = sprintf('histo_%s_%s', cfg.scale, cfg.model_set);
cfg.formats  = {'png', 'svg', 'fig'};
cfg.dpi      = 300;
cfg.export_main   = true;
cfg.export_rows   = true;
cfg.export_panels = true;
cfg.size_main  = [36 10 * (numel(cfg.stains) + 1)];
cfg.size_row   = [36 10];
cfg.size_panel = [9 7.5];
cfg.legend_on_every_panel = true;

% Supplemental figure: moderate stage (CAA17; occipital only)
cfg.export_supp       = true;
cfg.supp.region       = 'occip';
cfg.supp.stages       = {'control', 'moderate', 'severe'};
cfg.supp.stage_titles = struct('control', 'Control', 'moderate', 'Moderate (CAA17)', 'severe', 'Severe');
cfg.supp.contrast     = 'moderate - control';
cfg.supp.letters      = {'a', 'b', 'c', 'd'};      % one per stain row, then the difference row
cfg.supp.out_name     = sprintf('supp_moderate_CAA17_histo_%s_%s', cfg.scale, cfg.model_set);
cfg.supp.size_main    = [27 10 * (numel(cfg.stains) + 1)];
cfg.supp.size_row     = [27 10];

% Supplemental: subject-level fits (one joint EPVS + vessel model per subject)
cfg.export_subjects    = true;
cfg.subj.stage_order   = {'control', 'moderate', 'severe'};   % column order of subjects
cfg.subj.out_name      = sprintf('supp_subjects_histo_%s', cfg.scale);
cfg.subj.cm_per_plot   = [7 7];          % size of each plot in the multi-panel figures
cfg.subj.export_ev     = true;           % extra figure: EPVS - vessel for every subject
cfg.subj.ylim_ev       = [];             % [] = automatic
cfg.subj.stage_colors  = struct('control',  [0.35 0.35 0.35], ...
                                'moderate', [0.93 0.55 0.00], ...
                                'severe',   [0.70 0.13 0.13]);

%% Check that the brms results exist
root = fullfile(cfg.results_dir, cfg.scale);
if ~isfolder(root)
    error(['No brms results at:\n  %s\nRun fit_histo_shells_brms.R first with\n' ...
           '  input_filename <- "%s"\n  scales <- "%s"'], root, cfg.input_csv, cfg.scale);
end

%% Load every CSV once
data = struct();
for s = 1:numel(cfg.stains)
    for g = 1:numel(cfg.regions)
        st = cfg.stains{s}; rg = cfg.regions{g};
        data.fitted.(st).(rg) = load_results(root, cfg, st, rg, 'fitted');
        data.diff.(st).(rg)   = load_results(root, cfg, st, rg, 'diff');
        data.subj.(st).(rg)    = load_subject_results(root, cfg, st, rg, 'fitted');
        data.subj_ev.(st).(rg) = load_subject_results(root, cfg, st, rg, 'ev');
    end
end

%% Resolve y-axis limits (fixed / data / auto)
cfg = resolve_ylims(cfg, data);

% Row definitions: one curve row per stain, then the difference row
n_rows = numel(cfg.stains) + 1;
rows = struct('kind', {}, 'stain', {}, 'letter', {});
for s = 1:numel(cfg.stains)
    rows(end+1) = struct('kind', 'curve', 'stain', cfg.stains{s}, 'letter', cfg.letters{s}); %#ok<SAGROW>
end
rows(end+1) = struct('kind', 'diff', 'stain', '', 'letter', cfg.letters{n_rows});

if ~exist(cfg.out_dir, 'dir'), mkdir(cfg.out_dir); end

%% Main figure
if cfg.export_main
    fig   = new_figure(cfg.size_main);
    outer = tiledlayout(fig, n_rows, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
    for r = 1:n_rows
        rowT = tiledlayout(outer, 1, numel(cfg.regions), 'TileSpacing', 'compact', 'Padding', 'compact');
        rowT.Layout.Tile = r;
        draw_row(rowT, rows(r), data, cfg);
    end
    export_all(fig, fullfile(cfg.out_dir, [cfg.out_name '_main']), cfg);
end

%% One figure per row
if cfg.export_rows
    rdir = fullfile(cfg.out_dir, 'rows');
    if ~exist(rdir, 'dir'), mkdir(rdir); end
    for r = 1:n_rows
        fig  = new_figure(cfg.size_row);
        rowT = tiledlayout(fig, 1, numel(cfg.regions), 'TileSpacing', 'compact', 'Padding', 'compact');
        draw_row(rowT, rows(r), data, cfg);
        export_all(fig, fullfile(rdir, sprintf('%s_%s', cfg.out_name, rows(r).letter)), cfg);
        close(fig);
    end
end

%% One figure per plot
if cfg.export_panels
    pdir = fullfile(cfg.out_dir, 'panels');
    if ~exist(pdir, 'dir'), mkdir(pdir); end
    for r = 1:n_rows
        n_j = ternary(strcmp(rows(r).kind, 'curve'), numel(cfg.stages), numel(cfg.stains));
        for g = 1:numel(cfg.regions)
            for j = 1:n_j
                fig = new_figure(cfg.size_panel);
                ax  = axes(fig);
                [~, name] = draw_axes(ax, rows(r), g, j, data, cfg, cfg.legend_on_every_panel, true);
                export_all(fig, fullfile(pdir, sprintf('%s_%s_%s_%s', cfg.out_name, ...
                           rows(r).letter, cfg.regions{g}, name)), cfg);
                close(fig);
            end
        end
    end
end

%% Supplemental figure: moderate stage (CAA17), occipital
if cfg.export_supp
    supp_label = @(st, is_diff) y_label(st, cfg.scale, is_diff);
    supp_name  = @(st) st;
    run_supplemental(data, cfg, cfg.stains, cfg.ylim.curve, cfg.ylim.diff, supp_label, supp_name);
end

%% Supplemental figures: subject-level fits
if cfg.export_subjects
    run_subject_supplemental(data, cfg);
end

fprintf('Done. Figures in %s\n', cfg.out_dir);

%% ========================================================================
%% Drawing
%% ========================================================================
function draw_row(rowT, row, data, cfg)
n_j = ternary(strcmp(row.kind, 'curve'), numel(cfg.stages), numel(cfg.stains));
for g = 1:numel(cfg.regions)
    regT = tiledlayout(rowT, 1, n_j, 'TileSpacing', 'compact', 'Padding', 'compact');
    regT.Layout.Tile = g;
    if strcmp(row.kind, 'curve')
        ttl = sprintf('%s: %s', row.stain, cfg.region_names{g});
    else
        ttl = {cfg.region_names{g}, ['(' cfg.contrast ')']};
    end
    title(regT, ttl, 'FontName', cfg.font, 'FontSize', cfg.font_size + 2);
    for j = 1:n_j
        ax = nexttile(regT);
        draw_axes(ax, row, g, j, data, cfg, (g == numel(cfg.regions) && j == 1), false);
        if g == 1 && j == 1, panel_letter(ax, row.letter, cfg); end
    end
end
end

function [h, name] = draw_axes(ax, row, g, j, data, cfg, show_legend, standalone)
% curve row: j = stage index;  diff row: j = stain index
rgn = cfg.regions{g};
if strcmp(row.kind, 'curve')
    st  = row.stain;
    stg = cfg.stages{j};
    T   = data.fitted.(st).(rgn);
    if ~isempty(T), T = T(T.stage == stg, :); end
    yl   = cfg.ylim.curve.(st);
    ttl  = [upper(stg(1)) stg(2:end)];
    name = stg;
    if standalone, ttl = sprintf('%s, %s: %s', st, cfg.region_names{g}, ttl); end
else
    st  = cfg.stains{j};
    T   = data.diff.(st).(rgn);
    if ~isempty(T), T = T(T.contrast == cfg.contrast, :); end
    yl   = cfg.ylim.diff.(st);
    ttl  = st;
    name = st;
    if standalone, ttl = sprintf('%s (%s): %s', cfg.region_names{g}, cfg.contrast, st); end
end
h = [];
if ~isempty(T), h = plot_classes(ax, T, cfg); end
style_axes(ax, cfg, yl, y_label(st, cfg.scale, strcmp(row.kind, 'diff')));
title(ax, ttl, 'FontWeight', 'normal');
if show_legend, add_legend(ax, h); end
end

function s = y_label(stain, scale, is_diff)
switch scale
    case 'z',      u = 'z-score';
    case 'hmatch', u = 'hist.-matched';
    otherwise,     u = 'intensity';
end
if is_diff
    s = sprintf('\\Delta %s (%s)', stain, u);
else
    s = sprintf('%s (%s)', stain, u);
end
end

function h = plot_classes(ax, T, cfg)
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
    h(c) = plot(ax, x, m, '-', 'Color', col, 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 12);
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
function T = load_results(root, cfg, stain, rgn, kind)
% R layout: <root>/<stain>/<region>/group/<src>/all_stages/<scale>_<stain>_<region>_<src>_all_stages_*.csv
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
    stem = sprintf('%s_%s_%s_%s', cfg.scale, stain, rgn, srcs{k});
    f = fullfile(root, stain, rgn, 'group', srcs{k}, 'all_stages', [stem suffix]);
    if ~isfile(f)
        warning('Missing file, panel left empty: %s', f);
        continue
    end
    Ti = readtable(f, 'TextType', 'string');
    for v = setdiff(keep, Ti.Properties.VariableNames)
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
function run_subject_supplemental(data, cfg)
% Per stain: rows = regions, columns = subjects (ordered by stage).
% Each plot = that subject's joint-model EPVS and vessel curves.
% Optional extra figure: EPVS - vessel for every subject, rows = stains.
sp   = cfg.subj;
sdir = fullfile(cfg.out_dir, 'supplemental', 'subjects');
if ~exist(sdir, 'dir'), mkdir(sdir); end
nR = numel(cfg.regions);

for s = 1:numel(cfg.stains)
    st   = cfg.stains{s};
    subj = subject_order(data.subj.(st), cfg);
    if isempty(subj)
        warning('%s: no subject-level results found, skipped.', st);
        continue
    end
    nS   = height(subj);
    name = sprintf('%s_%s', sp.out_name, st);

    if cfg.export_main
        fig   = new_figure([sp.cm_per_plot(1) * nS + 2, sp.cm_per_plot(2) * nR + 2]);
        outer = tiledlayout(fig, nR, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
        title(outer, sprintf('%s: subject-level fits', st), 'FontName', cfg.font, ...
              'FontSize', cfg.font_size + 3);
        for g = 1:nR
            rowT = tiledlayout(outer, 1, nS, 'TileSpacing', 'compact', 'Padding', 'compact');
            rowT.Layout.Tile = g;
            subject_row(rowT, st, g, subj, data, cfg);
        end
        export_all(fig, fullfile(sdir, [name '_main']), cfg);
    end
    if cfg.export_rows
        for g = 1:nR
            fig  = new_figure([sp.cm_per_plot(1) * nS + 2, sp.cm_per_plot(2) + 2]);
            rowT = tiledlayout(fig, 1, nS, 'TileSpacing', 'compact', 'Padding', 'compact');
            subject_row(rowT, st, g, subj, data, cfg);
            export_all(fig, fullfile(sdir, sprintf('%s_%s', name, cfg.regions{g})), cfg);
            close(fig);
        end
    end
    if cfg.export_panels
        for g = 1:nR
            for k = 1:nS
                T = subject_rows(data.subj.(st).(cfg.regions{g}), subj.subjectID(k));
                if isempty(T), continue; end           % e.g. CAA17 has no frontal volume
                fig = new_figure(cfg.size_panel);
                ax  = axes(fig);
                subject_axes(ax, T, st, g, subj(k, :), cfg, cfg.legend_on_every_panel, true);
                export_all(fig, fullfile(sdir, sprintf('%s_%s_%s', name, ...
                           cfg.regions{g}, subj.subjectID(k))), cfg);
                close(fig);
            end
        end
    end
end

if sp.export_ev
    fig   = new_figure([sp.cm_per_plot(1) * nR * 1.4 + 2, sp.cm_per_plot(2) * numel(cfg.stains) + 2]);
    outer = tiledlayout(fig, numel(cfg.stains), nR, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(outer, 'EPVS - vessel by subject (subject-level joint models)', ...
          'FontName', cfg.font, 'FontSize', cfg.font_size + 3);
    for s = 1:numel(cfg.stains)
        for g = 1:nR
            ax = nexttile(outer);
            ev_axes(ax, data.subj_ev.(cfg.stains{s}).(cfg.regions{g}), cfg.stains{s}, g, cfg);
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

function subject_row(rowT, st, g, subj, data, cfg)
title(rowT, cfg.region_names{g}, 'FontName', cfg.font, 'FontSize', cfg.font_size + 2);
for k = 1:height(subj)
    ax = nexttile(rowT);
    T  = subject_rows(data.subj.(st).(cfg.regions{g}), subj.subjectID(k));
    if isempty(T)                                  % subject not imaged in this region
        axis(ax, 'off');
        text(ax, 0.5, 0.5, sprintf('%s\nno %s volume', upper(subj.subjectID(k)), ...
             lower(cfg.region_names{g})), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
             'FontName', cfg.font, 'FontSize', cfg.font_size, 'Color', [0.5 0.5 0.5]);
        continue
    end
    subject_axes(ax, T, st, g, subj(k, :), cfg, (k == 1), false);
end
end

function subject_axes(ax, T, st, g, subj_row, cfg, show_legend, standalone)
h = plot_classes(ax, T, cfg);
style_axes(ax, cfg, cfg.ylim.curve.(st), y_label(st, cfg.scale, false));
ttl = sprintf('%s (%s)', upper(subj_row.subjectID), subj_row.stage);
if standalone, ttl = sprintf('%s, %s: %s', st, cfg.region_names{g}, ttl); end
title(ax, ttl, 'FontWeight', 'normal');
if show_legend, add_legend(ax, h); end
end

function ev_axes(ax, T, st, g, cfg)
% One line + ribbon per subject; colour = stage, line style separates subjects
hold(ax, 'on');
styles = {'-', '--', ':', '-.'};
h = gobjects(0); lbl = {};
if ~isempty(T)
    ids = unique(T.subjectID, 'stable');
    stage_count = struct();
    for k = 1:numel(ids)
        S   = sortrows(T(T.subjectID == ids(k), :), 'distance');
        stg = char(S.stage(1));
        if isfield(stage_count, stg), stage_count.(stg) = stage_count.(stg) + 1;
        else, stage_count.(stg) = 1; end
        col = [0.5 0.5 0.5];
        if isfield(cfg.subj.stage_colors, stg), col = cfg.subj.stage_colors.(stg); end
        x = S.distance;
        fill(ax, [x; flipud(x)], [S.lower; flipud(S.upper)], col, 'FaceAlpha', 0.12, ...
             'EdgeColor', 'none', 'HandleVisibility', 'off');
        h(end+1) = plot(ax, x, S.mean, styles{mod(stage_count.(stg) - 1, numel(styles)) + 1}, ...
                        'Color', col, 'LineWidth', 1.5); %#ok<AGROW>
        lbl{end+1} = sprintf('%s (%s)', upper(ids(k)), stg); %#ok<AGROW>
    end
end
style_axes(ax, cfg, cfg.subj.ylim_ev, sprintf('EPVS - vessel, %s', y_label(st, cfg.scale, false)));
title(ax, sprintf('%s, %s', st, cfg.region_names{g}), 'FontWeight', 'normal');
if ~isempty(h), legend(ax, h, lbl, 'Location', 'best', 'Box', 'on', 'FontSize', cfg.font_size - 2); end
end

function T = load_subject_results(root, cfg, stain, rgn, kind)
% R layout: <root>/<stain>/<region>/subject/all_subjects/<scale>_<stain>_<region>_all_subjects_*.csv
switch kind
    case 'fitted'
        suffix = '_fitted.csv';
        keep   = {'subjectID', 'stage', 'class', 'distance', 'mean', 'lower', 'upper', 'low_support'};
    case 'ev'
        suffix = '_EPVS_minus_vessel.csv';
        keep   = {'subjectID', 'stage', 'distance', 'mean', 'lower', 'upper'};
end
f = fullfile(root, stain, rgn, 'subject', 'all_subjects', ...
             sprintf('%s_%s_%s_all_subjects%s', cfg.scale, stain, rgn, suffix));
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

%% ========================================================================
%% Y-axis limits
%% ========================================================================
function cfg = resolve_ylims(cfg, data)
% Fills cfg.ylim.curve.(stain) and cfg.ylim.diff.(stain) according to the
% chosen modes. Every plot (main, supplemental, subject-level) reads these.
yl = cfg.ylim;
fprintf('\nY-axis limits (curve: %s, diff: %s)\n', yl.curve_mode, yl.diff_mode);
for s = 1:numel(cfg.stains)
    st = cfg.stains{s};

    %%% Curves
    switch yl.curve_mode
        case 'fixed'
            curve = fixed_limit(yl.curve, st, 'curve');
        case 'data'
            vals = [];
            for g = 1:numel(cfg.regions)
                vals = [vals; band_values(data.fitted.(st).(cfg.regions{g}))]; %#ok<AGROW>
                if yl.include_subjects && isfield(data, 'subj')
                    vals = [vals; band_values(data.subj.(st).(cfg.regions{g}))]; %#ok<AGROW>
                end
            end
            curve = data_limit(vals, yl);
        case 'auto'
            curve = [];
        otherwise
            error('cfg.ylim.curve_mode must be ''fixed'', ''data'' or ''auto''');
    end

    %%% Differences
    switch yl.diff_mode
        case 'fixed'
            dlim = fixed_limit(yl.diff, st, 'diff');
        case 'data'
            vals = [];
            for g = 1:numel(cfg.regions)
                vals = [vals; band_values(data.diff.(st).(cfg.regions{g}))]; %#ok<AGROW>
            end
            dlim = data_limit([vals; 0], yl);           % always include zero
        case 'auto'
            dlim = [];
        otherwise
            error('cfg.ylim.diff_mode must be ''fixed'', ''data'' or ''auto''');
    end

    cfg.ylim.curve.(st) = curve;
    cfg.ylim.diff.(st)  = dlim;
    fprintf('  %-5s curve %-18s diff %s\n', st, lim_str(curve), lim_str(dlim));
end
fprintf('\n');
end

function lim = fixed_limit(lims, st, kind)
if isfield(lims, st) && numel(lims.(st)) == 2 && lims.(st)(1) < lims.(st)(2)
    lim = lims.(st);
else
    warning('No valid fixed %s limits for %s: using MATLAB default.', kind, st);
    lim = [];
end
end

function v = band_values(T)
if isempty(T), v = []; return; end
v = [T.lower; T.upper; T.mean];
v = v(isfinite(v));
end

function lim = data_limit(v, yl)
v = v(isfinite(v));
if isempty(v), lim = []; return; end
lo = min(v); hi = max(v);
pad = yl.pad_frac * max(hi - lo, eps);
lo = lo - pad; hi = hi + pad;
if ~isempty(yl.round_step) && yl.round_step > 0
    lo = floor(lo / yl.round_step) * yl.round_step;
    hi = ceil(hi  / yl.round_step) * yl.round_step;
end
if lo == hi, lo = lo - 1; hi = hi + 1; end
lim = [lo hi];
end

function s = lim_str(lim)
if isempty(lim), s = 'auto'; else, s = sprintf('[%.3g, %.3g]', lim(1), lim(2)); end
end