%% Joint-model EPVS - vessel curves (panel f) and their change with severity (panel g)
%{
Reads the joint-model outputs of fit_parenchyma_shells_brms.R
(group/joint/all_stages/) and plots:
  f  EPVS - vessel vs distance, one line per stage (control, severe),
     for µs and retardance in frontal and occipital
  g  (EPVS - vessel)_severe - (EPVS - vessel)_control vs distance

Both come from the joint model (one model per stage with EPVS and vessel
together), so the shared subject/stage offsets cancel and the intervals are
the model's 95% credible intervals.
%}

clear; clc; close all;

%% Settings (match plot_fig3_shells.m)
cfg.beta_root    = '/projectnb/npbssmic/ns/CAA/beta_stats';
cfg.base_name    = 'parenchyma_shells_2026-10-07_pct';
cfg.analysis_tag = 'main_vsub500';
cfg.scale        = 'pct';

cfg.regions      = {'front', 'occip'};
cfg.region_names = {'Frontal', 'Occipital'};
cfg.stages       = {'control', 'severe'};          % add 'moderate' for occipital CAA17
cfg.did_stage    = 'severe';                       % stage compared with control in panel g

cfg.min_distance = 0;       % e.g. 80 to hide the first (boundary) shell
cfg.max_distance = 520;     % e.g. 400 to hide the sparse far shells
cfg.xlim         = [0 540];
cfg.ylim.f       = [];      % [] = automatic
cfg.ylim.g       = [];

cfg.stage_colors.control  = [0.35 0.35 0.35];
cfg.stage_colors.moderate = [0.93 0.55 0.00];
cfg.stage_colors.severe   = [0.70 0.13 0.13];
cfg.did_color     = [0.12 0.42 0.75];
cfg.ribbon_alpha  = 0.2;
cfg.font      = 'Helvetica';
cfg.font_size = 10;

cfg.out_dir  = fullfile(cfg.beta_root, cfg.base_name, 'brms', cfg.analysis_tag, 'figures');
cfg.out_name = sprintf('Fig3_f-g_joint_EPVS_minus_vessel_%s', cfg.scale);

%% Paths and labels
root  = fullfile(cfg.beta_root, cfg.base_name, 'brms', cfg.analysis_tag, cfg.scale);
props = {'scattering', 'retardance'};
pname = struct('scattering', '\mu_s', 'retardance', 'Retardance');
did_contrast = sprintf('(EPVS - vessel) %s - control', cfg.did_stage);

%% Figure
fig   = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 36 20]);
outer = tiledlayout(fig, 2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
summary = {};

for row = 1:2                                   % 1 = f, 2 = g
    rowT = tiledlayout(outer, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    rowT.Layout.Tile = row;
    for g = 1:2
        rgn  = cfg.regions{g};
        regT = tiledlayout(rowT, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        regT.Layout.Tile = g;
        if row == 1
            ttl = cfg.region_names{g};
        else
            ttl = {cfg.region_names{g}, sprintf('change in EPVS - vessel (%s - control)', cfg.did_stage)};
        end
        title(regT, ttl, 'FontName', cfg.font, 'FontSize', cfg.font_size + 2);

        for p = 1:2
            prop = props{p};
            ax   = nexttile(regT); hold(ax, 'on');
            stem = sprintf('%s_%s_%s_joint_all_stages', cfg.scale, prop, rgn);
            adir = fullfile(root, prop, rgn, 'group', 'joint', 'all_stages');

            if row == 1
                %%% f: EPVS - vessel, one line per stage
                T = read_if_exists(fullfile(adir, [stem '_EPVS_minus_vessel.csv']));
                h = gobjects(0); lbl = {};
                for s = 1:numel(cfg.stages)
                    stg = cfg.stages{s};
                    if isempty(T), break; end
                    S = trim(T(T.stage == stg, :), cfg);
                    if isempty(S), continue; end
                    h(end+1) = ribbon_line(ax, S, cfg.stage_colors.(stg), cfg); %#ok<SAGROW>
                    lbl{end+1} = [upper(stg(1)) stg(2:end)]; %#ok<SAGROW>
                    summary{end+1} = summarize(S, rgn, prop, "EPVS - vessel, " + stg); %#ok<SAGROW>
                end
                style_axes(ax, cfg, cfg.ylim.f, ['EPVS - vessel, ' pname.(prop) ' (\Delta%)']);
                if g == 2 && p == 1 && ~isempty(h)
                    legend(ax, h, lbl, 'Location', 'best', 'Box', 'on');
                end
                if g == 1 && p == 1, panel_letter(ax, 'f', cfg); end
            else
                %%% g: change in EPVS - vessel with severity
                T = read_if_exists(fullfile(adir, [stem '_EPVS_minus_vessel_diff.csv']));
                if ~isempty(T)
                    S = trim(T(T.contrast == did_contrast, :), cfg);
                    if ~isempty(S)
                        ribbon_line(ax, S, cfg.did_color, cfg);
                        summary{end+1} = summarize(S, rgn, prop, did_contrast); %#ok<SAGROW>
                    end
                end
                style_axes(ax, cfg, cfg.ylim.g, ['\Delta(EPVS - vessel), ' pname.(prop) ' (\Delta%)']);
                if g == 1 && p == 1, panel_letter(ax, 'g', cfg); end
            end
            title(ax, pname.(prop), 'FontWeight', 'normal');
        end
    end
end

%% Export figure and summary table
if ~exist(cfg.out_dir, 'dir'), mkdir(cfg.out_dir); end
base = fullfile(cfg.out_dir, cfg.out_name);
exportgraphics(fig, [base '.png'], 'Resolution', 300);
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
savefig(fig, [base '.fig']);

if ~isempty(summary)
    Sum = vertcat(summary{:});
    writetable(Sum, [base '_summary.csv']);
    disp(Sum);
end
fprintf('Saved %s (.png, .pdf, .fig, _summary.csv)\n', base);

%% ========================================================================
%% Local functions
%% ========================================================================
function T = read_if_exists(f)
if ~isfile(f)
    warning('Missing file (was run_joint = TRUE?): %s', f);
    T = table(); return
end
T = readtable(f, 'TextType', 'string');
for v = intersect({'stage', 'contrast'}, T.Properties.VariableNames)
    T.(v{1}) = string(T.(v{1}));
end
end

function S = trim(S, cfg)
S = sortrows(S, 'distance');
S = S(S.distance >= cfg.min_distance & S.distance <= cfg.max_distance, :);
end

function h = ribbon_line(ax, S, col, cfg)
x = S.distance;
fill(ax, [x; flipud(x)], [S.lower; flipud(S.upper)], col, 'FaceAlpha', cfg.ribbon_alpha, ...
     'EdgeColor', 'none', 'HandleVisibility', 'off');
h = plot(ax, x, S.mean, '-', 'Color', col, 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 12);
end

function style_axes(ax, cfg, yl, ylabel_str)
yline(ax, 0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, cfg.xlim);
if ~isempty(yl), ylim(ax, yl); end
xticks(ax, 0:100:500);
xlabel(ax, 'Distance (\mum)');
ylabel(ax, ylabel_str);
set(ax, 'FontName', cfg.font, 'FontSize', cfg.font_size, 'Layer', 'top', ...
        'TickDir', 'out', 'LineWidth', 0.75);
end

function panel_letter(ax, letter, cfg)
text(ax, -0.32, 1.22, letter, 'Units', 'normalized', 'FontName', cfg.font, ...
     'FontSize', cfg.font_size + 6, 'FontWeight', 'bold', 'Clipping', 'off');
end

function R = summarize(S, rgn, prop, quantity)
% One row: mean over the plotted distances, and where the 95% CrI excludes 0
excl_pos = S.distance(S.lower > 0);
excl_neg = S.distance(S.upper < 0);
R = table(string(rgn), string(prop), string(quantity), ...
          mean(S.mean), min(S.mean), max(S.mean), ...
          string(strjoin(string(excl_pos), ',')), string(strjoin(string(excl_neg), ',')), ...
    'VariableNames', {'region','property','quantity','mean_over_distances', ...
                      'min_mean','max_mean','CrI_above_0_at_um','CrI_below_0_at_um'});
end