%% Compare the generalized additive models between severities
% Import GAM model results for epvs-swp and ves-swp
% Compare the slive curves
%   These are where the epvs-swp is held constant at a single value while
%   ves-swp is swept. The slices are discrete values, so are not entire
%   story. Figure for each pair of comparisons (ie mild vs severe)
% Create surface difference matrix
%   Displays where joint ves-swp/epvs-swp diverge by severity

%% Initialization
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
% Figure output directory (also where the GAM is stored)
swp_dir = '/projectnb/npbssmic/ns/CAA/figures/swp_gam_gmm/';
% Output directory for GAM comparison
fig_dir = '/projectnb/npbssmic/ns/CAA/figures/swp_gam_gmm/gam_comparison/';
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
%%% Import GAMs
fprintf('\n----Loading GAM struct----\n')
% Import GAM for comparing control vs. severe
gam_ctl_sev = 'GAM_compare_severe_control_24-Aug-2026.mat';
gam_ctl_sev = load(fullfile(swp_dir, gam_ctl_sev));
% Import GAM for comparing subjects/regions (not combined)
gam_sub_reg = 'GAM_subject_region_24-Aug-2026.mat';
gam_sub_reg = load(fullfile(swp_dir, gam_sub_reg));

%%% Unwrap the expected variables saved by the generation script
% Block 1 (compare_isosurfaces_severe_control) saved variable 'gam_cmp'
% Block 2 (fit_isosurface_by_dataset)          saved variable 'gam_iso'
% If your save() used different variable names, adjust the two lines below.
assert(isfield(gam_ctl_sev,'gam_cmp'), ...
    ['Expected variable ''gam_cmp'' inside the compare .mat. Found: %s. ' ...
     'Adjust to the actual saved variable name.'], ...
     strjoin(fieldnames(gam_ctl_sev),', '));
assert(isfield(gam_sub_reg,'gam_iso'), ...
    ['Expected variable ''gam_iso'' inside the subject/region .mat. Found: %s. ' ...
     'Adjust to the actual saved variable name.'], ...
     strjoin(fieldnames(gam_sub_reg),', '));
gam_cmp = gam_ctl_sev.gam_cmp;      % .pdif.(reg) -> compare results struct
gam_iso = gam_sub_reg.gam_iso;      % .pdif.(sub).(reg) -> single-dataset results

%%% Shared plotting configuration
% Data type to plot ('pdif' or 'raw'); the generation blocks used 'pdif'.
d = 'pdif';
% Close each figure after it is saved (appropriate for headless/cluster
% runs). Set false to keep figures open for interactive inspection.
close_after_save = true;
% Resamples for the severe - control difference band (option II). Use a
% small value while iterating; raise to 500 for manuscript figures. This
% REFITS both groups' GAMs per resample, so runtime scales with it.
ndiff_boot = 50;
diff_alpha = 0.05;

%%% Axis limits (top-level).
% --- Y-AXIS: dynamic and SHARED across control and severe, so the two can be
% compared side by side. For each metric, one y-limit is computed from the
% union of the control and severe curves + CI bands, evaluated over each
% figure's own visible x-window (so a wide severe x-range does not inflate a
% narrow control axis). It is applied to BOTH groups' figures of that metric.
% All other blocks autoscale y (dynamic per figure). The difference block is
% on the delta scale and autoscales.
%
% --- X-AXIS: programmable per block/region. Set any value to [] to autoscale.
% Severe uses per-region x-limits (a struct keyed by region name); the other
% blocks take a single [lo hi] applied to every region.
lim_xlim_subject = [];                                        % per-subject (autoscale)
lim_xlim_control = [0 40];                                    % control: both regions
lim_xlim_severe  = struct('front', [0 500], 'occip', [0 250]);% severe: per region
lim_xlim_overlay = [0 40];                                    % comparison A: as control
lim_xlim_diff    = [0 40];                                    % comparison B: as control

% Convenience handle to the per-subject/region results for the chosen dtype
G = gam_iso.(d);

% Shared control-vs-severe y-limits (per metric), computed over each figure's
% visible x-window. Returns [] for a metric if the data are unavailable, in
% which case that metric autoscales.
[shared_yl_scatter, shared_yl_retard] = ...
    shared_control_severe_ylim(G, lim_xlim_control, lim_xlim_severe);

%% Plot each subjects GAM across the full range
% Scattering and retardance on SEPARATE figures; frontal and occipital on
% separate figures. Source: per-subject/region fits (each on its own full
% EPVS-SWP range).
fprintf('\n----Plotting per-subject GAMs (full range)----\n')
% Real subjects only: exclude the combined 'ctl'/'sev' pseudo-subjects,
% which are plotted in their own labeled sections below.
real_subs = setdiff(fieldnames(G), {'ctl','sev'}, 'stable');
for ii = 1:numel(real_subs)
    sub = real_subs{ii};
    regs = fieldnames(G.(sub));
    for j = 1:numel(regs)
        reg = regs{j};
        xl_reg = resolve_region_xlim(lim_xlim_subject, reg);
        render_group_both_metrics(G.(sub).(reg), d, sub, reg, '', ...
                                  fig_dir, close_after_save, xl_reg, [], []);
    end
end

%% Plot the control GAMs across the full range
% Scattering and retardance on separate figures; frontal/occipital separate.
% Prefer the independent full-range control fit (G.ctl); if absent, fall
% back to the SHARED (intersection) grid from the compare struct, flagged.
fprintf('\n----Plotting combined-control GAMs (full range)----\n')
plot_combined_group(G, gam_cmp, d, 'ctl', 'control', 'control', ...
                    fig_dir, close_after_save, ...
                    lim_xlim_control, shared_yl_scatter, shared_yl_retard);

%% Plot the severe GAMs across the full range
% Scattering and retardance on separate figures; frontal/occipital separate.
% Prefer the independent full-range severe fit (G.sev); if absent, fall
% back to the SHARED (intersection) grid from the compare struct, flagged.
fprintf('\n----Plotting combined-severe GAMs (full range)----\n')
plot_combined_group(G, gam_cmp, d, 'sev', 'severe', 'severe', ...
                    fig_dir, close_after_save, ...
                    lim_xlim_severe, shared_yl_scatter, shared_yl_retard);

%% Comparison A: severe vs control OVERLAY (shared/joint grid)
% Both groups on the same axes over the joint EPVS-SWP range computed at
% generation time. Severe = solid, control = dashed; same color = same
% vessel-SWP percentile. Scattering and retardance on separate figures.
fprintf('\n----Plotting severe-vs-control OVERLAY (joint grid)----\n')
regs = fieldnames(gam_cmp.(d));
for j = 1:numel(regs)
    reg = regs{j};
    res = gam_cmp.(d).(reg);
    % Scattering
    render_overlay_metric(res.epv_grid, res.slice_quantiles, ...
        res.severe.ves_levels, ...
        res.severe.slice_scatter, res.severe.density_mask, res.severe.ci_scatter, ...
        res.control.slice_scatter, res.control.density_mask, res.control.ci_scatter, ...
        'Scattering coefficient', ...
        sprintf('%s severe vs control %s (scattering, joint grid)', d, reg), ...
        sprintf('Overlay_scatter_%s_%s', d, reg), fig_dir, close_after_save, lim_xlim_overlay, []);
    % Retardance
    render_overlay_metric(res.epv_grid, res.slice_quantiles, ...
        res.severe.ves_levels, ...
        res.severe.slice_retard, res.severe.density_mask, res.severe.ci_retard, ...
        res.control.slice_retard, res.control.density_mask, res.control.ci_retard, ...
        'Retardance', ...
        sprintf('%s severe vs control %s (retardance, joint grid)', d, reg), ...
        sprintf('Overlay_retard_%s_%s', d, reg), fig_dir, close_after_save, lim_xlim_overlay, []);
end

%% Comparison B: severe - control DIFFERENCE (shared/joint grid)
% Point-estimate difference = severe curve - control curve (from the stored
% full-fit GAMs). The band is a properly resampled difference band (option
% II): both groups are refit per bootstrap resample and the difference is
% banded directly - it is NOT derived from the stored per-group CIs.
% Scattering and retardance on separate figures. Zero line = no difference.
fprintf('\n----Plotting severe-MINUS-control DIFFERENCE (joint grid)----\n')
% Difference y autoscales (delta scale, distinct from the raw-metric axes).
diff_out = struct();
regs = fieldnames(gam_cmp.(d));
for j = 1:numel(regs)
    reg = regs{j};
    res = gam_cmp.(d).(reg);
    fprintf('\tResampling difference band for %s (%d resamples)...\n', reg, ndiff_boot);
    D = bootstrap_gam_difference(res, 'NBootstrap', ndiff_boot, ...
                                 'CIAlpha', diff_alpha, 'Verbose', true);
    diff_out.(d).(reg) = D;   % keep for the separate difference-surface work
    % Scattering difference
    render_diff_metric(D.epv_grid, D.slice_quantiles, D.ves_levels, ...
        D.diff_scatter, D.mask, D.ci_diff_scatter, ...
        '\Delta Scattering (severe - control)', ...
        sprintf('%s severe - control %s (scattering, joint grid)', d, reg), ...
        sprintf('Diff_scatter_%s_%s', d, reg), fig_dir, close_after_save, lim_xlim_diff, []);
    % Retardance difference
    render_diff_metric(D.epv_grid, D.slice_quantiles, D.ves_levels, ...
        D.diff_retard, D.mask, D.ci_diff_retard, ...
        '\Delta Retardance (severe - control)', ...
        sprintf('%s severe - control %s (retardance, joint grid)', d, reg), ...
        sprintf('Diff_retard_%s_%s', d, reg), fig_dir, close_after_save, lim_xlim_diff, []);
end
% Save the difference structs (curves + resampled bands) for downstream use
dt = datetime("now",'TimeZone','local','Format','dd-MMM-yyyy');
save(fullfile(fig_dir, sprintf('GAM_difference_severe_minus_control_%s.mat', string(dt))), ...
     'diff_out','ndiff_boot','diff_alpha','-v7.3');

fprintf('\n----Done. Figures written to %s----\n', fig_dir)


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function render_group_both_metrics(res, d, group_label, reg, extra_tag, ...
                                   dirout, close_after_save, xl, yl_scatter, yl_retard)
% Render one group's scattering and retardance slice curves as two SEPARATE
% figures (frontal/occipital already separated by caller). res must carry
% epv_grid, ves_levels, slice_quantiles, slice_scatter, slice_retard,
% density_mask, ci_scatter, ci_retard.
% xl / yl_scatter / yl_retard: axis limits ([] = autoscale that axis).
    if nargin < 8;  xl = [];         end
    if nargin < 9;  yl_scatter = []; end
    if nargin < 10; yl_retard = [];  end
    tag = '';
    if ~isempty(extra_tag); tag = [' ' extra_tag]; end
    % Scattering
    fig = render_group_metric(res.epv_grid, res.ves_levels, res.slice_quantiles, ...
        res.slice_scatter, res.density_mask, res.ci_scatter, ...
        'Scattering coefficient', ...
        sprintf('%s %s %s (scattering)%s', d, group_label, reg, tag), ...
        sprintf('Slice_scatter_%s_%s_%s', d, group_label, reg), dirout, xl, yl_scatter);
    if close_after_save; close(fig); end
    % Retardance
    fig = render_group_metric(res.epv_grid, res.ves_levels, res.slice_quantiles, ...
        res.slice_retard, res.density_mask, res.ci_retard, ...
        'Retardance', ...
        sprintf('%s %s %s (retardance)%s', d, group_label, reg, tag), ...
        sprintf('Slice_retard_%s_%s_%s', d, group_label, reg), dirout, xl, yl_retard);
    if close_after_save; close(fig); end
end


function plot_combined_group(G, gam_cmp, d, iso_key, cmp_key, label, ...
                             fig_dir, close_after_save, xl, yl_scatter, yl_retard)
% Plot a combined group (control or severe), scattering and retardance on
% separate figures, one region at a time.
% xl: x-limit spec, resolved PER REGION via resolve_region_xlim -- may be []
%     (autoscale), a single [lo hi] applied to every region, or a struct with
%     per-region fields (e.g. xl.front, xl.occip).
% yl_scatter / yl_retard: y-limits ([] = autoscale) applied to every region.
%
% Preference order:
%   1. G.(iso_key).(reg)  -- independent single-dataset fit on the group's
%      OWN full EPVS-SWP range (what "full range" means).
%   2. Fallback: gam_cmp.(d).(reg).(cmp_key) -- the group's substruct from
%      the comparison, which lives on the SHARED (intersection) grid, NOT
%      the group's full range. Title/filename tagged '(shared grid)' and a
%      warning printed so the two are never confused.
    if nargin < 9;  xl = [];         end
    if nargin < 10; yl_scatter = []; end
    if nargin < 11; yl_retard = [];  end
    if isfield(G, iso_key)
        src = G.(iso_key);
        regs = fieldnames(src);
        for j = 1:numel(regs)
            reg = regs{j};
            xl_reg = resolve_region_xlim(xl, reg);
            render_group_both_metrics(src.(reg), d, label, reg, '', ...
                                      fig_dir, close_after_save, xl_reg, yl_scatter, yl_retard);
        end
    else
        warning(['%s combined full-range fit not found (G.%s missing). ' ...
                 'Falling back to the SHARED (intersection) grid from the ' ...
                 'compare struct -- this is NOT the group''s full range. To get ' ...
                 'full-range figures, include the combined %s dataset in the ' ...
                 'per-subject/region GAM fitting block.'], label, iso_key, label);
        if ~isfield(gam_cmp, d); return; end
        regs = fieldnames(gam_cmp.(d));
        for j = 1:numel(regs)
            reg = regs{j};
            cmp = gam_cmp.(d).(reg);
            res = cmp.(cmp_key);
            res.epv_grid        = cmp.epv_grid;
            res.slice_quantiles = cmp.slice_quantiles;
            xl_reg = resolve_region_xlim(xl, reg);
            render_group_both_metrics(res, d, label, reg, '(shared grid)', ...
                                      fig_dir, close_after_save, xl_reg, yl_scatter, yl_retard);
        end
    end
end


function fig = render_group_metric(epv, ves_levels, squant, curves, mask, ci, ...
                                   metric_label, titlestr, fname, dirout, xl, yl)
% One group, one metric, own figure. Bold line = full model curve across
% the whole range; shaded band = stored bootstrap CI. Colors: blue (low
% ves-SWP) -> grey -> red. (mask is retained in the API but no longer used
% to break the line; the curve is drawn solid across the full range.)
% xl / yl: axis limits ([] = autoscale that axis).
    if nargin < 11; xl = []; end
    if nargin < 12; yl = []; end
    ns     = numel(ves_levels);
    colors = slice_colors(ns);
    has_ci = ~isempty(ci) && ~isempty(ci{1});

    fig = figure('Color','w','Position',[100 100 720 540]);
    ax  = axes(fig); hold(ax,'on');
    for sl = 1:ns
        col = colors(sl,:);
        y   = curves{sl};
        if has_ci && ~isempty(ci{sl}); fill_ci(ax, epv, ci{sl}, col); end
        plot(ax, epv, y, '-', 'Color', col, 'LineWidth', 2.2, ...
             'DisplayName', sprintf('Vessel SWP p%d', round(squant(sl)*100)));
    end
    hold(ax,'off');
    xlabel(ax,'EPVS SWP'); ylabel(ax, metric_label);
    title(ax, titlestr, 'Interpreter','none');
    % legend(ax,'Location','best');
    grid(ax,'on'); box(ax,'off');
    apply_limits(ax, xl, yl);
    save_one(fig, dirout, fname);
end


function render_overlay_metric(epv, squant, ves_levels, ...
        curves_sev, mask_sev, ci_sev, curves_ctl, mask_ctl, ci_ctl, ...
        metric_label, titlestr, fname, dirout, close_after_save, xl, yl)
% Severe (solid) vs control (dashed) on one axes, one metric, own figure.
% Same color = same vessel-SWP percentile.
% xl / yl: axis limits ([] = autoscale that axis).
    if nargin < 15; xl = []; end
    if nargin < 16; yl = []; end
    ns     = numel(ves_levels);
    colors = slice_colors(ns);
    has_ci = ~isempty(ci_sev) && ~isempty(ci_sev{1});

    fig = figure('Color','w','Position',[100 100 780 560]);
    ax  = axes(fig); hold(ax,'on');
    for sl = 1:ns
        col = colors(sl,:);
        p   = round(squant(sl)*100);
        % Severe (solid), full range
        if has_ci && ~isempty(ci_sev{sl}); fill_ci(ax, epv, ci_sev{sl}, col); end
        plot(ax, epv, curves_sev{sl}, '-', 'Color', col, 'LineWidth', 2.2, ...
             'DisplayName', sprintf('Severe p%d', p));
        % Control (dashed), full range
        if has_ci && ~isempty(ci_ctl{sl}); fill_ci(ax, epv, ci_ctl{sl}, col); end
        plot(ax, epv, curves_ctl{sl}, '--','Color', col, 'LineWidth', 2.2, ...
             'DisplayName', sprintf('Control p%d', p));
    end
    hold(ax,'off');
    xlabel(ax,'EPVS SWP'); ylabel(ax, metric_label);
    title(ax, titlestr, 'Interpreter','none');
    subtitle_safe(ax, 'Solid = severe, dashed = control; color = vessel-SWP percentile');
    % legend(ax,'Location','best','FontSize',7);
    grid(ax,'on'); box(ax,'off');
    apply_limits(ax, xl, yl);
    save_one(fig, dirout, fname);
    if close_after_save; close(fig); end
end


function render_diff_metric(epv, squant, ves_levels, diff_curves, mask, ci, ...
                            metric_label, titlestr, fname, dirout, close_after_save, xl, yl)
% Severe - control difference, one metric, own figure. Zero reference line
% marks "no difference". Band = resampled difference CI (option II).
% xl / yl: axis limits ([] = autoscale that axis). Note yl here is on the
% DELTA scale (severe - control), not the raw metric scale.
    if nargin < 12; xl = []; end
    if nargin < 13; yl = []; end
    ns     = numel(ves_levels);
    colors = slice_colors(ns);
    has_ci = ~isempty(ci) && ~isempty(ci{1});

    fig = figure('Color','w','Position',[100 100 780 560]);
    ax  = axes(fig); hold(ax,'on');
    % Zero reference
    plot(ax, [min(epv) max(epv)], [0 0], '-', 'Color', [.6 .6 .6], ...
         'LineWidth', 1.0, 'HandleVisibility','off');
    for sl = 1:ns
        col = colors(sl,:);
        y   = diff_curves{sl};
        if has_ci && ~isempty(ci{sl}); fill_ci(ax, epv, ci{sl}, col); end
        plot(ax, epv, y, '-', 'Color', col, 'LineWidth', 2.2, ...
             'DisplayName', sprintf('Severe - Control p%d', round(squant(sl)*100)));
    end
    hold(ax,'off');
    xlabel(ax,'EPVS SWP'); ylabel(ax, metric_label);
    title(ax, titlestr, 'Interpreter','none');
    subtitle_safe(ax, 'Band = resampled difference CI (full range)');
    % legend(ax,'Location','best','FontSize',8);
    grid(ax,'on'); box(ax,'off');
    apply_limits(ax, xl, yl);
    save_one(fig, dirout, fname);
    if close_after_save; close(fig); end
end


% ---- small shared rendering helpers (self-contained; mirror the toolkit) ----

function xlr = resolve_region_xlim(xl, reg)
% Resolve an x-limit spec for one region. xl may be:
%   []              -> [] (autoscale)
%   [lo hi]         -> same limits for every region
%   struct w/ .reg  -> per-region limits (missing region -> [] autoscale)
    if isempty(xl)
        xlr = [];
    elseif isstruct(xl)
        if isfield(xl, reg); xlr = xl.(reg); else; xlr = []; end
    else
        xlr = xl;
    end
end

function [yl_sc, yl_rt] = shared_control_severe_ylim(G, xl_ctl, xl_sev)
% Compute one shared y-limit per metric across the control and severe
% full-range figures, so the two groups can be compared side by side. Each
% figure contributes its curve+band y-range over its OWN visible x-window
% (control uses xl_ctl; severe uses per-region xl_sev). Returns [] for a
% metric when the source fits are unavailable (that metric then autoscales).
    specs_sc = {};   % scattering figure specs
    specs_rt = {};   % retardance figure specs
    if isfield(G, 'ctl')
        rr = fieldnames(G.ctl);
        for j = 1:numel(rr)
            reg = rr{j};
            specs_sc{end+1} = {G.ctl.(reg), 'slice_scatter', 'ci_scatter', ...
                               resolve_region_xlim(xl_ctl, reg)}; %#ok<AGROW>
            specs_rt{end+1} = {G.ctl.(reg), 'slice_retard',  'ci_retard',  ...
                               resolve_region_xlim(xl_ctl, reg)}; %#ok<AGROW>
        end
    end
    if isfield(G, 'sev')
        rr = fieldnames(G.sev);
        for j = 1:numel(rr)
            reg = rr{j};
            specs_sc{end+1} = {G.sev.(reg), 'slice_scatter', 'ci_scatter', ...
                               resolve_region_xlim(xl_sev, reg)}; %#ok<AGROW>
            specs_rt{end+1} = {G.sev.(reg), 'slice_retard',  'ci_retard',  ...
                               resolve_region_xlim(xl_sev, reg)}; %#ok<AGROW>
        end
    end
    yl_sc = shared_ylim_over_figures(specs_sc);
    yl_rt = shared_ylim_over_figures(specs_rt);
end

function yl = shared_ylim_over_figures(specs)
% Union of y-ranges (curves + CI bands) across a set of figures, each
% restricted to its own visible x-window. specs is a cell array; each entry
% is {res, curve_field, ci_field, xlim}. Returns [lo hi] with 5% padding, or
% [] if nothing usable.
    ymin = inf; ymax = -inf;
    for k = 1:numel(specs)
        res    = specs{k}{1};
        mfield = specs{k}{2};
        cfield = specs{k}{3};
        xl     = specs{k}{4};
        if ~isfield(res, mfield) || ~isfield(res, 'epv_grid'); continue; end
        epv = res.epv_grid(:);
        if isempty(xl)
            insel = true(size(epv));
        else
            insel = epv >= xl(1) & epv <= xl(2);
        end
        if ~any(insel); continue; end
        curves = res.(mfield);
        has_ci = isfield(res, cfield) && ~isempty(res.(cfield)) && ~isempty(res.(cfield){1});
        for sl = 1:numel(curves)
            y = curves{sl}(insel);
            ymin = min(ymin, min(y));  ymax = max(ymax, max(y));
            if has_ci && ~isempty(res.(cfield){sl})
                ci = res.(cfield){sl};
                ymin = min(ymin, min(ci(insel,1)));
                ymax = max(ymax, max(ci(insel,2)));
            end
        end
    end
    if ~isfinite(ymin) || ~isfinite(ymax) || ymin == ymax
        yl = [];
    else
        pad = 0.05 * (ymax - ymin);
        yl  = [ymin - pad, ymax + pad];
    end
end

function apply_limits(ax, xl, yl)
% Apply axis limits when provided; [] leaves that axis autoscaled.
    if ~isempty(xl); xlim(ax, xl); end
    if ~isempty(yl); ylim(ax, yl); end
end

function colors = slice_colors(ns)
    % base = [0.20 0.45 0.80; 0.50 0.50 0.50; 0.80 0.20 0.20];
    base = [hex2rgb('#1964B0');hex2rgb('#882D71');hex2rgb('#DB5829')];
    if ns == 3
        colors = base;
    else
        colors = interp1([0 0.5 1]', base, linspace(0,1,ns)', 'linear');
    end
end

function fill_ci(ax, x, ci, col)
    x  = x(:);
    xp = [x; flipud(x)];
    yp = [ci(:,1); flipud(ci(:,2))];
    fill(ax, xp, yp, col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
         'HandleVisibility','off');
end

function subtitle_safe(ax, str)
% subtitle() exists R2020b+; guard for older MATLAB.
    try
        subtitle(ax, str, 'FontSize', 8, 'Color', [.5 .5 .5]);
    catch
        % no-op on versions without subtitle()
    end
end

function save_one(fig, dirout, fname)
    if ~isfolder(dirout); mkdir(dirout); end
    fname = regexprep(fname, '[^\w\-]', '_');
    fname = regexprep(fname, '_+', '_');
    fname = regexprep(fname, '_+$', '');
    saveas(fig, fullfile(dirout, [fname '.fig']));
    saveas(fig, fullfile(dirout, [fname '.png']));
    exportgraphics(fig,fullfile(dirout,[fname '.svg']),'Resolution',600);
end