%% examine_parenchyma_shells: checks on the output of measure_parenchyma_shells_strat.m
%{
Reads the three CSVs of one run folder and writes figures + CSVs of every
plotted value to <run folder>/examine_<yyyyMMdd_HHmmss>_<EXPLORATORY|MANUSCRIPT>/.
Nothing in the run folder is modified and no image volumes are loaded.

Steps
  0  Baseline checks (table only): mean vs median, coverage, exclusion, etc.
  1  Shell curves under four settings, to separate what the baseline and the
     weighting each contribute:
        a  % from WM MEDIAN, objects weighted equally   (as in the current Fig. 5)
        b  % from WM MEAN,   objects weighted equally
        c  % from WM MEAN,   voxels weighted equally
        d  as c, with the largest EPVS and largest vessel component removed
     Group curves (subjects weighted equally) with a hierarchical bootstrap
     (objects resampled within subject), plus one figure of per-subject curves.
  2  Depth: near-boundary value (shells out to NEAR_UM) by depth band, EPVS
     and vessel side by side, so the two classes are compared at equal depth.
  3  Caliber (EPVS only): near-boundary value by within-volume caliber bin,
     and per-object scatter against local caliber with Spearman's rho.
     Vessel caliber is skipped: most vessels are one voxel wide.
  4  Near-boundary EPVS vs vessel per subject, with the EPVS - vessel
     difference and its bootstrap interval.

Normalization: mus_pct / ret_pct in the CSVs are % from the WM median.
mus_pct_mean / ret_pct_mean are computed here as % from the WM mean in
baseline.csv. Shell values are means of voxels, so the WM mean is the
baseline that makes "0" mean "same as white matter overall".

Requires Statistics and Machine Learning Toolbox (prctile, corr).
%}

clearvars -except -regexp ^caa\d*$
clc; close all; rng(1);

%% Configuration
cfg.out_base    = '/projectnb/npbssmic/ns/CAA/parenchyma_shells/';
cfg.run_dir     = 'latest';        % 'latest' or the full path of one run folder
cfg.exploratory = true;            % 100 bootstrap resamples; false = 500
cfg.near_um     = 120;             % shells with outer edge <= this = "near boundary"
cfg.n_min_bin   = 25;              % min voxels for an object x shell to count (object weighting)
cfg.n_min_obj   = 5;               % min objects for a stratified value to be shown
cfg.shell_um    = 40;
cfg.props       = ["mus", "ret"];
cfg.regions     = ["front", "occip"];
cfg.groups      = ["control", "severe"];
cfg.subjects    = ["caa26", "caa6", "caa22", "caa25", "caa17"];   % display order
cfg.classes     = ["EPVS", "vessel"];
cfg.colors      = [0.00 0.45 0.74; 0.85 0.33 0.10];               % EPVS blue, vessel orange
cfg.vis         = 'off';
cfg.n_boot      = 500;
if cfg.exploratory, cfg.n_boot = 100; end

%% Locate the run and read the tables
run_dir = cfg.run_dir;
if strcmpi(run_dir, 'latest')
    d = dir(cfg.out_base);
    d = d([d.isdir] & ~startsWith({d.name}, '.'));
    has = arrayfun(@(x) isfile(fullfile(x.folder, x.name, 'parenchyma_shells_shells.csv')), d);
    d = d(has);
    assert(~isempty(d), 'No run folder with parenchyma_shells_shells.csv under %s', cfg.out_base);
    [~, i] = max([d.datenum]);
    run_dir = fullfile(d(i).folder, d(i).name);
end
fprintf('Run folder: %s\n', run_dir);
rd = @(f) readtable(fullfile(run_dir, f), 'TextType', 'string', 'Delimiter', ',');
B  = rd('parenchyma_shells_baseline.csv');
Sh = rd('parenchyma_shells_shells.csv');
St = rd('parenchyma_shells_strat.csv');

tag = "MANUSCRIPT";
if cfg.exploratory, tag = "EXPLORATORY"; end
stamp   = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir = char(fullfile(run_dir, "examine_" + stamp + "_" + tag));
mkdir(out_dir);
fprintf('Writing to %s\n', out_dir);

cfg.n_shell = round(max(Sh.d_outer_um) / cfg.shell_um);
Sh = prepare(Sh, B, cfg);
St = prepare(St, B, cfg);

%% 0. Baseline checks
B.mus_mean_over_median_pct = 100 * (B.mus_mean ./ B.mus_median - 1);
B.ret_mean_over_median_pct = 100 * (B.ret_mean ./ B.ret_median - 1);
chk = B(:, {'subject','severity','region','mus_median','mus_mean','mus_mean_over_median_pct', ...
            'ret_median','ret_mean','ret_mean_over_median_pct','shell_coverage', ...
            'frac_owned_epvs','frac_owned_vessel','frac_owned_other', ...
            'ves_excluded_vox_frac','epvs_largest_frac','ves_largest_frac','n_epvs','n_vessel'});
disp(chk);
writetable(chk, fullfile(out_dir, 'baseline_checks.csv'));

%% 1. Shell curves: baseline and weighting
% {label, value-column suffix, weighting, drop largest component}
specs = { 'a  % from WM median, object-weighted',                   "_pct",      'object', false;
          'b  % from WM mean, object-weighted',                     "_pct_mean", 'object', false;
          'c  % from WM mean, voxel-weighted',                      "_pct_mean", 'voxel',  false;
          'd  % from WM mean, voxel-weighted, largest comp. removed', "_pct_mean", 'voxel',  true };
x = (1:cfg.n_shell) * cfg.shell_um;
rows = {};
for p = 1:numel(cfg.props)
    prop = cfg.props(p);
    fh = figure('Visible', cfg.vis, 'Color', 'w', 'Position', [50 50 1500 1150]);
    tl = tiledlayout(fh, size(specs, 1), 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    axs = gobjects(0);
    for k = 1:size(specs, 1)
        col = prop + specs{k, 2};
        for r = 1:2
            for g = 1:2
                ax = nexttile(tl); hold(ax, 'on'); axs(end+1) = ax; %#ok<SAGROW>
                yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
                for c = 1:2
                    T = Sh(Sh.region == cfg.regions(r) & Sh.severity == cfg.groups(g) & ...
                           Sh.class == cfg.classes(c), :);
                    if specs{k, 4}, T = T(~T.is_largest, :); end
                    R = group_curve(T, col, specs{k, 3}, cfg);
                    h = ribbon(ax, x, R, cfg.colors(c, :));
                    h.DisplayName = cfg.classes(c);
                    rows{end+1} = table(repmat(prop, cfg.n_shell, 1), repmat(string(specs{k, 1}), cfg.n_shell, 1), ...
                        repmat(cfg.regions(r), cfg.n_shell, 1), repmat(cfg.groups(g), cfg.n_shell, 1), ...
                        repmat(cfg.classes(c), cfg.n_shell, 1), x(:), R.mu(:), R.lo(:), R.hi(:), ...
                        'VariableNames', {'property','setting','region','group','class', ...
                                          'd_outer_um','mean','ci_lo','ci_hi'}); %#ok<SAGROW>
                end
                style(ax);
                title(ax, sprintf('%s, %s', region_name(cfg.regions(r)), cfg.groups(g)), 'FontWeight', 'normal');
                if r == 1 && g == 1, ylabel(ax, {specs{k, 1}, prop_label(prop)}); end
                if k == size(specs, 1), xlabel(ax, 'Distance (\mum)'); end
                if k == 1 && r == 1 && g == 1, legend(ax, 'Location', 'best', 'Box', 'off'); end
            end
        end
    end
    linkaxes(axs, 'y');     % one y-range for the whole figure, so the shift between rows is visible
    title(tl, sprintf('%s: effect of baseline and weighting (group mean, 95%% bootstrap interval, B = %d)', ...
          prop_label(prop), cfg.n_boot));
    save_fig(fh, out_dir, "step1_curves_" + prop);
end
writetable(vertcat(rows{:}), fullfile(out_dir, 'step1_curves.csv'));

%%% Per-subject curves (% from WM mean): solid = voxel-weighted, dashed = object-weighted
subs = cfg.subjects(ismember(cfg.subjects, unique(Sh.subject)));
for p = 1:numel(cfg.props)
    prop = cfg.props(p); col = prop + "_pct_mean";
    fh = figure('Visible', cfg.vis, 'Color', 'w', 'Position', [50 50 800 230 * numel(subs) + 80]);
    tl = tiledlayout(fh, numel(subs), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    axs = gobjects(0);
    for s = 1:numel(subs)
        for r = 1:2
            ax = nexttile(tl); hold(ax, 'on'); axs(end+1) = ax; %#ok<SAGROW>
            yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
            for c = 1:2
                T = Sh(Sh.subject == subs(s) & Sh.region == cfg.regions(r) & Sh.class == cfg.classes(c), :);
                if isempty(T), continue; end
                [V, N] = obj_shell(T, col, cfg.n_shell);
                plot(ax, x, curve(V, N, 'voxel',  cfg.n_min_bin), '-',  'Color', cfg.colors(c, :), ...
                     'LineWidth', 1.8, 'DisplayName', cfg.classes(c) + ", voxel-weighted");
                plot(ax, x, curve(V, N, 'object', cfg.n_min_bin), '--', 'Color', cfg.colors(c, :), ...
                     'LineWidth', 1.0, 'DisplayName', cfg.classes(c) + ", object-weighted");
            end
            style(ax);
            sev = unique(Sh.severity(Sh.subject == subs(s)));
            title(ax, sprintf('%s %s (%s)', upper(subs(s)), region_name(cfg.regions(r)), sev(1)), 'FontWeight', 'normal');
            if r == 1, ylabel(ax, [prop_label(prop) ' (% from WM mean)']); end
            if s == numel(subs), xlabel(ax, 'Distance (\mum)'); end
            if s == 1 && r == 1, legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 7); end
        end
    end
    linkaxes(axs, 'y');
    title(tl, sprintf('%s by subject', prop_label(prop)));
    save_fig(fh, out_dir, "step1_subjects_" + prop);
end

%% 2. Depth: near-boundary value by depth band, EPVS vs vessel at equal depth
dk  = unique(St.depth_idx(St.depth_idx > 0))';
dlab = strings(size(dk));
for i = 1:numel(dk), dlab(i) = St.depth_band(find(St.depth_idx == dk(i), 1)); end
rows = {};
for p = 1:numel(cfg.props)
    prop = cfg.props(p); col = prop + "_pct_mean";
    fh = figure('Visible', cfg.vis, 'Color', 'w', 'Position', [50 50 1500 420]);
    tl = tiledlayout(fh, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    axs = gobjects(0);
    for r = 1:2
        for g = 1:2
            ax = nexttile(tl); hold(ax, 'on'); axs(end+1) = ax; %#ok<SAGROW>
            yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
            gs = cfg.subjects(ismember(cfg.subjects, unique(St.subject(St.severity == cfg.groups(g) & St.region == cfg.regions(r)))));
            for c = 1:2
                Y = nan(numel(gs), numel(dk));
                for s = 1:numel(gs)
                    for i = 1:numel(dk)
                        T = St(St.subject == gs(s) & St.region == cfg.regions(r) & ...
                               St.class == cfg.classes(c) & St.depth_idx == dk(i), :);
                        [v, no, nv] = near_value(T, col, cfg);
                        if no < cfg.n_min_obj, v = NaN; end
                        Y(s, i) = v;
                        rows{end+1} = table(prop, gs(s), cfg.groups(g), cfg.regions(r), cfg.classes(c), ...
                            dlab(i), dk(i), v, no, nv, 'VariableNames', {'property','subject','group', ...
                            'region','class','depth_band','depth_idx','near_pct_from_mean','n_objects','n_vox'}); %#ok<SAGROW>
                    end
                    plot(ax, (1:numel(dk)) + 0.08 * (c - 1.5), Y(s, :), subj_marker(gs(s)), 'Color', cfg.colors(c, :), ...
                         'MarkerSize', 6, 'LineStyle', 'none', 'HandleVisibility', 'off');
                end
                plot(ax, (1:numel(dk)) + 0.08 * (c - 1.5), mean(Y, 1, 'omitnan'), '-', 'Color', cfg.colors(c, :), ...
                     'LineWidth', 1.8, 'DisplayName', cfg.classes(c));
            end
            set(ax, 'XLim', [0.5 numel(dk) + 0.5], 'XTick', 1:numel(dk), 'XTickLabel', dlab);
            style(ax);
            title(ax, sprintf('%s, %s', region_name(cfg.regions(r)), cfg.groups(g)), 'FontWeight', 'normal');
            xlabel(ax, 'Depth from gray-white interface');
            if r == 1 && g == 1
                ylabel(ax, sprintf('%s within %g \\mum (%% from WM mean)', prop_label(prop), cfg.near_um));
                legend(ax, 'Location', 'best', 'Box', 'off');
            end
        end
    end
    linkaxes(axs, 'y');
    title(tl, sprintf('%s near the boundary, by depth (markers = subjects, line = mean; voxel-weighted; >= %d objects)', ...
          prop_label(prop), cfg.n_min_obj));
    save_fig(fh, out_dir, "step2_depth_" + prop);
end
writetable(vertcat(rows{:}), fullfile(out_dir, 'step2_near_by_depth.csv'));

%% 3. Caliber (EPVS only)
rows = {}; rho_rows = {};
for p = 1:numel(cfg.props)
    prop = cfg.props(p); col = prop + "_pct_mean";
    fh = figure('Visible', cfg.vis, 'Color', 'w', 'Position', [50 50 1500 760]);
    tl = tiledlayout(fh, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = gobjects(0); ax2 = gobjects(0);
    for row = 1:2
        for r = 1:2
            for g = 1:2
                ax = nexttile(tl); hold(ax, 'on');
                yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
                gs = cfg.subjects(ismember(cfg.subjects, unique(Sh.subject(Sh.severity == cfg.groups(g) & Sh.region == cfg.regions(r)))));
                if row == 1, ax1(end+1) = ax; else, ax2(end+1) = ax; end %#ok<SAGROW>
                for s = 1:numel(gs)
                    if row == 1                                   % by within-volume caliber bin
                        cb = unique(St.cal_bin(St.cal_bin > 0 & St.class == "EPVS" & St.subject == gs(s) & St.region == cfg.regions(r)))';
                        y = nan(size(cb));
                        for i = 1:numel(cb)
                            T = St(St.subject == gs(s) & St.region == cfg.regions(r) & St.class == "EPVS" & St.cal_bin == cb(i), :);
                            [v, no, nv] = near_value(T, col, cfg);
                            if no < cfg.n_min_obj, v = NaN; end
                            y(i) = v;
                            rows{end+1} = table(prop, gs(s), cfg.groups(g), cfg.regions(r), cb(i), v, no, nv, ...
                                'VariableNames', {'property','subject','group','region','cal_bin', ...
                                                  'near_pct_from_mean','n_objects','n_vox'}); %#ok<SAGROW>
                        end
                        plot(ax, cb, y, ['-' subj_marker(gs(s))], 'Color', cfg.colors(1, :), 'MarkerSize', 6, ...
                             'LineWidth', 1.2, 'DisplayName', upper(gs(s)));
                        set(ax, 'XLim', [0.5 3.5], 'XTick', 1:3, 'XTickLabel', {'small','medium','large'});
                        xlabel(ax, 'EPVS local caliber (within-volume bin)');
                    else                                          % per object, continuous caliber
                        T = Sh(Sh.subject == gs(s) & Sh.region == cfg.regions(r) & Sh.class == "EPVS", :);
                        [V, N] = near_objects(T, col, cfg);
                        C = near_objects(T, "local_cal_um", cfg);
                        ok = N >= cfg.n_min_bin & isfinite(V) & isfinite(C);
                        rho = NaN;
                        if nnz(ok) >= cfg.n_min_obj
                            rho = corr(C(ok), V(ok), 'Type', 'Spearman');
                        end
                        plot(ax, C(ok), V(ok), subj_marker(gs(s)), 'Color', cfg.colors(1, :), 'MarkerSize', 4, ...
                             'LineStyle', 'none', 'DisplayName', sprintf('%s (\\rho = %.2f, n = %d)', upper(gs(s)), rho, nnz(ok)));
                        rho_rows{end+1} = table(prop, gs(s), cfg.groups(g), cfg.regions(r), rho, nnz(ok), ...
                            'VariableNames', {'property','subject','group','region','spearman_rho','n_objects'}); %#ok<SAGROW>
                        set(ax, 'XScale', 'log');
                        xlabel(ax, 'EPVS local caliber (\mum)');
                    end
                end
                style(ax);
                legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 7);
                title(ax, sprintf('%s, %s', region_name(cfg.regions(r)), cfg.groups(g)), 'FontWeight', 'normal');
                if r == 1 && g == 1
                    ylabel(ax, sprintf('%s within %g \\mum (%% from WM mean)', prop_label(prop), cfg.near_um));
                end
            end
        end
    end
    linkaxes(ax1, 'y'); linkaxes(ax2, 'y');
    title(tl, sprintf('%s near EPVS against EPVS caliber (top: binned, voxel-weighted; bottom: one point per EPVS)', prop_label(prop)));
    save_fig(fh, out_dir, "step3_caliber_" + prop);
end
writetable(vertcat(rows{:}),     fullfile(out_dir, 'step3_near_by_caliber.csv'));
writetable(vertcat(rho_rows{:}), fullfile(out_dir, 'step3_caliber_spearman.csv'));

%% 4. Near-boundary EPVS vs vessel, per subject, with the difference
rows = {};
for p = 1:numel(cfg.props)
    prop = cfg.props(p); col = prop + "_pct_mean";
    fh = figure('Visible', cfg.vis, 'Color', 'w', 'Position', [50 50 1000 700]);
    tl = tiledlayout(fh, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    axt = gobjects(1, 2); axb = gobjects(1, 2);
    for r = 1:2
        axt(r) = nexttile(tl, r);     hold(axt(r), 'on');
        axb(r) = nexttile(tl, 2 + r); hold(axb(r), 'on');
        yline(axt(r), 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        yline(axb(r), 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        rs = subs(ismember(subs, unique(Sh.subject(Sh.region == cfg.regions(r)))));
        for s = 1:numel(rs)
            est = nan(1, 2); bt = nan(cfg.n_boot, 2); raw = nan(1, 2); nobj = zeros(1, 2);
            for c = 1:2
                T = Sh(Sh.subject == rs(s) & Sh.region == cfg.regions(r) & Sh.class == cfg.classes(c), :);
                [V, N] = near_objects(T, col, cfg);
                Vr = near_objects(T, prop, cfg);
                nobj(c) = numel(V);
                if isempty(V), continue; end
                est(c) = wmean(V, N);
                raw(c) = wmean(Vr, N);
                for b = 1:cfg.n_boot
                    k = randi(numel(V), numel(V), 1);
                    bt(b, c) = wmean(V(k), N(k));
                end
                ci = prctile(bt(:, c), [2.5 97.5]);
                errorbar(axt(r), s + 0.12 * (c - 1.5), est(c), est(c) - ci(1), ci(2) - est(c), 'o', ...
                         'Color', cfg.colors(c, :), 'MarkerFaceColor', cfg.colors(c, :), 'MarkerSize', 6, ...
                         'LineWidth', 1.2, 'CapSize', 4, 'HandleVisibility', 'off');
            end
            dif = est(1) - est(2);
            dci = prctile(bt(:, 1) - bt(:, 2), [2.5 97.5]);   % classes resampled independently
            errorbar(axb(r), s, dif, dif - dci(1), dci(2) - dif, 'ko', 'MarkerFaceColor', 'k', ...
                     'MarkerSize', 6, 'LineWidth', 1.2, 'CapSize', 4);
            sev = unique(Sh.severity(Sh.subject == rs(s)));
            rows{end+1} = table(prop, rs(s), sev(1), cfg.regions(r), est(1), est(2), dif, dci(1), dci(2), ...
                raw(1), raw(2), nobj(1), nobj(2), 'VariableNames', {'property','subject','severity','region', ...
                'epvs_pct_from_mean','vessel_pct_from_mean','epvs_minus_vessel','diff_ci_lo','diff_ci_hi', ...
                'epvs_raw','vessel_raw','n_epvs','n_vessel'}); %#ok<SAGROW>
        end
        for c = 1:2
            plot(axt(r), NaN, NaN, 'o', 'Color', cfg.colors(c, :), 'MarkerFaceColor', cfg.colors(c, :), ...
                 'DisplayName', cfg.classes(c));
        end
        for ax = [axt(r) axb(r)]
            set(ax, 'XLim', [0.5 numel(rs) + 0.5], 'XTick', 1:numel(rs), 'XTickLabel', upper(rs));
            style(ax);
        end
        title(axt(r), region_name(cfg.regions(r)), 'FontWeight', 'normal');
        if r == 1
            ylabel(axt(r), sprintf('%s within %g \\mum (%% from WM mean)', prop_label(prop), cfg.near_um));
            ylabel(axb(r), 'EPVS - vessel (percentage points)');
            legend(axt(r), 'Location', 'best', 'Box', 'off');
        end
    end
    linkaxes(axt, 'y'); linkaxes(axb, 'y');
    title(tl, sprintf('%s near the boundary by subject (voxel-weighted; bars = 95%% bootstrap over objects, B = %d)', ...
          prop_label(prop), cfg.n_boot));
    save_fig(fh, out_dir, "step4_near_by_subject_" + prop);
end
writetable(vertcat(rows{:}), fullfile(out_dir, 'step4_near_epvs_vs_vessel.csv'));

%% Record settings
fid = fopen(fullfile(out_dir, 'examine_config.json'), 'w');
cfg.run_dir_used = run_dir;
fprintf(fid, '%s\n', jsonencode(cfg, 'PrettyPrint', true));
fclose(fid);
fprintf('Done. Outputs in %s\n', out_dir);


%% ======================= Local functions ===============================

function T = prepare(T, B, cfg)
% Row filter, shell index, % from WM mean, largest-component flag.
n0 = height(T);
T = T(T.mus_ok ~= 0 & T.ret_ok ~= 0 & T.n_vox > 0, :);
fprintf('  kept %d of %d rows (mus_ok & ret_ok)\n', height(T), n0);
T.sh = round(T.d_outer_um / cfg.shell_um);
T = T(T.sh >= 1 & T.sh <= cfg.n_shell, :);
[tf, ib] = ismember(T.subject + "_" + T.region, B.subject + "_" + B.region);
assert(all(tf), 'Rows without a baseline entry.');
T.mus_pct_mean = 100 * (T.mus - B.mus_mean(ib)) ./ B.mus_mean(ib);
T.ret_pct_mean = 100 * (T.ret - B.ret_mean(ib)) ./ B.ret_mean(ib);
G  = findgroups(T.subject, T.region, T.class);
mx = splitapply(@max, T.obj_vox, G);
T.is_largest = T.obj_vox == mx(G);       % the largest connected component of its class in its volume
end


function [V, N] = obj_shell(T, col, nS)
% Object x shell matrices: voxel-weighted mean of COL and voxel count.
% (Collapses the stratified table's several rows per object x shell.)
[~, ~, io] = unique(T.object_id);
nO = max([io; 0]);
V = nan(nO, nS); N = zeros(nO, nS);
if nO == 0, return; end
N = accumarray([io, T.sh], T.n_vox, [nO nS]);
V = accumarray([io, T.sh], T.(col) .* T.n_vox, [nO nS]) ./ N;    % NaN where N == 0
end


function c = curve(V, N, w, nmin)
% One subject's curve across shells from object x shell matrices.
if isempty(V), c = nan(1, size(V, 2)); return; end
switch w
    case 'object'                       % each object counts once
        V(N < nmin) = NaN;
        c = mean(V, 1, 'omitnan');
    case 'voxel'                        % each voxel counts once
        W = N; W(isnan(V)) = 0; V(isnan(V)) = 0;
        c = sum(V .* W, 1) ./ sum(W, 1);
end
end


function R = group_curve(T, col, w, cfg)
% Group curve: subjects weighted equally. 95% interval from a hierarchical
% bootstrap that resamples objects within each subject. The interval
% describes variability among the measured objects, not between subjects.
nS   = cfg.n_shell;
subs = unique(T.subject, 'stable');
M    = cell(numel(subs), 2);
R.subj = nan(numel(subs), nS);
for s = 1:numel(subs)
    [M{s, 1}, M{s, 2}] = obj_shell(T(T.subject == subs(s), :), col, nS);
    R.subj(s, :) = curve(M{s, 1}, M{s, 2}, w, cfg.n_min_bin);
end
R.subs = subs;
R.mu = nan(1, nS); R.lo = nan(1, nS); R.hi = nan(1, nS);
if isempty(subs), return; end
R.mu = mean(R.subj, 1, 'omitnan');
boot = nan(cfg.n_boot, nS);
for b = 1:cfg.n_boot
    cs = nan(numel(subs), nS);
    for s = 1:numel(subs)
        nO = size(M{s, 1}, 1);
        if nO == 0, continue; end
        k = randi(nO, nO, 1);
        cs(s, :) = curve(M{s, 1}(k, :), M{s, 2}(k, :), w, cfg.n_min_bin);
    end
    boot(b, :) = mean(cs, 1, 'omitnan');
end
R.lo = prctile(boot, 2.5, 1);
R.hi = prctile(boot, 97.5, 1);
end


function [V, N] = near_objects(T, col, cfg)
% Per object: voxel-weighted mean of COL over the near-boundary shells.
T = T(T.d_outer_um <= cfg.near_um, :);
[~, ~, io] = unique(T.object_id);
nO = max([io; 0]);
if nO == 0, V = zeros(0, 1); N = zeros(0, 1); return; end
N = accumarray(io, T.n_vox, [nO 1]);
V = accumarray(io, T.(col) .* T.n_vox, [nO 1]) ./ N;
end


function [v, n_obj, n_vox] = near_value(T, col, cfg)
% Voxel-weighted near-boundary value across all objects in T.
[V, N] = near_objects(T, col, cfg);
n_obj = numel(V); n_vox = sum(N);
v = wmean(V, N);
end


function m = wmean(V, N)
ok = isfinite(V) & N > 0;
if ~any(ok), m = NaN; else, m = sum(V(ok) .* N(ok)) / sum(N(ok)); end
end


function h = ribbon(ax, x, R, c)
ok = isfinite(R.lo) & isfinite(R.hi);
if any(ok)
    fill(ax, [x(ok) fliplr(x(ok))], [R.lo(ok) fliplr(R.hi(ok))], c, 'FaceAlpha', 0.2, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
end
h = plot(ax, x, R.mu, '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.5);
end


function style(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
grid(ax, 'on');
end


function m = subj_marker(sub)
switch char(sub)
    case {'caa26', 'caa22'}, m = 'd';
    case {'caa6',  'caa25'}, m = 's';
    otherwise,               m = '^';
end
end


function n = region_name(r)
if r == "front", n = 'Frontal'; else, n = 'Occipital'; end
end


function s = prop_label(p)
if p == "mus", s = '\mu_s'; else, s = 'Retardance'; end   % char, so it can sit in a cell label
end


function save_fig(fh, out_dir, base)
exportgraphics(fh, fullfile(out_dir, char(base + ".png")), 'Resolution', 300);
savefig(fh, fullfile(out_dir, char(base + ".fig")));
close(fh);
end