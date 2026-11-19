%% Shell profiles stratified by object size: EPVS vs vessel
%{
Reads the output of the parenchyma shell script (parenchyma_shells_*.mat,
or .csv) and plots mean shell profiles for small vs large objects.
Size variable is cfg.size_var: 'obj_vox' (object volume, voxels; default)
or 'obj_diam_um' (max inscribed diameter).

Size strata: thresholds on cfg.size_var are computed per region, POOLED
across EPVS and vessels and across all subjects/severities, so "small" and
"large" mean the same sizes for both classes and every group.
  ('epvs_median' / 'epvs_tertile' take thresholds from EPVS only.)
  split_mode 'median'  : small <= median, large > median
             'tertile' : small <= 33rd pct, large > 67th pct (middle dropped)
             'fixed'   : small <= fixed_thr, large > fixed_thr

Estimate per shell: mean over objects within subject, then mean over
subjects (subjects weighted equally; NaN if any subject lacks the shell).
CI: hierarchical bootstrap, objects resampled within subjects, subjects
fixed (2.5/97.5 percentiles). Describes uncertainty about the measured
subjects only. Thin lines = subject means (line style = subject).

Only rows with enough_voxels (and mus_ok / ret_ok when present) are used.
%}

clear; clc; close all;

%% Configuration
cfg.exploratory = true;
cfg.n_boot      = 100 + 400 * ~cfg.exploratory;   % 100 exploratory / 500 manuscript
cfg.seed        = 1;

cfg.data_dir = '/projectnb/npbssmic/ns/CAA/';
cfg.in_file  = fullfile(cfg.data_dir, 'parenchyma_shells_2026-10-07_pct.csv');  % '' = newest .mat

cfg.measures = {'ret_pct', 'mus_pct'};
cfg.labels   = {'Retardance (% from WM median)', '\mu_s (% from WM median)'};
cfg.labels_raw = {'Retardance (raw)', '\mu_s (raw)'};
cfg.flags    = {'ret_ok', 'mus_ok'};               % applied if the column exists
cfg.regions      = {'front', 'occip'};
cfg.region_names = {'Frontal', 'Occipital'};
cfg.severities   = {'control', 'severe', 'moderate'};   % moderate = CAA17, supplemental

cfg.split_mode = 'epvs_median';  % 'median' | 'tertile' | 'epvs_median' | 'epvs_tertile' | 'fixed'
cfg.fixed_thr  = 500;         % in units of size_var

cfg.size_var = 'obj_vox';     % 'obj_vox' | 'obj_diam_um'
cfg.vox_um3  = 20^3;          % voxel volume, for axis labels only
cfg.min_obj    = 5;           % min objects per subject x shell, else NaN

cfg.classes = ["EPVS", "vessel"];
cfg.colors  = [0.80 0.25 0.20;    % EPVS
               0.20 0.45 0.70];   % vessel
cfg.sub_ls  = {'--', ':', '-.'};

rng(cfg.seed, 'twister');
stamp = string(datetime('now', 'Format', 'yyyy-MM-dd_HHmm'));
tag = "";  if cfg.exploratory, tag = "_EXPLORATORY"; end
odir = fullfile(cfg.data_dir, 'figures', "shells_by_" + cfg.size_var + "_" + stamp + tag);
if ~exist(odir, 'dir'), mkdir(odir); end

%% Load
if isempty(cfg.in_file)
    f = dir(fullfile(cfg.data_dir, 'parenchyma_shells_*.mat'));
    f = f(~contains({f.name}, '_random'));
    assert(~isempty(f), 'No parenchyma_shells_*.mat found in %s', cfg.data_dir);
    [~, i] = max([f.datenum]);
    cfg.in_file = fullfile(f(i).folder, f(i).name);
end
fprintf('Reading %s\n', cfg.in_file);
[~, ~, ext] = fileparts(cfg.in_file);
if strcmpi(ext, '.mat')
    S = load(cfg.in_file, 'shells');  T = S.shells;  clear S
else
    T = readtable(cfg.in_file, 'TextType', 'string');
end
for v = ["subject", "severity", "region", "class", "object_id"]
    T.(v) = string(T.(v));
end

% Fall back to raw values if percentage-change columns are absent
vars = T.Properties.VariableNames;
for m = 1:numel(cfg.measures)
    if ~ismember(cfg.measures{m}, vars)
        raw = erase(cfg.measures{m}, '_pct');
        assert(ismember(raw, vars), 'Neither %s nor %s in file. Columns: %s', ...
               cfg.measures{m}, raw, strjoin(vars, ', '));
        warning('%s not in file; plotting raw %s instead', cfg.measures{m}, raw);
        cfg.measures{m} = raw;
        cfg.labels{m}   = cfg.labels_raw{m};
    end
end

T = T(tobool(T.enough_voxels) & ismember(T.class, cfg.classes), :);
T.d_mid = (T.d_inner_um + T.d_outer_um) / 2;
shells  = unique(T.d_mid)';

%% Size strata (shared thresholds per region)
[~, ia] = unique(T.object_id, 'stable');
O = T(ia, {'object_id', 'subject', 'severity', 'region', 'class', cfg.size_var});
O.objsize = O.(cfg.size_var);
T.objsize = T.(cfg.size_var);

T.stratum = strings(height(T), 1);
O.stratum = strings(height(O), 1);
thr = zeros(numel(cfg.regions), 2);
for r = 1:numel(cfg.regions)
    reg = string(cfg.regions{r});
    d  = O.objsize(O.region == reg);
    de = O.objsize(O.region == reg & O.class == "EPVS");
    if isempty(d), continue; end
    switch cfg.split_mode
        case 'median',       lo = median(d);           hi = lo;
        case 'tertile',      lo = prctile(d, 100/3);   hi = prctile(d, 200/3);
        case 'epvs_median',  lo = median(de);          hi = lo;
        case 'epvs_tertile', lo = prctile(de, 100/3);  hi = prctile(de, 200/3);
        case 'fixed',   lo = cfg.fixed_thr;         hi = lo;
        otherwise, error('Unknown split_mode: %s', cfg.split_mode);
    end
    thr(r, :) = [lo hi];
    fprintf('%s (%s): small <= %g, large > %g\n', cfg.region_names{r}, cfg.size_var, lo, hi);
    T = set_stratum(T, reg, lo, hi);
    O = set_stratum(O, reg, lo, hi);
end
strata = ["small", "large"];

% Object counts per subject x class x stratum
O2 = O(O.stratum ~= "", :);
counts = groupsummary(O2, {'severity', 'region', 'subject', 'class', 'stratum'});
writetable(counts, fullfile(odir, 'object_counts.csv'));
disp(counts);

%% Size distributions with thresholds
fig = figure('Color', 'w', 'Position', [100 100 900 380]);
tl = tiledlayout(1, numel(cfg.regions), 'TileSpacing', 'compact');
if strcmp(cfg.size_var, 'obj_vox')
    edges = logspace(log10(min(O.objsize)), log10(max(O.objsize)), 40);   % heavy-tailed
    xlab  = sprintf('Object volume (voxels; 1 voxel = %g µm^3)', cfg.vox_um3);
else
    edges = 0:20:(ceil(prctile(O.objsize, 99) / 20) * 20 + 20);
    xlab  = 'Max inscribed diameter (µm)';
end
for r = 1:numel(cfg.regions)
    nexttile; hold on
    reg = string(cfg.regions{r});
    for c = 1:numel(cfg.classes)
        d = O.objsize(O.region == reg & O.class == cfg.classes(c));
        histogram(d, edges, 'Normalization', 'probability', ...
                  'FaceColor', cfg.colors(c, :), 'FaceAlpha', 0.45, 'EdgeColor', 'none', ...
                  'DisplayName', sprintf('%s (n=%d)', cfg.classes(c), numel(d)));
    end
    xline(thr(r, 1), 'k--', 'HandleVisibility', 'off');
    if thr(r, 2) ~= thr(r, 1), xline(thr(r, 2), 'k--', 'HandleVisibility', 'off'); end
    if strcmp(cfg.size_var, 'obj_vox'), set(gca, 'XScale', 'log'); end
    title(cfg.region_names{r});  xlabel(xlab);  ylabel('Proportion');
    legend('Location', 'northeast', 'Box', 'off');
end
title(tl, sprintf('%s, all subjects (split: %s)', cfg.size_var, cfg.split_mode), 'Interpreter', 'none');
save_fig(fig, fullfile(odir, 'size_distribution'));

%% Shell profiles
summary = {};
figs = struct();
for m = 1:numel(cfg.measures)
    meas = cfg.measures{m};
    Tm = T(T.stratum ~= "" & ~isnan(T.(meas)), :);
    if ismember(cfg.flags{m}, Tm.Properties.VariableNames)
        Tm = Tm(tobool(Tm.(cfg.flags{m})), :);
    end
    all_ax = gobjects(0);

    for sv = 1:numel(cfg.severities)
        sev = string(cfg.severities{sv});
        Ts = Tm(Tm.severity == sev, :);
        if isempty(Ts), continue; end

        fig = figure('Color', 'w', 'Position', [100 100 950 720]);
        tl = tiledlayout(numel(cfg.regions), numel(strata), 'TileSpacing', 'compact');

        for r = 1:numel(cfg.regions)
            reg = string(cfg.regions{r});
            Tr = Ts(Ts.region == reg, :);
            subs = unique(Tr.subject);
            for q = 1:numel(strata)
                ax = nexttile; hold(ax, 'on');
                if isempty(Tr)
                    title(ax, sprintf('%s, %s: no data', cfg.region_names{r}, strata(q)));
                    axis(ax, 'off');  continue
                end
                n_obj = zeros(1, numel(cfg.classes));
                for c = 1:numel(cfg.classes)
                    Tc = Tr(Tr.class == cfg.classes(c) & Tr.stratum == strata(q), :);
                    n_obj(c) = numel(unique(Tc.object_id));
                    if isempty(Tc), continue; end
                    [est, lo, hi, subj] = hboot(Tc, subs, shells, meas, cfg);
                    col = cfg.colors(c, :);

                    band(ax, shells, lo, hi, col);
                    plot(ax, shells, est, '-', 'Color', col, 'LineWidth', 2, ...
                         'DisplayName', cfg.classes(c));
                    for k = 1:numel(subs)
                        plot(ax, shells, subj(k, :), cfg.sub_ls{k}, 'Color', [col 0.6], ...
                             'LineWidth', 0.8, 'HandleVisibility', 'off');
                    end

                    summary{end+1} = table( ...
                        repmat(string(meas), numel(shells), 1), repmat(sev, numel(shells), 1), ...
                        repmat(reg, numel(shells), 1), repmat(strata(q), numel(shells), 1), ...
                        repmat(cfg.classes(c), numel(shells), 1), shells(:), est(:), lo(:), hi(:), ...
                        repmat(n_obj(c), numel(shells), 1), ...
                        'VariableNames', {'measure', 'severity', 'region', 'stratum', 'class', ...
                                          'd_mid_um', 'est', 'ci_lo', 'ci_hi', 'n_obj'}); %#ok<SAGROW>
                end
                % Subject line-style key
                for k = 1:numel(subs)
                    plot(ax, NaN, NaN, cfg.sub_ls{k}, 'Color', [0.3 0.3 0.3], ...
                         'DisplayName', upper(subs(k)));
                end
                yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
                title(ax, sprintf('%s, %s (n_{EPVS}=%d, n_{ves}=%d)', ...
                      cfg.region_names{r}, strata(q), n_obj(1), n_obj(2)));
                xlabel(ax, 'Distance from structure (µm)');
                ylabel(ax, cfg.labels{m});
                xlim(ax, [0 max(cfg_edges(shells))]);
                legend(ax, 'Location', 'best', 'Box', 'off');
                all_ax(end+1) = ax; %#ok<SAGROW>
            end
        end
        title(tl, sprintf('%s: %s by %s (%s split, %d bootstraps%s)', ...
              sev, meas, cfg.size_var, cfg.split_mode, cfg.n_boot, strrep(tag, '_', ' ')), 'Interpreter', 'none');
        figs.(meas).(sev) = fig;
    end

    % Matched y-axes across all severities/regions/strata for this measure
    if ~isempty(all_ax), linkaxes(all_ax, 'y'); end
    for sv = 1:numel(cfg.severities)
        sev = cfg.severities{sv};
        if isfield(figs, meas) && isfield(figs.(meas), sev)
            save_fig(figs.(meas).(sev), fullfile(odir, sprintf('%s_%s_by_%s', meas, sev, cfg.size_var)));
        end
    end
end

if ~isempty(summary)
    writetable(vertcat(summary{:}), fullfile(odir, "shell_profiles_by_" + cfg.size_var + ".csv"));
end
save(fullfile(odir, 'cfg.mat'), 'cfg', 'thr');
fprintf('Saved to %s\n', odir);

%% ------------------------------------------------------------------------
function [est, lo, hi, subj] = hboot(Tc, subs, shells, meas, cfg)
% Objects resampled within subject; subjects fixed and weighted equally.
nS = numel(shells);  nSub = numel(subs);
subj = nan(nSub, nS);
boot = nan(cfg.n_boot, nS, nSub);
for k = 1:nSub
    Tk = Tc(Tc.subject == subs(k), :);
    if isempty(Tk), continue; end            % subject absent -> group NaN
    [ids, ~, oi] = unique(Tk.object_id);
    [~, si] = ismember(Tk.d_mid, shells);
    Y  = accumarray([oi si], Tk.(meas), [numel(ids) nS], @mean, NaN);
    nO = numel(ids);

    mk = mean(Y, 1, 'omitnan');
    mk(sum(~isnan(Y), 1) < cfg.min_obj) = NaN;
    subj(k, :) = mk;

    for b = 1:cfg.n_boot
        yb = Y(randi(nO, nO, 1), :);
        mb = mean(yb, 1, 'omitnan');
        mb(sum(~isnan(yb), 1) < cfg.min_obj) = NaN;
        boot(b, :, k) = mb;
    end
end
est = mean(subj, 1);                         % NaN if any subject missing
gb  = mean(boot, 3);
lo  = prctile(gb, 2.5, 1);
hi  = prctile(gb, 97.5, 1);
bad = mean(isnan(gb), 1) > 0.05 | isnan(est);  % unstable / missing shells
lo(bad) = NaN;  hi(bad) = NaN;
end

function Z = set_stratum(Z, reg, lo, hi)
inr = Z.region == reg;
Z.stratum(inr & Z.objsize <= lo) = "small";
Z.stratum(inr & Z.objsize >  hi) = "large";
end

function band(ax, x, lo, hi, col)
% Shaded CI, drawn per contiguous non-NaN run
ok = ~isnan(lo) & ~isnan(hi);
d  = diff([0 ok 0]);
st = find(d == 1);  en = find(d == -1) - 1;
for i = 1:numel(st)
    j = st(i):en(i);
    fill(ax, [x(j) fliplr(x(j))], [lo(j) fliplr(hi(j))], col, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end

function e = cfg_edges(shells)
% Outer edge of the last shell from midpoints (assumes uniform width)
w = median(diff(shells));
e = shells(end) + w / 2;
end

function b = tobool(x)
if islogical(x), b = x;
elseif isstring(x) || iscellstr(x), b = strcmpi(string(x), "true") | string(x) == "1";
else, b = logical(x);
end
end

function save_fig(fig, base)
exportgraphics(fig, base + ".png", 'Resolution', 300);
exportgraphics(fig, base + ".pdf", 'ContentType', 'vector');
end