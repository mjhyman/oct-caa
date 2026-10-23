%% Histology ribbon scatterplots: Bayesian fitted model + raw data
%{
Adapts the optical-property scatterplot script to histology (LHE, GFAP, CD68),
following the same logic: for each stain x region, overlay the EPVS and vessel
ribbons, and produce (1) each stage plotted separately and (2) the severe -
control difference plotted separately.

Sections:
  A) Bayesian fitted model, each stage separately (control, severe)
       source: R "*_stage{0,3}_fitted.csv"  (fit_mean + 95% CI ribbon)
  B) Bayesian fitted model, difference (severe - control)
       source: R "*_stage3_vs_stage0_model_diff.csv"  (diff_mean +/- diff_sd)
  C) Raw histology (from the ring .mat structs), each stage separately + diff
       source: average_within_region on the .mat structs (mean +/- SD)

Prerequisite for A/B: run the R histology Bayesian pipeline first so the
"*_fitted.csv" / "*_model_diff.csv" files exist under a
bayesian_stats_histo_<timestamp>/ run folder.
%}

clc; close all;

%% Shared setup: directories, iteration cells, styling
% --- Bayesian CSV location (R output) --------------------------------------
histo_stats = '/projectnb/npbssmic/ns/CAA/histology/stats/';
base_name   = 'histogram_matched_donuts_2026-May-05';   % R base_name (no ext)

% Pick the R run folder. Leave '' to auto-select the most recent run.
run_folder = '';   % e.g. 'bayesian_stats_histo_2026-05-05_12-00-00'
if isempty(run_folder)
    dd = dir(fullfile(histo_stats, 'bayesian_stats_histo_*'));
    dd = dd([dd.isdir]);
    assert(~isempty(dd), 'No bayesian_stats_histo_* run found in %s', histo_stats);
    [~, ix]    = max([dd.datenum]);
    run_folder = dd(ix).name;
    fprintf('Using Bayesian run: %s\n', run_folder);
end
stage_root = fullfile(histo_stats, run_folder, 'stage_vs_control'); % .../{region}/{src}/

% --- Figure output ---------------------------------------------------------
fig_out   = '/projectnb/npbssmic/ns/CAA/figures/histology_scatterplots/bayesian_fitted/';
save_flag = true;
if save_flag && ~exist(fig_out, 'dir'); mkdir(fig_out); end

% --- Iteration cells -------------------------------------------------------
stains       = {'lhe','gfap','cd68'};       % struct field / filename / output token
stain_title  = {'LHE','GFAP','CD68'};       % display

regs         = {'front','occip'};           % filename token + struct field
reg_title    = {'Frontal','Occipital'};
region_codes = {"_1","_7"};                 % raw .mat section codes (frontal/occipital)

srcs      = {'epvs','vessel'};              % R "src" token used in path/filenames
src_names = {'EPVS','Vessel'};
% EPVS = blue, Vessel = vermillion (orange)
epv_color = validatecolor('#1964B0');
ves_color = validatecolor('#DB5829');
src_color = {epv_color, ves_color};

% --- Axis units & per-stain y-limits (z-score) -----------------------------
% Histology ring radii are in voxels; relabel/convert here if your axis is um.
xlab = 'Distance (voxels)';

% Tune these to your data (start from the raw-plot ranges used previously).
ylims_fit  = struct('lhe',[-0.6 0.6], 'gfap',[-1.4 1.4], 'cd68',[-1.4 1.4]);
ylims_diff = struct('lhe',[-0.6 0.6], 'gfap',[-1.4 1.4], 'cd68',[-1.4 1.4]);

% --- Common style params ---------------------------------------------------
dot_size = 100;
lw       = 4;
fsize    = 30;

% Import raw data
fprintf('\nImporting raw data\n')
hdir = '/projectnb/npbssmic/ns/CAA/histology/';
lhe  = load(fullfile(hdir,'lhe_rings_21-Nov-2025.mat'));  lhe  = lhe.lhe;
gfap = load(fullfile(hdir,'gfap_rings_21-Nov-2025.mat')); gfap = gfap.gfap;
cd68 = load(fullfile(hdir,'cd68_rings_21-Nov-2025.mat')); cd68 = cd68.cd68;
raw_cells = {lhe, gfap, cd68};                 % aligned with `stains`

%% SECTION A - Bayesian fitted model: each stage separately (control, severe)
fit_stages = {'control','severe'};
fit_k      = [0 3];
fit_lab    = {'Control','Severe'};

% Load "*_stageK_fitted.csv" -> bf.(stage).(stain).(region).(src)
bf = struct();
for si = 1:numel(fit_stages)
    for j = 1:numel(stains)
        for rr = 1:numel(regs)
            for ss = 1:numel(srcs)
                fdir  = fullfile(stage_root, regs{rr}, srcs{ss});
                fname = sprintf('%s__%s_%s_%s_stage%d_fitted.csv', ...
                                base_name, stains{j}, regs{rr}, srcs{ss}, fit_k(si));
                try
                    tmp = readtable(fullfile(fdir, fname));
                catch
                    fprintf('\tSkipping (fitted) %s %s %s %s\n', ...
                            fit_stages{si}, stains{j}, regs{rr}, srcs{ss});
                    continue
                end
                bf.(fit_stages{si}).(stains{j}).(regs{rr}).(srcs{ss}).d     = tmp.distance(:);
                bf.(fit_stages{si}).(stains{j}).(regs{rr}).(srcs{ss}).mu    = tmp.fit_mean(:);
                bf.(fit_stages{si}).(stains{j}).(regs{rr}).(srcs{ss}).lower = tmp.lower(:);
                bf.(fit_stages{si}).(stains{j}).(regs{rr}).(srcs{ss}).upper = tmp.upper(:);
            end
        end
    end
end

% Plot: one figure per stain x region x stage
for si = 1:numel(fit_stages)
    for j = 1:numel(stains)
        for rr = 1:numel(regs)
            hs = gobjects(0); names = {}; have_any = false;

            figure; hold on;
            yline(0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

            for ss = 1:numel(srcs)
                try
                    S = bf.(fit_stages{si}).(stains{j}).(regs{rr}).(srcs{ss});
                catch
                    continue
                end
                h = add_ci_ribbon(S.d, S.mu, S.lower, S.upper, src_color{ss}, lw, dot_size);
                hs(end+1)    = h;             %#ok<SAGROW>
                names{end+1} = src_names{ss}; %#ok<SAGROW>
                have_any = true;
            end
            if ~have_any, close(gcf); continue; end

            ylim(ylims_fit.(stains{j}));
            xlabel(xlab);
            ylabel(sprintf('%s (z-score)', stain_title{j}));
            title(sprintf('Fitted Bayesian Model: %s %s %s (EPVS vs Vessel)', ...
                          fit_lab{si}, stain_title{j}, reg_title{rr}));
            legend(hs, names, 'Location', 'best');
            finalize_axes(fsize); drawnow;

            fname = fullfile(fig_out, sprintf('histo_fitted_%s_%s_%s', ...
                                              fit_stages{si}, stains{j}, regs{rr}));
            save_fig(fname, save_flag);
        end
    end
end


%% SECTION B - Bayesian fitted model: difference (severe - control)
diff_k = 3;   % severe vs control

% Load "*_stage3_vs_stage0_model_diff.csv" -> bd.(stain).(region).(src)
bd = struct();
for j = 1:numel(stains)
    for rr = 1:numel(regs)
        for ss = 1:numel(srcs)
            fdir  = fullfile(stage_root, regs{rr}, srcs{ss});
            fname = sprintf('%s__%s_%s_%s_stage%d_vs_stage0_model_diff.csv', ...
                            base_name, stains{j}, regs{rr}, srcs{ss}, diff_k);
            try
                tmp = readtable(fullfile(fdir, fname));
            catch
                fprintf('\tSkipping (diff) %s %s %s\n', stains{j}, regs{rr}, srcs{ss});
                continue
            end
            bd.(stains{j}).(regs{rr}).(srcs{ss}).d  = tmp.distance(:);
            bd.(stains{j}).(regs{rr}).(srcs{ss}).mu = tmp.diff_mean(:);
            bd.(stains{j}).(regs{rr}).(srcs{ss}).sd = tmp.diff_sd(:);
        end
    end
end

% Plot: one figure per stain x region
for j = 1:numel(stains)
    for rr = 1:numel(regs)
        hs = gobjects(0); names = {}; have_any = false;

        figure; hold on;
        yline(0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

        for ss = 1:numel(srcs)
            try
                S = bd.(stains{j}).(regs{rr}).(srcs{ss});
            catch
                continue
            end
            h = add_ribbon(S.d, S.mu, S.sd, src_color{ss}, lw, dot_size);
            hs(end+1)    = h;             %#ok<SAGROW>
            names{end+1} = src_names{ss}; %#ok<SAGROW>
            have_any = true;
        end
        if ~have_any, close(gcf); continue; end

        ylim(ylims_diff.(stains{j}));
        xlabel(xlab);
        ylabel(sprintf('%s: severe - control (z-score)', stain_title{j}));
        title(sprintf('Fitted Bayesian Model: %s %s (Severe - Control)', ...
                      stain_title{j}, reg_title{rr}));
        legend(hs, names, 'Location', 'best');
        finalize_axes(fsize); drawnow;

        fname = fullfile(fig_out, sprintf('histo_fitted_diff_%s_%s', stains{j}, regs{rr}));
        save_fig(fname, save_flag);
    end
end


%% SECTION C - Raw histology (ring structs): stages separately + difference
% Reuses the group-average logic from the raw histology script, but draws the
% EPVS/vessel overlay in the same ribbon style as the Bayesian sections above.

radii_sm = [40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486];
ctl_subs = {'CAA_26','CAA_6'};
sev_subs = {'CAA_22','CAA_25'};

% Manually set y-axis limits
ylims_raw.cd68 = [-0.5,1];
ylims_raw.gfap = [-0.5,1];
ylims_raw.lhe  = [-0.5,0.5];

for j = 1:numel(stains)
    ctl_group = extract_subjects(raw_cells{j}, ctl_subs);
    sev_group = extract_subjects(raw_cells{j}, sev_subs);

    for rr = 1:numel(regs)
        rcode = region_codes{rr};
        ctl = average_within_region(ctl_group, rcode, "rad40");
        sev = average_within_region(sev_group, rcode, "rad40");

        % --- each stage separately (mean +/- SD ribbon) ------------------
        grp     = {ctl, sev};
        grp_lab = {'Control','Severe'};
        grp_tag = {'control','severe'};
        for gi = 1:2
            G = grp{gi};
            figure; hold on;
            yline(0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
            h1 = add_ribbon(radii_sm, G.mean_epvs, G.std_epvs, src_color{1}, lw, dot_size);
            h2 = add_ribbon(radii_sm, G.mean_ves,  G.std_ves,  src_color{2}, lw, dot_size);

            ylim(ylims_raw.(stains{j}));
            xlabel(xlab);
            ylabel(sprintf('%s (z-score)', stain_title{j}));
            title(sprintf('Raw Histology: %s %s %s (EPVS vs Vessel, n=%d)', ...
                          grp_lab{gi}, stain_title{j}, reg_title{rr}, G.nsubs));
            legend([h1 h2], src_names, 'Location', 'best');
            finalize_axes(fsize); drawnow;

            fname = fullfile(fig_out, sprintf('histo_raw_%s_%s_%s', ...
                                              grp_tag{gi}, stains{j}, regs{rr}));
            save_fig(fname, save_flag);
        end

        % --- difference (severe - control) -------------------------------
        d = struct();
        [d.std_epvs, d.mean_epvs, d.std_ves, d.mean_ves] = stage_diff(ctl, sev);

        figure; hold on;
        yline(0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
        h1 = add_ribbon(radii_sm, d.mean_epvs, d.std_epvs, src_color{1}, lw, dot_size);
        h2 = add_ribbon(radii_sm, d.mean_ves,  d.std_ves,  src_color{2}, lw, dot_size);

        ylim(ylims_diff.(stains{j}));
        xlabel(xlab);
        ylabel(sprintf('%s: severe - control (z-score)', stain_title{j}));
        title(sprintf('Raw Histology: %s %s (Severe - Control)', ...
                      stain_title{j}, reg_title{rr}));
        legend([h1 h2], src_names, 'Location', 'best');
        finalize_axes(fsize); drawnow;

        fname = fullfile(fig_out, sprintf('histo_raw_diff_%s_%s', stains{j}, regs{rr}));
        save_fig(fname, save_flag);
    end
end


%% SECTION D - Model fit over raw data (shared axes)
% One figure per stain x region x stage. Shaded band + solid line = Bayesian
% fitted mean (95% CI); points + error bars = raw group mean +/- SD. Lets you
% see how closely the fitted model tracks the measured data for each source.
% Self-contained: loads raw structs and fitted CSVs itself (does not need A-C).

% Ensure raw ring structs are available (idempotent)
if ~exist('raw_cells', 'var')
    hdir = '/projectnb/npbssmic/ns/CAA/histology/';
    lhe  = load(fullfile(hdir,'lhe_rings_21-Nov-2025.mat'));  lhe  = lhe.lhe;
    gfap = load(fullfile(hdir,'gfap_rings_21-Nov-2025.mat')); gfap = gfap.gfap;
    cd68 = load(fullfile(hdir,'cd68_rings_21-Nov-2025.mat')); cd68 = cd68.cd68;
    raw_cells = {lhe, gfap, cd68};
    radii_sm  = [40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486];
    ctl_subs  = {'CAA_26','CAA_6'};
    sev_subs  = {'CAA_22','CAA_25'};
end

cmp_stages = {'control','severe'};
cmp_k      = [0 3];
cmp_lab    = {'Control','Severe'};

for si = 1:numel(cmp_stages)
    for j = 1:numel(stains)
        ctl_group = extract_subjects(raw_cells{j}, ctl_subs);
        sev_group = extract_subjects(raw_cells{j}, sev_subs);

        for rr = 1:numel(regs)
            rcode = region_codes{rr};
            if cmp_k(si) == 0
                R = average_within_region(ctl_group, rcode, "rad40");
            else
                R = average_within_region(sev_group, rcode, "rad40");
            end

            hs = gobjects(0); names = {}; have_any = false;
            figure; hold on;
            yline(0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

            for ss = 1:numel(srcs)
                % Fitted curve for this source (may be missing)
                F = load_fitted_csv(stage_root, base_name, ...
                                    stains{j}, regs{rr}, srcs{ss}, cmp_k(si));
                % Raw mean +/- SD for this source
                if strcmp(srcs{ss}, 'epvs')
                    rmu = R.mean_epvs; rsd = R.std_epvs;
                else
                    rmu = R.mean_ves;  rsd = R.std_ves;
                end

                [hFit, hRaw] = add_fit_over_raw(radii_sm, rmu, rsd, F, src_color{ss}, lw);
                if ~isempty(hFit)
                    hs(end+1) = hFit;                                %#ok<SAGROW>
                    names{end+1} = sprintf('%s (fit)', src_names{ss}); %#ok<SAGROW>
                    have_any = true;
                end
                if ~isempty(hRaw)
                    hs(end+1) = hRaw;                                %#ok<SAGROW>
                    names{end+1} = sprintf('%s (raw)', src_names{ss}); %#ok<SAGROW>
                    have_any = true;
                end
            end
            if ~have_any, close(gcf); continue; end

            ylim(ylims_fit.(stains{j}));
            xlabel(xlab);
            ylabel(sprintf('%s (z-score)', stain_title{j}));
            title(sprintf('Fit vs Raw: %s %s %s (n=%d)', ...
                          cmp_lab{si}, stain_title{j}, reg_title{rr}, R.nsubs));
            legend(hs, names, 'Location', 'best');
            finalize_axes(fsize); drawnow;

            fname = fullfile(fig_out, sprintf('histo_fit_vs_raw_%s_%s_%s', ...
                                              cmp_stages{si}, stains{j}, regs{rr}));
            save_fig(fname, save_flag);
        end
    end
end


%% ------------------------------------------------------------------------
%% Local function: SD ribbon + mean line + scatter
function h = add_ribbon(d, mu, sd, c, lw, dot_size)
    d = d(:); mu = mu(:); sd = sd(:);
    fill([d; flipud(d)], [mu - sd; flipud(mu + sd)], c, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h = plot(d, mu, '-', 'Color', c, 'LineWidth', lw);
    scatter(d, mu, dot_size, c, 'filled', 'HandleVisibility', 'off');
end

%% Local function: 95% CI ribbon + mean line + scatter
function h = add_ci_ribbon(d, mu, lower, upper, c, lw, dot_size)
    d = d(:); mu = mu(:); lower = lower(:); upper = upper(:);
    fill([d; flipud(d)], [lower; flipud(upper)], c, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h = plot(d, mu, '-', 'Color', c, 'LineWidth', lw);
    scatter(d, mu, dot_size, c, 'filled', 'HandleVisibility', 'off');
end

%% Local function: load one fitted CSV (returns [] if missing)
function S = load_fitted_csv(stage_root, base_name, stain, region, src, k)
    S = [];
    fdir  = fullfile(stage_root, region, src);
    fname = sprintf('%s__%s_%s_%s_stage%d_fitted.csv', base_name, stain, region, src, k);
    try
        tmp = readtable(fullfile(fdir, fname));
    catch
        return
    end
    S.d     = tmp.distance(:);
    S.mu    = tmp.fit_mean(:);
    S.lower = tmp.lower(:);
    S.upper = tmp.upper(:);
end

%% Local function: fitted CI band/line over raw mean +/- SD points
function [hFit, hRaw] = add_fit_over_raw(d_raw, mu_raw, sd_raw, F, c, lw)
    hFit = gobjects(0);
    hRaw = gobjects(0);

    % Fitted: shaded 95% CI band (hidden) + solid mean line (legend proxy)
    if ~isempty(F)
        dF = F.d(:); muF = F.mu(:); lo = F.lower(:); up = F.upper(:);
        fill([dF; flipud(dF)], [lo; flipud(up)], c, ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hFit = plot(dF, muF, '-', 'Color', c, 'LineWidth', lw);
    end

    % Raw: mean markers with +/- SD error bars (legend proxy)
    d_raw = d_raw(:); mu_raw = mu_raw(:); sd_raw = sd_raw(:);
    hRaw = errorbar(d_raw, mu_raw, sd_raw, 'o', 'Color', c, ...
                    'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 8, 'LineWidth', 1.5, 'LineStyle', 'none', ...
                    'CapSize', 8);
end

%% Local function: shared axis styling
function finalize_axes(fsize)
    grid on;
    set(gca, 'fontsize', fsize);
    set(get(gca, 'Title'), 'FontSize', 10);
    set(gca, 'fontname', 'Arial');
    set(gca, 'XColor', 'k', 'YColor', 'k');
    set(gca, 'TickLength', [0.04, 0.04]);
end

%% Local function: save png/fig/svg
function save_fig(fname, save_flag)
    if save_flag
        saveas(gcf, fname, 'png');
        saveas(gcf, fname, 'fig');
        exportgraphics(gcf, strcat(fname, '.svg'), 'Resolution', 600);
    end
end

%% Local function: extract a group's sections by subject ID
function sub = extract_subjects(stain, subs)
% Handles both baseName conventions: 'CAA_26_1' (LHE) and 'CAA26_1_GFAP_10x'.
names = string({stain.baseName});
keep  = false(size(names));
for i = 1:numel(subs)
    s_us = string(subs{i});
    s_no = erase(s_us, "_");
    keep = keep | contains(names, s_us) | contains(names, s_no);
end
sub = stain(keep);
end

%% Local function: average within one region across subjects
function [stats] = average_within_region(stain, region, rad_size)
stats = struct();
if isempty(stain)
    warning('average_within_region: empty group passed; returning empty stats.');
    return;
end

nrad     = length(fields(stain(1).(rad_size)));
rad_name = fields(stain(1).(rad_size));

names = {stain.baseName};
if endsWith(names{1}, "_10x")
    tf = contains(names, strcat(region, "_"));
else
    tf = endsWith(names, region);
end
stain_reg = stain(tf);
nsubs     = numel(stain_reg);

if nsubs == 0
    warning('average_within_region: no sections for region %s; returning NaNs.', region);
    stats.n_epvs = nan(nrad,1); stats.n_ves = nan(nrad,1);
    stats.mean_epvs = nan(nrad,1); stats.std_epvs = nan(nrad,1);
    stats.mean_ves  = nan(nrad,1); stats.std_ves  = nan(nrad,1);
    stats.nsubs = 0;
    return;
end

n_epvs    = zeros(nrad,nsubs);  n_ves    = zeros(nrad,nsubs);
mean_epvs = zeros(nrad,nsubs);  mean_ves = zeros(nrad,nsubs);

for s = 1:nsubs
    for ii = 1:nrad
        epvs = stain_reg(s).(rad_size).(rad_name{ii}).exp;
        ves  = stain_reg(s).(rad_size).(rad_name{ii}).ctl;
        n_epvs(ii,s)    = numel(epvs);
        n_ves(ii,s)     = numel(ves);
        mean_epvs(ii,s) = mean(epvs,'omitnan');
        mean_ves(ii,s)  = mean(ves, 'omitnan');
    end
end

stats.n_epvs    = sum(n_epvs, 2);
stats.n_ves     = sum(n_ves,  2);
stats.mean_epvs = mean(mean_epvs, 2, 'omitnan');
stats.std_epvs  = std(mean_epvs, 0, 2, 'omitnan');
stats.mean_ves  = mean(mean_ves,  2, 'omitnan');
stats.std_ves   = std(mean_ves,  0, 2, 'omitnan');
stats.nsubs     = nsubs;
end

%% Local function: difference between severities (severe - control)
function [std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(ctl, exp)
mean_epvs = exp.mean_epvs - ctl.mean_epvs;
std_epvs  = sqrt((exp.std_epvs).^2 ./ exp.n_epvs + (ctl.std_epvs).^2 ./ ctl.n_epvs);
mean_ves  = exp.mean_ves - ctl.mean_ves;
std_ves   = sqrt((exp.std_ves).^2 ./ exp.n_ves + (ctl.std_ves).^2 ./ ctl.n_ves);
end