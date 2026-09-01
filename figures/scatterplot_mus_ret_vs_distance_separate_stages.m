%% Ribbon plots overlaying stage-vs-control differences for EPVS and vessel
% For each optical property x region x nonzero stage (mild/moderate/severe),
% overlay two ribbons:
%   - EPVS   : (stage k - control) using the model posterior difference
%   - Vessel : (stage k - control) using the model posterior difference
% Mean and standard deviation are read from the "*_model_diff.csv" files
% produced by the R stage_vs_control pipeline:
%   columns: distance, diff_mean, diff_sd, lower, upper, n_draws
%
% TODO: Update this to plot the control and severe separately

clc; close all;

%% Directories and Data
beta_dir  = '/projectnb/npbssmic/ns/CAA/beta_stats/';
fig_out   = '/projectnb/npbssmic/ns/CAA/figures/fig3_op_vs_dist/';
base_name = 'caa_all_radii_percentage_diff_40um_donut_5-Mar-2026';
stage_root = fullfile(beta_dir, base_name, 'stage_vs_control');  % .../{region}/{src}/
save_flag = true;

if save_flag && ~exist(fig_out, 'dir'); mkdir(fig_out); end

%% Plot (stage k - control)

%%% Iteration cells
% Nonzero stages only (stage 0 is the control reference)
nz_stages  = {'stage1','stage2','stage3'};
stage_k    = [1 2 3];
stage_lab  = {'Mild','Moderate','Severe'};

op        = {'mus','ret'};                  % short names (struct fields / filenames out)
op_names  = {'scattering','retardance'};    % R "prop" names used in the CSV filenames
op_title  = {'\mu_s','Retardance'};

regs      = {'front','occip'};
reg_title = {'Frontal','Occipital'};

srcs      = {'epvs','vessel'};              % R "src" names used in the CSV path/filenames
src_names = {'EPVS','Vessel'};
% EPVS = blue, Vessel = vermillion (orange)
epv_color = validatecolor('#1964B0'); % blue 
ves_color = validatecolor('#DB5829'); % vermillion
src_color = {epv_color,ves_color};

%%% Import model posterior differences from "*_model_diff.csv"
% md.(stage).(op).(region).(src) = struct with fields d, mu, sd
md = struct();
for si = 1:numel(nz_stages)
    kk = stage_k(si);
    for j = 1:numel(op)
        for rr = 1:numel(regs)
            for ss = 1:numel(srcs)
                fdir  = fullfile(stage_root, regs{rr}, srcs{ss});
                fname = sprintf('%s__%s_%s_%s_stage%d_vs_stage0_model_diff.csv', ...
                                base_name, op_names{j}, regs{rr}, srcs{ss}, kk);
                fullpath = fullfile(fdir, fname);
                try
                    tmp = readtable(fullpath);
                catch
                    fprintf('\tSkipping %s %s %s %s\n', ...
                            nz_stages{si}, op{j}, regs{rr}, srcs{ss});
                    continue
                end
                md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss}).d  = tmp.distance(:);
                md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss}).mu = tmp.diff_mean(:);
                md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss}).sd = tmp.diff_sd(:);
            end
        end
    end
end

%%% Plot: one figure per optical property x region x nonzero stage
% Manually set ylims from data
ylims = [-20,10;-15,10];
dot_size = 100;
lw       = 4;
fsize    = 30;
% Iterate stages 
for si = 1:numel(nz_stages)
    for j = 1:numel(op)
        for rr = 1:numel(regs)

            % Collect whichever sources are available for this combination
            hs    = gobjects(0);
            names = {};
            have_any = false;

            figure; hold on;
            yline(0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

            for ss = 1:numel(srcs)
                try
                    S = md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss});
                catch
                    continue
                end
                h = add_ribbon(S.d, S.mu, S.sd, src_color{ss}, lw, dot_size);
                hs(end+1)    = h;            %#ok<SAGROW>
                names{end+1} = src_names{ss}; %#ok<SAGROW>
                have_any = true;
            end

            if ~have_any
                close(gcf);
                continue
            end

            % Axes / labels
            ylim(ylims(j, :));
            xlabel('Distance (\mum)');
            ylabel(sprintf('%s: %s - control', op_title{j}, lower(stage_lab{si})));
            tstr = sprintf('%s %s %s (EPVS vs Vessel)', ...
                           stage_lab{si}, op_title{j}, reg_title{rr});
            title(tstr);
            legend(hs, names, 'Location', 'best');

            % Figure properties (kept from original)
            grid on;
            set(gca, 'fontsize', fsize)
            set(get(gca, 'Title'), 'FontSize', 10);
            set(gca, 'fontname', 'Arial');
            set(gca, 'XColor', 'k', 'YColor', 'k');
            set(gca, 'TickLength', [0.04, 0.04]);
            pause(1)

            % Save as .pdf .png .fig
            if save_flag
                fname = fullfile(fig_out, ...
                                 sprintf('pdiff_overlay_%s_%s_%s', ...
                                         nz_stages{si}, op{j}, regs{rr}));
                saveas(gcf, fname, 'png');
                saveas(gcf, fname, 'fig');
                exportgraphics(gcf, strcat(fname,'.svg'),'Resolution',600);
            end
        end
    end
end

%% Plot Each stage separately
%%% Iteration cells
% Nonzero stages only (stage 0 is the control reference)
nz_stages  = {'control','severe'};
stage_k    = [0 3];
stage_lab  = {'Control','Severe'};

op        = {'mus','ret'};                  % short names (struct fields / filenames out)
op_names  = {'scattering','retardance'};    % R "prop" names used in the CSV filenames
op_title  = {'\mu_s','Retardance'};

regs      = {'front','occip'};
reg_title = {'Frontal','Occipital'};

srcs      = {'epvs','vessel'};              % R "src" names used in the CSV path/filenames
src_names = {'EPVS','Vessel'};
% EPVS = blue, Vessel = vermillion (orange)
epv_color = validatecolor('#1964B0'); % blue 
ves_color = validatecolor('#DB5829'); % vermillion
src_color = {epv_color,ves_color};

%%% Import model posterior differences from "*_model_diff.csv"
% md.(stage).(op).(region).(src) = struct with fields d, mu, sd
md = struct();
for si = 1:numel(nz_stages)
    kk = stage_k(si);
    for j = 1:numel(op)
        for rr = 1:numel(regs)
            for ss = 1:numel(srcs)
                fdir  = fullfile(stage_root, regs{rr}, srcs{ss});
                fname = sprintf('%s__%s_%s_%s_stage%d_fitted.csv', ...
                            base_name, op_names{j}, regs{rr}, srcs{ss}, kk);
                fullpath = fullfile(fdir, fname);
                try
                    tmp = readtable(fullpath);
                catch
                    fprintf('\tSkipping %s %s %s %s\n', ...
                            nz_stages{si}, op{j}, regs{rr}, srcs{ss});
                    continue
                end
                md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss}).d  = tmp.distance(:);
                md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss}).mu = tmp.fit_mean(:);
                md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss}).lower = tmp.lower(:);
                md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss}).upper = tmp.upper(:);
            end
        end
    end
end

%%% Plot: one figure per optical property x region x nonzero stage
% Manually set ylims from data
ylims = [-10,20;-20,5];
dot_size = 100;
lw       = 4;
fsize    = 30;
% Iterate stages 
for si = 1:numel(nz_stages)
    for j = 1:numel(op)
        for rr = 1:numel(regs)

            % Collect whichever sources are available for this combination
            hs    = gobjects(0);
            names = {};
            have_any = false;

            figure; hold on;
            yline(0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

            for ss = 1:numel(srcs)
                try
                    S = md.(nz_stages{si}).(op{j}).(regs{rr}).(srcs{ss});
                catch
                    continue
                end
                h = add_ci_ribbon(S.d, S.mu, S.lower, S.upper,...
                                  src_color{ss}, lw, dot_size);
                hs(end+1)    = h;            %#ok<SAGROW>
                names{end+1} = src_names{ss}; %#ok<SAGROW>
                have_any = true;
            end

            if ~have_any
                close(gcf);
                continue
            end

            % Axes / labels
            ylim(ylims(j, :));
            xlabel('Distance (\mum)');
            ylabel(sprintf('%s: %s', op_title{j}, lower(stage_lab{si})));
            tstr = sprintf('Fitted Bayesian Model: %s %s %s (EPVS vs Vessel)', ...
                           stage_lab{si}, op_title{j}, reg_title{rr});
            title(tstr);
            legend(hs, names, 'Location', 'best');

            % Figure properties (kept from original)
            grid on;
            set(gca, 'fontsize', fsize)
            set(get(gca, 'Title'), 'FontSize', 10);
            set(gca, 'fontname', 'Arial');
            set(gca, 'XColor', 'k', 'YColor', 'k');
            set(gca, 'TickLength', [0.04, 0.04]);
            pause(1)

            % Save as .pdf .png .fig
            if save_flag
                fname = fullfile(fig_out, ...
                                 sprintf('fitted_bayesian_scatterplot_%s_%s_%s', ...
                                         nz_stages{si}, op{j}, regs{rr}));
                saveas(gcf, fname, 'png');
                saveas(gcf, fname, 'fig');
                exportgraphics(gcf, strcat(fname,'.svg'),'Resolution',600);
            end
        end
    end
end

%% Local function: add standard deviation ribbon + mean line + scatter
function h = add_ribbon(d, mu, sd, c, lw, dot_size)
    d  = d(:); mu = mu(:); sd = sd(:);
    % Shaded mean +/- 1 SD ribbon (not shown in legend)
    fill([d; flipud(d)], [mu - sd; flipud(mu + sd)], c, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    % Mean line (this handle is used for the legend)
    h = plot(d, mu, '-', 'Color', c, 'LineWidth', lw);
    % Mean markers (not shown in legend)
    scatter(d, mu, dot_size, c, 'filled', 'HandleVisibility', 'off');
end

%% Local function: add confidence interval ribbon + mean line + scatter
function h = add_ci_ribbon(d, mu, lower, upper, c, lw, dot_size)
    d  = d(:); mu = mu(:); lower = lower(:); upper = upper(:);
    % Shaded mean w/ 95% confidence interval ribbon
    fill([d; flipud(d)], [lower; flipud(upper)], c, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    % Mean line (this handle is used for the legend)
    h = plot(d, mu, '-', 'Color', c, 'LineWidth', lw);
    % Mean markers (not shown in legend)
    scatter(d, mu, dot_size, c, 'filled', 'HandleVisibility', 'off');
end