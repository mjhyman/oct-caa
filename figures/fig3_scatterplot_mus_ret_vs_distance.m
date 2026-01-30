%% Scatter plot of mus + ret vs. distance
% Use the raw means of the data set
% Use Bayesian statistical model for standard deviations
clear; clc; close all;

%% Directories and Data
% Input directories
data_dir = '/projectnb/npbssmic/ns/CAA/';
beta_dir = '/projectnb/npbssmic/ns/CAA/beta_stats/';
% Directory to store csv
fig_out = '/projectnb/npbssmic/ns/CAA/figures/fig3_op_vs_dist/';
% Directoies w/ standard deviation tables for each method
raw_dir = 'caa_all_radii_40um_donut_14-01-2026';
med_dir = 'caa_all_radii_median_subtracted_40um_donut_15-01-2026';
per_dir = 'caa_all_radii_percentage_diff_40um_donut_13-01-2026';
% Raw data CSV names
raw = 'caa_all_radii_40um_donut_14-01-2026.xlsx';
med = 'caa_all_radii_median_subtracted_40um_donut_15-01-2026.xlsx';
per = 'caa_all_radii_percentage_diff_40um_donut_13-01-2026.xlsx';
% Save output
save_flag = true;

%% Import Std Dev. CSVs for all three methods
% standard deviation struct
sd = struct();
% cells for iterating
method = {'raw','med','per'};
op = {'mus','ret'};
regs = {'front','occip'};
% Directories
dirs.raw = raw_dir;
dirs.med = med_dir;
dirs.per = per_dir;
% Filenames
fnames.raw = raw(1:end-5);
fnames.med = med(1:end-5);
fnames.per = per(1:end-5);

% Iterate method
for ii = 1:3
    % Iterate optical property
    for j = 1:2
        % Iterate region
        for k = 1:2
            % Create filename
            fname = strcat(dirs.(method{ii}),'__',op{j},'_',regs{k},'_summary_stats.csv');
            sd.(method{ii}).(op{j}).(regs{k}) =...
                readtable(fullfile(beta_dir, dirs.(method{ii}), fname));
        end
    end
end


%% Import CSVs for all three methods

% import raw spreadsheet from data_dir
raw_mus = readtable(fullfile(data_dir, raw),"Sheet",'scattering');
raw_ret = readtable(fullfile(data_dir, raw),"Sheet",'retardance');

% Import median-subtracted CSV
med_mus = readtable(fullfile(data_dir, med),"Sheet",'scattering');
med_ret = readtable(fullfile(data_dir, med),"Sheet",'retardance');

% Import percentage difference CSV
per_mus = readtable(fullfile(data_dir, per),"Sheet",'scattering');
per_ret = readtable(fullfile(data_dir, per),"Sheet",'retardance');

%% Extract front and occipital lobes for each dataset

% Create structs for raw, med, per
data = struct();

% Raw data (unmodified)
raw_mus_front = raw_mus(strcmp(raw_mus.Region, 'front'), :);
raw_mus_occip = raw_mus(strcmp(raw_mus.Region, 'occip'), :);
raw_ret_front = raw_ret(strcmp(raw_ret.Region, 'front'), :);
raw_ret_occip = raw_ret(strcmp(raw_ret.Region, 'occip'), :);
data.raw.mus.front.exp = raw_mus_front(strcmp(raw_mus_front.Groups,'experimental'),:);
data.raw.mus.front.ctl = raw_mus_front(strcmp(raw_mus_front.Groups,'control'),:);
data.raw.mus.occip.exp = raw_mus_occip(strcmp(raw_mus_occip.Groups,'experimental'),:);
data.raw.mus.occip.ctl = raw_mus_occip(strcmp(raw_mus_occip.Groups,'control'),:);
data.raw.ret.front.exp = raw_ret_front(strcmp(raw_ret_front.Groups,'experimental'),:);
data.raw.ret.front.ctl = raw_ret_front(strcmp(raw_ret_front.Groups,'control'),:);
data.raw.ret.occip.exp = raw_ret_occip(strcmp(raw_ret_occip.Groups,'experimental'),:);
data.raw.ret.occip.ctl = raw_ret_occip(strcmp(raw_ret_occip.Groups,'control'),:);

% Median-offset data
med_mus_front = med_mus(strcmp(med_mus.Region, 'front'), :);
med_mus_occip = med_mus(strcmp(med_mus.Region, 'occip'), :);
med_ret_front = med_ret(strcmp(med_ret.Region, 'front'), :);
med_ret_occip = med_ret(strcmp(med_ret.Region, 'occip'), :);
data.med.mus.front.exp = med_mus_front(strcmp(med_mus_front.Groups,'experimental'),:);
data.med.mus.front.ctl = med_mus_front(strcmp(med_mus_front.Groups,'control'),:);
data.med.mus.occip.exp = med_mus_occip(strcmp(med_mus_occip.Groups,'experimental'),:);
data.med.mus.occip.ctl = med_mus_occip(strcmp(med_mus_occip.Groups,'control'),:);
data.med.ret.front.exp = med_ret_front(strcmp(med_ret_front.Groups,'experimental'),:);
data.med.ret.front.ctl = med_ret_front(strcmp(med_ret_front.Groups,'control'),:);
data.med.ret.occip.exp = med_ret_occip(strcmp(med_ret_occip.Groups,'experimental'),:);
data.med.ret.occip.ctl = med_ret_occip(strcmp(med_ret_occip.Groups,'control'),:);

% Percentage change
per_mus_front = per_mus(strcmp(per_mus.Region, 'front'), :);
per_mus_occip = per_mus(strcmp(per_mus.Region, 'occip'), :);
per_ret_front = per_ret(strcmp(per_ret.Region, 'front'), :);
per_ret_occip = per_ret(strcmp(per_ret.Region, 'occip'), :);
data.per.mus.front.exp = per_mus_front(strcmp(per_mus_front.Groups,'experimental'),:);
data.per.mus.front.ctl = per_mus_front(strcmp(per_mus_front.Groups,'control'),:);
data.per.mus.occip.exp = per_mus_occip(strcmp(per_mus_occip.Groups,'experimental'),:);
data.per.mus.occip.ctl = per_mus_occip(strcmp(per_mus_occip.Groups,'control'),:);
data.per.ret.front.exp = per_ret_front(strcmp(per_ret_front.Groups,'experimental'),:);
data.per.ret.front.ctl = per_ret_front(strcmp(per_ret_front.Groups,'control'),:);
data.per.ret.occip.exp = per_ret_occip(strcmp(per_ret_occip.Groups,'experimental'),:);
data.per.ret.occip.ctl = per_ret_occip(strcmp(per_ret_occip.Groups,'control'),:);

%% Take group average (experimental + control) at each distance
% Extract unique distances
distances = unique(raw_mus.distance);
% Create struct for storing averages
avg = struct();

f = fields(data);
op = {'mus','ret'};
regs = {'front','occip'};
grp = {'exp','ctl'};
% Iterate over datasets
for ii = 1:numel(f)
    % Iterate over optical properties
    for j = 1:numel(op)
        % Iterate regions
        for k = 1:numel(regs)
            % Iterate experimental + control
            for w = 1:numel(grp)
                % Extract table
                tmp = data.(f{ii}).(op{j}).(regs{k}).(grp{w});
                % Take mean across distances
                tmp_avg = arrayfun(@(d) mean(tmp.OpticalProperty( ...
                                            tmp.distance == d)), distances);
                avg.(f{ii}).(op{j}).(regs{k}).(grp{w}) = tmp_avg;
            end
        end
    end
end

%% Scatterplots (Just Percentage Change)

% Title Strings
op_title = {'\mu_s','Retardance'};
reg_title = {'Frontal','Occipital'};

% y-axis limits
ylims = [-8, 10; -15, 5];

% scatterplot and errorbar sizes
dot_size = 100;
lw = 4;
% Font size
fsize = 30;

% Iterate over optical properties
for j = 1:numel(op)
    % Iterate regions
    for k = 1:numel(regs)        
        % Extract experimental and control vectors
        exp = avg.per.(op{j}).(regs{k}).exp;
        ctl = avg.per.(op{j}).(regs{k}).ctl;
        
        % Extract standard deviations
        exp_sd = sd.per.(op{j}).(regs{k}).exp_sd;
        ctl_sd = sd.per.(op{j}).(regs{k}).ctrl_sd;
        
        % Ribbon scatter plot
        figure;
        tstr = sprintf('%s %s', op_title{j}, reg_title{k});
        mean_std_plot(distances, exp, exp_sd, ctl, ctl_sd, ylims(j,:),tstr)
        
        % Figure properties
        set(gca,'fontsize',fsize)
        set(gca,'fontname','Arial');
        set(gca,'XColor','k','YColor','k');
        set(gca,'TickLength',[0.04,0.04]);
        pause(1)

        % Save as .PDF and .FIG
        if save_flag
            saveas(gcf, fullfile(beta_dir, per_dir,...
                                sprintf('percentage_diff_%s_%s_scatterplot',...
                                op{j}, regs{k})), 'pdf');
            saveas(gcf, fullfile(beta_dir, per_dir,...
                        sprintf('percentage_diff_%s_%s_scatterplot',...
                        op{j}, regs{k})), 'fig');
        end

        %% Draw scatterplot w/ standard deviation bars
        %{
        p1 = scatter(distances, exp, dot_size, 'r', 'filled');
        hold on;
        errorbar(distances, exp, exp_sd, 'r', 'LineStyle', 'none','LineWidth',lw);
        p2 = scatter(distances, ctl, dot_size, 'k', 'filled');
        errorbar(distances, ctl, ctl_sd, 'k', 'LineStyle', 'none','LineWidth',lw);
        % axis labels + title
        xlabel('Distance (\mum)');
        ylabel('Percentage Change');
        title(sprintf('%s %s', op_title{j}, reg_title{k}));
        % Print legend on just one figure
        if j==1 && k == 1
            legend([p1,p2],{'Experimental', 'Control'},"Location",'northeast');
        end
        hold off;
        ylim(ylims(j,:));
        % Font size
        set(gca,'fontsize',fsize)
        set(gca,'fontname','Arial');
        set(gca,'XColor','k','YColor','k');
        set(gca,'TickLength',[0.04,0.04]);
        % Save as .PDF and .FIG
        if save_flag
            saveas(gcf, fullfile(beta_dir, per_dir,...
                                sprintf('percentage_diff_%s_%s_scatterplot',...
                                op{j}, regs{k})), 'pdf');
            saveas(gcf, fullfile(beta_dir, per_dir,...
                        sprintf('percentage_diff_%s_%s_scatterplot',...
                        op{j}, regs{k})), 'fig');
        end
        % Print the spearman's rho values to the console
        % Print Spearman's rho values to the console
        [rho_exp,p_exp] = corr(distances, exp, 'Type', 'Spearman');
        [rho_ctl,p_ctl] = corr(distances, ctl, 'Type', 'Spearman');
        fprintf(['Spearman''s rho for %s %s (Experimental): ' ...
            'rho = %.5f, p = %.5f\n'], op{j}, regs{k}, rho_exp, p_exp);
        fprintf(['Spearman''s rho for %s %s (Control): ' ...
            'rho = %.5f, p = %.5f\n'], op{j}, regs{k}, rho_ctl, p_ctl);
        %}
    end
end

%% Scatterplots (all 3 groups)

% Iterate over datasets
for ii = 1:numel(f)
    % Iterate over optical properties
    for j = 1:numel(op)
        % Iterate regions
        for k = 1:numel(regs)
            % Extract experimental and control vectors
            exp = avg.(f{ii}).(op{j}).(regs{k}).exp;
            ctl = avg.(f{ii}).(op{j}).(regs{k}).ctl;
            % Extract standard deviations
            exp_sd = sd.(f{ii}).(op{j}).(regs{k}).exp_sd;
            ctl_sd = sd.(f{ii}).(op{j}).(regs{k}).ctrl_sd;
            % Draw scatterplot
            figure;
            scatter(distances, exp, 'r', 'filled');
            hold on;
            errorbar(distances, exp, exp_sd, 'r', 'LineStyle', 'none');
            scatter(distances, ctl, 'b', 'filled');
            errorbar(distances, ctl, ctl_sd, 'b', 'LineStyle', 'none');
            xlabel('Distance');
            ylabel({op{j}});
            title(sprintf('%s %s %s', f{ii}, op{j}, regs{k}));
            legend('Experimental', 'Control');
            hold off;
        end
    end
end