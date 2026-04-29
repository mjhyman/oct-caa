%% Scatter plot of mus + ret vs. distance
% Use the raw means of the data set
% Use Bayesian statistical model for standard deviations
% clear;
clc; close all;

%% Directories and Data
% Input directories
data_dir = '/projectnb/npbssmic/ns/CAA/';
beta_dir = '/projectnb/npbssmic/ns/CAA/beta_stats/';
% Directory to store csv
fig_out = '/projectnb/npbssmic/ns/CAA/figures/fig3_op_vs_dist/';
% Directoies w/ standard deviation tables for each method
per_dir = 'caa_all_radii_percentage_diff_40um_donut_5-Mar-2026';
% Save output
save_flag = true;
% Current date/time for figure outputs
t = datetime("now",'TimeZone','local','Format','dd-MMM-yyyy');

%% Import Averages and Std Dev. from "summary_stats" CSVs
% summary statistics struct
ss = struct();
% posterior distribution struct
pd = struct();
% cells for iterating
stages = {'stage0','stage1','stage2','stage3'};
op = {'mus','ret'};
op_names = {'scattering','retardance'};
regs = {'front','occip'};

% Iterate stages (severity, 0-3)
for ii = 1:numel(stages)
    % Iterate optical property
    for j = 1:2
        % Iterate region
        for k = 1:2           
            % Generate filename and path
            fname = strcat(per_dir,'__',op_names{j},'_',stages{ii},'_',...
                           regs{k},'_summary_stats.csv');
            fullpath = fullfile(beta_dir,per_dir,stages{ii},...
                                strcat(stages{ii},'_',regs{k}),fname);
            % Import summary stats
            try
                tmp = readtable(fullpath);
            catch
                fprintf('\tSkipping %s %s %s\n',stages{ii},regs{k}, op{j})
                continue
            end
            % Extract the standard deviation for control and experimental
            ss.(stages{ii}).(op{j}).(regs{k}).ctl_mu = tmp.ctrl_mean;
            ss.(stages{ii}).(op{j}).(regs{k}).exp_mu = tmp.exp_mean;
            ss.(stages{ii}).(op{j}).(regs{k}).ctl_sd = tmp.ctrl_sd;
            ss.(stages{ii}).(op{j}).(regs{k}).exp_sd = tmp.exp_sd;
        end
    end
end
% Extract distances from last spreadsheet
distances = tmp.distance;

%% Scatterplots with ribbon for std error
% Title Strings
op_title = {'\mu_s','Retardance'};
reg_title = {'Frontal','Occipital'};
% scatterplot and errorbar sizes
dot_size = 100;
lw = 4;
% Font size
fsize = 30;
ylims = [-16, 30; -25, 10];

%%% Find Global Y-Limits for all scatterplots for consistency
% ylims = zeros(numel(op), 2); 
% for j = 1:numel(op)
%     min_val = inf;
%     max_val = -inf;
%     for ii = 1:numel(stages)
%         for k = 1:numel(regs)
%             % Extract data for this iteration
%             c_mu = ss.(stages{ii}).(op{j}).(regs{k}).ctl_mu;
%             e_mu = ss.(stages{ii}).(op{j}).(regs{k}).exp_mu;
%             c_sd = ss.(stages{ii}).(op{j}).(regs{k}).ctl_sd;
%             e_sd = ss.(stages{ii}).(op{j}).(regs{k}).exp_sd;
% 
%             % Calculate the floor and ceiling of the ribbons
%             % Use (:) to ensure we are looking at all elements in the arrays
%             local_min = min([c_mu - c_sd; e_mu - e_sd], [], 'all');
%             local_max = max([c_mu + c_sd; e_mu + e_sd], [], 'all');
% 
%             % Update global trackers
%             if local_min < min_val, min_val = local_min; end
%             if local_max > max_val, max_val = local_max; end
%         end
%     end
% 
%     % Apply a 5% buffer so the data doesn't touch the top/bottom axis
%     range_val = max_val - min_val;
%     ylims(j, :) = [min_val - 0.05*range_val, max_val + 0.05*range_val];
% end
% % Apply floor and ceiling for whole numbers for limits
% ylims(:,1) = floor(ylims(:,1));
% ylims(:,2) = ceil(ylims(:,2));


%%% Scatter / Ribbon for each
% Iterate stages
for ii = 1:numel(stages)
    % Iterate over optical properties
    for j = 1:numel(op)
        % Iterate regions
        for k = 1:numel(regs)                    
            % Extract experimental and control
            try
                ctl_mu = ss.(stages{ii}).(op{j}).(regs{k}).ctl_mu;
                exp_mu = ss.(stages{ii}).(op{j}).(regs{k}).exp_mu;
                ctl_sd = ss.(stages{ii}).(op{j}).(regs{k}).ctl_sd;
                exp_sd = ss.(stages{ii}).(op{j}).(regs{k}).exp_sd;        
            catch
                continue
            end
            
            % Ribbon scatter plot
            figure;
            tstr = sprintf('%s %s %s',stages{ii},op_title{j},reg_title{k});
            mean_std_plot(distances, exp_mu, exp_sd, ctl_mu, ctl_sd,...
                          ylims(j,:),tstr)
            
            % Figure properties
            set(gca,'fontsize',fsize)
            set(get(gca,'Title'),'FontSize',10);
            set(gca,'fontname','Arial');
            set(gca,'XColor','k','YColor','k');
            set(gca,'TickLength',[0.04,0.04]);
            pause(1)
    
            % Save as .PDF .PND and .FIG
            if save_flag
                % Define filename
                fname = fullfile(beta_dir, per_dir,...
                                sprintf('pdiff_%s_%s_%s_scatter_ribbon',...
                                        stages{ii}, op{j}, regs{k}));
                saveas(gcf, fname, 'pdf');
                saveas(gcf, fname, 'png');
                saveas(gcf, fname, 'fig');
            end
        end
    end
end