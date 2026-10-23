%% Plot the EPVS clustering + fit exponential decay
%{
The script "epvs_clustering_analyses" uses a threshold to classify EPVS as
seeds if they are above a certain size in total voxels. Then this script
measures the density of the surrounding small EPVS. There is an exponential
decay in the surrounding small EPVS.

This script will import all of the data, plot them, and fit an exponential
decay.

Grouping: control = CAA26 + CAA6 (averaged), moderate = CAA17,
severe = CAA22 + CAA25 (averaged). Pooled SEM = sqrt(sem1^2 + sem2^2)/2,
i.e., SEM of the mean of two independent estimates.
%}
% clear; clc; close all;

%% Import the results for each subejct
cluster_dir = '/projectnb/npbssmic/ns/CAA/cluster/';
% Iterate over subjects and regions and import cluster mat
subids = {'caa6','caa17','caa22','caa25','caa26'};
% Cluster filename
fname1 = 'epvs_cluster_';
fname2 = '_percentile_90_WatershedRadius_5_SkeletonDS_10';
% Initialize structure for storing data
clusters = struct();
for ii = 1:length(subids)
    % Define subfolder to be subject level
    subfolder = fullfile(cluster_dir, subids{ii});
    
    % Search for region subfolders
    front_dir = fullfile(subfolder,'front');
    occip_dir = fullfile(subfolder,'occip');

    % Load the cluster data for the current subject
    if isfolder(front_dir)
        fname_front = strcat(fname1,subids{ii},'_front',fname2);
        front = load(fullfile(front_dir,fname_front));
        clusters.(subids{ii}).front = front.results;
    end

    if isfolder(occip_dir)
        fname_occip = strcat(fname1,subids{ii},'_occip',fname2);
        occip = load(fullfile(occip_dir,fname_occip));
        clusters.(subids{ii}).occip = occip.results;
    end
end

% Output directory
figdir = '/projectnb/npbssmic/ns/CAA/figures/clustering/';

%% Plot each subject separately (frontal vs occipital)

% Region colours
rcol = struct('front','#DB5829','occip','#1964B0');
rlab = struct('front','Frontal','occip','Occipital');

for ii = 1:length(subids)
    sub = clusters.(subids{ii});
    regions = fields(sub);
    figure;
    hold on;
    for j = 1:numel(regions)
        dens = sub.(regions{j}).observed_mean_density;
        plot(dens, 'DisplayName', rlab.(regions{j}), ...
            'LineWidth', 1.5, 'Color', rcol.(regions{j}));
    end
    hold off;
    xlabel('Distance (mm)');
    ylabel('EPVS Density (count/mm^3)');
    title(['Subject: ' upper(subids{ii})]);
    legend('show', 'Location', 'northeast');
    set(gca,'FontSize',20)
    grid on;
    xlim([1,10]); ylim([0,11]);

    save_fig(gcf, figdir, ['epvs_cluster' fname2 '_' subids{ii}]);
end


%% Plot each severity (combined within region)

% Create struct for storing arrays for iterating
pstruct = struct();

%%% Average the control (caa26 and caa6)
% Average frontal
caa26f_avg = clusters.caa26.front.observed_mean_density;
caa6f_avg  = clusters.caa6.front.observed_mean_density;
caa26f_sem = clusters.caa26.front.observed_sem;
caa6f_sem  = clusters.caa6.front.observed_sem;
pstruct.ctl.front.avg = mean([caa26f_avg;caa6f_avg],1);
pstruct.ctl.front.sem = sqrt(caa26f_sem.^2 + caa6f_sem.^2) / 2;
% Average occipital
caa26o_avg = clusters.caa26.occip.observed_mean_density;
caa6o_avg  = clusters.caa6.occip.observed_mean_density;
caa26o_sem = clusters.caa26.occip.observed_sem;
caa6o_sem  = clusters.caa6.occip.observed_sem;
pstruct.ctl.occip.avg = mean([caa26o_avg;caa6o_avg],1);
pstruct.ctl.occip.sem = sqrt(caa26o_sem.^2 + caa6o_sem.^2) / 2;

% Assign the moderate
pstruct.mod.occip.avg = clusters.caa17.occip.observed_mean_density;
pstruct.mod.occip.sem = clusters.caa17.occip.observed_sem;

%%% Average the severe (Caa22 and caa25)
% Average frontal
caa22f_avg = clusters.caa22.front.observed_mean_density;
caa25f_avg = clusters.caa25.front.observed_mean_density;
caa22f_sem = clusters.caa22.front.observed_sem;
caa25f_sem = clusters.caa25.front.observed_sem;
pstruct.sev.front.avg = mean([caa22f_avg;caa25f_avg],1);
pstruct.sev.front.sem = sqrt(caa22f_sem.^2 + caa25f_sem.^2) / 2;
% Average occipital
caa22o_avg = clusters.caa22.occip.observed_mean_density;
caa25o_avg = clusters.caa25.occip.observed_mean_density;
caa22o_sem = clusters.caa22.occip.observed_sem;
caa25o_sem = clusters.caa25.occip.observed_sem;
pstruct.sev.occip.avg = mean([caa22o_avg;caa25o_avg],1);
pstruct.sev.occip.sem = sqrt(caa22o_sem.^2 + caa25o_sem.^2) / 2;


%% Plot each region with severity groups overlaid

% Group colours and labels (fixed across figures)
gcol = struct('ctl','#555555','mod','#882E72','sev','#DC050C');
glab = struct('ctl','Control','mod','Moderate','sev','Severe');
severity = fields(pstruct);
region_list = {'front','occip'};

for r = 1:numel(region_list)
    reg = region_list{r};
    figure;
    hold on;
    for ii = 1:numel(severity)
        grp = severity{ii};
        % Skip groups without this region (moderate has no frontal)
        if ~isfield(pstruct.(grp), reg)
            continue
        end
        avg = pstruct.(grp).(reg).avg;
        sem = pstruct.(grp).(reg).sem;
        errorbar(1:length(avg), avg, sem, 'DisplayName', glab.(grp), ...
            'LineWidth', 2, 'Color', gcol.(grp));
    end
    hold off;

    %%% Set figure properties
    xlim([1,10]); ylim([0,8]);
    xlabel('Distance (mm)');
    ylabel('EPVS Density (count/mm^3)');
    title(rlab.(reg));
    legend('show', 'Location', 'northeast');
    set(gca,'FontSize',30)

    save_fig(gcf, figdir, ['epvs_cluster' fname2 '_' reg]);
end


%% Local functions
function save_fig(fig, figdir, basename)
% Export figure as 600 dpi PNG and PDF
    drawnow; pause(1)
    exportgraphics(fig, fullfile(figdir, [basename '.png']), 'Resolution', 600)
    pause(1)
    exportgraphics(fig, fullfile(figdir, [basename '.pdf']), 'Resolution', 600)
    pause(1)
end