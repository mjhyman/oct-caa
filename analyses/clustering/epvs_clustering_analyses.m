%% Plot the EPVS clustering + fit exponential decay
%{
The script "epvs_clustering_analyses" uses a threshold to classify EPVS as
seeds if they are above a certain size in total voxels. Then this script
measures the density of the surrounding small EPVS. There is an exponential
decay in the surrounding small EPVS.

This script will import all of the data, plot them, and fit an exponential
decay.
%}
clear; clc; close all;

%% Import the results for each subejct
cluster_dir = '/projectnb/npbssmic/ns/CAA/cluster/';
% Iterate over subjects and regions and import cluster mat
subids = {'caa6','caa17','caa22','caa25','caa26'};
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
        front = load(fullfile(front_dir,'epvs_cluster.mat'));
        clusters.(subids{ii}).front = front.results;
    end

    if isfolder(occip_dir)
        occip = load(fullfile(occip_dir,'epvs_cluster.mat'));
        clusters.(subids{ii}).occip = occip.results;
    end
end

%% Plot subject/region separately

% Iterate subjects + region, plot observed_mean_density
for ii = 1:length(subids)
    % Extract subject
    sub = clusters.(subids{ii});
    regions = fields(sub);
    figure;
    % Iterate regions
    for j = 1:numel(regions)
        % Extract data for front and occipital regions
        dens = clusters.(subids{ii}).(regions{j}).observed_mean_density;
        % Plot observed mean density for front region
        plot(dens, 'DisplayName', regions{j});
        hold on;
    end
    hold off;
    xlabel('Distance (mm)');
    ylabel('EPVS Density (count / mm^3)');
    title(['Subject: ' subids{ii}]);
    legend show;
end

% Fit curve

%% Plot each severity (combined within region)

% Fit curve

