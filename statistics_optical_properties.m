%% Statstical analysis script for Optical Properties
%{
The purpose of this script is to analyze the optical properties in the
parenchymal tissue surrounding the EPVS and the vasculature. This covers
just the following cases:
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)

Outline:
- IMPORT struct containing optical properties (mus + retardance)
- LME model
-   epvs vs. non-epvs

TODO:
- verify labels for retardance measurements
- New statistical tests:
    - How does parenchyma change around EPVS compared to vessels?
        - Use only cases w/ EPVS
        - Vessles: EPVS vs. non-EPVS
    - Does CAA change optical properties globally (not just around EPVS)?
        - Include cases w/ and w/o 
        - Only analyze non-EPVS vessels
        - Compare non-EPVS cases vs. EPVS cases

%}

%% Prepare environment
clear; clc; close all;
% Add top-level directory
current_dir = pwd;
addpath(fullfile(current_dir));
% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/autofs/cluster/octdata2/users/mjhyman/oct_caa_analyses/optical_properties';
% load the parenchymal optical properties
load(fullfile(data_dir,'parenchyma_optical_properties.mat'));

%% Bar Charts
%{
EPVS subjects:
  - CAA17 occip, CAA22 front, CAA22 occip, CAA25 occip, CAA 26 occip
All subjects
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)
%}

%%% Import vessels
% Retrieve vessels (scattering)
ves_pmus_caa6f = parench.caa6.front.ves.pmus;
ves_pmus_caa6o = parench.caa6.occip.ves.pmus;
ves_pmus_17o = parench.caa17.occip.ves.pmus;
ves_pmus_22f = parench.caa22.front.ves.pmus;
ves_pmus_22o = parench.caa22.occip.ves.pmus;
ves_pmus_25f = parench.caa25.front.ves.pmus;
ves_pmus_25o = parench.caa25.occip.ves.pmus;
ves_pmus_26f = parench.caa26.front.ves.pmus;
ves_pmus_26o = parench.caa26.occip.ves.pmus;
ves_mus = vertcat(ves_pmus_caa6f(:),ves_pmus_caa6o(:),ves_pmus_17o(:),...
                ves_pmus_22f(:),ves_pmus_22o(:),...
                ves_pmus_25f(:),ves_pmus_25o(:),...
                ves_pmus_26f(:),ves_pmus_26o(:));
% Retrieve vessels (retardance)
ves_pret_caa6f = parench.caa6.front.ves.pret;
ves_pret_caa6o = parench.caa6.occip.ves.pret;
ves_pret_17o = parench.caa17.occip.ves.pret;
ves_pret_22f = parench.caa22.front.ves.pret;
ves_pret_22o = parench.caa22.occip.ves.pret;
ves_pret_25f = parench.caa25.front.ves.pret;
ves_pret_25o = parench.caa25.occip.ves.pret;
ves_pret_26f = parench.caa26.front.ves.pret;
ves_pret_26o = parench.caa26.occip.ves.pret;
ves_ret = vertcat(ves_pret_caa6f(:),ves_pret_caa6o(:),ves_pret_17o(:),...
                ves_pret_22f(:),ves_pret_22o(:),...
                ves_pret_25f(:),ves_pret_25o(:),...
                ves_pret_26f(:),ves_pret_26o(:));
% Count number of vessel measurements
n_ves = sum([length(ves_pmus_17o), length(ves_pmus_22f),...
    length(ves_pmus_22o), length(ves_pmus_25f),...
    length(ves_pmus_25o), length(ves_pmus_26f), length(ves_pmus_26o),...
    length(ves_pmus_caa6f), length(ves_pmus_caa6o)],'omitnan');

%%% Import EPVS
% Retrieve EPVS (scattering)
epvs_pmus_17 = parench.caa17.occip.epvs.pmus;
epvs_pmus_22f = parench.caa22.front.epvs.pmus;
epvs_pmus_22o = parench.caa22.occip.epvs.pmus;
epvs_pmus_25o = parench.caa25.occip.epvs.pmus;
epvs_pmus_26o = parench.caa26.occip.epvs.pmus;
epvs_mus = vertcat(epvs_pmus_17(:),epvs_pmus_22f(:),epvs_pmus_22o(:),...
                   epvs_pmus_25o(:),epvs_pmus_26o(:));
% Retrieve EPVS (retardance)
epvs_pret_17 = parench.caa17.occip.epvs.pret;
epvs_pret_22f = parench.caa22.front.epvs.pret;
epvs_pret_22o = parench.caa22.occip.epvs.pret;
epvs_pret_25o = parench.caa25.occip.epvs.pret;
epvs_pret_26o = parench.caa26.occip.epvs.pret;
epvs_ret = vertcat(epvs_pret_17(:),epvs_pret_22f(:),epvs_pret_22o(:),...
                   epvs_pret_25o(:),epvs_pret_26o(:));
% Counter number of EPVS measurements
n_epvs = sum([length(epvs_pmus_17), length(epvs_pmus_22f),...
    length(epvs_pmus_22o), length(epvs_pmus_25o),...
    length(epvs_pmus_26o)],'omitnan');

%%% Bar Chart - Scattering
% Subject mean EPVS scattering coefficient
mean_epvs_mus = [NaN,NaN,mean(epvs_pmus_17,'omitnan'),...
    mean(epvs_pmus_22f,'omitnan'),mean(epvs_pmus_22o,'omitnan'),...
    NaN, mean(epvs_pmus_25o,'omitnan'),...
    NaN, mean(epvs_pmus_26o,'omitnan')];
% Subject mean vessel scattering coefficient
mean_ves_mus = [mean(ves_pmus_caa6f,'omitnan'),mean(ves_pmus_caa6o,'omitnan'),...
                mean(ves_pmus_17o,'omitnan'),mean(ves_pmus_22f,'omitnan'),...
                mean(ves_pmus_22o,'omitnan'),mean(ves_pmus_25f,'omitnan'),...
                mean(ves_pmus_25o,'omitnan'),mean(ves_pmus_26f,'omitnan'),...
                mean(ves_pmus_26o,'omitnan')];
% Create barchart x-axis categories
x = categorical({'CAA 6 Front','CAA 6 Occip','CAA 17 Occip',...
    'CAA 22 Front','CAA 22 Occip','CAA 25 Front','CAA 25 Occip',...
    'CAA 26 Front','CAA 26 Occip'});
figure;
bar(x,[mean_ves_mus',mean_epvs_mus'])
title('Scattering Coefficient')
xtickangle(45)
set(gca,'FontSize',25)
legend({'Vessel','EPVS'})
% saveas(gcf,fullfile(data_dir,'mus_bar_chart.png'));

%%% Bar Chart - Retardance
% Subject mean EPVS scattering coefficient
mean_epvs_ret = [NaN,NaN,mean(epvs_pret_17,'omitnan'),...
    mean(epvs_pret_22f,'omitnan'),mean(epvs_pret_22o,'omitnan'),...
    NaN, mean(epvs_pret_25o,'omitnan'),...
    NaN, mean(epvs_pret_26o,'omitnan')];
% Subject mean vessel scattering coefficient
mean_ves_ret = [mean(ves_pret_caa6f,'omitnan'),mean(ves_pret_caa6o,'omitnan'),...
                mean(ves_pret_17o,'omitnan'),mean(ves_pret_22f,'omitnan'),...
                mean(ves_pret_22o,'omitnan'),mean(ves_pret_25f,'omitnan'),...
                mean(ves_pret_25o,'omitnan'),mean(ves_pret_26f,'omitnan'),...
                mean(ves_pret_26o,'omitnan')];
% Barchart
figure;
bar(x,[mean_ves_ret',mean_epvs_ret'])
title('Retardance')
xtickangle(45)
set(gca,'FontSize',25)
legend({'Vessel','EPVS'})
% saveas(gcf,fullfile(data_dir,'ret_bar_chart.png'));

%% T-test
%%% t-test scattering coefficient (EPVS vs. non-EPVS)
[h_mus,p_mus] = ttest2(ves_mus,epvs_mus,'Alpha',0.05);

%%% t-test retardance (EPVS vs. non-EPVS)
[h_ret,p_ret] = ttest2(ves_ret,epvs_ret,'Alpha',0.05);
%}

%% Generate a linear mixed-effects model
% The purpose is to incorporate a random effect into the analysis
% to account for a potential correlation between outcomes on the
% same person.
% TODO: iterate over brain regions (occipital vs. frontal)

% arrays for vessels (non-EPVS) and EPVS
ves_mus = zeros(n_ves,1);
ves_ret = zeros(n_ves,1);
epvs_mus = zeros(n_epvs,1);
epvs_ret = zeros(n_epvs,1);
% measurement array indices
vidx = 1;
eidx = 1;
% Subject ID vectors
ves_subid = zeros(n_ves,1);
epvs_subid = zeros(n_epvs,1);
% Counter for tissue ID
tid = 1;

% Iterate over subjects in parench
subs = fields(parench);
for ii = 1:length(subs)
    % retrieve subject
    sub = subs{ii};
    % retrieve regions for subject
    regions = fields(parench.(sub));
    % iterate over regions
    for j = 1:length(regions)
        % retrieve region
        reg = regions{j};
        % retrieve optical properties (vessels and epvs)
        segmentations = fields(parench.(sub).(reg));
        % Iterate over segmentations
        for k=1:length(segmentations)
            % retrieve segmentation
            seg = segmentations{k};
            op = parench.(sub).(reg).(seg);
            % Retrieve parenchyma scattering and retardance
            pmus = op.pmus;
            pret = op.pret;
            % count number of vessel measurements for subject
            n = length(pmus);
            % Add to vessel or epvs
            if strcmp(seg,'ves')
                % Add vessel measurements to array
                ves_mus(vidx:vidx+n-1) = pmus;
                ves_ret(vidx:vidx+n-1) = pret;
                % Update subject ID array
                ves_subid(vidx:vidx+n-1) = ones(n,1) .* tid;
                % Iterate counter
                vidx = vidx + n;
            elseif strcmp(seg,'epvs')
                % Add vessel measurements to array
                epvs_mus(eidx:eidx+n-1) = pmus;
                epvs_ret(eidx:eidx+n-1) = pret;
                % Update subject ID array
                epvs_subid(eidx:eidx+n-1) = ones(n,1) .* tid;
                % Iterate counter
                eidx = eidx + n;
            end
        end
        % Iterate tissue sample counter
        tid = tid + 1;
    end
end

%%% Create a table for fitting the LME model
% Column 1 = group label
% Column 2 = vascular metric value
% Create the group labels
g_epvs = repmat({'epvs'},[length(epvs_mus),1]);
g_ves = repmat({'vessel'},[length(ves_mus),1]);
g_epvs_ves = vertcat(g_epvs, g_ves);
% Scattering coefficient array
mus_epvs_ves = [epvs_mus; ves_mus];
% Retardance array
ret_epvs_ves = [epvs_ret; ves_ret];
% Combine the subject IDs into column vector
subids_epvs_ves = vertcat(epvs_subid, ves_subid);
% Table (group labels, subjectID, vascular metric values)
tbl_mus = table(g_epvs_ves,subids_epvs_ves,mus_epvs_ves,...
      'VariableNames',{'Groups','subID','OpticalProperty'});
tbl_ret = table(g_epvs_ves,subids_epvs_ves,ret_epvs_ves,...
      'VariableNames',{'Groups','subID','OpticalProperty'});
       
%%% Fit linear mixed-effects model (LME) (table)
% Define the model:
%   response = vascular metric array
%   random effect (intercept/subject): subID
%   fixed effect: Groups
fml = 'OpticalProperty ~ Groups + (1 | subID)';
% Fit the model for AD vs. HC
lme_mus = fitlme(tbl_mus,fml);
% Fit the model for CTE vs. HC
lme_ret = fitlme(tbl_ret,fml);

%%% Estimates of fixed effects 
[~,~,stats_mus] = fixedEffects(lme_mus);
[~,~,stats_ret] = fixedEffects(lme_ret);
% Store the p-value for the stats tests
p = zeros(2,1);
p(1) = stats_mus{2,6};
p(2) = stats_ret{2,6};

%% Histograms
% Scattering
figure; hold on;
histogram(ves_mus,'BinWidth',1);
histogram(epvs_mus,'BinWidth',1);
xlim([0,40]);
% ylim([0,200]);
legend({'Vessel','EPVS'})
title('\mu_s'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('cm^-^1')

% Retardance
figure; hold on;
histogram(ves_ret,'BinWidth',1);
histogram(epvs_ret,'BinWidth',1);
xlim([5,50]);
% ylim([0,200]);
legend({'Vessel','EPVS'})
title('Retardance'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Degrees')
