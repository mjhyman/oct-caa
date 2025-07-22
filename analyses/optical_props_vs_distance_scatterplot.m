%% Analyze the Optical Properties
%{
The purpose of this script is to measure the optical properties in the
parenchymal tissue surrounding the EPVS and the vasculature. This covers
just the following cases:
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)

Outline:
- IMPORT struct containing measurements at various distances
- Create scatterplot (property vs. distance)
- Create box/whisker plot
%}

%% Prepare environment
clc; close all;
% Add top-level directory
current_dir = pwd;
addpath(fullfile(current_dir));
% Directory for loading seg, mus, ret, mask, epvs
data_dir = ['/autofs/cluster/octdata3/users/mjhyman/' ...
    'oct_caa_analyses/optical_properties'];
scat_out = fullfile(data_dir,'/scatter_plots');
bw_out = fullfile(data_dir,'/bw_plots');
% Load parenchyam struct
load(fullfile(data_dir,"parenchyma_optical_properties_40um_thick_26Jun2025.mat"));
% Scatter plot dot size
scat_size = 100;
% Isotropic voxel size (microns)
vox = 20;
% Filename substring
substr = '_40um_donut_26Jun2025';

%% Scatterplot of optical property vs. distance (most & least severe)
% ONLY most severe (CAA22 occip) + least severe (CAA26 occip) EPVS cases
% y-axis = average optical property  for epvs or vessel
% x-axis = distance of ring from edge of epvs or vessel
% Repeat this for all 3 properties (mus, ret, ori)
% This will only use the "outer" ring, which is disjoint from the annotated
% vessel/epvs

%%% Extract number of radii
radii = fieldnames(parench.caa22.occip);
x = zeros(length(radii),1);
% Retrieve inner radii (this will be x-axis)
for ii = 1:length(radii)
    x(ii) = str2num(radii{ii}(4:end));
end
% Convert x-axis to microns
x = x.*vox;

%%% Scatterplots of most & least severe
%{
%%% Most severe subject
% set subject ID and region
sub = 'caa22';
reg = 'occip';
% subfolder for most severe subject
fout = fullfile(scat_out,'/most_severe_epvs/');
% subtitle
sub_tit = 'CAA22 Occipital';
% Call function to measure average for this subject at each radius
[ves_mus, ves_ret, ves_ori, epvs_mus, epvs_ret, epvs_ori] =...
    meas_sub(parench,sub,reg);
% Call scatterplot function
subject_scatter(x,ves_mus, ves_ret, ves_ori,...
                epvs_mus, epvs_ret, epvs_ori,fout,sub_tit)

%%% Least severe subject
% set subject ID and region
sub = 'caa26';
reg = 'occip';
% subfolder for most severe subject
fout = fullfile(scat_out,'/least_severe_epvs/');
% subtitle
sub_tit = 'CAA26 Occipital';
% Call function to measure average for this subject at each radius
[ves_mus, ves_ret, ves_ori, epvs_mus, epvs_ret, epvs_ori] =...
    meas_sub(parench,sub,reg);
% Call scatterplot function
subject_scatter(x,ves_mus, ves_ret, ves_ori,...
                epvs_mus, epvs_ret, epvs_ori,fout,sub_tit)
%}

%% Scatterplot of optical property vs. distance (disjoint ring)
% y-axis = average optical property  for epvs or vessel
% x-axis = distance of ring from edge of epvs or vessel
% repeat this for all 3 properties (mus, ret, ori)
% This will only use the "outer" ring, which is disjoint from the annotated
% vessel/epvs

%%% Initialization
% optical properties list
pnames = {'pmus','pret','pori'};
% retrieve subject IDs
subjects = fieldnames(parench);
% Retrieve radii within struct
radii = fieldnames(parench.(subjects{1}).front);

%%% Combined array for pairs of optical proprety vs. distance
% Mean
comb_mean_ves_mus = zeros(length(radii),1);
comb_mean_ves_ret = zeros(length(radii),1);
comb_mean_ves_ori = zeros(length(radii),1);
comb_mean_epvs_mus = zeros(length(radii),1);
comb_mean_epvs_ret = zeros(length(radii),1);
comb_mean_epvs_ori = zeros(length(radii),1);
% standard error of measurement (SEM)
comb_sem_ves_mus = zeros(length(radii),1);
comb_sem_ves_ret = zeros(length(radii),1);
comb_sem_ves_ori = zeros(length(radii),1);
comb_sem_epvs_mus = zeros(length(radii),1);
comb_sem_epvs_ret = zeros(length(radii),1);
comb_sem_epvs_ori = zeros(length(radii),1);

%%% Combined array for pairs of optical proprety vs. distance
% Mean
front_mean_ves_mus = zeros(length(radii),1);
front_mean_ves_ret = zeros(length(radii),1);
front_mean_ves_ori = zeros(length(radii),1);
front_mean_epvs_mus = zeros(length(radii),1);
front_mean_epvs_ret = zeros(length(radii),1);
front_mean_epvs_ori = zeros(length(radii),1);
% standard error of measurement (SEM)
front_sem_ves_mus = zeros(length(radii),1);
front_sem_ves_ret = zeros(length(radii),1);
front_sem_ves_ori = zeros(length(radii),1);
front_sem_epvs_mus = zeros(length(radii),1);
front_sem_epvs_ret = zeros(length(radii),1);
front_sem_epvs_ori = zeros(length(radii),1);

%%% Combined array for pairs of optical proprety vs. distance
% Mean
occip_mean_ves_mus = zeros(length(radii),1);
occip_mean_ves_ret = zeros(length(radii),1);
occip_mean_ves_ori = zeros(length(radii),1);
occip_mean_epvs_mus = zeros(length(radii),1);
occip_mean_epvs_ret = zeros(length(radii),1);
occip_mean_epvs_ori = zeros(length(radii),1);
% standard error of measurement (SEM)
occip_sem_ves_mus = zeros(length(radii),1);
occip_sem_ves_ret = zeros(length(radii),1);
occip_sem_ves_ori = zeros(length(radii),1);
occip_sem_epvs_mus = zeros(length(radii),1);
occip_sem_epvs_ret = zeros(length(radii),1);
occip_sem_epvs_ori = zeros(length(radii),1);

%%% Iterate distance, subjects, optical property
% iterate distance from edge of epvs/vessel
for ii = 1:length(radii)
    % local radii name
    rname = string(radii(ii));
    %%% initialize arrays for combined
    comb_ves_mus = [];
    comb_ves_ret = [];
    comb_ves_ori = [];
    comb_epvs_mus = [];
    comb_epvs_ret = [];
    comb_epvs_ori = [];
    %%% initialize arrays for frontal
    front_ves_mus = [];
    front_ves_ret = [];
    front_ves_ori = [];
    front_epvs_mus = [];
    front_epvs_ret = [];
    front_epvs_ori = [];
    %%% initialize arrays for occiptial
    occip_ves_mus = [];
    occip_ves_ret = [];
    occip_ves_ori = [];
    occip_epvs_mus = [];
    occip_epvs_ret = [];
    occip_epvs_ori = [];

    %%% iterate subjects
    for j = 1:length(subjects)
        % local subject ID
        sub = subjects{j};
        % retrieve regions for this subject
        regions = fieldnames(parench.(sub));
        % iterate regions
        for k = 1:length(regions)
            % region for this iteration
            reg = regions{k};
            % Retrieve vessel measurements (mus, ret, ori)
            ves = parench.(sub).(reg).(rname).outter.ves;
            
            %%% Add data to combined vectors
            % Take average of mus, ret, ori
            comb_ves_mus = [comb_ves_mus, mean(ves.pmus,'omitnan')];
            comb_ves_ret = [comb_ves_ret, mean(ves.pret,'omitnan')];
            comb_ves_ori = [comb_ves_ori, mean(ves.pori,'omitnan')];
            %%% Add to frontal or occipital
            if strcmp(reg,'front')
                front_ves_mus = [front_ves_mus, mean(ves.pmus,'omitnan')];
                front_ves_ret = [front_ves_ret, mean(ves.pret,'omitnan')];
                front_ves_ori = [front_ves_ori, mean(ves.pori,'omitnan')];
            elseif strcmp(reg,'occip')
                occip_ves_mus = [occip_ves_mus, mean(ves.pmus,'omitnan')];
                occip_ves_ret = [occip_ves_ret, mean(ves.pret,'omitnan')];
                occip_ves_ori = [occip_ves_ori, mean(ves.pori,'omitnan')];
            end

            %%% retrieve EPVS measurements (if exist)
            if isfield(parench.(sub).(reg).(rname).outter, 'epvs')
                % Create local
                epvs = parench.(sub).(reg).(rname).outter.epvs;
                %%% Add to combined
                comb_epvs_mus = [comb_epvs_mus, mean(epvs.pmus,'omitnan')];
                comb_epvs_ret = [comb_epvs_ret, mean(epvs.pret,'omitnan')];
                comb_epvs_ori = [comb_epvs_ori, mean(epvs.pori,'omitnan')];
                %%% Add to Frontl or occipital
                if strcmp(reg,'front')
                    front_epvs_mus = [front_epvs_mus, mean(epvs.pmus,'omitnan')];
                    front_epvs_ret = [front_epvs_ret, mean(epvs.pret,'omitnan')];
                    front_epvs_ori = [front_epvs_ori, mean(epvs.pori,'omitnan')];
                elseif strcmp(reg,'occip')
                    occip_epvs_mus = [occip_epvs_mus, mean(epvs.pmus,'omitnan')];
                    occip_epvs_ret = [occip_epvs_ret, mean(epvs.pret,'omitnan')];
                    occip_epvs_ori = [occip_epvs_ori, mean(epvs.pori,'omitnan')];
                end
            end
        end
    end

    %% Take average across subjects at distance
    %%% Combined - take average across subjects
    % Add average optical property to main array
    comb_mean_ves_mus(ii) = mean(comb_ves_mus);
    comb_mean_ves_ret(ii) = mean(comb_ves_ret);
    comb_mean_ves_ori(ii) = real(mean(comb_ves_ori));
    comb_mean_epvs_mus(ii) = mean(comb_epvs_mus);
    comb_mean_epvs_ret(ii) = mean(comb_epvs_ret);
    comb_mean_epvs_ori(ii) = real(mean(comb_epvs_ori));
    % Measure standard error of mean of array
    comb_sem_ves_mus(ii) = std(comb_ves_mus)./length(comb_ves_mus);
    comb_sem_ves_ret(ii) = std(comb_ves_ret)./length(comb_ves_ret);
    comb_sem_ves_ori(ii) = std(real(comb_ves_ori))./length(comb_ves_ori);
    comb_sem_epvs_mus(ii) = std(comb_epvs_mus)./length(comb_epvs_mus);
    comb_sem_epvs_ret(ii) = std(comb_epvs_ret)./length(comb_epvs_ret);
    comb_sem_epvs_ori(ii) = std(real(comb_epvs_ori))./length(comb_epvs_ret);

    %%% Front - take average across subjects
    % Add average optical property to main array
    front_mean_ves_mus(ii) = mean(front_ves_mus);
    front_mean_ves_ret(ii) = mean(front_ves_ret);
    front_mean_ves_ori(ii) = real(mean(front_ves_ori));
    front_mean_epvs_mus(ii) = mean(front_epvs_mus);
    front_mean_epvs_ret(ii) = mean(front_epvs_ret);
    front_mean_epvs_ori(ii) = real(mean(front_epvs_ori));
    % Measure standard error of mean of array
    front_sem_ves_mus(ii) = std(front_ves_mus)./length(front_ves_mus);
    front_sem_ves_ret(ii) = std(front_ves_ret)./length(front_ves_ret);
    front_sem_ves_ori(ii) = std(real(front_ves_ori))./length(front_ves_ori);
    front_sem_epvs_mus(ii) = std(front_epvs_mus)./length(front_epvs_mus);
    front_sem_epvs_ret(ii) = std(front_epvs_ret)./length(front_epvs_ret);
    front_sem_epvs_ori(ii) = std(real(front_epvs_ori))./length(front_epvs_ret);

    %%% Occip - take average across subjects
    % Add average optical property to main array
    occip_mean_ves_mus(ii) = mean(occip_ves_mus);
    occip_mean_ves_ret(ii) = mean(occip_ves_ret);
    occip_mean_ves_ori(ii) = real(mean(occip_ves_ori));
    occip_mean_epvs_mus(ii) = mean(occip_epvs_mus);
    occip_mean_epvs_ret(ii) = mean(occip_epvs_ret);
    occip_mean_epvs_ori(ii) = real(mean(occip_epvs_ori));
    % Measure standard error of mean of array
    occip_sem_ves_mus(ii) = std(occip_ves_mus)./length(occip_ves_mus);
    occip_sem_ves_ret(ii) = std(occip_ves_ret)./length(occip_ves_ret);
    occip_sem_ves_ori(ii) = std(real(occip_ves_ori))./length(occip_ves_ori);
    occip_sem_epvs_mus(ii) = std(occip_epvs_mus)./length(occip_epvs_mus);
    occip_sem_epvs_ret(ii) = std(occip_epvs_ret)./length(occip_epvs_ret);
    occip_sem_epvs_ori(ii) = std(real(occip_epvs_ori))./length(occip_epvs_ret);

    %% Wilcoxon Signed Rank Test - EPVS vs. Vessel
    % TODO: create a new vessel vector for combined that only keeps the
    % regions containing EPVS

    %%% Combined
    % tail = 'both';
    % [comb_mus] = wilco(comb_epvs_mus,comb_ves_mus,tail);
    % [comb_ret] = wilco(comb_epvs_ret,comb_ves_ret,tail);
    % [comb_ori] = wilco(comb_epvs_ori,comb_ves_ori,tail);
    % pause(0.1)
end

%% Scatterplot limits for each optical property
mus_yl = [9.5, 13];
ret_yl = [24, 30.5];
ori_yl = [0.2, 0.65];

%%% COMBINED scatterplots (EPVS and ves)
% Mus
xlab = 'Distance (\mum)';
ylab = '\mus (cm^-^1)';
tit = 'Combined: \mus vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('COMBINED_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_mus, comb_sem_ves_mus,...
                    comb_mean_epvs_mus, comb_sem_epvs_mus,...
                    xlab, ylab, mus_yl, tit, dir_out, fname)
% Retardance
xlab = 'Distance (\mum)';
ylab = 'Retardance (degrees)';
tit = 'Combined: Retardance vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('COMBINED_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_ret, comb_sem_ves_ret,...
                    comb_mean_epvs_ret, comb_sem_epvs_ret,...
                    xlab, ylab, ret_yl, tit, dir_out, fname)
% Orientation
xlab = 'Distance (\mum)';
ylab = '\sigma_{orientation} (radians)';
tit = 'Combined: CSDO vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('COMBINED_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_ori, comb_sem_ves_ori,...
                    comb_mean_epvs_ori, comb_sem_epvs_ori,...
                    xlab, ylab, ori_yl, tit, dir_out, fname)

%%% FRONTAL scatterplots (EPVS and ves)
% Mus
xlab = 'Distance (\mum)';
ylab = '\mus (cm^-^1)';
tit = 'Frontal: \mus vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('FRONT_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_mus, front_sem_ves_mus,...
                    front_mean_epvs_mus, front_sem_epvs_mus,...
                    xlab, ylab, mus_yl, tit, dir_out, fname)
% Retardance
xlab = 'Distance (\mum)';
ylab = 'Retardance (degrees)';
tit = 'Frontal: Retardance vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('FRONT_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_ret, front_sem_ves_ret,...
                    front_mean_epvs_ret, front_sem_epvs_ret,...
                    xlab, ylab, ret_yl, tit, dir_out, fname)
% Orientation
xlab = 'Distance (\mum)';
ylab = '\sigma_{orientation} (radians)';
tit = 'Frontal: CSDO vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('FRONT_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_ori, front_sem_ves_ori,...
                    front_mean_epvs_ori, front_sem_epvs_ori,...
                    xlab, ylab, ori_yl, tit, dir_out, fname)

%%% OCCIPITAL scatterplots (EPVS and ves)
% Mus
xlab = 'Distance (\mum)';
ylab = '\mus (cm^-^1)';
tit = 'Occipital: \mus vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('OCCIP_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_mus, occip_sem_ves_mus,...
                    occip_mean_epvs_mus, occip_sem_epvs_mus,...
                    xlab, ylab, mus_yl, tit, dir_out, fname)
% Retardance
xlab = 'Distance (\mum)';
ylab = 'Retardance (degrees)';
tit = 'Occipital: Retardance vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('OCCIP_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_ret, occip_sem_ves_ret,...
                    occip_mean_epvs_ret, occip_sem_epvs_ret,...
                    xlab, ylab, ret_yl, tit, dir_out, fname)
% Orientation
xlab = 'Distance (\mum)';
ylab = '\sigma_{orientation} (radians)';
tit = 'Occipital: CSDO vs. Distance';
dir_out = fullfile(scat_out,'/all_subjects/');
fname = strcat('OCCIP_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_ori, occip_sem_ves_ori,...
                    occip_mean_epvs_ori, occip_sem_epvs_ori,...
                    xlab, ylab, ori_yl, tit, dir_out, fname)

%% Box/whisker Plots

%{
% Reorder the subjects list in terms of severity
subjects = {'caa6','caa26','caa17','caa25','caa22'};

%%% Iterate distance, subjects, optical property
% iterate distance from edge of epvs/vessel
for ii = 1:length(radii)
    % local radii name
    rname = string(radii(ii));
    % cell arrays to store the values grouped by subject
    ves_mus = {}; ves_ret = {}; ves_ori = {}; ves_name = {};
    epvs_mus = {}; epvs_ret = {}; epvs_ori = {}; epvs_name = {};
    % cell arrays to store names for each observation
    ves_mus_name = {}; ves_ret_name = {}; ves_ori_name = {};
    epvs_mus_name = {}; epvs_ret_name = {}; epvs_ori_name = {};
    % counter for tracking volume index and epvs index
    vol_idx = 1;
    epvs_idx = 1;
    % iterate subjects
    for j = 1:length(subjects)
        % local subject ID
        sub = subjects{j};
        % retrieve regions for this subject
        regions = fieldnames(parench.(sub));
        % iterate regions
        for k = 1:length(regions)
            % region for this iteration
            reg = regions{k};
            % Retrieve vessel measurements (mus, ret, ori)
            ves = parench.(sub).(reg).(rname).outter.ves;
            
            %%% Retrieve mus, ret, ori
            ves_mus{vol_idx} = rmmissing(ves.pmus);
            ves_ret{vol_idx} = rmmissing(ves.pret);
            ves_ori{vol_idx} = real(rmmissing(ves.pori));
            
            %%% Create x-axis names for each datapoint
            % Create x-axis name for boxplot
            xname = [sub,'_',reg];
            % # observations for each optical property
            nmus = length(ves_mus{vol_idx});
            nret = length(ves_ret{vol_idx});
            nori = length(ves_ori{vol_idx});
            % replicate x-name for each observation
            mus_name = repmat({xname},nmus,1);
            ret_name = repmat({xname},nret,1);
            ori_name = repmat({xname},nori,1);
            % Add vector of names to main cell array
            ves_mus_name{vol_idx} = mus_name;
            ves_ret_name{vol_idx} = ret_name;
            ves_ori_name{vol_idx} = ori_name;
            % Iterate counter
            vol_idx = vol_idx + 1;
            % retrieve EPVS measurements (if exist)
            if isfield(parench.(sub).(reg).(rname).outter, 'epvs')
                epvs = parench.(sub).(reg).(rname).outter.epvs;
                % Create x-axis name for boxplot
                epvs_name{vol_idx} = [sub,'_',reg];
                % Take average of mus, ret, ori
                epvs_mus{epvs_idx} = rmmissing(epvs.pmus);
                epvs_ret{epvs_idx} = rmmissing(epvs.pret);
                epvs_ori{epvs_idx} = real(rmmissing(epvs.pori));
               
                %%% Create x-axis names for each datapoint
                % Create x-axis name for boxplot
                xname = [sub,'_',reg];
                % # observations for each optical property
                nmus = length(epvs_mus{epvs_idx});
                nret = length(epvs_ret{epvs_idx});
                nori = length(epvs_ori{epvs_idx});
                % replicate x-name for each observation
                mus_name = repmat({xname},nmus,1);
                ret_name = repmat({xname},nret,1);
                ori_name = repmat({xname},nori,1);
                % Add vector of names to main cell array
                epvs_mus_name{epvs_idx} = mus_name;
                epvs_ret_name{epvs_idx} = ret_name;
                epvs_ori_name{epvs_idx} = ori_name;
                % Iterate counter
                epvs_idx = epvs_idx + 1;
            end
        end
    end
    
    %%% Box/Whisker Plot (mus)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(ves_mus_name{:});
    % Create y-axis of all concatenated observations
    y = [ves_mus{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,50]); set(bx,'LineWidth',2);
    ylabel('\mus (cm^-^1)');
    title({'Vessels - \mus',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('ves_mus_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;

    %%% Box/Whisker Plot (ret)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(ves_ret_name{:});
    % Create y-axis of all concatenated observations
    y = [ves_ret{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,50]); set(bx,'LineWidth',2);
    ylabel('Retardance (degrees)');
    title({'Vessels - Retardance',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('ves_ret_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;

    %%% Box/Whisker Plot (ori)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(ves_ori_name{:});
    % Create y-axis of all concatenated observations
    y = [ves_ori{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,3]); set(bx,'LineWidth',2);
    ylabel('\sigma_{orientation} (radians)');
    title({'Vessels - Orientation',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('ves_ori_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;
    
    %%% EPVS Box/Whisker Plot (mus)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(epvs_mus_name{:});
    % Create y-axis of all concatenated observations
    y = [epvs_mus{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,30]); set(bx,'LineWidth',2);
    ylabel('\mus (cm^-^1)');
    title({'EPVS - \mus',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('epvs_mus_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;

    %%% EPVS Box/Whisker Plot (ret)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(epvs_ret_name{:});
    % Create y-axis of all concatenated observations
    y = [epvs_ret{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,50]); set(bx,'LineWidth',2);
    ylabel('Retardance (degrees)');
    title({'EPVS - Retardance',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('epvs_ret_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;

    %%% EPVS Box/Whisker Plot (ori)
    % Create group labels for each observation (for unbalanced b/w plots)
    g = vertcat(epvs_ori_name{:});
    % Create y-axis of all concatenated observations
    y = [epvs_ori{:}];
    % create boxplot
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    bx = boxplot(y,g); ylim([0,3]); set(bx,'LineWidth',2);
    ylabel('\sigma_{orientation} (radians)');
    title({'EPVS - Orientation',['Inner radius = ', num2str(x(ii)),'\mum']});
    set(gca,'fontsize',24);
    % Save output
    fname = strcat('epvs_ori_radius_',num2str(x(ii)),substr,'.png');
    fout = fullfile(bw_out,fname);
    pause(1); saveas(gcf,fout); pause(0.5); close all;
end

%}

%% Function to scatterplot
function scatter_op_vs_dist(x, ves_op, ves_sem, epvs_op, epvs_sem,...
    xlab, ylab, yl, tit, dir_out, fname)
% Scatter plot with error bars for optical property vs distance
% INPUTS:
%   - x (vector): x-axis positions
%   - ves_op (vector): optical propery for vessels
%   - ves_sem (vector): optical propery standard error for vessels
%   - epvs_op (vector): optical propery for EPVS
%   - epvs_sem (vector): optical propery standard error for EPVS
%   - xlab (str): x-axis label
%   - ylab (str): y-axis label
%   - yl (vector): y-axis limits
%   - tit (str): figure title
%   - dir_out (str): output directory
%   - fname (str): output filename

    % Init figure
    figure('Position', [100, 100, 900, 900],'Resize', 'off');
    
    % Scatterplot with error bars
    errorbar(x,ves_op,ves_sem,"-s",'Color','k',...
        "MarkerFaceColor",'k','MarkerSize',20,'LineWidth',4);
    hold on;
    errorbar(x,epvs_op,epvs_sem,"-s",'Color','r',...
        "MarkerFaceColor",'r','MarkerSize',20,'LineWidth',4);
    
    % Limits, Labels, and legend
    ylim(yl);
    xlabel(xlab); ylabel(ylab); title(tit);
    % legend({'Vessels','EPVS'});
    
    % Font and fontsize
    fontname("SansSerif")
    set(gca,'fontsize',24);
    
    % Save output
    fout = fullfile(dir_out,fname);
    pause(1); saveas(gcf,fout);
end

%% Function to measure the mean optical property across a single subject
function [ves_mus, ves_ret, ves_ori, epvs_mus, epvs_ret, epvs_ori] =...
    meas_sub(parench,sub,reg)

%%% Initialize arrays for storing optical properties measurements
% Retrieve radii within struct
radii = fieldnames(parench.(sub).(reg));
% Initialize array for pairs of optical proprety vs. distance
ves_mus = zeros(length(radii),1);
ves_ret = zeros(length(radii),1);
ves_ori = zeros(length(radii),1);
epvs_mus = zeros(length(radii),1);
epvs_ret = zeros(length(radii),1);
epvs_ori = zeros(length(radii),1);

%%% Iterate distance & optical property
% iterate distance from edge of epvs/vessel
for ii = 1:length(radii)
    % local radii name
    rname = string(radii(ii));
    
    %%% Vessel measurements
    % Retrieve vessel measurements (mus, ret, ori)
    ves = parench.(sub).(reg).(rname).outter.ves;
    % Take average of mus, ret, ori
    ves_mus(ii) = mean(ves.pmus,'omitnan');
    ves_ret(ii)  = mean(ves.pret,'omitnan');
    ves_ori(ii) = mean(ves.pori,'omitnan');

    %%% EPVS measurements (if exist)
    if isfield(parench.(sub).(reg).(rname).outter, 'epvs')
        epvs = parench.(sub).(reg).(rname).outter.epvs;
        % Take average of mus, ret, ori
        epvs_mus(ii) = mean(epvs.pmus,'omitnan');
        epvs_ret(ii) = mean(epvs.pret,'omitnan');
        epvs_ori(ii) = real(mean(epvs.pori,'omitnan'));
    end
end
end

%% Scatterplot for individual subject
function subject_scatter(x,ves_mus, ves_ret, ves_ori,...
    epvs_mus, epvs_ret, epvs_ori,scat_out,sub_tit)
% SUBJECT_SCATTER scatter plot for an individual subject
% INPUTS:
%   x (cell array): x-axis subject labels
%   ves_mus (vector): vessel scattering coefficient vector
%   ves_ret (vector): vessel retardance vector
%   ves_ori (vector): vessel circular std. dev. orientation
%   epvs_mus (vector): EPVS scattering coefficient vector
%   epvs_ret (vector): EPVS retardance vector
%   epvs_ori (vector): EPVS circular std. dev. orientation
%   scat_out (string): filepath for saving
%   sub_tit (string): subtitle string

% vessel - mus
plot_and_save(x,ves_mus,'\mus (mm^-^1)',sub_tit,...
              'Vessels \mus vs. Distance',scat_out,...
              'vessels_mus_vs_distance.png')
% vessel - retardance
plot_and_save(x,ves_ret,'Retardance (degrees)',sub_tit,...
              'Vessels Retardance vs. Distance',scat_out,...
              'vessels_ret_vs_distance.png')
% vessel - orientation
plot_and_save(x,ves_ori,'\sigma_{orientation} (radians)',sub_tit,...
              'Vessels \sigma_{orientation} vs. Distance',scat_out,...
              'vessels_ori_vs_distance.png')
% EPVS - mus
plot_and_save(x,epvs_mus,'\mus (mm^-^1)',sub_tit,...
              'EPVS \mus vs. Distance',scat_out,...
              'epvs_mus_vs_distance.png')
% EPVS - retardance
plot_and_save(x,epvs_ret,'Retardance (degrees)',sub_tit,...
              'EPVS Retardance vs. Distance',scat_out,...
              'epvs_ret_vs_distance.png')
% EPVS - orientation
plot_and_save(x,epvs_ori,'\sigma_{orientation} (radians)',sub_tit,...
              'EPVS \sigma_{orientation} vs. Distance',scat_out,...
              'epvs_ori_vs_distance.png')

    function plot_and_save(x,y,ylab,tstr1,tstr2,scat_out,fname)
        set(gcf,'Units','Normalized','OuterPosition',[0, 0.04, 1, 0.96]);
        scatter(x,y,100,'filled','k');
        xlabel('Distance (\mum)'); ylabel(ylab);
        title({tstr1,tstr2}); set(gca,'fontsize',24);
        scat_out = fullfile(scat_out,fname);
        pause(1); saveas(gcf,scat_out); pause(0.5);
    end

end