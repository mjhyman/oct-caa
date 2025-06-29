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
load(fullfile(data_dir,"parenchyma_optical_properties_100um_thick.mat"));
% Scatter plot dot size
scat_size = 100;
% Isotropic voxel size (microns)
vox = 20;
% Filename substring
substr = '_100um_donut';

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

% optical properties list
pnames = {'pmus','pret','pori'};

% retrieve subject IDs
subjects = fieldnames(parench);

% Retrieve radii within struct
radii = fieldnames(parench.(subjects{1}).front);

% Initialize array for pairs of optical proprety vs. distance
ves_mus = zeros(length(radii),1);
ves_ret = zeros(length(radii),1);
ves_ori = zeros(length(radii),1);
epvs_mus = zeros(length(radii),1);
epvs_ret = zeros(length(radii),1);
epvs_ori = zeros(length(radii),1);

% Count the total number of tissue volumes
nvol = 0;
nepvs = 0;
for ii = 1:length(subjects)
    % retrieve subject
    sub = subjects{ii};
    % count number of tissue volumes
    n = length(fieldnames(parench.(sub)));
    % Iterate tissue volume counter
    nvol = nvol + n;
    % counter number of volumes with epvs
    regions = fieldnames(parench.(sub));
    for j=1:length(regions)
        if isfield(parench.(sub).(regions{j}).(radii{1}).outter,'epvs')
            nepvs = nepvs + 1;
        end
    end
end

%%% Iterate distance, subjects, optical property
% iterate distance from edge of epvs/vessel
for ii = 1:length(radii)
    % local radii name
    rname = string(radii(ii));
    % initialize array for adding each subject's optical properties
    ves_mus_tmp = zeros(nvol,1);
    ves_ret_tmp = zeros(nvol,1);
    ves_ori_tmp = zeros(nvol,1);
    epvs_mus_tmp = zeros(nepvs,1);
    epvs_ret_tmp = zeros(nepvs,1);
    epvs_ori_tmp = zeros(nepvs,1);
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
            % Take average of mus, ret, ori
            ves_mus_tmp(vol_idx) = mean(ves.pmus,'omitnan');
            ves_ret_tmp(vol_idx) = mean(ves.pret,'omitnan');
            ves_ori_tmp(vol_idx) = mean(ves.pori,'omitnan');
            % Iterate counter
            vol_idx = vol_idx + 1;
            % retrieve EPVS measurements (if exist)
            if isfield(parench.(sub).(reg).(rname).outter, 'epvs')
                epvs = parench.(sub).(reg).(rname).outter.epvs;
                % Take average of mus, ret, ori
                epvs_mus_tmp(epvs_idx) = mean(epvs.pmus,'omitnan');
                epvs_ret_tmp(epvs_idx) = mean(epvs.pret,'omitnan');
                epvs_ori_tmp(epvs_idx) = mean(epvs.pori,'omitnan');
                % Iterate counter
                epvs_idx = epvs_idx + 1;
            end
        end
    end
    % Add average optical property to main array
    ves_mus(ii) = mean(ves_mus_tmp);
    ves_ret(ii) = mean(ves_ret_tmp);
    ves_ori(ii) = real(mean(ves_ori_tmp));
    epvs_mus(ii) = mean(epvs_mus_tmp);
    epvs_ret(ii) = mean(epvs_ret_tmp);
    epvs_ori(ii) = real(mean(epvs_ori_tmp));
end

%%% Create scatterplot
% vessel - mus
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
scatter(x,ves_mus,scat_size,'filled','k');
xlabel('Distance (\mum)'); ylabel('\mus (cm^-^1)');
title('Vessels \mus vs. Distance'); set(gca,'fontsize',24);
fname = strcat('vessels_mus_vs_distance',substr,'.png');
fout = fullfile(scat_out,'/all_subjects/',fname);
pause(1); saveas(gcf,fout);
% vessel - retardance
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
scatter(x,ves_ret,scat_size,'filled','k');
xlabel('Distance (\mum)'); ylabel('Retardance (degrees)');
title('Vessels Retardance vs. Distance'); set(gca,'fontsize',24);
fname = strcat('vessels_ret_vs_distance',substr,'.png');
fout = fullfile(scat_out,'/all_subjects/',fname);
pause(1); saveas(gcf,fout);
% vessel - orientation
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
scatter(x,ves_ori,scat_size,'filled','k');
xlabel('Distance (\mum)'); ylabel('\sigma_{orientation} (radians)');
title('Vessels \sigma_{orientation} vs. Distance'); set(gca,'fontsize',24);
fname = strcat('vessels_ori_vs_distance',substr,'.png');
fout = fullfile(scat_out,'/all_subjects/',fname);
pause(1); saveas(gcf,fout);
% EPVS - mus
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
scatter(x,epvs_mus,scat_size,'filled','k');
xlabel('Distance (\mum)'); ylabel('\mus (cm^-^1)');
title('EPVS \mus vs. Distance'); set(gca,'fontsize',24);
fname = strcat('epvs_mus_vs_distance',substr,'.png');
fout = fullfile(scat_out,'/all_subjects/',fname);
pause(1); saveas(gcf,fout);
% EPVS - retardance
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
scatter(x,epvs_ret,scat_size,'filled','k');
xlabel('Distance (\mum)'); ylabel('Retardance (degrees)');
title('EPVS Retardance vs. Distance'); set(gca,'fontsize',24);
fname = strcat('epvs_ret_vs_distance',substr,'.png');
fout = fullfile(scat_out,'/all_subjects/',fname);
pause(1); saveas(gcf,fout);
% EPVS - orientation
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
scatter(x,epvs_ori,scat_size,'filled','k');
xlabel('Distance (\mum)'); ylabel('\sigma_{orientation} (radians)');
title('EPVS \sigma_{orientation} vs. Distance'); set(gca,'fontsize',24);
fname = strcat('epvs_ori_vs_distance',substr,'.png');
fout = fullfile(scat_out,'/all_subjects/',fname);
pause(1); saveas(gcf,fout);

%% Box/whisker Plots

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