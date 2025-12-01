%% Test linear mixed effects model
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
- Linearity assumptions for LMM
- LME model
%}

%% Prepare environment
clc; close all;
addpath(genpath('/projectnb/npbssmic/s/mhyman/oct-caa'));

% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/projectnb/npbssmic/ns/CAA/';
fig_out = '/projectnb/npbssmic/ns/CAA/figures/';

%%% 40 um donut
% load the parenchymal optical properties
load(fullfile(data_dir, ...
    'parenchyma_optical_properties_40um_thick_03Nov2025.mat'));
% String indicating EPVS ring radius (in voxels) to access structure
rad = 'rad2';

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
ves_pmus_caa6f = rmmissing(parench.caa6.front.(rad).outter.ves.pmus);
ves_pmus_caa6o = rmmissing(parench.caa6.occip.(rad).outter.ves.pmus);
ves_pmus_17o = rmmissing(parench.caa17.occip.(rad).outter.ves.pmus);
ves_pmus_22f = rmmissing(parench.caa22.front.(rad).outter.ves.pmus);
ves_pmus_22o = rmmissing(parench.caa22.occip.(rad).outter.ves.pmus);
ves_pmus_25f = rmmissing(parench.caa25.front.(rad).outter.ves.pmus);
ves_pmus_25o = rmmissing(parench.caa25.occip.(rad).outter.ves.pmus);
ves_pmus_26f = rmmissing(parench.caa26.front.(rad).outter.ves.pmus);
ves_pmus_26o = rmmissing(parench.caa26.occip.(rad).outter.ves.pmus);
ves_mus = vertcat(ves_pmus_caa6f(:),ves_pmus_caa6o(:),ves_pmus_17o(:),...
                ves_pmus_22f(:),ves_pmus_22o(:),...
                ves_pmus_25f(:),ves_pmus_25o(:),...
                ves_pmus_26f(:),ves_pmus_26o(:));
% Retrieve vessels (retardance)
ves_pret_caa6f = rmmissing(parench.caa6.front.(rad).outter.ves.pret);
ves_pret_caa6o = rmmissing(parench.caa6.occip.(rad).outter.ves.pret);
ves_pret_17o = rmmissing(parench.caa17.occip.(rad).outter.ves.pret);
ves_pret_22f = rmmissing(parench.caa22.front.(rad).outter.ves.pret);
ves_pret_22o = rmmissing(parench.caa22.occip.(rad).outter.ves.pret);
ves_pret_25f = rmmissing(parench.caa25.front.(rad).outter.ves.pret);
ves_pret_25o = rmmissing(parench.caa25.occip.(rad).outter.ves.pret);
ves_pret_26f = rmmissing(parench.caa26.front.(rad).outter.ves.pret);
ves_pret_26o = rmmissing(parench.caa26.occip.(rad).outter.ves.pret);
ves_ret = vertcat(ves_pret_caa6f(:),ves_pret_caa6o(:),ves_pret_17o(:),...
                ves_pret_22f(:),ves_pret_22o(:),...
                ves_pret_25f(:),ves_pret_25o(:),...
                ves_pret_26f(:),ves_pret_26o(:));
% Retrieve vessel circular mean of standard deviation 
ves_pori_std_caa6f = rmmissing(parench.caa6.front.(rad).outter.ves.pori);
ves_pori_std_caa6o = rmmissing(parench.caa6.occip.(rad).outter.ves.pori);
ves_pori_std_17o = rmmissing(parench.caa17.occip.(rad).outter.ves.pori);
ves_pori_std_22f = rmmissing(parench.caa22.front.(rad).outter.ves.pori);
ves_pori_std_22o = rmmissing(parench.caa22.occip.(rad).outter.ves.pori);
ves_pori_std_25f = rmmissing(parench.caa25.front.(rad).outter.ves.pori);
ves_pori_std_25o = rmmissing(parench.caa25.occip.(rad).outter.ves.pori);
ves_pori_std_26f = rmmissing(parench.caa26.front.(rad).outter.ves.pori);
ves_pori_std_26o = rmmissing(parench.caa26.occip.(rad).outter.ves.pori);
ves_pori_std = vertcat(ves_pori_std_caa6f(:),ves_pori_std_caa6o(:),...
                ves_pori_std_17o(:),...
                ves_pori_std_22f(:),ves_pori_std_22o(:),...
                ves_pori_std_25f(:),ves_pori_std_25o(:),...
                ves_pori_std_26f(:),ves_pori_std_26o(:));

%%% Import EPVS
% Retrieve EPVS (scattering)
epvs_pmus_17 = rmmissing(parench.caa17.occip.(rad).outter.epvs.pmus);
epvs_pmus_22f = rmmissing(parench.caa22.front.(rad).outter.epvs.pmus);
epvs_pmus_22o = rmmissing(parench.caa22.occip.(rad).outter.epvs.pmus);
epvs_pmus_25o = rmmissing(parench.caa25.occip.(rad).outter.epvs.pmus);
epvs_pmus_26o = rmmissing(parench.caa26.occip.(rad).outter.epvs.pmus);
epvs_mus = vertcat(epvs_pmus_17(:),epvs_pmus_22f(:),epvs_pmus_22o(:),...
                   epvs_pmus_25o(:),epvs_pmus_26o(:));
% Retrieve EPVS (retardance)
epvs_pret_17 = rmmissing(parench.caa17.occip.(rad).outter.epvs.pret);
epvs_pret_22f = rmmissing(parench.caa22.front.(rad).outter.epvs.pret);
epvs_pret_22o = rmmissing(parench.caa22.occip.(rad).outter.epvs.pret);
epvs_pret_25o = rmmissing(parench.caa25.occip.(rad).outter.epvs.pret);
epvs_pret_26o = rmmissing(parench.caa26.occip.(rad).outter.epvs.pret);
epvs_ret = vertcat(epvs_pret_17(:),epvs_pret_22f(:),epvs_pret_22o(:),...
                   epvs_pret_25o(:),epvs_pret_26o(:));
% Retrieve EPVS circular standard deviation
epvs_pori_std_17 = rmmissing(parench.caa17.occip.(rad).outter.epvs.pret);
epvs_pori_std_22f = rmmissing(parench.caa22.front.(rad).outter.epvs.pret);
epvs_pori_std_22o = rmmissing(parench.caa22.occip.(rad).outter.epvs.pret);
epvs_pori_std_25o = rmmissing(parench.caa25.occip.(rad).outter.epvs.pret);
epvs_pori_std_26o = rmmissing(parench.caa26.occip.(rad).outter.epvs.pret);
epvs_pori_std = vertcat(epvs_pori_std_17(:),epvs_pori_std_22f(:),...
                        epvs_pori_std_22o(:),...
                        epvs_pori_std_25o(:),epvs_pori_std_26o(:));

%% Bar Charts - Scattering + Retardance of all subjects
%{
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
%}

%% Estimates mean + std. dev. for mu_intercept and beta_condition

%%% mu_intercept = average optical property around vessel (control)
% scattering
mus_mean = mean(ves_mus(:));
mus_std = std(ves_mus(:));
% Retardance
ret_mean = mean(ves_ret(:));
ret_std = std(ves_ret(:));
% orientation
ori_mean = mean(ves_pori_std(:));
ori_std = std(ves_pori_std(:));

%%% beta_condition = estimated diff. between control vs. exp
% beta_cond < 0 <--> experimental > control
% beta_cond > 0 <--> experimental < control
% scattering
mus_beta_mean = mean(ves_mus) - mean(epvs_mus);
mus_beta_std = sqrt( ...
                    var(ves_mus, 1) / length(ves_mus) + ...
                    var(epvs_mus, 1) / length(epvs_mus));
% ret
ret_beta_mean = mean(ves_ret) - mean(epvs_ret);
ret_beta_std = sqrt( ...
                    var(ves_ret, 1) / length(ves_ret) + ...
                    var(epvs_ret, 1) / length(epvs_ret));
% ori
ori_beta_mean = mean(ves_pori_std) - mean(epvs_pori_std);
ori_beta_std = sqrt( ...
                    var(ves_pori_std, 1) / length(ves_pori_std) + ...
                    var(epvs_pori_std, 1) / length(epvs_pori_std));