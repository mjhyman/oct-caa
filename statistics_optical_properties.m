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

%}

%% Prepare environment
% clear;
clc; close all;
% Add top-level directory
current_dir = pwd;
addpath(fullfile(current_dir));
% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/optical_properties';
fig_out = '/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/figures/statistics';
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
ves_pmus_caa6f = rmmissing(parench.caa6.front.ves.pmus);
ves_pmus_caa6o = rmmissing(parench.caa6.occip.ves.pmus);
ves_pmus_17o = rmmissing(parench.caa17.occip.ves.pmus);
ves_pmus_22f = rmmissing(parench.caa22.front.ves.pmus);
ves_pmus_22o = rmmissing(parench.caa22.occip.ves.pmus);
ves_pmus_25f = rmmissing(parench.caa25.front.ves.pmus);
ves_pmus_25o = rmmissing(parench.caa25.occip.ves.pmus);
ves_pmus_26f = rmmissing(parench.caa26.front.ves.pmus);
ves_pmus_26o = rmmissing(parench.caa26.occip.ves.pmus);
ves_mus = vertcat(ves_pmus_caa6f(:),ves_pmus_caa6o(:),ves_pmus_17o(:),...
                ves_pmus_22f(:),ves_pmus_22o(:),...
                ves_pmus_25f(:),ves_pmus_25o(:),...
                ves_pmus_26f(:),ves_pmus_26o(:));
% Retrieve vessels (retardance)
ves_pret_caa6f = rmmissing(parench.caa6.front.ves.pret);
ves_pret_caa6o = rmmissing(parench.caa6.occip.ves.pret);
ves_pret_17o = rmmissing(parench.caa17.occip.ves.pret);
ves_pret_22f = rmmissing(parench.caa22.front.ves.pret);
ves_pret_22o = rmmissing(parench.caa22.occip.ves.pret);
ves_pret_25f = rmmissing(parench.caa25.front.ves.pret);
ves_pret_25o = rmmissing(parench.caa25.occip.ves.pret);
ves_pret_26f = rmmissing(parench.caa26.front.ves.pret);
ves_pret_26o = rmmissing(parench.caa26.occip.ves.pret);
ves_ret = vertcat(ves_pret_caa6f(:),ves_pret_caa6o(:),ves_pret_17o(:),...
                ves_pret_22f(:),ves_pret_22o(:),...
                ves_pret_25f(:),ves_pret_25o(:),...
                ves_pret_26f(:),ves_pret_26o(:));
% Retrieve vessel circular mean of orientation 
ves_pori_mean_caa6f = rmmissing(parench.caa6.front.ves.pori_mean);
ves_pori_mean_caa6o = rmmissing(parench.caa6.occip.ves.pori_mean);
ves_pori_mean_17o = rmmissing(parench.caa17.occip.ves.pori_mean);
ves_pori_mean_22f = rmmissing(parench.caa22.front.ves.pori_mean);
ves_pori_mean_22o = rmmissing(parench.caa22.occip.ves.pori_mean);
ves_pori_mean_25f = rmmissing(parench.caa25.front.ves.pori_mean);
ves_pori_mean_25o = rmmissing(parench.caa25.occip.ves.pori_mean);
ves_pori_mean_26f = rmmissing(parench.caa26.front.ves.pori_mean);
ves_pori_mean_26o = rmmissing(parench.caa26.occip.ves.pori_mean);
ves_pori_mean = vertcat(ves_pori_mean_caa6f(:),ves_pori_mean_caa6o(:),...
                ves_pori_mean_17o(:),...
                ves_pori_mean_22f(:),ves_pori_mean_22o(:),...
                ves_pori_mean_25f(:),ves_pori_mean_25o(:),...
                ves_pori_mean_26f(:),ves_pori_mean_26o(:));
% Retrieve vessel circular mean of standard deviation 
ves_pori_std_caa6f = rmmissing(parench.caa6.front.ves.pori_std);
ves_pori_std_caa6o = rmmissing(parench.caa6.occip.ves.pori_std);
ves_pori_std_17o = rmmissing(parench.caa17.occip.ves.pori_std);
ves_pori_std_22f = rmmissing(parench.caa22.front.ves.pori_std);
ves_pori_std_22o = rmmissing(parench.caa22.occip.ves.pori_std);
ves_pori_std_25f = rmmissing(parench.caa25.front.ves.pori_std);
ves_pori_std_25o = rmmissing(parench.caa25.occip.ves.pori_std);
ves_pori_std_26f = rmmissing(parench.caa26.front.ves.pori_std);
ves_pori_std_26o = rmmissing(parench.caa26.occip.ves.pori_std);
ves_pori_std = vertcat(ves_pori_std_caa6f(:),ves_pori_std_caa6o(:),...
                ves_pori_std_17o(:),...
                ves_pori_std_22f(:),ves_pori_std_22o(:),...
                ves_pori_std_25f(:),ves_pori_std_25o(:),...
                ves_pori_std_26f(:),ves_pori_std_26o(:));

%%% Import EPVS
% Retrieve EPVS (scattering)
epvs_pmus_17 = rmmissing(parench.caa17.occip.epvs.pmus);
epvs_pmus_22f = rmmissing(parench.caa22.front.epvs.pmus);
epvs_pmus_22o = rmmissing(parench.caa22.occip.epvs.pmus);
epvs_pmus_25o = rmmissing(parench.caa25.occip.epvs.pmus);
epvs_pmus_26o = rmmissing(parench.caa26.occip.epvs.pmus);
epvs_mus = vertcat(epvs_pmus_17(:),epvs_pmus_22f(:),epvs_pmus_22o(:),...
                   epvs_pmus_25o(:),epvs_pmus_26o(:));
% Retrieve EPVS (retardance)
epvs_pret_17 = rmmissing(parench.caa17.occip.epvs.pret);
epvs_pret_22f = rmmissing(parench.caa22.front.epvs.pret);
epvs_pret_22o = rmmissing(parench.caa22.occip.epvs.pret);
epvs_pret_25o = rmmissing(parench.caa25.occip.epvs.pret);
epvs_pret_26o = rmmissing(parench.caa26.occip.epvs.pret);
epvs_ret = vertcat(epvs_pret_17(:),epvs_pret_22f(:),epvs_pret_22o(:),...
                   epvs_pret_25o(:),epvs_pret_26o(:));
% Retrieve EPVS circular mean
epvs_pori_mean_17 = rmmissing(parench.caa17.occip.epvs.pori_mean);
epvs_pori_mean_22f = rmmissing(parench.caa22.front.epvs.pori_mean);
epvs_pori_mean_22o = rmmissing(parench.caa22.occip.epvs.pori_mean);
epvs_pori_mean_25o = rmmissing(parench.caa25.occip.epvs.pori_mean);
epvs_pori_mean_26o = rmmissing(parench.caa26.occip.epvs.pori_mean);
epvs_pori_mean = vertcat(epvs_pori_mean_17(:),epvs_pori_mean_22f(:),...
                        epvs_pori_mean_22o(:),...
                        epvs_pori_mean_25o(:),epvs_pori_mean_26o(:));
% Retrieve EPVS circular standard deviation
epvs_pori_std_17 = rmmissing(parench.caa17.occip.epvs.pori_std);
epvs_pori_std_22f = rmmissing(parench.caa22.front.epvs.pori_std);
epvs_pori_std_22o = rmmissing(parench.caa22.occip.epvs.pori_std);
epvs_pori_std_25o = rmmissing(parench.caa25.occip.epvs.pori_std);
epvs_pori_std_26o = rmmissing(parench.caa26.occip.epvs.pori_std);
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

%% Linear Mixed Effects Models
% The purpose is to incorporate a random effect into the analysis
% to account for a potential correlation between outcomes on the
% same person.
%{
- test 1: all subjects (done)
    - vessels w/ EPVS vs. vessels w/o EPVS 
- test 2: only subjects w/ EPVS (TBD)
    - EPVS vs. non-EPVS vessels
- test 3: all subjects (TBD)
    - exclude EPVS parenchyma. only analyze vessels
    - EPVS cases vs. non-EPVS cases 
    - vessels vs. vessels
%}

%%% Test 1
% Count number of vessel measurements across all volumes
n_ves = sum([length(ves_pmus_17o), length(ves_pmus_22f),...
    length(ves_pmus_22o), length(ves_pmus_25f),...
    length(ves_pmus_25o), length(ves_pmus_26f), length(ves_pmus_26o),...
    length(ves_pmus_caa6f), length(ves_pmus_caa6o)],'omitnan');
% Count number of EPVS measurements
n_epvs = sum([length(epvs_pmus_17), length(epvs_pmus_22f),...
    length(epvs_pmus_22o), length(epvs_pmus_25o),...
    length(epvs_pmus_26o)],'omitnan');
% Test index
test_idx = 1;
% Perform LMM Test 1 (same subID for exp/cont within volume)
[stats1, p1, exp_mus1, exp_ret1, exp_sori1,...
    ctl_mus1, ctl_ret1, ctl_sori1,...
    tbl_mus, tbl_ret, tbl_sori] =...
        lmm_test(test_idx, parench, n_ves, n_epvs);
% Export the tables to CSV
fout = fullfile(data_dir,'lmm_test1_same_id.xlsx');
writetable(tbl_mus,fout,'Sheet',1);
writetable(tbl_ret,fout,'Sheet',2);
writetable(tbl_sori,fout,'Sheet',3);
%%% Perform LMM Test 1 (same subID for exp/cont within volume)
[stats1_diff_subid, ~, ~, ~, ~, ~,~, ~,...
    tbl_mus_diff, tbl_ret_diff, tbl_sori_diff] =...
        lmm_test_diff_subid(test_idx, parench, n_ves, n_epvs);
% Export the tables to CSV
fout = fullfile(data_dir,'lmm_test1_diff_id.xlsx');
writetable(tbl_mus_diff,fout,'Sheet',1);
writetable(tbl_ret_diff,fout,'Sheet',2);
writetable(tbl_sori_diff,fout,'Sheet',3);


%% Test 2 only subjects w/ EPVS
% CAA 17 occip, CAA 22 frontal, CAA 22 occip, CAA 25 occip CAA 26 occip
% # vessel measurements across subjects w/ EPVS
n_ctl = sum([length(ves_pmus_17o), length(ves_pmus_22f),...
                length(ves_pmus_22o), length(ves_pmus_25o),...
                length(ves_pmus_26o)],'omitnan');
% # EPVS measurements across subjects w/ EPVS
n_exp = sum([length(epvs_pmus_17), length(epvs_pmus_22f),...
                length(epvs_pmus_22o), length(epvs_pmus_25o),...
                length(epvs_pmus_26o)],'omitnan');
% Test index
test_idx = 2;
% Perform LMM Test 1 (same subID for exp/cont within volume)
[stats2, p2, exp_mus2, exp_ret2, exp_sori2,...
    ctl_mus2, ctl_ret2, ctl_sori2,...
    tbl_mus, tbl_ret, tbl_sori] =...
        lmm_test(test_idx, parench, n_ctl, n_exp);
% Export the tables to CSV
fout = fullfile(data_dir,'lmm_test2_same_id.xlsx');
writetable(tbl_mus,fout,'Sheet',1);
writetable(tbl_ret,fout,'Sheet',2);
writetable(tbl_sori,fout,'Sheet',3);
%%% Perform LMM Test 2 (same subID for exp/cont within volume)
[stats2_diff_subid, ~, ~, ~, ~, ~,~, ~,...
    tbl_mus_diff, tbl_ret_diff, tbl_sori_diff] =...
        lmm_test_diff_subid(test_idx, parench, n_ctl, n_exp);
% Export the tables to CSV
fout = fullfile(data_dir,'lmm_test2_diff_id.xlsx');
writetable(tbl_mus_diff,fout,'Sheet',1);
writetable(tbl_ret_diff,fout,'Sheet',2);
writetable(tbl_sori_diff,fout,'Sheet',3);


%%% Test 3 (All subjects, only vessels)
% control: CAA6 front, CAA25 front, CAA26 front, CAA 6 occip
% exp: CAA17 occip, CAA22 front, CAA22 occip, CAA25 occip, CAA26 occip
% Count number of vessel measurements across controls (no EPVS cases)
n_ctl = sum([length(ves_pmus_caa6f), length(ves_pmus_25f),...
                length(ves_pmus_26f), length(ves_pmus_caa6o)],'omitnan');
% Count number of vessel measurements across exp (cases w/ EPVS)
n_exp = sum([length(ves_pmus_17o), length(ves_pmus_22f),...
                length(ves_pmus_22o), length(ves_pmus_25o),...
                length(ves_pmus_26o)],'omitnan');
% Test index
test_idx = 3;
% Perform LMM Test 3 (same subID for exp/cont within volume)
[stats3, p3, exp_mus3, exp_ret3, exp_sori3,...
        ctl_mus3, ctl_ret3, ctl_sori3,...
        tbl_mus, tbl_ret, tbl_sori] =...
    lmm_test(test_idx, parench, n_ctl, n_exp);
% Export the tables to CSV
fout = fullfile(data_dir,'lmm_test3_same_id.xlsx');
writetable(tbl_mus,fout,'Sheet',1);
writetable(tbl_ret,fout,'Sheet',2);
writetable(tbl_sori,fout,'Sheet',3);
%%% Perform LMM Test 3 (same subID for exp/cont within volume)
[stats3_diff_subid, ~, ~, ~, ~, ~,~, ~,...
    tbl_mus_diff, tbl_ret_diff, tbl_sori_diff] =...
        lmm_test_diff_subid(test_idx, parench, n_ctl, n_exp);
% Export the tables to CSV
fout = fullfile(data_dir,'lmm_test3_diff_id.xlsx');
writetable(tbl_mus_diff,fout,'Sheet',1);
writetable(tbl_ret_diff,fout,'Sheet',2);
writetable(tbl_sori_diff,fout,'Sheet',3);

%% Histograms - Optical Properties of all subjects

%%% Test 1
% Scattering
figure; hold on;
histogram(ctl_mus1,'BinWidth',1);
histogram(exp_mus1,'BinWidth',1);
xlim([0,40]); % ylim([0,200]);
legend({'Vessel','EPVS'})
title('Test 1 - \mu_s'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('cm^-^1')
% Retardance
figure; hold on;
histogram(ctl_ret1,'BinWidth',1);
histogram(exp_ret1,'BinWidth',1);
xlim([5,50]); % ylim([0,200]);
legend({'Vessel','EPVS'})
title('Test 1 - Retardance'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Degrees')
% Circular mean of orientation
figure; hold on;
histogram(real(ctl_mori1));
histogram(real(exp_mori1));
legend({'Vessel','EPVS'})
title('Test 1 - Mean Orientation'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Radians')
% Circular variance of orientation
figure; hold on;
histogram(real(ctl_sori1));
histogram(real(exp_sori1));
legend({'Vessel','EPVS'})
title('Test 1 - Std. Dev. Orientation'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Radians')

%%% Test 2
% Scattering
figure; hold on;
histogram(ctl_mus2,'BinWidth',1);
histogram(exp_mus2,'BinWidth',1);
xlim([0,40]); % ylim([0,200]);
legend({'Vessel','EPVS'})
title('Test 2 - \mu_s'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('cm^-^1')
% Retardance
figure; hold on;
histogram(ctl_ret2,'BinWidth',1);
histogram(exp_ret2,'BinWidth',1);
xlim([5,50]); % ylim([0,200]);
legend({'Vessel','EPVS'})
title('Test 2 - Retardance'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Degrees')
% Circular mean of orientation
figure; hold on;
histogram(real(ctl_mori2));
histogram(real(exp_mori2));
legend({'Vessel','EPVS'})
title('Test 2 - Mean Orientation'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Radians')
% Circular variance of orientation
figure; hold on;
histogram(real(ctl_sori2));
histogram(real(exp_sori2));
legend({'Vessel','EPVS'})
title('Test 2 - Std. Dev. Orientation'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Radians')

%%% Test 3
% Scattering
figure; hold on;
histogram(ctl_mus3,'BinWidth',1);
histogram(exp_mus3,'BinWidth',1);
xlim([0,40]); % ylim([0,200]);
legend({'Vessel','EPVS'})
title('Test 3 - \mu_s'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('cm^-^1')
% Retardance
figure; hold on;
histogram(ctl_ret3,'BinWidth',1);
histogram(exp_ret3,'BinWidth',1);
xlim([5,50]); % ylim([0,200]);
legend({'Vessel','EPVS'})
title('Test 3 - Retardance'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Degrees')
% Circular mean of orientation
figure; hold on;
histogram(real(ctl_mori3));
histogram(real(exp_mori3));
legend({'Vessel','EPVS'})
title('Test 3 - Mean Orientation'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Radians')
% Circular variance of orientation
figure; hold on;
histogram(real(ctl_sori3));
histogram(real(exp_sori3));
legend({'Vessel','EPVS'})
title('Test 3 - Std. Dev. Orientation'); set(gca,'fontsize',30)
ylabel('Counts'); xlabel('Radians')

%% Box/Whisker: Optical Properties vs. N EPVS

%%% Compute N EPVS for each region
% Create a structure array with the desired fields
top = struct();
top.caa6 = caa6;
top.caa17 = caa17;
top.caa22 = caa22;
top.caa25 = caa25;
top.caa26 = caa26;

% List of structure names
sub_ids = fieldnames(top);

% Struct to store the number of EPVS
n_epvs = struct();

% Structuring element for dilating epvs
se = strel('sphere',5);

% Iterate over each structure in the list
for i = 1:length(sub_ids)
	% Retrieve top-level names of struct
    sub = sub_ids{i};
    tmp = top.(sub);
    
	% Retrieve regions of this sample
	regions = fieldnames(top.(sub));

    % Access the regions directly by field names
    for ii = 1:length(regions)
		% Check if it contains EPVS
        if isfield(tmp.(regions{ii}),'epvs')
            % Extract EPVS
            epvs = tmp.(regions{ii}).epvs;
            % Dilate EPVS
            epvs = imdilate(epvs,se);
            % Label connected components
            [~, n] = bwlabeln(epvs);
		    % Add # EPVS to struct
		    n_epvs.(sub).(regions{ii}) = n;
        else
            n_epvs.(sub).(regions{ii}) = 0;
        end
    end
end

%%  Box/Whisker plots for vessel measurements
% X-labels for box/whiskers
% caa6f, caa6o, caa25f, caa26f, caa17o, caa22f, caa22o, caa25o, caa26o
xlab = {'CAA6 Front - 0','CAA6 Occip - 0','CAA25 Front - 0',...
    'CAA26 Front - 0','CAA 26 Occip - 15','CAA 25 Occip - 36',...
    'CAA 22 Occip - 198','CAA17 Occip - 241','CAA22 Front - 620'};
% Retrieve number of observations from each sample
n = zeros(9,1);
n(1) = length(ves_pmus_caa6f);
n(2) = length(ves_pmus_caa6o);
n(3) = length(ves_pmus_25f);
n(4) = length(ves_pmus_26f);
n(5) = length(ves_pmus_26o);
n(6) = length(ves_pmus_25o);
n(7) = length(ves_pmus_22o);
n(8) = length(ves_pmus_17o);
n(9) = length(ves_pmus_22f);
% Create grouping element for box/whisker
grp = [];
for ii=1:length(n)
    grp = [grp, ones(1,n(ii)).*ii];
end

%%% Scattering vs. N EPVS
x = [ves_pmus_caa6f,ves_pmus_caa6o,ves_pmus_25f,ves_pmus_26f,...
    ves_pmus_26o,ves_pmus_25o,ves_pmus_22o,ves_pmus_17o,ves_pmus_22f];
figure('units','normalized','outerposition',[0 0 1 1])
h = boxplot(x,grp,'Notch','on','Labels',xlab);
title('Scattering around Vessels')
xlabel('# EPVS'); ylabel('cm^-^1');
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_vessel_mus_vs_N_epvs.png');
saveas(gcf,fout);

%%% Retardance vs. N EPVS
x = [ves_pret_caa6f,ves_pret_caa6o,ves_pret_25f,ves_pret_26f,...
    ves_pret_26o,ves_pret_25o,ves_pret_22o,ves_pret_17o,ves_pret_22f];
figure('units','normalized','outerposition',[0 0 1 1])
h = boxplot(x,grp,'Notch','on','Labels',xlab);
title('Retardance around Vessels')
xlabel('# EPVS'); ylabel('Degrees')
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_vessel_ret_vs_N_epvs.png');
saveas(gcf,fout);

%%% Circular std of orientation vs. N EPVS
x = [ves_pori_std_caa6f,ves_pori_std_caa6o,ves_pori_std_25f,...
    ves_pori_std_26f,ves_pori_std_26o,ves_pori_std_25o,...
    ves_pori_std_22o,ves_pori_std_17o,ves_pori_std_22f];
figure('units','normalized','outerposition',[0 0 1 1])
h = boxplot(x,grp,'Notch','on','Labels',xlab);
title('Standard Deviation of Orientation around Vessels')
xlabel('# EPVS'); ylabel('Radians')
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_vessel_std_orientation_vs_N_epvs.png');
saveas(gcf,fout);

%% Histograms for optical property vs EPVS w/ just EPVS

%%% X-labels for box/whiskers
% caa6f, caa6o, caa25f, caa26f, caa17o, caa22f, caa22o, caa25o, caa26o
xlab = {'CAA 26 Occip - 15','CAA 25 Occip - 36',...
    'CAA 22 Occip - 198','CAA17 Occip - 241','CAA22 Front - 620'};
% Retrieve number of observations from each sample
n = zeros(9,1);
n(1) = length(ves_pmus_26o);
n(2) = length(ves_pmus_25o);
n(3) = length(ves_pmus_22o);
n(4) = length(ves_pmus_17o);
n(5) = length(ves_pmus_22f);
% Create grouping element for box/whisker
grp = [];
for ii=1:length(n)
    grp = [grp, ones(1,n(ii)).*ii];
end

%%% Scattering vs. N EPVS
x = [ves_pmus_26o,ves_pmus_25o,ves_pmus_22o,ves_pmus_17o,ves_pmus_22f];
figure('units','normalized','outerposition',[0 0 1 1])
h=boxplot(x,grp,'Notch','on','Labels',xlab);
title('Scattering around EPVS')
xlabel('# EPVS'); ylabel('cm^-^1')
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_EPVS_mus_vs_N_epvs.png');
saveas(gcf,fout);

%%% Retardance vs. N EPVS
x = [ves_pret_26o,ves_pret_25o,ves_pret_22o,ves_pret_17o,ves_pret_22f];
figure('units','normalized','outerposition',[0 0 1 1])
h=boxplot(x,grp,'Notch','on','Labels',xlab);
title('Retardance around EPVS')
xlabel('# EPVS'); ylabel('Degrees')
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_EPVS_ret_vs_N_epvs.png');
saveas(gcf,fout);

%%% Circular mean of orientation vs. N EPVS
x = [ves_pori_mean_26o,ves_pori_mean_25o,...
    ves_pori_mean_22o,ves_pori_mean_17o,ves_pori_mean_22f];
figure('units','normalized','outerposition',[0 0 1 1])
h=boxplot(x,grp,'Notch','on','Labels',xlab);
title('Circular Mean of Orientation around EPVS')
xlabel('# EPVS'); ylabel('Radians')
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_EPVS_mean_orientation_vs_N_epvs.png');
saveas(gcf,fout);

%%% Circular std of orientation vs. N EPVS
x = [ves_pori_std_26o,ves_pori_std_25o,...
    ves_pori_std_22o,ves_pori_std_17o,ves_pori_std_22f];
figure('units','normalized','outerposition',[0 0 1 1])
h=boxplot(x,grp,'Notch','on','Labels',xlab);
title('Standard Deviation of Orientation around EPVS')
xlabel('# EPVS'); ylabel('Radians')
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_EPVS_std_orientation_vs_N_epvs.png');
saveas(gcf,fout);
