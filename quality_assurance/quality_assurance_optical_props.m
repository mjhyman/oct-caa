%% Measure white matter parenchymal Optical Properties
%{
The purpose of this script is to measure the optical properties in the
parenchymal white matter tissue. The average of these values will be
compared against the values in the NIFTI volumes. This is a sanity check.

The full volumes are here:
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)

Outline:
- IMPORT struct containing:
    - tissue masks (exclude background signal)
    - EPVS
    - segmented vasculature (from Etienne)
    - optical properties (mus + retardance)
- MASK
    - Apply WM tissue mask to vessels and optical properties
- Measure average WM value of retardance, scattering, orientation
%}

%% Prepare environment
clear; clc; close all;
% Add top-level directory + subdirectories
addpath(genpath(fullfile(pwd, '..')))
% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/projectnb/npbssmic/ns/CAA/';
% Output figure directory
fig_out = '/projectnb/npbssmic/ns/CAA/bw_plots/';
% Stats output directory
stats_out = '/projectnb/npbssmic/ns/CAA/metrics/';
% Voxel dimensions (microns) for all runs
res = [20,20,20]; % resolution in microns

%%% Flag for loading .MAT struct and creating WM masks
% flag for reloading the .MAT struct for each subject
flag_load_caa_structs = true;

% down sampling factor for the boxplot arrays
ds = 4;

%% Load each subject's .MAT struct and create WM mask
%%% Load the .MAT structs
if flag_load_caa_structs
    % CAA 6
    fprintf('Loading CAA6\n')
    caa6 = load(fullfile(data_dir,'/caa6/caa6.mat'));
    caa6 = caa6.caa6;
    fprintf('Finished Loading CAA6\n')
    % CAA 17
    fprintf('Loading CAA17\n')
    caa17 = load(fullfile(data_dir,'/caa17/occip/caa17.mat'));
    caa17 = caa17.caa17;
    fprintf('Finished Loading CAA17\n')
    % CAA 22
    fprintf('Loading CAA22\n')
    caa22 = load(fullfile(data_dir,'/caa22/caa22.mat'));
    caa22 = caa22.caa22;
    fprintf('Finished Loading CAA22\n')
    % CAA 25
    fprintf('Loading CAA25\n')
    caa25 = load(fullfile(data_dir,'/caa25/caa25.mat'));
    caa25 = caa25.caa25;
    fprintf('Finished Loading CAA25\n')
    % CAA 26
    fprintf('Loading CAA26\n')
    caa26 = load(fullfile(data_dir,'/caa26/caa26.mat'));
    caa26 = caa26.caa26;
    fprintf('Finished Loading CAA26\n')
end

%% Measure mean WM + GM optical props

% Add all subjects to a struct
subjects = struct();
subjects.caa6 = caa6;
subjects.caa17 = caa17;
subjects.caa22 = caa22;
subjects.caa25 = caa25;
subjects.caa26 = caa26;

% Struct for storing WM and GM values
op = struct();

% Iterate over subjects and regions
subs = fields(subjects);
for ii = 1:length(subs)
    % Count number of regions
    reg = fields(subjects.(subs{ii}));
    % iterate over regions
    for j=1:length(reg)
        % Retrieve mus, ret, wm, mask, seg, epvs
        mus = subjects.(subs{ii}).(reg{j}).mus;
        ret = subjects.(subs{ii}).(reg{j}).ret_full;
        wm = subjects.(subs{ii}).(reg{j}).mask_wm;
        mask = subjects.(subs{ii}).(reg{j}).mask;
        seg =  subjects.(subs{ii}).(reg{j}).seg;
        epvs =  subjects.(subs{ii}).(reg{j}).epvs;
        % Set the vessels and EPVS to NaN in mus, ret
        mus(seg) = NaN;
        ret(seg) = NaN;
        mus(epvs) = NaN;
        ret(epvs) = NaN;
        % Measure WM values
        op.(subs{ii}).(reg{j}).wm.mus = single(mus(wm));
        op.(subs{ii}).(reg{j}).wm.ret = single(ret(wm));
        % Measure GM values
        gm = logical(mask .* ~wm);
        op.(subs{ii}).(reg{j}).gm.mus = single(mus(gm));
        op.(subs{ii}).(reg{j}).gm.ret = single(ret(gm));
    end
end

%% Create Table of mean values

rows = {}; % for collecting subject/region rows

subs = fields(op);
for i = 1:length(subs)
    sub = subs{i};
    regs = fields(op.(sub));
    for j = 1:length(regs)
        reg = regs{j};
        % Default values as NaN (in case fields missing)
        entry = {sub, reg, NaN, NaN, NaN, NaN}; % Subject, Region, ret_wm, ret_gm, mus_wm, mus_gm
        % ret_wm
        if isfield(op.(sub).(reg), 'wm') && isfield(op.(sub).(reg).wm, 'ret')
            entry{3} = mean(op.(sub).(reg).wm.ret, 'omitnan');
        end
        % ret_gm
        if isfield(op.(sub).(reg), 'gm') && isfield(op.(sub).(reg).gm, 'ret')
            entry{4} = mean(op.(sub).(reg).gm.ret, 'omitnan');
        end
        % mus_wm
        if isfield(op.(sub).(reg), 'wm') && isfield(op.(sub).(reg).wm, 'mus')
            entry{5} = mean(op.(sub).(reg).wm.mus, 'omitnan');
        end
        % mus_gm
        if isfield(op.(sub).(reg), 'gm') && isfield(op.(sub).(reg).gm, 'mus')
            entry{6} = mean(op.(sub).(reg).gm.mus, 'omitnan');
        end

        rows = [rows; entry];
    end
end

% Convert to table and write to CSV
T = cell2table(rows, 'VariableNames', {'Subject','Region','ret_wm','ret_gm','mus_wm','mus_gm'});
tout = fullfile(stats_out,'mean_op.csv');
writetable(T, tout);
disp('Export complete: all_op_measurements_columns.csv');


%% Box/Whisker White Matter

%%% Mus 
vec = struct();
vec(1).mus = op.caa26.front.wm.mus;
vec(2).mus = op.caa26.occip.wm.mus;
vec(3).mus = op.caa6.front.wm.mus;
vec(4).mus = op.caa6.occip.wm.mus;
vec(5).mus = op.caa17.occip.wm.mus;
vec(6).mus = op.caa25.front.wm.mus;
vec(7).mus = op.caa25.occip.wm.mus;
vec(8).mus = op.caa22.front.wm.mus;
vec(9).mus = op.caa22.occip.wm.mus;

%%% Retardance
vec(1).ret = op.caa26.front.wm.ret;
vec(2).ret = op.caa26.occip.wm.ret;
vec(3).ret = op.caa6.front.wm.ret;
vec(4).ret = op.caa6.occip.wm.ret;
vec(5).ret = op.caa17.occip.wm.ret;
vec(6).ret = op.caa25.front.wm.ret;
vec(7).ret = op.caa25.occip.wm.ret;
vec(8).ret = op.caa22.front.wm.ret;
vec(9).ret = op.caa22.occip.wm.ret;

%%% Create box-whisker vectors for mus and retardance
fprintf('Creating Vectors for WM\n')
[mus_x, mus_y, ret_x, ret_y] = create_vectors(vec,ds);
fprintf('Finished creating Vectors for WM\n')

%%% WM box/whisker plots
% x-axis label
xlab = {'CAA 26 Front', 'CAA 26 Occip', 'CAA 6 Front','CAA 6 Occip',...
        'CAA 17 Occip', 'CAA 25 Front', 'CAA 25 Occip',...
        'CAA 22 Front', 'CAA 22 Occip'};

% mus - WM
figure('units','normalized','outerposition',[0 0 1 1])
h = boxchart(single(mus_x),mus_y);
xticks(1:9); xticklabels(xlab);
title('White Matter \mu_s')
xlabel('Subject / Region'); ylabel('cm^-^1');
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_mus_wm.png');
saveas(gcf,fout); pause(0.5); close;

% ret - WM
figure('units','normalized','outerposition',[0 0 1 1])
h = boxchart(single(ret_x),ret_y);
xticks(1:9); xticklabels(xlab);
title('White Matter Retardance')
xlabel('Subject / Region'); ylabel('Degrees');
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_ret_wm.png');
saveas(gcf,fout); pause(0.5); close;

%% Box/Whisker Gray Matter

%%% Mus 
vec = struct();
vec(1).mus = op.caa26.front.gm.mus;
vec(2).mus = op.caa26.occip.gm.mus;
vec(3).mus = op.caa6.front.gm.mus;
vec(4).mus = op.caa6.occip.gm.mus;
vec(5).mus = op.caa17.occip.gm.mus;
vec(6).mus = op.caa25.front.gm.mus;
vec(7).mus = op.caa25.occip.gm.mus;
vec(8).mus = op.caa22.front.gm.mus;
vec(9).mus = op.caa22.occip.gm.mus;

%%% Retardance
vec(1).ret = op.caa26.front.gm.ret;
vec(2).ret = op.caa26.occip.gm.ret;
vec(3).ret = op.caa6.front.gm.ret;
vec(4).ret = op.caa6.occip.gm.ret;
vec(5).ret = op.caa17.occip.gm.ret;
vec(6).ret = op.caa25.front.gm.ret;
vec(7).ret = op.caa25.occip.gm.ret;
vec(8).ret = op.caa22.front.gm.ret;
vec(9).ret = op.caa22.occip.gm.ret;

%%% Create box-whisker vectors for mus and retardance
fprintf('Creating Vectors for GM\n')
[mus_x, mus_y, ret_x, ret_y] = create_vectors(vec,ds);
fprintf('Finished creating Vectors for GM\n')

%%% WM box/whisker plots
% x-axis label
xlab = {'CAA 26 Front', 'CAA 26 Occip', 'CAA 6 Front','CAA 6 Occip',...
        'CAA 17 Occip', 'CAA 25 Front', 'CAA 25 Occip',...
        'CAA 22 Front', 'CAA 22 Occip'};

% mus - WM
figure('units','normalized','outerposition',[0 0 1 1])
h = boxchart(single(mus_x),mus_y);
xticks(1:9); xticklabels(xlab);
title('Gray Matter \mu_s')
xlabel('Subject / Region'); ylabel('cm^-^1');
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_mus_gm.png');
saveas(gcf,fout); pause(0.5); close;

% ret - WM
figure('units','normalized','outerposition',[0 0 1 1])
h = boxchart(single(ret_x),ret_y);
xticks(1:9); xticklabels(xlab);
title('Gray Matter Retardance')
xlabel('Subject / Region'); ylabel('Degrees');
set(h,{'linew'},{2}); set(gca,'FontSize',20)
fout = fullfile(fig_out,'box_whisker_ret_gm.png');
saveas(gcf,fout); pause(0.5); close;

%% Create vectors for Box/Whisker 
function [mus_x,mus_y,ret_x,ret_y] = create_vectors(vec, ds)
% INPUTS
%   vec (struct): optical properties
%   ds (int): down sampling factor
% OUTPUTS:
%   mus_x (vector): mus x-axis vectors
%   mus_y (vector): mus y-axis vectors
%   ret_x (vector): ret x-axis vectors
%   ret_y (vector): ret y-axis vectors

% Initialize x-axis offset and values
mus_x = [];
mus_y = [];
ret_x = [];
ret_y = [];

% Iterate over vectors
for ii = 1:length(vec)
    %%% mus
    n = length(vec(ii).mus);
    % Indexing for b/w plots
    mus_x = [mus_x, uint8(ones(1,n).*ii)];
    % Vector of values
    mus_y = [mus_y, [vec(ii).mus]'];
    
    %%% retardance
    n = length(vec(ii).ret);
    % Indexing for b/w plots
    ret_x = [ret_x, uint8(ones(1,n).*ii)];
    % Vector of values
    ret_y = [ret_y, [vec(ii).ret]'];
end

mus_x = uint8(mus_x);
ret_x = uint8(ret_x);
mus_y = single(mus_y);
ret_y = single(ret_y);

% Downsample each vector for speed
mus_x = mus_x(1:ds:end);
mus_y = mus_y(1:ds:end);
ret_x = ret_x(1:ds:end);
ret_y = ret_y(1:ds:end);

%%% Remove outliers for efficiency
% mus
out_idx = isoutlier(mus_y);
mus_y = mus_y(~out_idx);
mus_x = mus_x(~out_idx);
% ret
out_idx = isoutlier(ret_y);
ret_y = ret_y(~out_idx);
ret_x = ret_x(~out_idx);


end