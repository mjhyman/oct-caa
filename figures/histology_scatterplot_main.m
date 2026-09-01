%% Histopathology: absolute mean +/- SD per distance (control & severe) + difference
%{
For each stain (LHE, GFAP, CD68) and region (frontal, occipital) this script
plots, at each annular-ring distance:
  - the ABSOLUTE mean +/- SD for the control group (separate figure)
  - the ABSOLUTE mean +/- SD for the severe group  (separate figure)
  - the severe - control difference               (separate figure)

Intensity units = z-score w.r.t. the section's median stain intensity.
Subject is the unit of analysis: the per-distance value is averaged over
subjects and the SD is the spread across subject means.

Corrected 2x2 design (mild/moderate removed as erroneous leftovers):
  control = {CAA_26, CAA_6}   severe = {CAA_22, CAA_25}
%}

%% Top-level settings
clearvars -except caa6 caa17 caa22 caa25 caa26 swp_struct
clc; close all;
% Get the current folder
currentFolder = pwd;
% Move one directory up
parentFolder = fileparts(currentFolder);
% Add the parent folder and all its subfolders to the MATLAB search path
addpath(genpath(parentFolder));

%%% Directories (SCC)
% Input directory
hdir = '/projectnb/npbssmic/ns/CAA/histology/';
% Output directories: absolute curves vs. severe-control differences
absdir  = ['/projectnb/npbssmic/ns/CAA/figures/' ...
           'histology_scatterplots/histo_matched/absolute/'];
diffdir = ['/projectnb/npbssmic/ns/CAA/figures/' ...
           'histology_scatterplots/histo_matched/stage_diff/'];
if ~exist(absdir,'dir'),  mkdir(absdir);  end
if ~exist(diffdir,'dir'), mkdir(diffdir); end

%%% Load the stain and segmentation structs
lhe  = load(fullfile(hdir,'lhe_rings_21-Nov-2025.mat'));  lhe  = lhe.lhe;
gfap = load(fullfile(hdir,'gfap_rings_21-Nov-2025.mat')); gfap = gfap.gfap;
cd68 = load(fullfile(hdir,'cd68_rings_21-Nov-2025.mat')); cd68 = cd68.cd68;

%%% measurement outer radii (units = voxels)
radii_sm = [40, 81, 121, 162, 202, 243, 283, 324, 364, 405, 445, 486];

%% Subject identifiers for each stage (corrected 2x2 design)
ctl_subs = {'CAA_26','CAA_6'};
sev_subs = {'CAA_22','CAA_25'};

%% Per-stain settings and the main loop
% One driver replaces the previous copy-paste blocks. Absolute and difference
% plots reuse the same per-stain y-limits; widen them if a curve clips.
stain_names = {'lhe','gfap','cd68'};
stain_data  = {lhe,   gfap,   cd68};
stain_ylim  = {[-0.5 0.5], [-1.2 1.2], [-1.2 1.2]};
stain_ystep = {0.1, 0.2, 0.2};

region_tags  = {'frontal','occipital'};
region_codes = {"_1","_7"};             % frontal = _1, occipital = _7
region_names = {'Frontal','Occipital'};

close all;
for si = 1:numel(stain_names)
    sname = stain_names{si};
    sdata = stain_data{si};
    ylims = stain_ylim{si};
    ytick = ylims(1):stain_ystep{si}:ylims(2);

    % Group substructs by subject ID (fills the old "find indices" stubs and
    % works for both baseName conventions: CAA_26_1 vs CAA26_1_GFAP_10x)
    ctl_group = extract_subjects(sdata, ctl_subs);
    sev_group = extract_subjects(sdata, sev_subs);

    for ri = 1:numel(region_tags)
        rtag  = region_tags{ri};
        rcode = region_codes{ri};
        rname = region_names{ri};

        % Group averages: mean + between-subject SD at each distance
        ctl = average_within_region(ctl_group, rcode, "rad40");
        sev = average_within_region(sev_group, rcode, "rad40");

        % --- Absolute: control -------------------------------------------
        tstr  = sprintf('%s %s (Control, n=%d)', upper(sname), rname, ctl.nsubs);
        fname = fullfile(absdir, sprintf('%s_%s_control_ribbon', sname, rtag));
        plot_and_save(radii_sm, ctl, ylims, ytick, tstr, fname);

        % --- Absolute: severe --------------------------------------------
        tstr  = sprintf('%s %s (Severe, n=%d)', upper(sname), rname, sev.nsubs);
        fname = fullfile(absdir, sprintf('%s_%s_severe_ribbon', sname, rtag));
        plot_and_save(radii_sm, sev, ylims, ytick, tstr, fname);

        % --- Difference: severe - control --------------------------------
        d = struct();
        [d.std_epvs, d.mean_epvs, d.std_ves, d.mean_ves] = stage_diff(ctl, sev);
        tstr  = sprintf('%s %s (Severe - Control)', upper(sname), rname);
        fname = fullfile(diffdir, sprintf('%s_%s_severe-ctl_ribbon', sname, rtag));
        plot_and_save(radii_sm, d, ylims, ytick, tstr, fname);
    end
end


%% ------------------------------------------------------------------------
%% Helper: extract a group's sections by subject ID
function sub = extract_subjects(stain, subs)
% Return the sections of `stain` whose baseName belongs to any subject in
% `subs`. Handles both naming conventions in this dataset:
%   LHE:        'CAA_26_1'          (underscore after CAA)
%   GFAP/CD68:  'CAA26_1_GFAP_10x'  (no underscore, stain suffix)
names = string({stain.baseName});
keep  = false(size(names));
for i = 1:numel(subs)
    s_us = string(subs{i});      % e.g. "CAA_26"
    s_no = erase(s_us, "_");     % e.g. "CAA26"
    keep = keep | contains(names, s_us) | contains(names, s_no);
end
sub = stain(keep);
end


%% Helper: average within one region across subjects
function [stats] = average_within_region(stain, region, rad_size)
% Mean stain intensity and between-subject SD across subjects within one brain
% region, at each annular-ring distance. Subject is the unit of analysis: the
% per-distance mean is averaged over subjects and the SD is the spread across
% those subject means (so n_subjects = 1 -> SD = 0).
%
% INPUTS
%   stain (struct array): sections for ONE group (e.g. all control sections)
%   region (string): "_1" (frontal) or "_7" (occipital)
%   rad_size (string): "rad40"
% OUTPUT
%   stats (struct): mean_epvs/std_epvs, mean_ves/std_ves (nrad x 1),
%                   n_epvs/n_ves (total voxels per distance), nsubs

stats = struct();
if isempty(stain)
    warning('average_within_region: empty group passed; returning empty stats.');
    return;
end

% Number of radii (from the group's first section; region subset may be empty)
nrad     = length(fields(stain(1).(rad_size)));
rad_name = fields(stain(1).(rad_size));

% Sections matching the requested region (both naming conventions)
names = {stain.baseName};
if endsWith(names{1},"_10x")
    tf = contains(names, strcat(region,"_"));
else
    tf = endsWith(names, region);
end
stain_reg = stain(tf);
nsubs     = numel(stain_reg);

% No section for this region -> return NaN-filled stats rather than erroring
if nsubs == 0
    warning('average_within_region: no sections for region %s; returning NaNs.', region);
    stats.n_epvs = nan(nrad,1); stats.n_ves = nan(nrad,1);
    stats.mean_epvs = nan(nrad,1); stats.std_epvs = nan(nrad,1);
    stats.mean_ves  = nan(nrad,1); stats.std_ves  = nan(nrad,1);
    stats.nsubs = 0;
    return;
end

% Per-subject, per-distance mean and voxel count
n_epvs    = zeros(nrad,nsubs);  n_ves    = zeros(nrad,nsubs);
mean_epvs = zeros(nrad,nsubs);  mean_ves = zeros(nrad,nsubs);

for s = 1:nsubs
    for ii = 1:nrad
        epvs = stain_reg(s).(rad_size).(rad_name{ii}).exp;
        ves  = stain_reg(s).(rad_size).(rad_name{ii}).ctl;
        n_epvs(ii,s)    = numel(epvs);
        n_ves(ii,s)     = numel(ves);
        mean_epvs(ii,s) = mean(epvs,'omitnan');
        mean_ves(ii,s)  = mean(ves, 'omitnan');
    end
end

% Aggregate across subjects (generalized from the old hard-coded 2-subject sum)
stats.n_epvs    = sum(n_epvs, 2);
stats.n_ves     = sum(n_ves,  2);
stats.mean_epvs = mean(mean_epvs, 2, 'omitnan');
stats.std_epvs  = std(mean_epvs, 0, 2, 'omitnan');
stats.mean_ves  = mean(mean_ves,  2, 'omitnan');
stats.std_ves   = std(mean_ves,  0, 2, 'omitnan');
stats.nsubs     = nsubs;
end


%% Helper: difference between severities (severe - control)
function [std_epvs, mean_epvs, std_ves, mean_ves] = stage_diff(ctl, exp)
% EPVS
mean_epvs = exp.mean_epvs - ctl.mean_epvs;
std_epvs  = sqrt((exp.std_epvs).^2 ./ exp.n_epvs + ...
                 (ctl.std_epvs).^2 ./ ctl.n_epvs);
% Vessel
mean_ves  = exp.mean_ves - ctl.mean_ves;
std_ves   = sqrt((exp.std_ves).^2 ./ exp.n_ves + ...
                 (ctl.std_ves).^2 ./ ctl.n_ves);
end


%% Helper: draw one ribbon figure and save png/fig/svg
function plot_and_save(radii, stats, ylims, ytick, tstr, fname)
figure('Position',[100,100,1000,1000]);
histo_mean_std_plot(radii, stats.mean_epvs, stats.std_epvs, ...
                    stats.mean_ves, stats.std_ves, ylims, ytick, tstr);
saveas(gcf, fname, 'png');
saveas(gcf, fname, 'fig');
exportgraphics(gcf, strcat(fname,'.svg'), 'Resolution', 600);
end


%% Helper (retained for the subject-level TODO; not used in the main flow)
function [stats] = iterate_subject_radii(stain, rad_size)
% Per-section mean + within-section SD at each distance, keyed by baseName.
% Kept for future subject-level scatterplots; the group-level flow above uses
% average_within_region instead.
nrad     = length(fields(stain(1).(rad_size)));
rad_name = fields(stain(1).(rad_size));
subs     = {stain(:).baseName};
stats    = struct();

for s = 1:length(subs)
    n_epvs = zeros(nrad,1); n_ves = zeros(nrad,1);
    mean_epvs = zeros(nrad,1); std_epvs = zeros(nrad,1);
    mean_ves  = zeros(nrad,1); std_ves  = zeros(nrad,1);
    for ii = 1:nrad
        epvs = stain(s).(rad_size).(rad_name{ii}).exp;
        ves  = stain(s).(rad_size).(rad_name{ii}).ctl;
        n_epvs(ii)    = numel(epvs);
        n_ves(ii)     = numel(ves);
        mean_epvs(ii) = mean(epvs,'omitnan');
        std_epvs(ii)  = std(epvs,'omitnan');
        mean_ves(ii)  = mean(ves,'omitnan');
        std_ves(ii)   = std(ves,'omitnan');
    end
    stats.(subs{s}).n_epvs    = n_epvs;
    stats.(subs{s}).n_ves     = n_ves;
    stats.(subs{s}).mean_epvs = mean_epvs;
    stats.(subs{s}).std_epvs  = std_epvs;
    stats.(subs{s}).mean_ves  = mean_ves;
    stats.(subs{s}).std_ves   = std_ves;
end
end