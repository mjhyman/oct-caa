%% Plot control & severe stages of the Histology GAM (per stain, per region)
% Overview:
%{
- import the main histology GAM struct produced by create_histo_swp_gam.m,
  whose fields are the combined 2x2 groups:
      histo_gam.(stain).{control_front, control_occip, severe_front, severe_occip}
- for each stain, overlay the control and severe GAM slices, with frontal and
  occipital on SEPARATE plots.
%}

%% Initialization
clc; close all;
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
ncores = feature('numcores');
maxNumCompThreads(ncores);

%% Top-level settings
%%% Directories (SCC)
% Input directory
ddir   = '/projectnb/npbssmic/ns/CAA/histology/gam';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/figures/histology_gam/GAM_histo_manuscript_18-Oct-2026/';

%%% Import GAM
% Must be the histo_gam struct saved by create_histo_swp_gam.m
% (file name pattern: GAM_histo_<exploratory|manuscript>_<date>.mat)
gam_file = 'GAM_histo_manuscript_18-Oct-2026.mat';
fprintf('\nImporting histology GAM struct: %s\n', gam_file)
gam = load(fullfile(ddir, gam_file));
gam = gam.histo_gam;

%%% Output figure settings
% Date for file output
t = datetime('today');
% Font size for figures
fsize = 30;

%% Map combined-group fields to region -> stage
% The GAM struct fields are '<stage>_<regioncode>' (e.g. 'control_front',
% 'severe_occip'). Reorganize into gam_remapped.(stain).(region).(stage) so the
% control and severe slices for one region can be overlaid together.
%   region codes: front -> frontal, occip -> occipital
%   stages present in the corrected 2x2 design: control, severe
fprintf('\nRemapping GAM fields to region/stage\n')
region_map = struct('front','frontal', 'occip','occipital');

stains       = fieldnames(gam);
gam_remapped = struct();
for i = 1:numel(stains)
    st = stains{i};
    % Pre-initialize so downstream isfield checks are safe
    gam_remapped.(st) = struct('frontal', struct(), 'occipital', struct());

    grps = fieldnames(gam.(st));
    for j = 1:numel(grps)
        g     = grps{j};              % e.g. 'control_front'
        parts = split(g, '_');        % {'control','front'}
        if numel(parts) < 2,        continue; end
        stage = parts{1};             % control / severe
        rcode = parts{2};             % front / occip
        if ~isfield(region_map, rcode), continue; end
        region = region_map.(rcode);

        gam_remapped.(st).(region).(stage) = gam.(st).(g);
    end
end

%% Plot control and severe for each stain + region (separate plots)
fprintf('\nCreating histo vs SWP GAM figures (control vs severe)\n')
stains  = fieldnames(gam_remapped);
regions = {'frontal','occipital'};   % frontal and occipital on separate plots

% Per-stain axis limits (SWP on x, stain z-score on y). Tune per stain; the
% fallback covers any stain not listed in the lookup.
plot_lims = struct( ...
    'lhe',  struct('xlim',[0 10], 'ylim',[-0.5 0.5]), ...
    'gfap', struct('xlim',[0 10], 'ylim',[-1 1]), ...
    'cd68', struct('xlim',[0 10], 'ylim',[-0.25 0.25]));
default_lims = struct('xlim',[0 10], 'ylim',[-1.5 1.5]);

for ii = 1:numel(stains)
    st = stains{ii};

    % Resolve this stain's axis limits (both regions share them)
    if isfield(plot_lims, st)
        lims = plot_lims.(st);
    else
        lims = default_lims;
    end
    xlims = lims.xlim;
    ylims = lims.ylim;

    for j = 1:numel(regions)
        reg = regions{j};

        % Skip if this stain/region has no fitted stages
        if ~isfield(gam_remapped.(st), reg) || ...
                isempty(fieldnames(gam_remapped.(st).(reg)))
            fprintf('  Skipping %s %s (no control/severe data)\n', st, reg);
            continue
        end

        fprintf('\n----%s %s----\n', st, reg)
        dirout = fullfile(figdir, st, reg);

        % Overlay the absolute control/severe slice curves (NOT the difference)
        histo_swp_plot_gam_stages( ...
            gam_remapped.(st).(reg), ...
            st, ...
            reg, ...
            xlims, ...
            ylims, dirout);
    end
end

%% Plot severe - control difference surfaces (from the separate diff .MAT)
% Loads the histo_gam_diff struct saved by create_histo_swp_gam.m into
% .../gam/severe_minus_control/ and renders each stain x region difference
% surface (severe - control over the epv_swp x ves_swp plane).
fprintf('\nPlotting severe - control difference surfaces\n')
diffdir_in  = fullfile(ddir,   'severe_minus_control');   % where the .MAT lives
diffdir_out = fullfile(figdir, 'severe_minus_control');   % where figures go
regttl_map  = struct('front','Frontal', 'occip','Occipital');

% Filename written by create_histo_swp_gam.m (GAM_histo_diff_<tag>_<date>.mat).
% Leave '' to auto-pick the most recent diff file.
diff_file = '';
if isempty(diff_file)
    dfl = dir(fullfile(diffdir_in, 'GAM_histo_diff_*.mat'));
    if isempty(dfl)
        warning('No GAM_histo_diff_*.mat found in %s; skipping difference plots.', diffdir_in);
    else
        [~, ix]   = max([dfl.datenum]);
        diff_file = dfl(ix).name;
    end
end

if ~isempty(diff_file)
    fprintf('  Using difference file: %s\n', diff_file)
    gd = load(fullfile(diffdir_in, diff_file));
    gd = gd.histo_gam_diff;

    stains_d = fieldnames(gd);
    for ii = 1:numel(stains_d)
        st     = stains_d{ii};
        regs_d = fieldnames(gd.(st));

        % Per-stain symmetric color limit so both regions share one scale
        allabs = [];
        for j = 1:numel(regs_d)
            D = gd.(st).(regs_d{j});
            if isfield(D,'diff'), allabs = [allabs; abs(D.diff(:))]; end %#ok<AGROW>
        end
        cval = max(allabs, [], 'omitnan');
        if isempty(cval) || ~isfinite(cval) || cval == 0, cval = 1; end

        for j = 1:numel(regs_d)
            reg = regs_d{j};
            D   = gd.(st).(reg);
            if ~isfield(D,'diff') || isempty(D.diff)
                fprintf('  Skipping %s %s (empty diff)\n', st, reg);
                continue
            end
            if isfield(regttl_map, reg), regttl = regttl_map.(reg); else, regttl = reg; end

            dirout = fullfile(diffdir_out, st);
            if ~exist(dirout,'dir'), mkdir(dirout); end
            outstem = fullfile(dirout, sprintf('GAM_diff_%s_%s_severe-control', st, reg));
            plot_diff_surface(D, sprintf('%s %s: severe - control', upper(st), regttl), ...
                              outstem, cval);
            fprintf('  Saved diff surface: %s %s\n', st, reg);
        end
    end
end

fprintf('\nDone, go graduate and get a real job.\n')


%% ------------------------------------------------------------------------
%% Local function: plot and save a severe - control difference surface
function plot_diff_surface(D, tstr, outstem, clim_val)
if ~isfield(D, 'diff') || isempty(D.diff), return; end
if nargin < 4 || isempty(clim_val)
    clim_val = max(abs(D.diff(:)), [], 'omitnan');
    if isempty(clim_val) || ~isfinite(clim_val) || clim_val == 0, clim_val = 1; end
end
figure('Position', [100,100,900,750]);
imagesc(D.epv_vec, D.ves_vec, D.diff, 'AlphaData', ~isnan(D.diff));
set(gca, 'YDir', 'normal', 'Color', [0.92 0.92 0.92]); % grey = outside support
axis tight; hold on;
caxis([-clim_val clim_val]);
colormap(bluewhitered_local());
cb = colorbar; cb.Label.String = 'severe - control (z-score)';

% Outline where the bootstrap CI excludes zero (only if bootstrap was run)
if isfield(D, 'do_bootstrap') && D.do_bootstrap && isfield(D, 'ci_excludes_zero')
    [~, hc] = contour(D.epv_vec, D.ves_vec, double(D.ci_excludes_zero), [0.5 0.5]);
    set(hc, 'LineColor', 'k', 'LineWidth', 1.2);
end

xlabel('epv\_swp'); ylabel('ves\_swp'); title(tstr);
set(gca, 'fontname', 'Arial');
saveas(gcf, outstem, 'png');
saveas(gcf, outstem, 'fig');
exportgraphics(gcf, strcat(outstem, '.svg'), 'Resolution', 600);
close(gcf);
end


%% Local function: diverging blue-white-red colormap (self-contained)
function cmap = bluewhitered_local(n)
if nargin < 1, n = 256; end
b = [0.23 0.30 0.75]; w = [1 1 1]; r = [0.75 0.15 0.15];
half = floor(n/2);
top  = [linspace(b(1),w(1),half)', linspace(b(2),w(2),half)', linspace(b(3),w(3),half)'];
bot  = [linspace(w(1),r(1),n-half)', linspace(w(2),r(2),n-half)', linspace(w(3),r(3),n-half)'];
cmap = [top; bot];
end


