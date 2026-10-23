function histo_swp_plot_gam_stages(gam, stain_name, reg, xlims, ylims, dirout)
% HISTO_SWP_PLOT_GAM_STAGES  Plot the ABSOLUTE GAM slice curves for each stage
%   (e.g. control, severe) on SEPARATE figures, one figure per stage. This is
%   the stage-level counterpart to histo_swp_compare_gam_slices, which plots
%   pairwise differences (severe - control) rather than the stages themselves.
%
%   Within each stage figure, vessel-SWP slice level -> color (p10/p50/p90).
%
%   gam        : struct where each field is a stage result struct from
%                fit_swp_histo_gam (fields: epv_swp_grid, slice_ves_levels,
%                slice_epv{:}, ci_slice{:})
%   stain_name : stain name (lhe, cd68, gfap)
%   reg        : brain region string
%   xlims      : x-axis limits
%   ylims      : y-axis limits
%   dirout     : output directory for saved figures

stages   = fields(gam);
n_stages = numel(stages);
if n_stages == 0
    fprintf('  No stages to plot for %s %s\n', stain_name, reg);
    return
end

% Vessel-SWP slice labels + colors (10% dark blue, 50% dark purple, 90% vermillion)
slice_labels = {'p10','p50','p90'};
colors = validatecolor({'#1964b0','#882d71','#DB5829'}, 'multiple');
ylab   = 'Stain Z Score';

if ~exist(dirout,'dir'); mkdir(dirout); end

% ==================================================================
% One figure per stage
% ==================================================================
for gi = 1:n_stages
    stage = stages{gi};
    r     = gam.(stage);

    fig = figure('Position',[100,100,1000,1000]);
    yline(0,'k--','LineWidth',1,'Alpha',0.4,'HandleVisibility','off');
    hold on;

    h_leg = gobjects(0);
    t_leg = {};

    swp      = r.epv_swp_grid(:);
    n_slices = numel(r.slice_ves_levels);

    % Overlay this stage's vessel-SWP slices
    for sl = 1:n_slices
        col = colors(mod(sl-1, size(colors,1)) + 1, :);

        % Absolute slice curve + CI, truncated to a common length defensively
        c  = r.slice_epv{sl}(:);
        ci = r.ci_slice{sl};
        m  = min([numel(swp), numel(c), size(ci,1)]);
        swp_i = swp(1:m);
        c_i   = c(1:m);
        ci_i  = ci(1:m, :);

        % Shaded CI band (not shown in legend)
        fill([swp_i; flipud(swp_i)], [ci_i(:,1); flipud(ci_i(:,2))], col, ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        % Mean slice curve (legend proxy)
        if sl <= numel(slice_labels), slab = slice_labels{sl};
        else,                         slab = sprintf('ves%d', sl); end
        h = plot(swp_i, c_i, '-', 'Color', col, 'LineWidth', 2);
        h_leg(end+1) = h;                          %#ok<AGROW>
        t_leg{end+1} = sprintf('EPVS %s', slab);   %#ok<AGROW>
    end
    hold off;

    ylim(ylims); xlim(xlims);
    legend(h_leg, t_leg, 'Location', 'best', 'FontSize', 8);
    xlabel('EPVS SWP');
    ylabel(ylab, 'FontSize', 8);
    title(sprintf('%s %s %s: GAM (EPVS slices)', ...
                  stain_name, reg, capitalize(stage)), 'FontSize', 10);
    grid on; box on;
    set(gca, 'FontSize', 30);
    pause(0.5)

    % Save output (one file set per stage)
    fname = sprintf('GAM_stages_%s_%s_%s', stain_name, reg, stage);
    saveas(fig, fullfile(dirout, [fname, '.fig']));
    exportgraphics(fig, fullfile(dirout, [fname, '.png']), 'Resolution', 600);
    exportgraphics(fig, fullfile(dirout, [fname, '.svg']), 'Resolution', 600);
    fprintf('Saved: %s\n', fname);
    pause(0.5)
    close all;
end
end


% =========================================================================
% LOCAL HELPER: capitalise first letter of a string
% =========================================================================
function s = capitalize(s)
    if ~isempty(s)
        s(1) = upper(s(1));
    end
end