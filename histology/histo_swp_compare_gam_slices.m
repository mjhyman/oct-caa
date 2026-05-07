function histo_swp_compare_gam_slices(gam, stain_name, reg,...
                                      xlims, ylims, dirout)
% COMPARE_GAM_ALL_SEVERITIES  Pairwise pointwise difference plots for all
%                              severity combinations in a single figure.
%                              The values of the stain pixels are z scores
%
%   gam            : struct where each field is a severity group results
%                    struct from fit_swp_psoct_gam
%   stain_name     : the stain name (lhe, cd68, gfap)
%   reg (str)      : brain region string
%   ylims (vector) : y-axis limits
%   dirout         : output directory for saved figures

% Retrieve the number of severities within the GAM
severities = fields(gam);
n_groups = numel(fields(gam));
pairs    = nchoosek(1:n_groups, 2);   % all unique pairs, e.g. [1,2],[1,3]...
n_pairs  = size(pairs, 1);            % 6 for 4 groups
% Labels for each vessel-swp slice
slice_labels = {'p10', 'p50', 'p90'};
% Colors for each vessel-swp slice
% 10% = dark blue, 50% = dark teal, 90% = vermillion
colors = {'#1964b0','#882d71','#DB5829'};
colors = validatecolor(colors,"multiple");
ylab = 'Stain Z Score';

% ==================================================================
% Iterate the severity combinations
% ==================================================================
% Iterate pairs
for p = 1:n_pairs
    fig = figure('Units','Normalized','Position',[0, 0, 0.5, 0.9],...
                 'Resize','off');
    % Extract pairs
    g1  = pairs(p, 1);
    g2  = pairs(p, 2);

    % If g2 is control, then assign control to g1
    if strcmpi(severities{g2},'control')
        tmp = g1;
        g1 = g2;
        g2 = tmp;
    end

    % Assign r1 (baseline) and r2 (comparison)
    r1  = gam.(severities{g1});
    r2  = gam.(severities{g2});
    
    %%% Find grid with fewest number of elements
    n1 = numel(r1.epv_swp_grid);
    n2 = numel(r2.epv_swp_grid);
    % If # elements is equivalent, choose smallest range
    if n1 == n2
        % Set the number of points (either n1, n2 work)
        n_pts = n1;
        n1_max = max([r1.epv_swp_grid]);
        n2_max = max([r2.epv_swp_grid]);
        if n1_max < n2_max
            idx = 1;
        else
            idx = 2;
        end
    % If not equivalent, then find array with fewer elements
    else
        % Find minimum # elements
        [n_pts, idx] = min([n1, n2]);
    end
    % Select the appropriate SWP grid and # slices
    if idx == 1
        swp = r1.epv_swp_grid(1:n_pts);
        n_slices = numel(r1.slice_ves_levels);
    else
        swp = r2.epv_swp_grid(1:n_pts);
        n_slices = numel(r2.slice_ves_levels);
    end        

    % Create baseline dashed line (control)
    yline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.4);
    hold on;
    
    % Overaly all slices
    for sl = 1:n_slices
        % Retrieve the slices + CIs
        c1  = r1.slice_epv{sl}(1:n_pts);
        c2  = r2.slice_epv{sl}(1:n_pts);
        ci1 = r1.ci_slice{sl}(1:n_pts, :);
        ci2 = r2.ci_slice{sl}(1:n_pts, :);
                    
        % Create difference curve and confidence interval            
        diff_curve = c2 - c1;
        diff_lo    = ci2(:,1) - ci1(:,2);
        diff_hi    = ci2(:,2) - ci1(:,1);

        %%% Shaded CI band
        % Select color for ploting
        col        = colors(sl, :);
        x_patch = [swp;     flipud(swp)];
        y_patch = [diff_lo; flipud(diff_hi)];
        fill(x_patch, y_patch, col, ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none');

        % Mark points where CI excludes zero
        sig = diff_lo > 0 | diff_hi < 0;
        if any(sig)
            plot(swp(sig), diff_curve(sig), 's', ...
                 'Color', col, 'MarkerSize', 3, ...
                 'MarkerFaceColor', col, 'HandleVisibility', 'off');
        end

        %%% Main difference curve
        plot(swp, diff_curve, '-', ...
             'Color',       col,      ...
             'LineWidth',   2,        ...
             'DisplayName', sprintf('EPVS %s', slice_labels{sl}));
    end
    
    % Disable Hold and set y-axis limits
    hold off;
    ylim(ylims);
    xlim(xlims);

    % Only show legend on first panel to save space
    lstr = {'Baseline','','10th Percentile','','50th Percentile',...
             '','90th Percentile'};
    if p == 1
        legend(lstr,'Location', 'best', 'FontSize', 8);
    end
    
    % x-axis + y-axis labels
    xlabel('EPVS SWP');
    ylabel(ylab,'FontSize', 8);
    tstr = sprintf('%s %s: %s  minus  %s',...
                    stain_name,...
                    reg,...
                    capitalize(severities{g2}),...
                    capitalize(severities{g1}));
    title(tstr,'FontSize', 10);

    % Shade zero-crossing reference
    grid on; box on;
    set(gca,'Fontsize',30)
    pause(0.5)

    % Save output
    fname = sprintf('GAM_pairwise_%s_%s_minus_%s',...
                    stain_name,severities{g2},severities{g1});
    saveas(fig, fullfile(dirout, [fname, '.fig']));
    saveas(fig, fullfile(dirout, [fname, '.png']));
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