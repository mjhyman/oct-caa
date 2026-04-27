function swp_compare_gam_slices(results_struct, severities, dirout)
% COMPARE_GAM_ALL_SEVERITIES  Pairwise pointwise difference plots for all
%                              severity combinations in a single figure.
%
%   results_struct : struct where each field is a severity group results
%                    struct from fit_swp_psoct_gam
%   severities     : cell array of group names e.g.
%                    {'control','mild','moderate','severe'}
%   dirout         : output directory for saved figures
%
% EXAMPLE
%   compare_gam_all_severities(results, ...
%       {'control','mild','moderate','severe'}, '/path/to/output');

    n_groups = numel(severities);
    pairs    = nchoosek(1:n_groups, 2);   % all unique pairs, e.g. [1,2],[1,3]...
    n_pairs  = size(pairs, 1);            % 6 for 4 groups

    slice_labels = {'p10', 'p50', 'p90'};
    colors       = [0.271, 0.459, 0.706;       % dark blue  — low EPVS
                    0.878, 0.953, 0.973;       % beige      — mid EPVS
                    0.957, 0.427, 0.263];      % orange/red — high EPVS
    outcomes     = {'scatter', 'retard'};
    ylabels      = {'Scattering difference', 'Retardance difference'};

    % ==================================================================
    % One figure per outcome (scattering, retardance)
    % Layout: n_pairs rows x 1 column of slice curve difference panels
    % ==================================================================
    for o = 1:2

        fig = figure('Position', [60, 60, 700, 280 * n_pairs], 'Color', 'w');

        for p = 1:n_pairs

            g1  = pairs(p, 1);
            g2  = pairs(p, 2);
            r1  = results_struct.(severities{g1});
            r2  = results_struct.(severities{g2});

            n_pts   = min(numel(r1.ves_swp_grid), numel(r2.ves_swp_grid));
            n_slices = numel(r1.slice_epvs_levels);
            swp     = r1.ves_swp_grid(1:n_pts);

            ax = subplot(n_pairs, 1, p);
            hold(ax, 'on');
            yline(ax, 0, 'k--', 'LineWidth', 1, 'Alpha', 0.4);

            for sl = 1:n_slices
                if o == 1
                    c1  = r1.slice_scatter_swp{sl}(1:n_pts);
                    c2  = r2.slice_scatter_swp{sl}(1:n_pts);
                    ci1 = r1.ci_slice_scatter{sl}(1:n_pts, :);
                    ci2 = r2.ci_slice_scatter{sl}(1:n_pts, :);
                else
                    c1  = r1.slice_retard_swp{sl}(1:n_pts);
                    c2  = r2.slice_retard_swp{sl}(1:n_pts);
                    ci1 = r1.ci_slice_retard{sl}(1:n_pts, :);
                    ci2 = r2.ci_slice_retard{sl}(1:n_pts, :);
                end

                col        = colors(sl, :);
                diff_curve = c2 - c1;
                diff_lo    = ci2(:,1) - ci1(:,2);
                diff_hi    = ci2(:,2) - ci1(:,1);

                % Shaded CI band
                x_patch = [swp;     flipud(swp)];
                y_patch = [diff_lo; flipud(diff_hi)];
                fill(ax, x_patch, y_patch, col, ...
                     'FaceAlpha', 0.15, 'EdgeColor', 'none');

                % Mark points where CI excludes zero
                sig = diff_lo > 0 | diff_hi < 0;
                if any(sig)
                    plot(ax, swp(sig), diff_curve(sig), 's', ...
                         'Color', col, 'MarkerSize', 3, ...
                         'MarkerFaceColor', col, 'HandleVisibility', 'off');
                end

                % Main difference curve
                plot(ax, swp, diff_curve, '-', ...
                     'Color',       col,      ...
                     'LineWidth',   2,        ...
                     'DisplayName', sprintf('EPVS %s', slice_labels{sl}));
            end

            hold(ax, 'off');

            % Only show legend on first panel to save space
            if p == 1
                legend(ax, 'Location', 'best', 'FontSize', 8);
            end

            % Only show x-label on bottom panel
            if p == n_pairs
                xlabel(ax, 'Vessel SWP');
            else
                set(ax, 'XTickLabel', []);
            end

            ylabel(ax, ylabels{o}, 'FontSize', 8);
            title(ax, sprintf('%s  minus  %s', ...
                              capitalize(severities{g2}), ...
                              capitalize(severities{g1})), ...
                  'FontSize', 10);

            % Shade zero-crossing reference
            grid(ax, 'on'); box(ax, 'off');
        end

        sgtitle(fig, sprintf('%s: pairwise GAM differences across severities', ...
                             ylabels{o}), ...
                'FontSize', 12, 'FontWeight', 'bold');

        fname = sprintf('GAM_pairwise_%s', outcomes{o});
        saveas(fig, fullfile(dirout, [fname, '.fig']));
        saveas(fig, fullfile(dirout, [fname, '.png']));
        fprintf('Saved: %s\n', fname);
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