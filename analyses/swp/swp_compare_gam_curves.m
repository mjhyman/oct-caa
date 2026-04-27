function swp_compare_gam_curves(gam_struct, severities, dirout)
% COMPARE_GAM_SURFACES_MATRIX  Matrix of pairwise 2D surface differences.
%   Upper triangle shows severe-minus-milder differences.
%   Diagonal is blank. Lower triangle mirrors upper (negated).
%   One figure per outcome.

    n_groups = numel(severities);
    outcomes = {'scatter', 'retard'};
    titles   = {'Scattering', 'Retardance'};

    % Build common grid from union of all group grids
    all_swp_min  = min(cellfun(@(g) min(gam_struct.(g).ves_swp_grid),  severities));
    all_swp_max  = max(cellfun(@(g) max(gam_struct.(g).ves_swp_grid),  severities));
    all_epvs_min = min(cellfun(@(g) min(gam_struct.(g).epv_ves_grid),  severities));
    all_epvs_max = max(cellfun(@(g) max(gam_struct.(g).epv_ves_grid),  severities));

    n_pts    = 100;
    swp_com  = linspace(all_swp_min,  all_swp_max,  n_pts);
    epvs_com = linspace(all_epvs_min, all_epvs_max, n_pts);

    [SWP_com, EPVS_com] = meshgrid(swp_com, epvs_com);
    T_common = table(SWP_com(:), EPVS_com(:), SWP_com(:) .* EPVS_com(:), ...
        'VariableNames', {'ves_swp', 'epv_swp', 'ves_epv_swp'});

    % Pre-predict all groups on the common grid
    pred = struct();
    for g = 1:n_groups
        grp = severities{g};
        pred.(grp).scatter = reshape( ...
            predict(gam_struct.(grp).gam_scatter, T_common), n_pts, n_pts);
        pred.(grp).retard  = reshape( ...
            predict(gam_struct.(grp).gam_retard,  T_common), n_pts, n_pts);
    end

    for o = 1:2
        fig = figure('Position', [60, 60, 320*n_groups, 280*n_groups], 'Color', 'w');
        outcome = outcomes{o};
        
        % 1. Set Figure-level Colormap immediately so all children inherit it
        cmap_256 = redblue_cmap(256);
        colormap(fig, cmap_256);

        % 2. Compute all pairwise differences to set shared colour limits
        all_diffs = [];
        pairs = nchoosek(1:n_groups, 2);
        for p = 1:size(pairs,1)
            g1 = pairs(p,1); g2 = pairs(p,2);
            d  = pred.(severities{g2}).(outcome) - ...
                 pred.(severities{g1}).(outcome);
            all_diffs = [all_diffs; d(:)]; %#ok<AGROW>
        end
        
        % 3. Force absolute symmetry around zero
        max_abs_diff = max(abs(all_diffs(:)), [], 'omitnan');
        if isempty(max_abs_diff) || max_abs_diff == 0, max_abs_diff = 1e-5; end
        clim_diff = [-max_abs_diff, max_abs_diff];

        panel = 0;
        last_ax = []; % To store a handle for the colorbar reference
        for row = 1:n_groups
            for col = 1:n_groups
                panel = panel + 1;
                ax = subplot(n_groups, n_groups, panel);

            % Upper triangle    
                if col > row
                    % Off-diagonal — show difference surface
                    g1   = severities{row};
                    g2   = severities{col};
                    diff_data = pred.(g2).(outcome) - pred.(g1).(outcome);

                    imagesc(ax, swp_com, epvs_com, diff_data);
                    % CRITICAL: Match axes limits to the global symmetric scale
                    set(ax, 'YDir', 'normal', 'CLim', clim_diff);
                    last_ax = ax; 

                    % Formatting
                    if row == n_groups
                        xlabel(ax, 'Vessel SWP', 'FontSize', 7);
                    else
                        set(ax, 'XTickLabel', []);
                    end
                    if col == 1
                        ylabel(ax, 'EPVS SWP', 'FontSize', 7);
                    else
                        set(ax, 'YTickLabel', []);
                    end

                    title(ax, sprintf('%s - %s', upper(severities{col}),...
                        upper(severities{row})), ...
                          'FontSize', 8);
                else
                    % Diagonal ��� show group label
                    text(ax, 0.5, 0.5, upper(severities{row}), ...
                         'HorizontalAlignment', 'center', ...
                         'FontSize', 12, 'FontWeight', 'bold', ...
                         'Units', 'normalized');
                    axis(ax, 'off');
                end
            end
        end

        % 4. Create single shared colorbar linked to the symmetric limits
        if ~isempty(last_ax)
            cb = colorbar(last_ax, 'Position', [0.93, 0.1, 0.02, 0.8]);
            cb.Label.String = sprintf('%s difference', titles{o});
            % Force the colorbar to use the same symmetric limits
            set(cb, 'Limits', clim_diff);
        end

        sgtitle(fig, sprintf('%s GAM surface differences', titles{o}), ...
                'FontSize', 13, 'FontWeight', 'bold');

        fname = sprintf('GAM_surface_matrix_%s', outcome);
        saveas(fig, fullfile(dirout, [fname, '.fig']));
        saveas(fig, fullfile(dirout, [fname, '.png']));
        fprintf('Saved: %s\n', fname);
    end
end

% =========================================================================
% LOCAL HELPER: red-blue diverging colormap centred at zero
% =========================================================================
function cmap = redblue_cmap(n)
    % Create a diverging Blue-White-Red colormap
    % Blue (Negative) -> White (Zero) -> Red (Positive)
    m = floor(n/2);
    
    % Blue side: Red and Green climb to 1, Blue stays 1
    blue_to_white = [linspace(0, 1, m)', linspace(0, 1, m)', ones(m, 1)];
    
    % Red side: Green and Blue drop to 0, Red stays 1
    white_to_red = [ones(m, 1), linspace(1, 0, m)', linspace(1, 0, m)'];
    
    if mod(n, 2) == 0
        cmap = [blue_to_white; white_to_red];
    else
        % Ensure exactly one white row in the middle for odd n
        cmap = [blue_to_white; [1 1 1]; white_to_red];
    end
end