function toolkit = gam_slice_toolkit()
% GAM_SLICE_TOOLKIT  Gateway function returning handles to shared local
%                     helpers used by fit_isosurface_by_dataset.m and
%                     compare_isosurfaces_severe_control.m, so the two
%                     top-level functions do not duplicate (and risk
%                     diverging on) fitting/prediction/plotting logic.
%
% USAGE
%   toolkit = gam_slice_toolkit();
%   T_work  = toolkit.clean_and_subsample(T, 1e5, true);
%
% This file must be on the MATLAB path alongside tune_swp_psoct_gam.m,
% fit_isosurface_by_dataset.m, and compare_isosurfaces_severe_control.m.

    toolkit.clean_and_subsample = @clean_and_subsample;
    toolkit.fit_gam_pair        = @fit_gam_pair;
    toolkit.predict_slices      = @predict_slices;
    toolkit.slice_density_mask  = @slice_density_mask;
    toolkit.bootstrap_slices    = @bootstrap_slices;
    toolkit.plot_slices_single  = @plot_slices_single;
    toolkit.plot_slices_compare = @plot_slices_compare;
end


% =========================================================================
% Clean required columns, drop NaN/Inf rows, stratified-subsample by
% ves_swp if needed, and engineer the ves_epv_swp interaction column.
% =========================================================================
function T_work = clean_and_subsample(T, MaxSample, Verbose)
    required_cols = {'scattering', 'retardance', 'ves_swp', 'epv_swp'};
    missing = required_cols(~ismember(required_cols, T.Properties.VariableNames));
    if ~isempty(missing)
        error('gam_slice_toolkit:missingColumns', ...
              'Table is missing required column(s): %s', strjoin(missing, ', '));
    end

    T_work = T(:, required_cols);
    arr    = table2array(T_work);
    bad    = any(isnan(arr) | isinf(arr), 2);
    if Verbose && any(bad)
        fprintf('[gam_slice_toolkit] Removed %d rows with NaN/Inf.\n', sum(bad));
    end
    T_work(bad, :) = [];

    if height(T_work) < 10
        error('gam_slice_toolkit:insufficientData', ...
              'Too few valid rows (%d) after removing NaN/Inf.', height(T_work));
    end

    N = height(T_work);
    if N > MaxSample
        if Verbose
            fprintf('[gam_slice_toolkit] Subsampling %d -> %d voxels (stratified by ves_swp).\n', ...
                    N, round(MaxSample));
        end
        n_strata = 10;
        edges    = quantile(T_work.ves_swp, linspace(0, 1, n_strata + 1));
        edges(1) = -Inf; edges(end) = Inf;
        n_per    = floor(MaxSample / n_strata);
        keep_idx = [];
        for s = 1:n_strata
            in_s = find(T_work.ves_swp > edges(s) & T_work.ves_swp <= edges(s+1));
            if isempty(in_s), continue; end
            n_draw   = min(n_per, numel(in_s));
            keep_idx = [keep_idx; in_s(randperm(numel(in_s), n_draw))]; %#ok<AGROW>
        end
        T_work = T_work(keep_idx, :);
    end

    T_work.ves_epv_swp = T_work.ves_swp .* T_work.epv_swp;
end


% =========================================================================
% Fit scattering + retardance GAMs at a fixed complexity.
% =========================================================================
function gam_pair = fit_gam_pair(T_work, MaxSplits, NumTrees, LearnRate)
    gam_pair.gam_scatter = fitrgam(T_work, 'scattering ~ ves_swp + epv_swp + ves_epv_swp', ...
        'NumTreesPerPredictor',          NumTrees,  ...
        'MaxNumSplitsPerPredictor',      MaxSplits, ...
        'InitialLearnRateForPredictors', LearnRate);
    gam_pair.gam_retard = fitrgam(T_work, 'retardance ~ ves_swp + epv_swp + ves_epv_swp', ...
        'NumTreesPerPredictor',          NumTrees,  ...
        'MaxNumSplitsPerPredictor',      MaxSplits, ...
        'InitialLearnRateForPredictors', LearnRate);
    gam_pair.MaxSplits = MaxSplits;
    gam_pair.NumTrees  = NumTrees;
    gam_pair.LearnRate = LearnRate;
end


% =========================================================================
% Predict scattering/retardance ~ epv_swp at each fixed ves_swp level.
% =========================================================================
function [slice_scatter, slice_retard] = predict_slices(gam_pair, epv_grid, ves_levels)
    epv_grid = epv_grid(:);
    n_pts    = numel(epv_grid);
    n_slices = numel(ves_levels);
    slice_scatter = cell(n_slices, 1);
    slice_retard  = cell(n_slices, 1);
    for sl = 1:n_slices
        ves_fixed = ves_levels(sl) * ones(n_pts, 1);
        T_sl = table(ves_fixed, epv_grid, ves_fixed .* epv_grid, ...
            'VariableNames', {'ves_swp', 'epv_swp', 'ves_epv_swp'});
        slice_scatter{sl} = predict(gam_pair.gam_scatter, T_sl);
        slice_retard{sl}  = predict(gam_pair.gam_retard,  T_sl);
    end
end


% =========================================================================
% Data-density mask per slice: true = fewer than MinDensity voxels near
% (ves_level, epv point), i.e. that segment of the curve is not well
% supported by data. ves bin edges are midpoints between consecutive
% ves_levels (requires ves_levels sorted ascending, as quantiles are).
% =========================================================================
function mask = slice_density_mask(T_work, epv_grid, ves_levels, MinDensity)
    epv_grid = epv_grid(:);
    n_pts    = numel(epv_grid);
    n_slices = numel(ves_levels);

    ves_levels_sorted = ves_levels(:)';
    ves_edges = zeros(1, n_slices + 1);
    ves_edges(1)   = -Inf;
    ves_edges(end) = Inf;
    for i = 2:n_slices
        ves_edges(i) = (ves_levels_sorted(i-1) + ves_levels_sorted(i)) / 2;
    end

    epv_edges = [-Inf; (epv_grid(1:end-1) + epv_grid(2:end)) / 2; Inf];

    counts = histcounts2(T_work.ves_swp, T_work.epv_swp, ves_edges, epv_edges); % n_slices x n_pts
    mask   = (counts < MinDensity)';  % n_pts x n_slices
    if size(mask, 2) ~= n_slices || size(mask, 1) ~= n_pts
        error('gam_slice_toolkit:densityMaskShape', ...
              'Internal error: density mask shape mismatch.');
    end
end


% =========================================================================
% Bootstrap CIs on the slice curves only (refit both GAMs B times on
% row-resampled data). Uses parfor if Parallel Computing Toolbox is
% available. NOTE: resampling is voxel-level -- see caveats in the
% top-level function headers regarding pseudoreplication across subjects.
% =========================================================================
function [ci_scatter, ci_retard] = bootstrap_slices(T_work, MaxSplits, NumTrees, ...
        LearnRate, epv_grid, ves_levels, NBootstrap, CIAlpha, Verbose)

    epv_grid = epv_grid(:);
    n_pts    = numel(epv_grid);
    n_slices = numel(ves_levels);
    n_fit    = height(T_work);

    T_slices = cell(n_slices, 1);
    for sl = 1:n_slices
        ves_fixed    = ves_levels(sl) * ones(n_pts, 1);
        T_slices{sl} = table(ves_fixed, epv_grid, ves_fixed .* epv_grid, ...
            'VariableNames', {'ves_swp', 'epv_swp', 'ves_epv_swp'});
    end

    boot_scatter = zeros(n_pts, NBootstrap, n_slices);
    boot_retard  = zeros(n_pts, NBootstrap, n_slices);

    use_par = ~isempty(ver('parallel'));
    if Verbose
        if use_par
            fprintf('[gam_slice_toolkit] Bootstrapping %d resamples (parallel)...\n', NBootstrap);
        else
            fprintf('[gam_slice_toolkit] Bootstrapping %d resamples...\n', NBootstrap);
        end
    end

    if use_par
        parfor b = 1:NBootstrap
            idx_b = randi(n_fit, n_fit, 1);
            T_b   = T_work(idx_b, :);
            mdl_s = fitrgam(T_b, 'scattering ~ ves_swp + epv_swp + ves_epv_swp', ...
                'NumTreesPerPredictor', NumTrees, 'MaxNumSplitsPerPredictor', MaxSplits, ...
                'InitialLearnRateForPredictors', LearnRate);
            mdl_r = fitrgam(T_b, 'retardance ~ ves_swp + epv_swp + ves_epv_swp', ...
                'NumTreesPerPredictor', NumTrees, 'MaxNumSplitsPerPredictor', MaxSplits, ...
                'InitialLearnRateForPredictors', LearnRate);
            tmp_s = zeros(n_pts, n_slices);
            tmp_r = zeros(n_pts, n_slices);
            for sl = 1:n_slices
                tmp_s(:, sl) = predict(mdl_s, T_slices{sl});
                tmp_r(:, sl) = predict(mdl_r, T_slices{sl});
            end
            boot_scatter(:, b, :) = tmp_s;
            boot_retard(:,  b, :) = tmp_r;
        end
    else
        for b = 1:NBootstrap
            if Verbose && mod(b, 10) == 0
                fprintf('[gam_slice_toolkit]   resample %d / %d\n', b, NBootstrap);
            end
            idx_b = randi(n_fit, n_fit, 1);
            T_b   = T_work(idx_b, :);
            mdl_s = fitrgam(T_b, 'scattering ~ ves_swp + epv_swp + ves_epv_swp', ...
                'NumTreesPerPredictor', NumTrees, 'MaxNumSplitsPerPredictor', MaxSplits, ...
                'InitialLearnRateForPredictors', LearnRate);
            mdl_r = fitrgam(T_b, 'retardance ~ ves_swp + epv_swp + ves_epv_swp', ...
                'NumTreesPerPredictor', NumTrees, 'MaxNumSplitsPerPredictor', MaxSplits, ...
                'InitialLearnRateForPredictors', LearnRate);
            for sl = 1:n_slices
                boot_scatter(:, b, sl) = predict(mdl_s, T_slices{sl});
                boot_retard(:,  b, sl) = predict(mdl_r, T_slices{sl});
            end
        end
    end

    alpha = CIAlpha;
    lo    = 100 *  alpha / 2;
    hi    = 100 * (1 - alpha / 2);

    ci_scatter = cell(n_slices, 1);
    ci_retard  = cell(n_slices, 1);
    for sl = 1:n_slices
        ci_scatter{sl} = prctile(boot_scatter(:, :, sl), [lo, hi], 2);
        ci_retard{sl}  = prctile(boot_retard(:,  :, sl), [lo, hi], 2);
    end
end


% =========================================================================
% Lighten an RGB triple toward white by `frac` (0 = unchanged, 1 = white).
% Used as a portable stand-in for line transparency: passing a 4-element
% [R G B alpha] vector directly into plot(...,'Color',...) is not reliably
% supported at call time (it errors under Octave and is not the officially
% documented MATLAB pattern -- true line-alpha requires setting it on the
% line handle post-creation). A lightened solid color avoids the issue
% entirely and renders identically in MATLAB and Octave.
% =========================================================================
function c = lighten_color(col, frac)
    c = col + (1 - col) * frac;
end


% =========================================================================
% Slice color palette: blue (low ves_swp) -> grey (mid) -> red (high),
% interpolated for n_slices != 3.
% =========================================================================
function colors = local_slice_colors(n_slices)
    base_colors = [0.20 0.45 0.80;
                   0.50 0.50 0.50;
                   0.80 0.20 0.20];
    if n_slices == 3
        colors = base_colors;
    else
        t      = linspace(0, 1, n_slices)';
        t_base = [0, 0.5, 1]';
        colors = interp1(t_base, base_colors, t, 'linear');
    end
end


% =========================================================================
% Draw one group's slice curves onto an existing (held) axes.
%   - Faded dotted line: full model curve across the whole grid.
%   - Solid/dashed emphasized line: NaN'd out wherever slice_density_mask
%     flagged insufficient data, so under-supported segments visually
%     drop out rather than looking equally trustworthy.
%   group_label: '' for a single-dataset plot; a group name (e.g.
%     'Severe') to disambiguate legend entries in a comparison plot.
% =========================================================================
function plot_group_slices(ax, epv_grid, slice_vals, mask, ci_cell, colors, ...
        ves_levels, slice_quantiles, has_ci, linestyle, group_label)

    n_slices = numel(ves_levels);
    for sl = 1:n_slices
        col = colors(sl, :);
        y   = slice_vals{sl};

        y_valid = y;
        if ~isempty(mask)
            y_valid(mask(:, sl)) = NaN;
        end

        if has_ci && ~isempty(ci_cell) && ~isempty(ci_cell{sl})
            fill_ci_ax(ax, epv_grid, ci_cell{sl}, col);
        end

        % Faded full curve (includes low-density / extrapolated segments)
        plot(ax, epv_grid, y, ':', 'Color', lighten_color(col, 0.55), ...
             'LineWidth', 1.0, 'HandleVisibility', 'off');

        if isempty(group_label)
            lbl = sprintf('Vessel SWP = %.3f (p%d)', ...
                          ves_levels(sl), round(slice_quantiles(sl) * 100));
        else
            lbl = sprintf('%s, Vessel SWP = %.3f (p%d)', group_label, ...
                          ves_levels(sl), round(slice_quantiles(sl) * 100));
        end

        plot(ax, epv_grid, y_valid, linestyle, ...
             'Color', col, 'LineWidth', 2.2, 'DisplayName', lbl);
    end
end


% =========================================================================
% Single-dataset slice plot: 1 figure, 2 panels (scattering, retardance).
% =========================================================================
function fig = plot_slices_single(epv_grid, ves_levels, slice_quantiles, ...
        slice_scatter, slice_retard, mask, ci_scatter, ci_retard, TitleStr, dirout)

    n_slices     = numel(ves_levels);
    slice_colors = local_slice_colors(n_slices);
    has_ci       = ~isempty(ci_scatter) && ~isempty(ci_scatter{1});

    fig = figure('Name', sprintf('Isosurface slices - %s', TitleStr), ...
                 'Position', [80, 100, 1200, 480], 'Color', 'w');

    ax1 = subplot(1, 2, 1); hold(ax1, 'on');
    plot_group_slices(ax1, epv_grid, slice_scatter, mask, ci_scatter, ...
                       slice_colors, ves_levels, slice_quantiles, has_ci, '-', '');
    hold(ax1, 'off');
    xlabel(ax1, 'EPVS SWP'); ylabel(ax1, 'Scattering coefficient');
    title(ax1, 'Scattering'); legend(ax1, 'Location', 'best', 'FontSize', 8);
    grid(ax1, 'on'); box(ax1, 'off');

    ax2 = subplot(1, 2, 2); hold(ax2, 'on');
    plot_group_slices(ax2, epv_grid, slice_retard, mask, ci_retard, ...
                       slice_colors, ves_levels, slice_quantiles, has_ci, '-', '');
    hold(ax2, 'off');
    xlabel(ax2, 'EPVS SWP'); ylabel(ax2, 'Retardance');
    title(ax2, 'Retardance'); legend(ax2, 'Location', 'best', 'FontSize', 8);
    grid(ax2, 'on'); box(ax2, 'off');

    plabel = strjoin(arrayfun(@(q) sprintf('p%d', round(q * 100)), ...
                     slice_quantiles, 'UniformOutput', false), '/');
    sub = 'Dotted = below MinDensity (model extrapolation)';
    if has_ci
        sub = [sub, sprintf('  |  bootstrap CI')];
    end
    sgtitle(fig, sprintf('%s -- GAM slice curves (fanning at vessel SWP %s)', TitleStr, plabel), ...
            'FontSize', 13, 'FontWeight', 'bold');

    if ~isempty(dirout)
        save_fig(fig, dirout, sprintf('IsoSlices_%s', strrep(TitleStr, ' ', '_')));
    end
end


% =========================================================================
% Two-group comparison slice plot: 1 figure, 2 panels. Group A solid,
% group B dashed; same color = same vessel-SWP slice level.
% =========================================================================
function fig = plot_slices_compare(epv_grid, ...
        ves_levels_A, slice_scatter_A, slice_retard_A, mask_A, ci_scatter_A, ci_retard_A, label_A, ...
        ves_levels_B, slice_scatter_B, slice_retard_B, mask_B, ci_scatter_B, ci_retard_B, label_B, ...
        slice_quantiles, TitleStr, dirout)

    n_slices     = numel(ves_levels_A);
    slice_colors = local_slice_colors(n_slices);
    has_ci       = ~isempty(ci_scatter_A) && ~isempty(ci_scatter_A{1});

    fig = figure('Name', sprintf('Isosurface comparison - %s', TitleStr), ...
                 'Position', [60, 80, 1300, 520], 'Color', 'w');

    ax1 = subplot(1, 2, 1); hold(ax1, 'on');
    plot_group_slices(ax1, epv_grid, slice_scatter_A, mask_A, ci_scatter_A, ...
                       slice_colors, ves_levels_A, slice_quantiles, has_ci, '-',  label_A);
    plot_group_slices(ax1, epv_grid, slice_scatter_B, mask_B, ci_scatter_B, ...
                       slice_colors, ves_levels_B, slice_quantiles, has_ci, '--', label_B);
    hold(ax1, 'off');
    xlabel(ax1, 'EPVS SWP'); ylabel(ax1, 'Scattering coefficient');
    title(ax1, 'Scattering'); legend(ax1, 'Location', 'best', 'FontSize', 7);
    grid(ax1, 'on'); box(ax1, 'off');

    ax2 = subplot(1, 2, 2); hold(ax2, 'on');
    plot_group_slices(ax2, epv_grid, slice_retard_A, mask_A, ci_retard_A, ...
                       slice_colors, ves_levels_A, slice_quantiles, has_ci, '-',  label_A);
    plot_group_slices(ax2, epv_grid, slice_retard_B, mask_B, ci_retard_B, ...
                       slice_colors, ves_levels_B, slice_quantiles, has_ci, '--', label_B);
    hold(ax2, 'off');
    xlabel(ax2, 'EPVS SWP'); ylabel(ax2, 'Retardance');
    title(ax2, 'Retardance'); legend(ax2, 'Location', 'best', 'FontSize', 7);
    grid(ax2, 'on'); box(ax2, 'off');

    plabel = strjoin(arrayfun(@(q) sprintf('p%d', round(q * 100)), ...
                     slice_quantiles, 'UniformOutput', false), '/');
    sgtitle(fig, sprintf('%s: %s (solid) vs %s (dashed) -- fanning at vessel SWP %s', ...
            TitleStr, label_A, label_B, plabel), ...
            'FontSize', 13, 'FontWeight', 'bold');

    if ~isempty(dirout)
        save_fig(fig, dirout, sprintf('IsoCompare_%s', strrep(TitleStr, ' ', '_')));
    end
end


% =========================================================================
% Shaded CI band. Caller must hold(ax,'on')/hold(ax,'off') around calls.
% =========================================================================
function fill_ci_ax(ax, x, ci, col)
    x = x(:);
    x_patch = [x;       flipud(x)];
    y_patch = [ci(:,1); flipud(ci(:,2))];
    fill(ax, x_patch, y_patch, col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    plot(ax, x, ci(:,1), '-', 'Color', lighten_color(col, 0.55), 'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(ax, x, ci(:,2), '-', 'Color', lighten_color(col, 0.55), 'LineWidth', 0.8, 'HandleVisibility', 'off');
end


% =========================================================================
% Save figure as .fig and .png.
% =========================================================================
function save_fig(fig, dirout, fname)
    if ~isfolder(dirout)
        mkdir(dirout);
    end
    fname = regexprep(fname, '_+$', '');
    saveas(fig, fullfile(dirout, [fname, '.fig']));
    saveas(fig, fullfile(dirout, [fname, '.png']));
end
