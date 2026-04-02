function results = fit_swp_psoct_gam(T, varargin)
% FIT_SWP_PSOCT_GAM  Fit GAMs for scattering coefficient and retardance as
%                    functions of vessel SWP and EPVS SWP, with a full 2D
%                    interaction surface using MATLAB's fitrgam.
%
%   fitrgam uses gradient boosted regression trees - NOT splines.
%   An explicit product interaction term (SWP * EPVS_SWP) is included so
%   that the effect of each proximity measure can vary continuously with
%   the level of the other.  Predictions are evaluated on a full 2D
%   meshgrid so the complete joint response surface is captured.
%
% USAGE
%   results = fit_swp_psoct_gam(T, 'TitleStr', 'MyData', 'dirout', '/path')
%   results = fit_swp_psoct_gam(T, 'TitleStr', 'MyData', 'dirout', '/path', ...
%                               'NumTrees', 200, 'PlotResults', false)
%
% REQUIRED INPUT
%   T : table with columns:
%         'scattering'  - scattering coefficient (voxel-level)
%         'retardance'  - optical retardance      (voxel-level)
%         'SWP'         - vessel size-weighted proximity
%         'EPVS_SWP'    - EPVS size-weighted proximity
%
% REQUIRED NAME-VALUE PAIRS
%   'TitleStr'      String used in figure title and output filename
%   'dirout'        Output directory for saved figures
%
% OPTIONAL NAME-VALUE PAIRS
%   'NumTrees'      Trees per predictor (controls fit complexity)     (default: 100)
%   'MaxSplits'     Max splits per tree (controls term smoothness)    (default: 10)
%   'LearnRate'     Shrinkage / learning rate per tree                (default: 0.1)
%   'NumGridPts'    Points along each axis of the 2D prediction grid  (default: 100)
%   'MaxSample'     Max voxels used for fitting (subsampling)         (default: 1e5)
%   'PlotResults'   Show smooth effect plots                          (default: true)
%   'Verbose'       Print model summaries to command window           (default: true)
%   'BootstrapCI'   Compute bootstrap CIs on marginal slice curves   (default: true)
%   'NBootstrap'    Number of bootstrap resamples                     (default: 100)
%   'CIAlpha'       Alpha level for CI (0.05 = 95% CI)               (default: 0.05)
%   'MinDensity'    Min data count per cell to show surface (masking) (default: 5)
%
% OUTPUT  results struct with fields:
%   .gam_scatter         - fitted GAM object for scattering
%   .gam_retard          - fitted GAM object for retardance
%   .swp_grid            - SWP axis vector (length NumGridPts)
%   .epvs_grid           - EPVS_SWP axis vector (length NumGridPts)
%   .SWP_mesh            - NumGridPts x NumGridPts meshgrid of SWP values
%   .EPVS_mesh           - NumGridPts x NumGridPts meshgrid of EPVS values
%   .pred_scatter_2d     - NumGridPts x NumGridPts predicted scattering surface
%   .pred_retard_2d      - NumGridPts x NumGridPts predicted retardance surface
%   .density_mask        - logical mask: true where data are too sparse
%   .slice_epvs_levels   - EPVS quantile values used for marginal slices
%   .slice_scatter_swp   - cell array of scattering curves at each EPVS slice
%   .slice_retard_swp    - cell array of retardance curves at each EPVS slice
%   .ci_slice_scatter    - cell array of Nx2 bootstrap CIs per EPVS slice
%   .ci_slice_retard     - cell array of Nx2 bootstrap CIs per EPVS slice
%   .T_fit               - (possibly subsampled) table used for fitting
%
% NOTES
%   * The interaction is captured via an engineered product feature:
%       SWP_EPVS = SWP .* EPVS_SWP
%     The formula becomes:  outcome ~ SWP + EPVS_SWP + SWP_EPVS
%   * Prediction tables always set SWP_EPVS = SWP .* EPVS_SWP so the
%     interaction column is internally consistent with the main effects.
%   * The 2D surface evaluates predictions across the full joint space.
%     Cells with fewer than MinDensity data points are masked (set to NaN)
%     because predictions there are extrapolations.
%   * Marginal slice plots show the SWP effect at low / median / high
%     EPVS_SWP (10th, 50th, 90th percentiles) so interaction fanning is
%     immediately visible.
%   * Bootstrap CIs are computed on the marginal slices (not the full 2D
%     surface, which would be prohibitively expensive).

    % =====================================================================
    % 1. Parse and validate inputs
    % =====================================================================
    p = inputParser();
    p.FunctionName = 'fit_swp_psoct_gam';

    addRequired(p,  'T',             @istable);
    addParameter(p, 'NumTrees',      100,   @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'MaxSplits',     10,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'LearnRate',     0.1,   @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    addParameter(p, 'NumGridPts',    100,   @(x) isnumeric(x) && isscalar(x) && x > 1);
    addParameter(p, 'MaxSample',     1e5,   @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'PlotResults',   true,  @islogical);
    addParameter(p, 'Verbose',       true,  @islogical);
    addParameter(p, 'BootstrapCI',   true,  @islogical);
    addParameter(p, 'NBootstrap',    100,   @(x) isnumeric(x) && isscalar(x) && x >= 10);
    addParameter(p, 'CIAlpha',       0.05,  @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'MinDensity',    5,     @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'TitleStr',      '',    @ischar);
    addParameter(p, 'dirout',        pwd,   @ischar);

    parse(p, T, varargin{:});
    opts = p.Results;

    % Check required columns exist
    required_cols = {'scattering', 'retardance', 'SWP', 'EPVS_SWP'};
    missing = required_cols(~ismember(required_cols, T.Properties.VariableNames));
    if ~isempty(missing)
        error('fit_swp_psoct_gam:missingColumns', ...
              'Table is missing required column(s): %s', strjoin(missing, ', '));
    end

    % Check columns are numeric and not all-NaN
    for c = required_cols
        col = T.(c{1});
        if ~isnumeric(col)
            error('fit_swp_psoct_gam:nonNumeric', 'Column "%s" must be numeric.', c{1});
        end
        if all(isnan(col))
            error('fit_swp_psoct_gam:allNaN', 'Column "%s" is all NaN.', c{1});
        end
    end

    % Check toolbox availability
    if ~license('test', 'Statistics_Toolbox')
        error('fit_swp_psoct_gam:noToolbox', ...
              'Statistics and Machine Learning Toolbox is required for fitrgam.');
    end

    % =====================================================================
    % 2. Remove rows with NaN or Inf in any required column
    % =====================================================================
    T_work   = T(:, required_cols);
    arr      = table2array(T_work);
    bad_rows = any(isnan(arr) | isinf(arr), 2);
    n_removed = sum(bad_rows);
    T_work(bad_rows, :) = [];

    if opts.Verbose && n_removed > 0
        fprintf('[fit_swp_psoct_gam] Removed %d rows containing NaN or Inf.\n', n_removed);
    end

    if height(T_work) < 10
        error('fit_swp_psoct_gam:insufficientData', ...
              'Too few valid rows (%d) after removing NaN/Inf.', height(T_work));
    end

    % =====================================================================
    % 3. Subsample if necessary (stratified by EPVS_SWP quantile)
    % =====================================================================
    N = height(T_work);
    if N > opts.MaxSample
        if opts.Verbose
            fprintf('[fit_swp_psoct_gam] Subsampling %d -> %d voxels (stratified by EPVS SWP).\n', ...
                    N, round(opts.MaxSample));
        end
        n_strata   = 10;
        edges      = quantile(T_work.EPVS_SWP, linspace(0, 1, n_strata + 1));
        edges(1)   = -Inf;
        edges(end) = Inf;
        n_per      = floor(opts.MaxSample / n_strata);
        keep_idx   = [];
        for s = 1:n_strata
            in_stratum = find(T_work.EPVS_SWP > edges(s) & T_work.EPVS_SWP <= edges(s+1));
            if isempty(in_stratum), continue; end
            n_draw   = min(n_per, numel(in_stratum));
            selected = in_stratum(randperm(numel(in_stratum), n_draw));
            keep_idx = [keep_idx; selected]; %#ok<AGROW>
        end
        T_work = T_work(keep_idx, :);
    end

    % =====================================================================
    % 4. Engineer interaction feature
    %    SWP_EPVS = SWP .* EPVS_SWP
    %    This product term allows the effect of SWP to vary continuously
    %    with EPVS_SWP level (and vice versa) within the additive framework
    %    of fitrgam.  The prediction table must always set:
    %        SWP_EPVS = SWP_column .* EPVS_SWP_column
    %    to remain internally consistent.
    % =====================================================================
    T_work.SWP_EPVS = T_work.SWP .* T_work.EPVS_SWP;

    formula_scatter = 'scattering ~ SWP + EPVS_SWP + SWP_EPVS';
    formula_retard  = 'retardance ~ SWP + EPVS_SWP + SWP_EPVS';

    % =====================================================================
    % 5. Fit GAMs
    % =====================================================================
    if opts.Verbose
        fprintf('[fit_swp_psoct_gam] Fitting GAM for scattering (N=%d voxels)...\n', height(T_work));
    end

    gam_scatter = fitrgam(T_work, formula_scatter, ...
        'NumTreesPerPredictor',          opts.NumTrees,  ...
        'MaxNumSplitsPerPredictor',      opts.MaxSplits, ...
        'InitialLearnRateForPredictors', opts.LearnRate);

    if opts.Verbose
        fprintf('[fit_swp_psoct_gam] Fitting GAM for retardance...\n');
    end

    gam_retard = fitrgam(T_work, formula_retard, ...
        'NumTreesPerPredictor',          opts.NumTrees,  ...
        'MaxNumSplitsPerPredictor',      opts.MaxSplits, ...
        'InitialLearnRateForPredictors', opts.LearnRate);

    if opts.Verbose
        yhat_s = resubPredict(gam_scatter);
        yhat_r = resubPredict(gam_retard);
        r2_s   = compute_r2(T_work.scattering, yhat_s);
        r2_r   = compute_r2(T_work.retardance,  yhat_r);
        fprintf('[fit_swp_psoct_gam] Scattering GAM resubstitution R^2 = %.4f\n', r2_s);
        fprintf('[fit_swp_psoct_gam] Retardance  GAM resubstitution R^2 = %.4f\n', r2_r);
        fprintf('  (Resubstitution R^2 is optimistic; use cross-validation for generalization.)\n');
    end

    % =====================================================================
    % 6. Build full 2D prediction grid
    %
    %    meshgrid produces two n_pts x n_pts matrices:
    %      SWP_mesh(i,j)  = swp_grid(j)    [SWP varies along columns]
    %      EPVS_mesh(i,j) = epvs_grid(i)   [EPVS varies along rows]
    %    Predictions are reshaped back to n_pts x n_pts after predict().
    %
    %    The interaction column is set to SWP_mesh(:) .* EPVS_mesh(:) so
    %    it is always the exact product of the two main effect columns.
    % =====================================================================
    n_pts     = opts.NumGridPts;
    swp_grid  = linspace(min(T_work.SWP),      max(T_work.SWP),      n_pts)';
    epvs_grid = linspace(min(T_work.EPVS_SWP), max(T_work.EPVS_SWP), n_pts)';

    [SWP_mesh, EPVS_mesh] = meshgrid(swp_grid, epvs_grid);

    T_pred_2d = table(SWP_mesh(:), EPVS_mesh(:), SWP_mesh(:) .* EPVS_mesh(:), ...
        'VariableNames', {'SWP', 'EPVS_SWP', 'SWP_EPVS'});

    if opts.Verbose
        fprintf('[fit_swp_psoct_gam] Predicting over %d x %d = %d grid points...\n', ...
                n_pts, n_pts, n_pts^2);
    end

    pred_scatter_2d = reshape(predict(gam_scatter, T_pred_2d), n_pts, n_pts);
    pred_retard_2d  = reshape(predict(gam_retard,  T_pred_2d), n_pts, n_pts);

    % =====================================================================
    % 7. Data density mask
    %    Count observed data points in each grid cell and mask cells that
    %    have fewer than MinDensity points (extrapolation regions).
    % =====================================================================
    % Build bin edges that exactly bracket the grid axes
    swp_edges  = [-Inf; (swp_grid(1:end-1)  + swp_grid(2:end))  / 2; Inf];
    epvs_edges = [-Inf; (epvs_grid(1:end-1) + epvs_grid(2:end)) / 2; Inf];

    % histcounts2: rows = EPVS bins, cols = SWP bins (matches meshgrid orientation)
    density = histcounts2(T_work.EPVS_SWP, T_work.SWP, epvs_edges, swp_edges);

    density_mask = density < opts.MinDensity;   % true where too sparse

    pred_scatter_masked = pred_scatter_2d;
    pred_retard_masked  = pred_retard_2d;
    pred_scatter_masked(density_mask) = NaN;
    pred_retard_masked(density_mask)  = NaN;

    if opts.Verbose
        pct_masked = 100 * mean(density_mask(:));
        fprintf('[fit_swp_psoct_gam] %.1f%% of grid cells masked (< %d data points).\n', ...
                pct_masked, opts.MinDensity);
    end

    % =====================================================================
    % 8. Marginal slice predictions
    %    Evaluate the SWP effect at three EPVS_SWP levels:
    %    10th, 50th, and 90th percentiles.
    %    This shows whether the SWP curve fans out (interaction) or stays
    %    parallel (no interaction) across EPVS levels.
    % =====================================================================
    slice_quantiles  = [0.10, 0.50, 0.90];
    epvs_levels      = quantile(T_work.EPVS_SWP, slice_quantiles);
    n_slices         = numel(epvs_levels);

    slice_scatter_swp = cell(n_slices, 1);
    slice_retard_swp  = cell(n_slices, 1);

    for sl = 1:n_slices
        epvs_fixed = epvs_levels(sl) * ones(n_pts, 1);
        T_sl = table(swp_grid, epvs_fixed, swp_grid .* epvs_fixed, ...
            'VariableNames', {'SWP', 'EPVS_SWP', 'SWP_EPVS'});
        slice_scatter_swp{sl} = predict(gam_scatter, T_sl);
        slice_retard_swp{sl}  = predict(gam_retard,  T_sl);
    end

    % =====================================================================
    % 9. Bootstrap CIs on marginal slices
    %    Resamples rows with replacement, refits both GAMs, collects
    %    predictions on the fixed SWP grid at each EPVS slice level.
    %    Uses parfor if Parallel Computing Toolbox is available.
    % =====================================================================
    ci_slice_scatter = cell(n_slices, 1);
    ci_slice_retard  = cell(n_slices, 1);

    if opts.BootstrapCI
        if opts.Verbose
            fprintf('[fit_swp_psoct_gam] Running %d bootstrap resamples for slice CIs...\n', ...
                    opts.NBootstrap);
        end

        B       = opts.NBootstrap;
        n_fit   = height(T_work);
        b_trees  = opts.NumTrees;
        b_splits = opts.MaxSplits;
        b_lr     = opts.LearnRate;

        % boot arrays: rows = grid points, cols = bootstrap resamples,
        % pages = EPVS slices
        boot_scatter = zeros(n_pts, B, n_slices);
        boot_retard  = zeros(n_pts, B, n_slices);

        use_par = ~isempty(ver('parallel'));
        if use_par && opts.Verbose
            fprintf('[fit_swp_psoct_gam] Parallel Computing Toolbox detected - using parfor.\n');
        end

        % Pre-build slice prediction tables for speed inside loop
        T_slices = cell(n_slices, 1);
        for sl = 1:n_slices
            epvs_fixed = epvs_levels(sl) * ones(n_pts, 1);
            T_slices{sl} = table(swp_grid, epvs_fixed, swp_grid .* epvs_fixed, ...
                'VariableNames', {'SWP', 'EPVS_SWP', 'SWP_EPVS'});
        end

        if use_par
            parfor b = 1:B
                idx_b = randi(n_fit, n_fit, 1);
                T_b   = T_work(idx_b, :);

                mdl_s = fitrgam(T_b, 'scattering ~ SWP + EPVS_SWP + SWP_EPVS', ...
                    'NumTreesPerPredictor',          b_trees,  ...
                    'MaxNumSplitsPerPredictor',      b_splits, ...
                    'InitialLearnRateForPredictors', b_lr);
                mdl_r = fitrgam(T_b, 'retardance ~ SWP + EPVS_SWP + SWP_EPVS',  ...
                    'NumTreesPerPredictor',          b_trees,  ...
                    'MaxNumSplitsPerPredictor',      b_splits, ...
                    'InitialLearnRateForPredictors', b_lr);

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
            for b = 1:B
                if opts.Verbose && mod(b, 10) == 0
                    fprintf('[fit_swp_psoct_gam]   Bootstrap resample %d / %d\n', b, B);
                end
                idx_b = randi(n_fit, n_fit, 1);
                T_b   = T_work(idx_b, :);

                mdl_s = fitrgam(T_b, 'scattering ~ SWP + EPVS_SWP + SWP_EPVS', ...
                    'NumTreesPerPredictor',          b_trees,  ...
                    'MaxNumSplitsPerPredictor',      b_splits, ...
                    'InitialLearnRateForPredictors', b_lr);
                mdl_r = fitrgam(T_b, 'retardance ~ SWP + EPVS_SWP + SWP_EPVS',  ...
                    'NumTreesPerPredictor',          b_trees,  ...
                    'MaxNumSplitsPerPredictor',      b_splits, ...
                    'InitialLearnRateForPredictors', b_lr);

                for sl = 1:n_slices
                    boot_scatter(:, b, sl) = predict(mdl_s, T_slices{sl});
                    boot_retard(:,  b, sl) = predict(mdl_r, T_slices{sl});
                end
            end
        end

        % Percentile-based CIs per slice
        alpha = opts.CIAlpha;
        lo    = 100 *  alpha / 2;
        hi    = 100 * (1 - alpha / 2);

        for sl = 1:n_slices
            % Use row vector [lo, hi] for prctile — column vector [lo; hi]
            % gives ambiguous output shape in some MATLAB versions.
            ci_slice_scatter{sl} = prctile(boot_scatter(:, :, sl), [lo, hi], 2);  % Nx2
            ci_slice_retard{sl}  = prctile(boot_retard(:,  :, sl), [lo, hi], 2);
        end

        if opts.Verbose
            fprintf('[fit_swp_psoct_gam] Bootstrap CIs computed (alpha = %.2f).\n', alpha);
        end
    end

    % =====================================================================
    % 10. Plot
    % =====================================================================
    if opts.PlotResults
        plot_2d_effects(swp_grid, epvs_grid, SWP_mesh, EPVS_mesh,       ...
                        pred_scatter_masked,  pred_retard_masked,         ...
                        pred_scatter_2d,      pred_retard_2d,             ...
                        slice_scatter_swp,    slice_retard_swp,           ...
                        ci_slice_scatter,     ci_slice_retard,            ...
                        epvs_levels,          slice_quantiles,            ...
                        opts.CIAlpha,         T_work,                    ...
                        opts.TitleStr,        opts.dirout);
    end

    % =====================================================================
    % 11. Package results
    % =====================================================================
    results.gam_scatter        = gam_scatter;
    results.gam_retard         = gam_retard;
    results.swp_grid           = swp_grid;
    results.epvs_grid          = epvs_grid;
    results.SWP_mesh           = SWP_mesh;
    results.EPVS_mesh          = EPVS_mesh;
    results.pred_scatter_2d    = pred_scatter_2d;       % unmasked
    results.pred_retard_2d     = pred_retard_2d;        % unmasked
    results.pred_scatter_masked = pred_scatter_masked;  % NaN in sparse cells
    results.pred_retard_masked  = pred_retard_masked;
    results.density_mask       = density_mask;
    results.density            = density;
    results.slice_epvs_levels  = epvs_levels;
    results.slice_scatter_swp  = slice_scatter_swp;
    results.slice_retard_swp   = slice_retard_swp;
    results.ci_slice_scatter   = ci_slice_scatter;
    results.ci_slice_retard    = ci_slice_retard;
    results.T_fit              = T_work;

    if opts.Verbose
        fprintf('[fit_swp_psoct_gam] Done.\n');
    end
end


% =========================================================================
% LOCAL FUNCTION: R-squared
% =========================================================================
function r2 = compute_r2(y, yhat)
    ss_res = sum((y - yhat).^2);
    ss_tot = sum((y - mean(y)).^2);
    if ss_tot == 0
        r2 = NaN;
    else
        r2 = 1 - ss_res / ss_tot;
    end
end


% =========================================================================
% LOCAL FUNCTION: plot 2D surface effects + marginal slice curves
%
%   Layout (2 outcomes x 2 plot types = 4 panels per figure):
%
%   Figure 1 — Scattering
%     Panel 1 (left):  heatmap of the 2D surface (SWP x EPVS)
%     Panel 2 (right): marginal slice curves at 3 EPVS levels
%
%   Figure 2 — Retardance
%     Panel 1 (left):  heatmap of the 2D surface
%     Panel 2 (right): marginal slice curves
% =========================================================================
function plot_2d_effects(swp_grid, epvs_grid, SWP_mesh, EPVS_mesh,   ...
                         pred_scatter_masked, pred_retard_masked,       ...
                         pred_scatter_2d,     pred_retard_2d,           ...
                         slice_scatter_swp,   slice_retard_swp,         ...
                         ci_slice_scatter,    ci_slice_retard,          ...
                         epvs_levels,         slice_quantiles,          ...
                         ci_alpha,            T_work,                   ...
                         tstr,                dirout)

    % has_ci must check the *content* of the first cell, not the cell array
    % itself.  When BootstrapCI=false the cell array is non-empty but every
    % cell contains [].  ~isempty(ci_slice_scatter) would incorrectly return
    % true in that case and pass [] into fill_ci_ax, causing an index error.
    has_ci  = ~isempty(ci_slice_scatter) && ~isempty(ci_slice_scatter{1});
    ci_pct  = 100 * (1 - ci_alpha);
    n_slices = numel(epvs_levels);

    % Colours for the three EPVS slices: low=blue, mid=grey, high=red
    slice_colors = [0.20 0.45 0.80;   % low EPVS (10th pct)
                    0.50 0.50 0.50;   % mid EPVS (50th pct)
                    0.80 0.20 0.20];  % high EPVS (90th pct)

    % Rug subsample for data density overlay on heatmaps
    n_rug   = min(3000, height(T_work));
    rug_idx = randperm(height(T_work), n_rug);
    rug_swp  = T_work.SWP(rug_idx);
    rug_epvs = T_work.EPVS_SWP(rug_idx);

    % ------------------------------------------------------------------
    % Shared colour limits: use the unmasked surface so masked NaNs do
    % not distort the colour axis.
    % ------------------------------------------------------------------
    clim_scatter = [min(pred_scatter_2d(:)), max(pred_scatter_2d(:))];
    clim_retard  = [min(pred_retard_2d(:)),  max(pred_retard_2d(:))];

    % ==================================================================
    % Figure 1: Scattering
    % ==================================================================
    fig1 = figure('Name',     'GAM 2D Surface - Scattering', ...
                  'Position', [60, 120, 1300, 520],           ...
                  'Color',    'w');

    % --- Panel 1: 2D heatmap ---
    ax1 = subplot(1, 2, 1);
    imagesc(ax1, swp_grid, epvs_grid, pred_scatter_masked);
    set(ax1, 'YDir', 'normal');
    set(ax1, 'CLim', clim_scatter);   % clim() syntax requires R2022a+; set() is cross-version safe
    colormap(ax1, parula);
    cb1 = colorbar(ax1);
    cb1.Label.String = 'Predicted scattering';
    hold(ax1, 'on');
    % Overlay rug to show data density
    scatter(ax1, rug_swp, rug_epvs, 2, 'w', 'filled', ...
            'MarkerFaceAlpha', 0.08);
    % Mark the three slice levels as horizontal dashed lines
    for sl = 1:n_slices
        yline(ax1, epvs_levels(sl), '--', ...
              'Color', slice_colors(sl,:), 'LineWidth', 1.2, ...
              'Label', sprintf('p%d', round(slice_quantiles(sl)*100)), ...
              'LabelHorizontalAlignment', 'left', ...
              'FontSize', 7);
    end
    hold(ax1, 'off');
    xlabel(ax1, 'Vessel SWP');
    ylabel(ax1, 'EPVS SWP');
    title(ax1, 'Scattering: joint response surface');
    subtitle(ax1, 'NaN (grey) = sparse data region', ...
             'FontSize', 8, 'Color', [.5 .5 .5]);
    grid(ax1, 'off'); box(ax1, 'on');

    % --- Panel 2: Marginal slice curves ---
    ax2 = subplot(1, 2, 2);
    hold(ax2, 'on');
    for sl = 1:n_slices
        col = slice_colors(sl, :);
        lbl = sprintf('EPVS SWP = %.3f (p%d)', ...
                      epvs_levels(sl), round(slice_quantiles(sl)*100));
        if has_ci
            fill_ci_ax(ax2, swp_grid, ci_slice_scatter{sl}, col);
        end
        plot(ax2, swp_grid, slice_scatter_swp{sl}, '-', ...
             'Color', col, 'LineWidth', 2.2, 'DisplayName', lbl);
    end
    hold(ax2, 'off');
    legend(ax2, 'Location', 'best', 'FontSize', 8);
    xlabel(ax2, 'Vessel SWP');
    ylabel(ax2, 'Predicted scattering');
    title(ax2, 'Scattering ~ f(SWP) at 3 EPVS levels');
    sub2 = 'Fanning = interaction present';
    if has_ci
        sub2 = sprintf('%s  |  %g%% CI', sub2, ci_pct);
    end
    subtitle(ax2, sub2, 'FontSize', 8, 'Color', [.5 .5 .5]);
    grid(ax2, 'on'); box(ax2, 'off');

    sgtitle(fig1, sprintf('GAM Scattering — %s', tstr), ...
            'FontSize', 13, 'FontWeight', 'bold');

    save_figure(fig1, dirout, sprintf('GAM_Scatter_2D_%s', strrep(tstr,' ','_')));

    % ==================================================================
    % Figure 2: Retardance
    % ==================================================================
    fig2 = figure('Name',     'GAM 2D Surface - Retardance', ...
                  'Position', [100, 80, 1300, 520],           ...
                  'Color',    'w');

    % --- Panel 1: 2D heatmap ---
    ax3 = subplot(1, 2, 1);
    imagesc(ax3, swp_grid, epvs_grid, pred_retard_masked);
    set(ax3, 'YDir', 'normal');
    set(ax3, 'CLim', clim_retard);    % clim() syntax requires R2022a+; set() is cross-version safe
    colormap(ax3, parula);
    cb3 = colorbar(ax3);
    cb3.Label.String = 'Predicted retardance';
    hold(ax3, 'on');
    scatter(ax3, rug_swp, rug_epvs, 2, 'w', 'filled', ...
            'MarkerFaceAlpha', 0.08);
    for sl = 1:n_slices
        yline(ax3, epvs_levels(sl), '--', ...
              'Color', slice_colors(sl,:), 'LineWidth', 1.2, ...
              'Label', sprintf('p%d', round(slice_quantiles(sl)*100)), ...
              'LabelHorizontalAlignment', 'left', ...
              'FontSize', 7);
    end
    hold(ax3, 'off');
    xlabel(ax3, 'Vessel SWP');
    ylabel(ax3, 'EPVS SWP');
    title(ax3, 'Retardance: joint response surface');
    subtitle(ax3, 'NaN (grey) = sparse data region', ...
             'FontSize', 8, 'Color', [.5 .5 .5]);
    grid(ax3, 'off'); box(ax3, 'on');

    % --- Panel 2: Marginal slice curves ---
    ax4 = subplot(1, 2, 2);
    hold(ax4, 'on');
    for sl = 1:n_slices
        col = slice_colors(sl, :);
        lbl = sprintf('EPVS SWP = %.3f (p%d)', ...
                      epvs_levels(sl), round(slice_quantiles(sl)*100));
        if has_ci
            fill_ci_ax(ax4, swp_grid, ci_slice_retard{sl}, col);
        end
        plot(ax4, swp_grid, slice_retard_swp{sl}, '-', ...
             'Color', col, 'LineWidth', 2.2, 'DisplayName', lbl);
    end
    hold(ax4, 'off');
    legend(ax4, 'Location', 'best', 'FontSize', 8);
    xlabel(ax4, 'Vessel SWP');
    ylabel(ax4, 'Predicted retardance');
    title(ax4, 'Retardance ~ f(SWP) at 3 EPVS levels');
    sub4 = 'Fanning = interaction present';
    if has_ci
        sub4 = sprintf('%s  |  %g%% CI', sub4, ci_pct);
    end
    subtitle(ax4, sub4, 'FontSize', 8, 'Color', [.5 .5 .5]);
    grid(ax4, 'on'); box(ax4, 'off');

    sgtitle(fig2, sprintf('GAM Retardance — %s', tstr), ...
            'FontSize', 13, 'FontWeight', 'bold');

    save_figure(fig2, dirout, sprintf('GAM_Retard_2D_%s', strrep(tstr,' ','_')));

    % ==================================================================
    % Figure 3: 3D surface plots (optional visual for presentations)
    % ==================================================================
    fig3 = figure('Name',     'GAM 3D Surfaces', ...
                  'Position', [120, 60, 1300, 520], ...
                  'Color',    'w');

    ax5 = subplot(1, 2, 1);
    surf(ax5, SWP_mesh, EPVS_mesh, pred_scatter_masked, ...
         'EdgeAlpha', 0.0, 'FaceAlpha', 0.9);
    colormap(ax5, parula);
    colorbar(ax5);
    xlabel(ax5, 'Vessel SWP');
    ylabel(ax5, 'EPVS SWP');
    zlabel(ax5, 'Scattering');
    title(ax5, 'Scattering: 3D response surface');
    shading(ax5, 'flat');    % 'interp' would interpolate across NaN-masked cells, filling the gaps
    view(ax5, -35, 30);
    grid(ax5, 'on');

    ax6 = subplot(1, 2, 2);
    surf(ax6, SWP_mesh, EPVS_mesh, pred_retard_masked, ...
         'EdgeAlpha', 0.0, 'FaceAlpha', 0.9);
    colormap(ax6, parula);
    colorbar(ax6);
    xlabel(ax6, 'Vessel SWP');
    ylabel(ax6, 'EPVS SWP');
    zlabel(ax6, 'Retardance');
    title(ax6, 'Retardance: 3D response surface');
    shading(ax6, 'flat');    % 'interp' would interpolate across NaN-masked cells, filling the gaps
    view(ax6, -35, 30);
    grid(ax6, 'on');

    sgtitle(fig3, sprintf('GAM 3D Surfaces — %s', tstr), ...
            'FontSize', 13, 'FontWeight', 'bold');

    save_figure(fig3, dirout, sprintf('GAM_3D_%s', strrep(tstr,' ','_')));
end


% =========================================================================
% LOCAL FUNCTION: shaded CI band on a specified axes
%   Caller is responsible for hold(ax,'on') before calling this function
%   and hold(ax,'off') after the plotting loop completes.  This function
%   deliberately does NOT call hold() itself so that an error mid-loop
%   cannot leave hold permanently enabled on the axes.
% =========================================================================
function fill_ci_ax(ax, x, ci, col)
    x_patch = [x;       flipud(x)];
    y_patch = [ci(:,1); flipud(ci(:,2))];
    fill(ax, x_patch, y_patch, col, ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(ax, x, ci(:,1), '-', 'Color', [col, 0.30], 'LineWidth', 0.8);
    plot(ax, x, ci(:,2), '-', 'Color', [col, 0.30], 'LineWidth', 0.8);
end


% =========================================================================
% LOCAL FUNCTION: save figure as .fig and .png
% =========================================================================
function save_figure(fig, dirout, fname)
    if ~isfolder(dirout)
        mkdir(dirout);
    end
    saveas(fig, fullfile(dirout, [fname, '.fig']));
    saveas(fig, fullfile(dirout, [fname, '.png']));
end