function results = fit_swp_histo_gam(T, varargin)
% FIT_SWP_HISTO_GAM  Fit Genrralized Additive Model (GAM) for
%   histopathology stains as
%   functions of ves_swp and epv_swp, with a
%   full 2D interaction surface using MATLAB's fitrgam.
%
%   fitrgam uses gradient boosted regression trees - NOT splines.
%   An explicit product interaction term (ves_swp * epv_swp) is included so
%   that the effect of each proximity measure can vary continuously with
%   the level of the other.  Predictions are evaluated on a full 2D
%   meshgrid so the complete joint response surface is captured.
%
%   Marginal slice curves sweep EPV SWP at three vessel SWP quantile
%   levels (10th, 50th, 90th percentile). This allows for 2D plots of the
%   stain vs. EPVS-SWP at three vessel-SWP percentiles
%
% USAGE
%   results = fit_swp_histo_gam(T, 'TitleStr', 'MyData', 'dirout', '/path')
%   results = fit_swp_histo_gam(T, 'TitleStr', 'MyData', 'dirout', '/path', ...
%                               'NumTrees', 200, 'PlotResults', false)
%
% REQUIRED INPUT
%   T : table with columns:
%         'lhe'  - luxol-fast blue H&E (myelin)
%         'gfap'  - glial fibrillary acidic protein (activated glia)
%         'cd68' - cluster of differentiation 68 (monocytes)
%         'ves_swp'     - vessel size-weighted proximity
%         'epv_swp'     - EPVS size-weighted proximity
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
%   'Verbose'       Print model summaries to command window           (default: true)
%   'BootstrapCI'   Compute bootstrap CIs on marginal slice curves    (default: true)
%   'NBootstrap'    Number of bootstrap resamples                     (default: 100)
%   'CIAlpha'       Alpha level for CI (0.05 = 95% CI)               (default: 0.05)
%   'MinDensity'    Min data count per cell to show surface (masking) (default: 5)
%
% OUTPUT  results struct with fields:
%   .gam                  - fitted GAM object for stain
%   .ves_swp_grid         - ves_swp axis vector (length NumGridPts)
%   .epv_swp_grid         - epv_swp axis vector (length NumGridPts)
%   .ves_swp_mesh         - NumGridPts x NumGridPts meshgrid of ves_swp values
%   .epv_swp_mesh         - NumGridPts x NumGridPts meshgrid of epv_swp values
%   .pred_2d      - NumGridPts x NumGridPts predicted scattering surface
%   .pred_masked  - same but NaN in sparse data regions
%   .density_mask         - logical mask: true where data are too sparse
%   .density              - NumGridPts x NumGridPts data point counts per cell
%   .slice_ves_levels     - ves_swp quantile values used for marginal slices
%   .slice_quantiles      - quantile fractions used [0.10, 0.50, 0.90]
%   .slice_epv            - cell {3x1} scattering ~ epv_swp at each ves slice
%   .ci_slice             - cell {3x1} Nx2 [lower, upper] bootstrap CI
%   .T_fit                - (possibly subsampled) table used for fitting
%
% NOTES
%   * The interaction is captured via an engineered product feature:
%       ves_epv_swp = ves_swp .* epv_swp
%     The formula becomes:  outcome ~ ves_swp + epv_swp + ves_epv_swp
%   * Prediction tables always set ves_epv_swp = ves_swp .* epv_swp so
%     the interaction column is internally consistent with the main effects.
%   * The 2D surface evaluates predictions across the full joint space.
%     Cells with fewer than MinDensity data points are masked (NaN).
%   * Marginal slice plots show the EPV SWP effect at low / median / high
%     vessel SWP (10th, 50th, 90th percentiles).  Fanning of the three
%     curves indicates a meaningful interaction between the two predictors.
%   * Bootstrap CIs are computed on the marginal slices only (not the full
%     2D surface, which would be prohibitively expensive).
%

% =====================================================================
% 1. Parse and validate inputs
% =====================================================================
p = inputParser();

addRequired(p,  'T',             @istable);
addParameter(p, 'NumTrees',      100,   @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'MaxSplits',     10,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'LearnRate',     0.1,   @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'NumGridPts',    100,   @(x) isnumeric(x) && isscalar(x) && x > 1);
addParameter(p, 'MaxSample',     1e5,   @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Verbose',       true,  @islogical);
addParameter(p, 'BootstrapCI',   true,  @islogical);
addParameter(p, 'NBootstrap',    100,   @(x) isnumeric(x) && isscalar(x) && x >= 10);
addParameter(p, 'CIAlpha',       0.05,  @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'MinDensity',    5,     @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'TitleStr',      '');
addParameter(p, 'dirout',        pwd);
addParameter(p, 'stain_name','');
addParameter(p, 'section_name','');

parse(p, T, varargin{:});
opts = p.Results;

% Print the current iteration
fprintf('\n----Stain %s, section %s----\n', opts.stain_name, opts.section_name)

% Check required columns exist
required_cols = {'epv_swp','ves_swp','stain'};
missing = required_cols(~ismember(required_cols, T.Properties.VariableNames));
if ~isempty(missing)
    error('fit_swp_histo_gam:missingColumns', ...
          'Table is missing required column(s): %s', strjoin(missing, ', '));
end

% Check columns are numeric and not all-NaN
for c = required_cols
    col = T.(c{1});
    if ~isnumeric(col)
        error('fit_swp_histo_gam:nonNumeric', 'Column "%s" must be numeric.', c{1});
    end
    if all(isnan(col))
        error('fit_swp_histo_gam:allNaN', 'Column "%s" is all NaN.', c{1});
    end
end

% Check toolbox availability
if ~license('test', 'Statistics_Toolbox')
    error('fit_swp_histo_gam:noToolbox', ...
          'Statistics and Machine Learning Toolbox is required for fitrgam.');
end

% =====================================================================
% 2. Remove rows with NaN or Inf in any required column
% =====================================================================
T_work    = T(:, required_cols);
arr       = table2array(T_work);
bad_rows  = any(isnan(arr) | isinf(arr), 2);
n_removed = sum(bad_rows);
T_work(bad_rows, :) = [];

if opts.Verbose && n_removed > 0
    fprintf('[fit_swp_histo_gam] Removed %d rows containing NaN or Inf.\n', n_removed);
end

if height(T_work) < 10
    error('fit_swp_histo_gam:insufficientData', ...
          'Too few valid rows (%d) after removing NaN/Inf.', height(T_work));
end

% =====================================================================
% 3. Subsample if necessary (stratified by ves_swp quantile)
%    Stratify on ves_swp because it is the held predictor in the
%    marginal slices - preserving its distribution matters most.
% =====================================================================
N = height(T_work);
if N > opts.MaxSample
    if opts.Verbose
        fprintf('[fit_swp_histo_gam] Subsampling %d -> %d voxels (stratified by ves_swp).\n', ...
                N, round(opts.MaxSample));
    end
    n_strata   = 10;
    edges      = quantile(T_work.ves_swp, linspace(0, 1, n_strata + 1));
    edges(1)   = -Inf;
    edges(end) = Inf;
    n_per      = floor(opts.MaxSample / n_strata);
    keep_idx   = [];
    for s = 1:n_strata
        in_stratum = find(T_work.ves_swp > edges(s) & T_work.ves_swp <= edges(s+1));
        if isempty(in_stratum), continue; end
        n_draw   = min(n_per, numel(in_stratum));
        selected = in_stratum(randperm(numel(in_stratum), n_draw));
        keep_idx = [keep_idx; selected]; %#ok<AGROW>
    end
    T_work = T_work(keep_idx, :);
end

% =====================================================================
% 4. Engineer interaction feature
%    ves_epv_swp = ves_swp .* epv_swp
%    Allows the effect of epv_swp to vary continuously with ves_swp
%    level (and vice versa) within fitrgam's additive framework.
%    Every prediction table must set ves_epv_swp = ves_swp .* epv_swp
%    to remain internally consistent with the training formula.
% =====================================================================
T_work.ves_epv_swp = T_work.ves_swp .* T_work.epv_swp;

formula = 'stain ~ ves_swp + epv_swp + ves_epv_swp';

% =====================================================================
% 5. Fit GAMs
% =====================================================================
if opts.Verbose
    fprintf('[fit_swp_histo_gam] Fitting GAM (N=%d voxels)...\n', height(T_work));
end

gam = fitrgam(T_work, formula, ...
    'NumTreesPerPredictor',          opts.NumTrees,  ...
    'MaxNumSplitsPerPredictor',      opts.MaxSplits, ...
    'InitialLearnRateForPredictors', opts.LearnRate);

if opts.Verbose
    fprintf('[fit_swp_histo_gam] Finished fitting GAM...\n');
end

if opts.Verbose
    yhat = resubPredict(gam);
    r2   = compute_r2(T_work.stain, yhat);
    fprintf('[fit_swp_histo_gam] GAM resubstitution R^2 = %.4f\n', r2);
    fprintf('  (Resubstitution R^2 is optimistic; use cross-validation for generalization.)\n');
end

% =====================================================================
% 6. Build full 2D prediction grid
%
%    meshgrid(ves_swp_grid, epv_swp_grid) produces:
%      ves_swp_mesh(i,j) = ves_swp_grid(j)   [ves_swp varies along columns]
%      epv_swp_mesh(i,j) = epv_swp_grid(i)   [epv_swp varies along rows]
%    Predictions are reshaped back to n_pts x n_pts after predict().
%    The interaction column is always ves_swp_mesh(:) .* epv_swp_mesh(:).
% =====================================================================
n_pts        = opts.NumGridPts;
ves_swp_grid = linspace(min(T_work.ves_swp), max(T_work.ves_swp), n_pts)';
epv_swp_grid = linspace(min(T_work.epv_swp), max(T_work.epv_swp), n_pts)';

[ves_swp_mesh, epv_swp_mesh] = meshgrid(ves_swp_grid, epv_swp_grid);

T_pred_2d = table(ves_swp_mesh(:), epv_swp_mesh(:), ...
                  ves_swp_mesh(:) .* epv_swp_mesh(:), ...
    'VariableNames', {'ves_swp', 'epv_swp', 'ves_epv_swp'});

if opts.Verbose
    fprintf('[fit_swp_histo_gam] Predicting over %d x %d = %d grid points...\n', ...
            n_pts, n_pts, n_pts^2);
end

pred_2d = reshape(predict(gam, T_pred_2d), n_pts, n_pts);

% =====================================================================
% 7. Data density mask
%    Count data points per grid cell; mask cells below MinDensity.
%    histcounts2(row_var, col_var, row_edges, col_edges):
%      rows = epv_swp, cols = ves_swp  ->  matches meshgrid orientation.
% =====================================================================
swp_edges  = [-Inf; (ves_swp_grid(1:end-1) + ves_swp_grid(2:end)) / 2; Inf];
epvs_edges = [-Inf; (epv_swp_grid(1:end-1) + epv_swp_grid(2:end)) / 2; Inf];

density      = histcounts2(T_work.epv_swp, T_work.ves_swp, epvs_edges, swp_edges);
density_mask = density < opts.MinDensity;

pred_masked = pred_2d;
pred_masked(density_mask) = NaN;

if opts.Verbose
    pct_masked = 100 * mean(density_mask(:));
    fprintf('[fit_swp_histo_gam] %.1f%% of grid cells masked (< %d data points).\n', ...
            pct_masked, opts.MinDensity);
end

% =====================================================================
% 8. Marginal slice predictions
%
%    Sweep epv_swp across its full range while holding ves_swp fixed
%    at its 10th, 50th, and 90th percentiles.  This directly answers:
%    "How do scattering and retardance change with EPVS proximity,
%     for low / medium / high vessel proximity?"
%
%    The interaction column must equal ves_fixed .* epv_swp_grid at
%    every point to stay consistent with the training formula.
% =====================================================================
slice_quantiles = [0.10, 0.50, 0.90];
ves_levels      = quantile(T_work.ves_swp, slice_quantiles);
n_slices        = numel(ves_levels);

slice_epv = cell(n_slices, 1);

for sl = 1:n_slices
    ves_fixed = ves_levels(sl) * ones(n_pts, 1);   % ves_swp held constant per slice
    T_sl = table(ves_fixed, epv_swp_grid, ves_fixed .* epv_swp_grid, ...
        'VariableNames', {'ves_swp', 'epv_swp', 'ves_epv_swp'});
    slice_epv{sl} = predict(gam, T_sl);
end

bootstrap_cols = {'stain', 'ves_swp', 'epv_swp', 'ves_epv_swp'};
T_boot_base = T_work(:, bootstrap_cols);

% =====================================================================
% 9. Bootstrap CIs on marginal slices
%    Resample rows with replacement, refit both GAMs, collect
%    predictions on the fixed epv_swp grid at each ves_swp slice.
%    Uses parfor if Parallel Computing Toolbox is available.
% =====================================================================
ci_slice = cell(n_slices, 1);

if opts.BootstrapCI
    if opts.Verbose
        fprintf('[fit_swp_histo_gam] Running %d bootstrap resamples for slice CIs...\n', ...
                opts.NBootstrap);
    end

    B        = opts.NBootstrap;
    n_fit    = height(T_work);
    b_trees  = opts.NumTrees;
    b_splits = opts.MaxSplits;
    b_lr     = opts.LearnRate;

    % Pre-allocate: rows = grid points, cols = bootstrap resamples,
    % pages = ves_swp slices
    boot = zeros(n_pts, B, n_slices);

    % Pre-build slice prediction tables outside the bootstrap loop.
    % Each table sweeps epv_swp with ves_swp held at one quantile level.
    T_slices = cell(n_slices, 1);
    for sl = 1:n_slices
        ves_fixed    = ves_levels(sl) * ones(n_pts, 1);
        T_slices{sl} = table(ves_fixed, epv_swp_grid, ves_fixed .* epv_swp_grid, ...
            'VariableNames', {'ves_swp', 'epv_swp', 'ves_epv_swp'});
    end

    use_par = ~isempty(ver('parallel'));
    if use_par && opts.Verbose
        fprintf('[fit_swp_histo_gam] Parallel Computing Toolbox detected - using parfor.\n');
    end

    if use_par
        parfor b = 1:B
            idx_b = randi(n_fit, n_fit, 1);
            T_b   = T_boot_base(idx_b, :);

            mdl = fitrgam(T_b, 'stain ~ ves_swp + epv_swp + ves_epv_swp', ...
                'NumTreesPerPredictor',          b_trees,  ...
                'MaxNumSplitsPerPredictor',      b_splits, ...
                'InitialLearnRateForPredictors', b_lr);

            tmp = zeros(n_pts, n_slices);
            for sl = 1:n_slices
                tmp(:, sl) = predict(mdl, T_slices{sl});
            end
            boot(:, b, :) = tmp;
        end
    else
        for b = 1:B
            if opts.Verbose && mod(b, 10) == 0
                fprintf('[fit_swp_histo_gam]   Bootstrap resample %d / %d\n', b, B);
            end
            idx_b = randi(n_fit, n_fit, 1);
            T_b   = T_boot_base(idx_b, :);

            mdl = fitrgam(T_b, 'stain ~ ves_swp + epv_swp + ves_epv_swp', ...
                'NumTreesPerPredictor',          b_trees,  ...
                'MaxNumSplitsPerPredictor',      b_splits, ...
                'InitialLearnRateForPredictors', b_lr);

            for sl = 1:n_slices
                boot(:, b, sl) = predict(mdl, T_slices{sl});
            end
        end
    end

    % Percentile-based CIs per slice.
    % Row vector [lo, hi] avoids ambiguous output shape in older MATLAB.
    alpha = opts.CIAlpha;
    lo    = 100 *  alpha / 2;
    hi    = 100 * (1 - alpha / 2);

    for sl = 1:n_slices
        ci_slice{sl} = prctile(boot(:, :, sl), [lo, hi], 2);  % Nx2
    end

    if opts.Verbose
        fprintf('[fit_swp_histo_gam] Bootstrap CIs computed (alpha = %.2f).\n', alpha);
    end
end

% =====================================================================
% 10. Package results
% =====================================================================
results.gam                 = gam;
results.ves_swp_grid        = ves_swp_grid;
results.epv_swp_grid        = epv_swp_grid;
results.ves_swp_mesh        = ves_swp_mesh;
results.epv_swp_mesh        = epv_swp_mesh;
results.pred_2d             = pred_2d;      % unmasked
results.pred_masked         = pred_masked;  % NaN in sparse cells
results.density_mask        = density_mask;
results.density             = density;
results.slice_ves_levels    = ves_levels;           % ves_swp quantiles used
results.slice_quantiles     = slice_quantiles;      % [0.10, 0.50, 0.90]
results.slice_epv           = slice_epv;    % stain ~ epv_swp curves
results.ci_slice            = ci_slice;     % Nx2 CIs per ves slice
results.T_fit               = T_work;

if opts.Verbose
    fprintf('[fit_swp_histo_gam] Done.\n');
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


