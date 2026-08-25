function tuning = tune_swp_psoct_gam(T, varargin)
% TUNE_SWP_PSOCT_GAM  Sweep fitrgam complexity parameters (MaxSplits,
%                     NumTrees) for scattering/retardance GAMs and report
%                     resubstitution vs. k-fold cross-validated R^2 for
%                     each combination, plus a separate MinDensity /
%                     grid-masking sweep -- so per-dataset tuning can be
%                     done from printed diagnostics instead of by eye.
%
%   Model complexity (MaxSplits, NumTrees) requires refitting for each
%   combination. MinDensity only affects which grid cells get masked in
%   the *display* of an already-fitted surface, so it is swept separately
%   and cheaply, without refitting, using one representative fit.
%
% USAGE
%   tuning = tune_swp_psoct_gam(T)
%   tuning = tune_swp_psoct_gam(T, 'MaxSplitsGrid', [5 10 20 40], ...
%                                   'KFold', 5, 'Verbose', true)
%
% REQUIRED INPUT
%   T : table with columns 'scattering', 'retardance', 'ves_swp', 'epv_swp'
%       (same requirements as fit_swp_psoct_gam)
%
% OPTIONAL NAME-VALUE PAIRS
%   'MaxSplitsGrid'   Values of MaxNumSplitsPerPredictor to try   (default: [5 10 20 40])
%   'NumTreesGrid'    Values of NumTreesPerPredictor to try       (default: [50 100 200])
%   'LearnRate'       Fixed learning rate used for all fits       (default: 0.1)
%   'KFold'           Number of CV folds                          (default: 5)
%   'MaxSample'       Max voxels used for fitting (subsampling)   (default: 1e5)
%   'MinDensityGrid'  MinDensity thresholds to report masking for (default: [3 5 10 20 40])
%   'NumGridPts'      Grid resolution used for the density sweep  (default: 100)
%   'Verbose'         Print progress and a summary table          (default: true)
%   'TitleStr'        Label used in printed output / saved table  (default: '')
%   'dirout'          If non-empty, writes tuning.summary_table   (default: '')
%                      to <dirout>/GAM_tuning_<TitleStr>.csv
%
% OUTPUT  tuning struct with fields:
%   .complexity_table   table, one row per (MaxSplits, NumTrees) combo, with
%                        columns: MaxSplits, NumTrees, R2_scatter_resub,
%                        R2_scatter_cv, gap_scatter, R2_retard_resub,
%                        R2_retard_cv, gap_retard, n_fit
%   .density_table       table, one row per MinDensity threshold, with
%                        columns: MinDensity, pct_masked
%   .recommended         struct with suggested MaxSplits/NumTrees: the
%                        combo with the smallest max(gap_scatter, gap_retard)
%                        among combos whose CV R^2 is within 90% of the
%                        best observed CV R^2 (i.e. "simplest model that
%                        isn't leaving much CV performance on the table")
%   .T_fit               (possibly subsampled) table used for fitting
%
% NOTES
%   * "gap" = R2_resub - R2_cv. A large gap at fixed N means the model is
%     fitting noise rather than signal -- typical with small/unbalanced
%     per-stage sample sizes. Prefer the simplest combo (fewest splits/
%     trees) whose CV R^2 is close to the best CV R^2 observed, rather
%     than the combo with the single highest resubstitution R^2.
%   * CV R^2 here is computed via RegressionGAM's built-in crossval/
%     kfoldLoss, which uses out-of-fold predictions -- this is honest
%     held-out performance, unlike resubstitution R^2.
%   * This function does NOT fix the pseudoreplication issue (voxels
%     nested within subjects/volumes): k-fold splits are voxel-level,
%     the same as in fit_swp_psoct_gam's bootstrap. CV R^2 here answers
%     "does this model generalize to held-out voxels from the same
%     pool," not "does it generalize to a new subject." Keep both
%     caveats explicit in the manuscript if this diagnostic informs a
%     reported figure.

    % =====================================================================
    % 1. Parse inputs
    % =====================================================================
    p = inputParser();
    p.FunctionName = 'tune_swp_psoct_gam';

    addRequired(p,  'T',               @istable);
    addParameter(p, 'MaxSplitsGrid',   [5 10 20 40],   @isnumeric);
    addParameter(p, 'NumTreesGrid',    [50 100 200],   @isnumeric);
    addParameter(p, 'LearnRate',       0.1,            @(x) isnumeric(x) && x > 0 && x <= 1);
    addParameter(p, 'KFold',           5,              @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'MaxSample',       1e5,            @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MinDensityGrid',  [3 5 10 20 40], @isnumeric);
    addParameter(p, 'NumGridPts',      100,            @(x) isnumeric(x) && isscalar(x) && x > 1);
    addParameter(p, 'Verbose',         true,           @islogical);
    addParameter(p, 'TitleStr',        '');
    addParameter(p, 'dirout',          '');

    parse(p, T, varargin{:});
    opts = p.Results;

    required_cols = {'scattering', 'retardance', 'ves_swp', 'epv_swp'};
    missing = required_cols(~ismember(required_cols, T.Properties.VariableNames));
    if ~isempty(missing)
        error('tune_swp_psoct_gam:missingColumns', ...
              'Table is missing required column(s): %s', strjoin(missing, ', '));
    end
    if ~license('test', 'Statistics_Toolbox')
        error('tune_swp_psoct_gam:noToolbox', ...
              'Statistics and Machine Learning Toolbox is required for fitrgam.');
    end

    % =====================================================================
    % 2. Clean + subsample (same logic as fit_swp_psoct_gam, kept local
    %    and self-contained so this helper can be run standalone before
    %    ever calling the main fitting function)
    % =====================================================================
    T_work = T(:, required_cols);
    arr    = table2array(T_work);
    T_work(any(isnan(arr) | isinf(arr), 2), :) = [];

    if height(T_work) < 10
        error('tune_swp_psoct_gam:insufficientData', ...
              'Too few valid rows (%d) after removing NaN/Inf.', height(T_work));
    end

    N = height(T_work);
    if N > opts.MaxSample
        if opts.Verbose
            fprintf('[tune_swp_psoct_gam] Subsampling %d -> %d voxels (stratified by ves_swp).\n', ...
                    N, round(opts.MaxSample));
        end
        n_strata = 10;
        edges    = quantile(T_work.ves_swp, linspace(0, 1, n_strata + 1));
        edges(1) = -Inf; edges(end) = Inf;
        n_per    = floor(opts.MaxSample / n_strata);
        keep_idx = [];
        for s = 1:n_strata
            in_stratum = find(T_work.ves_swp > edges(s) & T_work.ves_swp <= edges(s+1));
            if isempty(in_stratum), continue; end
            n_draw   = min(n_per, numel(in_stratum));
            selected = in_stratum(randperm(numel(in_stratum), n_draw));
            keep_idx = [keep_idx; selected]; %#ok<AGROW>
        end
        T_work = T_work(keep_idx, :);
    end

    T_work.ves_epv_swp = T_work.ves_swp .* T_work.epv_swp;
    n_fit = height(T_work);

    var_scatter = var(T_work.scattering);
    var_retard  = var(T_work.retardance);

    % =====================================================================
    % 3. Sweep MaxSplits x NumTrees, refitting both GAMs each time
    % =====================================================================
    grid_splits = opts.MaxSplitsGrid(:)';
    grid_trees  = opts.NumTreesGrid(:)';
    n_combo     = numel(grid_splits) * numel(grid_trees);

    MaxSplits        = zeros(n_combo, 1);
    NumTrees          = zeros(n_combo, 1);
    R2_scatter_resub  = zeros(n_combo, 1);
    R2_scatter_cv     = zeros(n_combo, 1);
    R2_retard_resub   = zeros(n_combo, 1);
    R2_retard_cv      = zeros(n_combo, 1);

    if opts.Verbose
        fprintf('[tune_swp_psoct_gam] Sweeping %d (MaxSplits x NumTrees) combinations, N=%d voxels, %d-fold CV...\n', ...
                n_combo, n_fit, opts.KFold);
    end

    combo_i = 0;
    for ns = grid_splits
        for nt = grid_trees
            combo_i = combo_i + 1;
            MaxSplits(combo_i) = ns;
            NumTrees(combo_i)  = nt;

            gam_s = fitrgam(T_work, 'scattering ~ ves_swp + epv_swp + ves_epv_swp', ...
                'NumTreesPerPredictor',          nt, ...
                'MaxNumSplitsPerPredictor',      ns, ...
                'InitialLearnRateForPredictors', opts.LearnRate);
            gam_r = fitrgam(T_work, 'retardance ~ ves_swp + epv_swp + ves_epv_swp', ...
                'NumTreesPerPredictor',          nt, ...
                'MaxNumSplitsPerPredictor',      ns, ...
                'InitialLearnRateForPredictors', opts.LearnRate);

            yhat_s = resubPredict(gam_s);
            yhat_r = resubPredict(gam_r);
            R2_scatter_resub(combo_i) = local_r2(T_work.scattering, yhat_s);
            R2_retard_resub(combo_i)  = local_r2(T_work.retardance,  yhat_r);

            cv_s = crossval(gam_s, 'KFold', opts.KFold);
            cv_r = crossval(gam_r, 'KFold', opts.KFold);
            R2_scatter_cv(combo_i) = 1 - kfoldLoss(cv_s) / var_scatter;
            R2_retard_cv(combo_i)  = 1 - kfoldLoss(cv_r) / var_retard;

            if opts.Verbose
                fprintf(['  [%2d/%2d] MaxSplits=%3d NumTrees=%3d | ' ...
                         'scatter R2 resub=%.3f cv=%.3f (gap=%.3f) | ' ...
                         'retard R2 resub=%.3f cv=%.3f (gap=%.3f)\n'], ...
                        combo_i, n_combo, ns, nt, ...
                        R2_scatter_resub(combo_i), R2_scatter_cv(combo_i), ...
                        R2_scatter_resub(combo_i) - R2_scatter_cv(combo_i), ...
                        R2_retard_resub(combo_i),  R2_retard_cv(combo_i), ...
                        R2_retard_resub(combo_i)  - R2_retard_cv(combo_i));
            end
        end
    end

    gap_scatter = R2_scatter_resub - R2_scatter_cv;
    gap_retard  = R2_retard_resub  - R2_retard_cv;
    n_fit_col   = repmat(n_fit, n_combo, 1);

    complexity_table = table(MaxSplits, NumTrees, ...
        R2_scatter_resub, R2_scatter_cv, gap_scatter, ...
        R2_retard_resub,  R2_retard_cv,  gap_retard, ...
        n_fit_col, ...
        'VariableNames', {'MaxSplits', 'NumTrees', ...
        'R2_scatter_resub', 'R2_scatter_cv', 'gap_scatter', ...
        'R2_retard_resub', 'R2_retard_cv', 'gap_retard', 'n_fit'});

    % =====================================================================
    % 4. Recommend simplest combo within 90% of best observed CV R^2
    %    "Best" combined CV performance = mean of the two CV R^2 values.
    %    "Simplest" = smallest MaxSplits*NumTrees among qualifying combos.
    % =====================================================================
    combined_cv = (R2_scatter_cv + R2_retard_cv) / 2;
    best_cv     = max(combined_cv);
    threshold   = 0.90 * best_cv;
    qualifying  = combined_cv >= threshold;

    complexity_score = MaxSplits .* NumTrees;
    idx_pool = find(qualifying);
    [~, ix]  = min(complexity_score(idx_pool));
    rec_idx  = idx_pool(ix);

    recommended = struct( ...
        'MaxSplits',  MaxSplits(rec_idx), ...
        'NumTrees',   NumTrees(rec_idx), ...
        'R2_scatter_cv', R2_scatter_cv(rec_idx), ...
        'R2_retard_cv',  R2_retard_cv(rec_idx), ...
        'gap_scatter',   gap_scatter(rec_idx), ...
        'gap_retard',    gap_retard(rec_idx), ...
        'note', ['Simplest (MaxSplits x NumTrees) combo whose mean CV R^2 ', ...
                 'is within 90% of the best mean CV R^2 observed in the sweep.']);

    if opts.Verbose
        fprintf('\n[tune_swp_psoct_gam] Recommended: MaxSplits=%d, NumTrees=%d\n', ...
                recommended.MaxSplits, recommended.NumTrees);
        fprintf('  CV R^2: scatter=%.3f, retard=%.3f | resub-CV gap: scatter=%.3f, retard=%.3f\n', ...
                recommended.R2_scatter_cv, recommended.R2_retard_cv, ...
                recommended.gap_scatter, recommended.gap_retard);
        fprintf('  (Recommendation favors the simplest model within 90%% of peak CV performance --\n');
        fprintf('   inspect the full table before trusting this on a very small/unbalanced subset.)\n');
    end

    % =====================================================================
    % 5. MinDensity / masking sweep -- single representative fit, no
    %    refitting needed since MinDensity only changes which grid cells
    %    are masked, not the model or its predictions.
    % =====================================================================
    n_pts        = opts.NumGridPts;
    ves_swp_grid = linspace(min(T_work.ves_swp), max(T_work.ves_swp), n_pts)';
    epv_swp_grid = linspace(min(T_work.epv_swp), max(T_work.epv_swp), n_pts)';

    swp_edges  = [-Inf; (ves_swp_grid(1:end-1) + ves_swp_grid(2:end)) / 2; Inf];
    epvs_edges = [-Inf; (epv_swp_grid(1:end-1) + epv_swp_grid(2:end)) / 2; Inf];
    density    = histcounts2(T_work.epv_swp, T_work.ves_swp, epvs_edges, swp_edges);

    MinDensity = opts.MinDensityGrid(:);
    pct_masked = zeros(numel(MinDensity), 1);
    for i = 1:numel(MinDensity)
        pct_masked(i) = 100 * mean(density(:) < MinDensity(i));
    end
    density_table = table(MinDensity, pct_masked, ...
        'VariableNames', {'MinDensity', 'pct_masked'});

    if opts.Verbose
        fprintf('\n[tune_swp_psoct_gam] Grid masking at NumGridPts=%d (%d x %d = %d cells):\n', ...
                n_pts, n_pts, n_pts, n_pts^2);
        for i = 1:numel(MinDensity)
            fprintf('  MinDensity=%3d -> %.1f%% of cells masked\n', MinDensity(i), pct_masked(i));
        end
        fprintf(['  (If most cells are masked at your intended MinDensity, lower it, lower\n' ...
                 '   NumGridPts, or accept that this dataset only supports a coarser surface.)\n']);
    end

    % =====================================================================
    % 6. Package + optional CSV export
    % =====================================================================
    tuning.complexity_table = complexity_table;
    tuning.density_table    = density_table;
    tuning.recommended      = recommended;
    tuning.T_fit            = T_work;

    if ~isempty(opts.dirout)
        if ~isfolder(opts.dirout)
            mkdir(opts.dirout);
        end
        fname = sprintf('GAM_tuning_%s.csv', strrep(opts.TitleStr, ' ', '_'));
        fname = regexprep(fname, '_+\.csv$', '.csv');
        writetable(complexity_table, fullfile(opts.dirout, fname));
        if opts.Verbose
            fprintf('\n[tune_swp_psoct_gam] Wrote complexity table to %s\n', ...
                    fullfile(opts.dirout, fname));
        end
    end
end


% =========================================================================
% LOCAL FUNCTION: R-squared (resubstitution)
% =========================================================================
function r2 = local_r2(y, yhat)
    ss_res = sum((y - yhat).^2);
    ss_tot = sum((y - mean(y)).^2);
    if ss_tot == 0
        r2 = NaN;
    else
        r2 = 1 - ss_res / ss_tot;
    end
end
