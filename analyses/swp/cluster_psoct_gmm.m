function results = cluster_psoct_gmm(T, varargin)
% CLUSTER_PSOCT_GMM  GMM clustering of PS-OCT voxels on scattering and
%                    retardance, with post-hoc stratification by EPVS SWP
%                    and Vessel SWP quantile bins.
%
% DESIGN
%   Analysis 1 — Stratify by EPVS SWP quantile bins (10 bins):
%     Features: scattering, retardance
%     Post-hoc: cluster composition (proportion per cluster) across bins
%
%   Analysis 2 — Stratify by Vessel SWP quantile bins (10 bins):
%     Features: scattering, retardance
%     Post-hoc: cluster composition across bins
%
%   A single GMM is fit on the pooled data (all voxels combined). Cluster
%   labels are globally consistent, so prevalence shifts across proximity
%   bins are directly comparable. BIC is used to select the optimal number
%   of components k up to KMax.
%
% USAGE
%   results = cluster_psoct_gmm(T)
%   results = cluster_psoct_gmm(T, 'KMax', 8, 'NQuantileBins', 10)
%
% REQUIRED INPUT
%   T : table with columns:
%         'scattering'  - scattering coefficient (voxel-level)
%         'retardance'  - optical retardance      (voxel-level)
%         'SWP'         - vessel size-weighted proximity
%         'EPVS_SWP'    - EPVS size-weighted proximity
%
% OPTIONAL NAME-VALUE PAIRS
%   'KMax'          Maximum number of GMM components to evaluate    (default: 8)
%   'KFixed'        Fix k to this value, skipping BIC search        (default: [], use BIC)
%   'NQuantileBins' Number of quantile bins for stratification      (default: 10)
%   'MaxSample'     Max voxels for GMM fitting (subsampled)         (default: 1e5)
%   'CovType'       GMM covariance type: 'full','diagonal','spherical' (default: 'full')
%   'NReplicates'   GMM fitting replicates (guards against local minima) (default: 5)
%   'ZScore'        Z-score scattering and retardance before GMM fitting  (default: true)
%                   Set false when features are already on a common scale
%                   (e.g. percentage difference values for both features)
%   'PlotResults'   Produce output figures                          (default: true)
%   'Verbose'       Print progress to command window                (default: true)
%
% OUTPUT  results struct with fields:
%   .gmm              - fitted GMM object (gmdistribution)
%   .k_optimal        - selected number of components
%   .bic_values       - BIC for each k tested (length KMax)
%   .cluster_labels   - Nx1 cluster assignment for T_fit rows
%   .cluster_probs    - NxK soft assignment probabilities
%   .T_fit            - subsampled table used for GMM fitting
%   .bin_edges_epvs   - bin edges used for EPVS SWP stratification
%   .bin_edges_swp    - bin edges used for Vessel SWP stratification
%   .comp_epvs        - (NQuantileBins x K) cluster composition per EPVS bin
%   .comp_swp         - (NQuantileBins x K) cluster composition per Vessel SWP bin
%   .cluster_stats    - struct: mean/std of scattering & retardance per cluster
%
% EXAMPLE
%   results = cluster_psoct_gmm(T, 'KMax', 6, 'NReplicates', 3);

    % -------------------------------------------------------------------------
    % 1. Parse and validate inputs
    % -------------------------------------------------------------------------
    p = inputParser();
    p.FunctionName = 'cluster_psoct_gmm';

    addRequired(p,  'T',              @istable);
    addParameter(p, 'KMax',           8,       @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'KFixed',         [],      @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    addParameter(p, 'NQuantileBins',  10,      @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'MaxSample',      1e5,     @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'CovType',        'full',  @(x) ismember(x, {'full','diagonal','spherical'}));
    addParameter(p, 'NReplicates',    5,       @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'ZScore',         true,    @islogical);
    addParameter(p, 'PlotResults',    true,    @islogical);
    addParameter(p, 'Verbose',        true,    @islogical);
    addRequired(p,  'TitleStr');
    addRequired(p,  'dirout');

    parse(p, T, varargin{:});
    opts = p.Results;

    required_cols = {'scattering', 'retardance', 'SWP', 'EPVS_SWP'};
    missing = required_cols(~ismember(required_cols, T.Properties.VariableNames));
    if ~isempty(missing)
        error('cluster_psoct_gmm:missingColumns', ...
              'Table is missing required column(s): %s', strjoin(missing, ', '));
    end
    for c = required_cols
        if ~isnumeric(T.(c{1}))
            error('cluster_psoct_gmm:nonNumeric', 'Column "%s" must be numeric.', c{1});
        end
    end
    if ~license('test', 'Statistics_Toolbox')
        error('cluster_psoct_gmm:noToolbox', ...
              'Statistics and Machine Learning Toolbox is required.');
    end

    % -------------------------------------------------------------------------
    % 2. Clean data — remove NaN / Inf rows
    % -------------------------------------------------------------------------
    T_work   = T(:, required_cols);
    arr      = table2array(T_work);
    bad      = any(isnan(arr) | isinf(arr), 2);
    T_work(bad, :) = [];
    if opts.Verbose
        fprintf('[cluster_psoct_gmm] %d valid voxels after removing %d NaN/Inf rows.\n', ...
                height(T_work), sum(bad));
    end

    % -------------------------------------------------------------------------
    % 3. Uniform subsampling across the joint (SWP, EPVS_SWP) distribution
    %
    %    Problem with 1D quantile stratification: SWP variables are highly
    %    right-skewed — most voxels cluster at low values, leaving the upper
    %    tail severely under-represented after downsampling.
    %
    %    Solution: partition the 2D (SWP x EPVS_SWP) space into a uniform
    %    grid of cells (each axis divided into n_grid_bins equal-width bins
    %    spanning the observed range), then draw at most n_per_cell voxels
    %    from each occupied cell. This gives equal representation per unit
    %    area of the joint distribution, including sparse upper-tail regions.
    %
    %    n_per_cell is chosen so that the total sample is <= MaxSample.
    %    Empty cells are skipped; cells with fewer than n_per_cell voxels
    %    contribute all their voxels (no upsampling).
    % -------------------------------------------------------------------------
    N = height(T_work);
    if N > opts.MaxSample
        n_grid_bins = 20;   % bins per axis -> up to 400 cells

        % Equal-width bin edges across the observed range of each variable
        edges_swp  = linspace(min(T_work.SWP),      max(T_work.SWP),      n_grid_bins + 1);
        edges_epvs = linspace(min(T_work.EPVS_SWP), max(T_work.EPVS_SWP), n_grid_bins + 1);

        % Extend outer edges to capture boundary voxels exactly
        edges_swp(1)    = -Inf;   edges_swp(end)    = Inf;
        edges_epvs(1)   = -Inf;   edges_epvs(end)   = Inf;

        % Assign each voxel to a 2D cell
        bin_swp  = discretize(T_work.SWP,      edges_swp);
        bin_epvs = discretize(T_work.EPVS_SWP, edges_epvs);

        % Count occupied cells to set n_per_cell
        cell_ids       = (bin_swp - 1) * n_grid_bins + bin_epvs;
        occupied_cells = unique(cell_ids);
        n_occupied     = numel(occupied_cells);
        n_per_cell     = max(1, floor(opts.MaxSample / n_occupied));

        if opts.Verbose
            fprintf(['[cluster_psoct_gmm] 2D uniform subsampling: ' ...
                     '%d -> ~%d voxels across %d/%d occupied cells ' ...
                     '(%d per cell).\n'], ...
                    N, min(N, n_occupied * n_per_cell), ...
                    n_occupied, n_grid_bins^2, n_per_cell);
        end

        keep_idx = [];
        for c = 1:n_occupied
            in_cell  = find(cell_ids == occupied_cells(c));
            n_draw   = min(n_per_cell, numel(in_cell));
            keep_idx = [keep_idx; ...
                        in_cell(randperm(numel(in_cell), n_draw))]; %#ok<AGROW>
        end
        T_fit = T_work(keep_idx, :);
    else
        T_fit = T_work;
    end

    % -------------------------------------------------------------------------
    % 4. Feature matrix construction
    %    If ZScore = true  : z-score each feature using mean/std from T_fit.
    %                        Required when features have different units/ranges
    %                        (e.g. raw mus and retardance).
    %    If ZScore = false : use features as-is. Appropriate when features are
    %                        already on a common scale (e.g. percentage difference
    %                        values for both mus and retardance).
    %    In both cases the same transform is applied to T_work (all voxels)
    %    using parameters estimated from T_fit, avoiding data leakage.
    % -------------------------------------------------------------------------
    mu_s  = mean(T_fit.scattering);  sig_s = std(T_fit.scattering);
    mu_r  = mean(T_fit.retardance);  sig_r = std(T_fit.retardance);

    if sig_s == 0 || sig_r == 0
        error('cluster_psoct_gmm:zeroVariance', ...
              'Scattering or retardance has zero variance after cleaning.');
    end

    if opts.ZScore
        if opts.Verbose
            fprintf('[cluster_psoct_gmm] Z-scoring features (mu_s=%.4f, sig_s=%.4f, mu_r=%.4f, sig_r=%.4f).\n', ...
                    mu_s, sig_s, mu_r, sig_r);
        end
        X     = [(T_fit.scattering  - mu_s) / sig_s, ...
                 (T_fit.retardance  - mu_r) / sig_r];
        X_all = [(T_work.scattering - mu_s) / sig_s, ...
                 (T_work.retardance - mu_r) / sig_r];
    else
        if opts.Verbose
            fprintf('[cluster_psoct_gmm] Skipping z-score — using features on original scale.\n');
        end
        X     = [T_fit.scattering,  T_fit.retardance];
        X_all = [T_work.scattering, T_work.retardance];
        % Set mu/sig to identity transform so norm_params is always populated
        mu_s = 0;  sig_s = 1;
        mu_r = 0;  sig_r = 1;
    end

    norm_params = struct('mu_s', mu_s, 'sig_s', sig_s, ...
                         'mu_r', mu_r, 'sig_r', sig_r, ...
                         'z_scored', opts.ZScore);

    % -------------------------------------------------------------------------
    % 5. BIC-optimal GMM selection
    %    Fit GMMs for k = 1 .. KMax, pick k with lowest BIC.
    % -------------------------------------------------------------------------
    if opts.Verbose
        fprintf('[cluster_psoct_gmm] Fitting GMMs for k = 1 to %d (N_fit = %d)...\n', ...
                opts.KMax, size(X, 1));
    end

    if ~isempty(opts.KFixed)
        % --- Fixed k: fit a single GMM at the requested k ----------------
        if opts.Verbose
            fprintf('[cluster_psoct_gmm] KFixed = %d specified — skipping BIC search.\n', ...
                    opts.KFixed);
        end
        bic_values = nan(opts.KFixed, 1);   % populated only at KFixed index
        try
            gmm_best = fitgmdist(X, opts.KFixed, ...
                'CovarianceType',    opts.CovType,      ...
                'RegularizationValue', 1e-5,            ...
                'Replicates',        opts.NReplicates,  ...
                'Options',           statset('MaxIter', 500, 'TolFun', 1e-6));
            bic_values(opts.KFixed) = gmm_best.BIC;
            k_optimal = opts.KFixed;
        catch ME
            error('cluster_psoct_gmm:fitFailed', ...
                  'fitgmdist failed for KFixed=%d: %s', opts.KFixed, ME.message);
        end
        if opts.Verbose
            fprintf('[cluster_psoct_gmm] k = %d  BIC = %.2f\n', ...
                    k_optimal, gmm_best.BIC);
        end

    else
        % --- BIC search: fit GMMs for k = 1 .. KMax ----------------------
        bic_values = nan(opts.KMax, 1);
        gmm_models = cell(opts.KMax, 1);

        for k = 1:opts.KMax
            if opts.Verbose
                fprintf('[cluster_psoct_gmm]   k = %d ...', k);
            end
            try
                gmm_k = fitgmdist(X, k, ...
                    'CovarianceType',    opts.CovType,      ...
                    'RegularizationValue', 1e-5,            ...
                    'Replicates',        opts.NReplicates,  ...
                    'Options',           statset('MaxIter', 500, 'TolFun', 1e-6));
                bic_values(k)  = gmm_k.BIC;
                gmm_models{k}  = gmm_k;
                if opts.Verbose
                    fprintf('  BIC = %.2f\n', gmm_k.BIC);
                end
            catch ME
                if opts.Verbose
                    fprintf('  FAILED (%s)\n', ME.message);
                end
            end
        end

        % Select k with lowest BIC (ignoring failed fits)
        [~, k_optimal] = min(bic_values);
        gmm_best       = gmm_models{k_optimal};

        if opts.Verbose
            fprintf('[cluster_psoct_gmm] Optimal k = %d (BIC = %.2f)\n', ...
                    k_optimal, bic_values(k_optimal));
        end
    end

    % -------------------------------------------------------------------------
    % 6. Assign cluster labels to ALL voxels (T_work) using fitted GMM
    %    This ensures stratification uses the full voxel population.
    % -------------------------------------------------------------------------
    [labels_all, ~, probs_all] = cluster(gmm_best, X_all);

    % -------------------------------------------------------------------------
    % 7. Compute cluster summary statistics (original scale)
    % -------------------------------------------------------------------------
    K = k_optimal;
    cluster_stats.mean_scatter = zeros(K, 1);
    cluster_stats.std_scatter  = zeros(K, 1);
    cluster_stats.mean_retard  = zeros(K, 1);
    cluster_stats.std_retard   = zeros(K, 1);
    cluster_stats.n_voxels     = zeros(K, 1);

    for k = 1:K
        idx_k = labels_all == k;
        cluster_stats.mean_scatter(k) = mean(T_work.scattering(idx_k));
        cluster_stats.std_scatter(k)  = std(T_work.scattering(idx_k));
        cluster_stats.mean_retard(k)  = mean(T_work.retardance(idx_k));
        cluster_stats.std_retard(k)   = std(T_work.retardance(idx_k));
        cluster_stats.n_voxels(k)     = sum(idx_k);
    end

    if opts.Verbose
        fprintf('[cluster_psoct_gmm] Cluster summary (original scale):\n');
        fprintf('  %6s  %12s  %12s  %10s\n', ...
                'Cluster', 'Mean Scatter', 'Mean Retard', 'N voxels');
        for k = 1:K
            fprintf('  %6d  %12.4f  %12.4f  %10d\n', k, ...
                    cluster_stats.mean_scatter(k), ...
                    cluster_stats.mean_retard(k),  ...
                    cluster_stats.n_voxels(k));
        end
    end

    % -------------------------------------------------------------------------
    % 8. Post-hoc stratification: cluster composition per quantile bin
    %
    %    Analysis 1: bin by EPVS_SWP quantiles
    %    Analysis 2: bin by Vessel SWP quantiles
    %
    %    comp(b, k) = fraction of voxels in bin b assigned to cluster k
    % -------------------------------------------------------------------------
    Q = opts.NQuantileBins;

    [comp_epvs, bin_edges_epvs] = compute_composition( ...
        T_work.EPVS_SWP, labels_all, K, Q);

    [comp_swp, bin_edges_swp] = compute_composition( ...
        T_work.SWP, labels_all, K, Q);

    % -------------------------------------------------------------------------
    % 9. Plot
    % -------------------------------------------------------------------------
    if opts.PlotResults
        if isempty(opts.KFixed)
            plot_bic(bic_values, k_optimal, opts.KMax);
        end
        plot_feature_space(X, labels_all, K, norm_params, ...
                           cluster_stats, gmm_best, opts.dirout, opts.TitleStr);
        plot_composition(comp_epvs, bin_edges_epvs, K, ...
                         'EPVS SWP', cluster_stats, opts.dirout, opts.TitleStr);
        plot_composition(comp_swp,  bin_edges_swp,  K, ...
                         'Vessel SWP', cluster_stats, opts.dirout, opts.TitleStr);
    end

    % -------------------------------------------------------------------------
    % 10. Package results
    % -------------------------------------------------------------------------
    results.gmm             = gmm_best;
    results.k_optimal       = k_optimal;
    results.bic_values      = bic_values;
    results.cluster_labels  = labels_all;
    results.cluster_probs   = probs_all;
    results.T_fit           = T_fit;
    results.T_work          = T_work;
    results.norm_params     = norm_params;
    results.bin_edges_epvs  = bin_edges_epvs;
    results.bin_edges_swp   = bin_edges_swp;
    results.comp_epvs       = comp_epvs;
    results.comp_swp        = comp_swp;
    results.cluster_stats   = cluster_stats;

    if opts.Verbose
        fprintf('[cluster_psoct_gmm] Done.\n');
    end
end


% =========================================================================
% LOCAL: compute cluster composition across quantile bins
%   x       : Nx1 proximity values (EPVS_SWP or SWP)
%   labels  : Nx1 cluster assignments (integers 1..K)
%   K       : number of clusters
%   Q       : number of quantile bins
% Returns:
%   comp        : QxK matrix of cluster proportions per bin (rows sum to 1)
%   bin_edges   : (Q+1)x1 vector of bin edges (original scale)
% =========================================================================
function [comp, bin_edges] = compute_composition(x, labels, K, Q)
    % Compute quantile-based edges on the proximity variable
    probs      = linspace(0, 1, Q + 1);
    bin_edges  = quantile(x, probs);
    bin_edges(1)   = -Inf;
    bin_edges(end) = Inf;

    comp = zeros(Q, K);
    for q = 1:Q
        in_bin = x > bin_edges(q) & x <= bin_edges(q + 1);
        if ~any(in_bin), continue; end
        bin_labels = labels(in_bin);
        for k = 1:K
            comp(q, k) = sum(bin_labels == k) / sum(in_bin);
        end
    end
end


% =========================================================================
% LOCAL: BIC curve plot
% =========================================================================
function plot_bic(bic_values, k_optimal, k_max)
    figure('Name', 'GMM BIC Model Selection', ...
           'Position', [100 100 520 380], 'Color', 'w');

    k_range = 1:k_max;
    plot(k_range, bic_values, 'o-', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 2, 'MarkerSize', 7, ...
         'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on;
    plot(k_optimal, bic_values(k_optimal), 'o', ...
         'MarkerSize', 12, 'MarkerFaceColor', [0.85 0.20 0.20], ...
         'MarkerEdgeColor', [0.6 0.1 0.1], 'LineWidth', 1.5);
    hold off;

    xlabel('Number of GMM components (k)');
    ylabel('BIC');
    title('GMM Model Selection via BIC');
    subtitle(sprintf('Optimal k = %d  (lowest BIC)', k_optimal), ...
             'FontSize', 9, 'Color', [.4 .4 .4]);
    xticks(k_range);
    grid on; box off;
end


% =========================================================================
% LOCAL: 2D feature space scatter coloured by cluster
% =========================================================================
function plot_feature_space(X, labels, K, norm_params, cluster_stats, gmm,...
                          dirout, tstr)
    figure('Name', 'GMM Clusters in Feature Space', ...
           'Position', [150 150 620 500], 'Color', 'w');

    cmap = cluster_colormap(K);

    % Subsample for scatter (avoid overplotting millions of points)
    n_show  = min(15000, size(X, 1));
    idx_show = randperm(size(X, 1), n_show);

    hold on;
    for k = 1:K
        in_k   = labels(idx_show) == k;
        scatter(X(idx_show(in_k), 1), X(idx_show(in_k), 2), 4, ...
                cmap(k,:), 'filled', 'MarkerFaceAlpha', 0.25);
    end

    % Draw GMM ellipses (1 SD contour) for each component
    theta = linspace(0, 2*pi, 200);
    for k = 1:K
        mu_k  = gmm.mu(k, :);
        if ndims(gmm.Sigma) == 3
            S = gmm.Sigma(:,:,k);
        else
            S = diag(gmm.Sigma(:,k));
        end
        [V, D]  = eig(S);
        ellipse = V * sqrt(D) * [cos(theta); sin(theta)];
        plot(mu_k(1) + ellipse(1,:), mu_k(2) + ellipse(2,:), '-', ...
             'Color', cmap(k,:), 'LineWidth', 2);
        plot(mu_k(1), mu_k(2), '+', ...
             'Color', cmap(k,:) * 0.6, 'MarkerSize', 10, 'LineWidth', 2);
    end
    hold off;

    % Build legend labels with original-scale cluster means
    leg_entries = cell(K, 1);
    for k = 1:K
        leg_entries{k} = sprintf('C%d  (μ_s=%.3f, μ_r=%.3f)', k, ...
            cluster_stats.mean_scatter(k), cluster_stats.mean_retard(k));
    end
    legend(leg_entries, 'Location', 'best', 'FontSize', 8);

    xlabel(sprintf('Scattering (z-scored,  μ=%.3f, σ=%.3f)', ...
                   norm_params.mu_s, norm_params.sig_s));
    ylabel(sprintf('Retardance  (z-scored,  μ=%.3f, σ=%.3f)', ...
                   norm_params.mu_r, norm_params.sig_r));
    title('GMM Cluster Assignments — Feature Space');
    subtitle('Ellipses = 1 SD contour per component', ...
             'FontSize', 9, 'Color', [.4 .4 .4]);
    grid on; box off;
    
    % Top Title
    sgtitle(sprintf('GMM Feature Space %s',tstr),'FontSize', 12,...
            'FontWeight', 'bold');
    % Save figure as .FIG and .PNG
    fname = strcat('GMM_feature_space_',strrep(tstr,' ','_'));
    fout = fullfile(dirout, strcat(fname,'.fig'));
    saveas(gcf,fout);
    fout = fullfile(dirout, strcat(fname,'.png'));
    saveas(gcf,fout);
end

% =========================================================================
% LOCAL: stacked bar of cluster composition per quantile bin
% =========================================================================
function plot_composition(comp, bin_edges, K, prox_label, cluster_stats,...
                          dirout, tstr)
    Q    = size(comp, 1);
    cmap = cluster_colormap(K);

    % Midpoints of bins (original scale, ignoring ±Inf sentinels)
    finite_edges        = bin_edges;
    finite_edges(1)     = bin_edges(2) - (bin_edges(3) - bin_edges(2));
    finite_edges(end)   = bin_edges(end-1) + (bin_edges(end-1) - bin_edges(end-2));
    bin_mids = (finite_edges(1:end-1) + finite_edges(2:end)) / 2;

    figure('Name', sprintf('Cluster Composition by %s Quantile', prox_label), ...
           'Position', [200 200 800 420], 'Color', 'w');

    % --- Left panel: stacked bar ---
    ax1 = subplot(1, 2, 1);
    b   = bar(1:Q, comp * 100, 'stacked', 'EdgeColor', 'none');
    for k = 1:K
        b(k).FaceColor = cmap(k, :);
    end
    xticks(1:Q);
    xticklabels(arrayfun(@(v) sprintf('%.2f', v), bin_mids, ...
                         'UniformOutput', false));
    xtickangle(45);
    xlabel(sprintf('%s (bin midpoint)', prox_label));
    ylabel('Voxels in cluster (%)');
    title(sprintf('Cluster composition\nacross %s bins', prox_label));
    ylim([0 100]);
    grid on; box off;

    % Build legend
    leg_entries = cell(K, 1);
    for k = 1:K
        leg_entries{k} = sprintf('C%d (μ_s=%.3f, μ_r=%.3f)', k, ...
            cluster_stats.mean_scatter(k), cluster_stats.mean_retard(k));
    end
    legend(leg_entries, 'Location', 'eastoutside', 'FontSize', 8);

    % --- Right panel: line plot of each cluster's proportion ---
    ax2 = subplot(1, 2, 2);
    hold on;
    for k = 1:K
        plot(1:Q, comp(:, k) * 100, 'o-', ...
             'Color', cmap(k,:), 'LineWidth', 2, 'MarkerSize', 5, ...
             'MarkerFaceColor', cmap(k,:));
    end
    hold off;
    xticks(1:Q);
    xticklabels(arrayfun(@(v) sprintf('%.2f', v), bin_mids, ...
                         'UniformOutput', false));
    xtickangle(45);
    xlabel(sprintf('%s (bin midpoint)', prox_label));
    ylabel('Proportion (%)');
    title(sprintf('Cluster prevalence\nacross %s bins', prox_label));
    legend(leg_entries, 'Location', 'eastoutside', 'FontSize', 8);
    grid on; box off;

    linkaxes([ax1, ax2], 'x');

    sgtitle(sprintf(['Post-hoc Stratification by %s Quantile Bins' ...
                    '  (k=%d GMM). %s'], ...
                    prox_label, K, tstr), ...
            'FontSize', 12, 'FontWeight', 'bold');
    % Save figure as .FIG and .PNG
    fout = fullfile(dirout, strcat('GMM_',strrep(tstr,' ','_'),'.fig'));
    saveas(gcf,fout);
    fout = fullfile(dirout, strcat('GMM_',strrep(tstr,' ','_'),'.png'));
    saveas(gcf,fout);
end


% =========================================================================
% LOCAL: generate a distinguishable colormap for K clusters
% =========================================================================
function cmap = cluster_colormap(K)
    % Uses a hand-picked palette for K<=8, falls back to hsv for larger K
    base = [0.22 0.49 0.72;   % blue
            0.82 0.24 0.24;   % red
            0.20 0.63 0.37;   % green
            0.95 0.60 0.10;   % orange
            0.58 0.35 0.68;   % purple
            0.30 0.75 0.75;   % teal
            0.90 0.40 0.60;   % pink
            0.50 0.35 0.20];  % brown
    if K <= size(base, 1)
        cmap = base(1:K, :);
    else
        cmap = hsv(K);
    end
end