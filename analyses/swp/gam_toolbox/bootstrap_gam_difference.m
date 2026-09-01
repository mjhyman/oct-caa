function D = bootstrap_gam_difference(res, varargin)
% BOOTSTRAP_GAM_DIFFERENCE  Compute the severe-minus-control slice-curve
%   difference and a properly resampled confidence band for it, from a
%   compare_isosurfaces_severe_control results struct (one region).
%
%   The stored per-group CIs (res.severe.ci_scatter etc.) are independent
%   percentile bands and CANNOT be combined post hoc into a difference
%   band, because only their percentiles were kept, not the bootstrap
%   draws. This function therefore recomputes the difference band directly:
%   on each resample it draws severe and control voxels independently (they
%   are independent samples), refits both groups' GAMs at the SAME
%   complexity used originally, predicts each group's slice curve on the
%   shared epv_grid, and takes the difference. Percentiles across resamples
%   give a band that correctly propagates BOTH groups' sampling variability
%   into severe - control.
%
%   The difference is only meaningful when the two groups' slices are held
%   at the SAME absolute vessel-SWP values, i.e. the comparison was run
%   with SharedPercentiles = true (the default). This is asserted.
%
% REQUIRED
%   res : one region's struct from compare_isosurfaces_severe_control, i.e.
%         gam_cmp.(d).(reg). Must contain:
%           .epv_grid, .slice_quantiles, .shared_percentiles_used
%           .severe.T_fit, .severe.gam_pair (.MaxSplits/.NumTrees/.LearnRate),
%           .severe.ves_levels, .severe.slice_scatter, .severe.slice_retard,
%           .severe.density_mask   (and the same .control fields)
%
% NAME-VALUE
%   'NBootstrap' resamples for the difference band          (default: 500)
%   'CIAlpha'    alpha for the percentile band              (default: 0.05)
%   'Verbose'    print progress                             (default: true)
%
% OUTPUT  struct D:
%   .epv_grid           shared EPVS-SWP axis (from res)
%   .ves_levels         vessel-SWP slice levels (shared/pooled)
%   .slice_quantiles    quantile fractions
%   .diff_scatter       cell{n_slices,1}: severe - control scattering curve
%   .diff_retard        cell{n_slices,1}: severe - control retardance curve
%   .ci_diff_scatter    cell{n_slices,1}: Nx2 resampled difference band
%   .ci_diff_retard     cell{n_slices,1}: Nx2 resampled difference band
%   .mask               n_pts x n_slices logical: true where EITHER group
%                       was below MinDensity (difference unreliable there)
%   .nbootstrap, .ci_alpha
%
% NOTE
%   Point estimates reuse the stored full-fit curves
%   (res.severe.slice_* - res.control.slice_*), so the plotted difference
%   line is exactly severe_curve - control_curve from the original GAMs;
%   only the band is (re)bootstrapped.

    p = inputParser();
    p.FunctionName = 'bootstrap_gam_difference';
    addRequired(p,  'res',        @isstruct);
    addParameter(p, 'NBootstrap', 500,  @(x) isnumeric(x) && isscalar(x) && x >= 10);
    addParameter(p, 'CIAlpha',    0.05, @(x) isnumeric(x) && x > 0 && x < 1);
    addParameter(p, 'Verbose',    true, @islogical);
    parse(p, res, varargin{:});
    opts = p.Results;

    % -----------------------------------------------------------------
    % Validate that a difference is meaningful for this comparison
    % -----------------------------------------------------------------
    if isfield(res, 'shared_percentiles_used') && ~res.shared_percentiles_used
        error('bootstrap_gam_difference:notShared', ...
            ['This comparison was run with SharedPercentiles = false, so ' ...
             'severe and control slice k sit at DIFFERENT absolute vessel-SWP ' ...
             'values. A severe - control difference per slice is not meaningful. ' ...
             'Re-run the comparison with SharedPercentiles = true.']);
    end

    ves_levels = res.severe.ves_levels(:);
    if ~isequal(numel(ves_levels), numel(res.control.ves_levels)) || ...
       max(abs(ves_levels - res.control.ves_levels(:))) > 1e-9 * max(1, max(abs(ves_levels)))
        error('bootstrap_gam_difference:vesMismatch', ...
            ['Severe and control ves_levels differ; cannot difference slice-by-slice. ' ...
             'This should not happen when SharedPercentiles = true.']);
    end

    epv_grid = res.epv_grid(:);
    n_pts    = numel(epv_grid);
    n_slices = numel(ves_levels);
    squant   = res.slice_quantiles;

    % -----------------------------------------------------------------
    % Point-estimate difference curves from the stored full-fit curves
    % -----------------------------------------------------------------
    diff_scatter = cell(n_slices, 1);
    diff_retard  = cell(n_slices, 1);
    for sl = 1:n_slices
        diff_scatter{sl} = res.severe.slice_scatter{sl} - res.control.slice_scatter{sl};
        diff_retard{sl}  = res.severe.slice_retard{sl}  - res.control.slice_retard{sl};
    end

    % Combined density mask: unreliable where EITHER group is sparse
    mask = res.severe.density_mask | res.control.density_mask;

    % -----------------------------------------------------------------
    % Refit ingredients
    % -----------------------------------------------------------------
    Ts = res.severe.T_fit;   ns = height(Ts);
    Tc = res.control.T_fit;  nc = height(Tc);

    sp = res.severe.gam_pair;   % complexity for severe
    cp = res.control.gam_pair;  % complexity for control

    f_scatter = 'scattering ~ ves_swp + epv_swp + ves_epv_swp';
    f_retard  = 'retardance ~ ves_swp + epv_swp + ves_epv_swp';

    % Pre-build slice prediction tables on the shared grid (same for both
    % groups, since ves_levels is shared)
    T_slices = cell(n_slices, 1);
    for sl = 1:n_slices
        ves_fixed    = ves_levels(sl) * ones(n_pts, 1);
        T_slices{sl} = table(ves_fixed, epv_grid, ves_fixed .* epv_grid, ...
            'VariableNames', {'ves_swp', 'epv_swp', 'ves_epv_swp'});
    end

    B = opts.NBootstrap;
    boot_diff_scatter = zeros(n_pts, B, n_slices);
    boot_diff_retard  = zeros(n_pts, B, n_slices);

    use_par = ~isempty(ver('parallel'));
    if opts.Verbose
        if use_par
            fprintf('[bootstrap_gam_difference] %d paired resamples (parallel)...\n', B);
        else
            fprintf('[bootstrap_gam_difference] %d paired resamples...\n', B);
        end
    end

    if use_par
        parfor b = 1:B
            [ds, dr] = one_diff_resample(Ts, ns, Tc, nc, sp, cp, ...
                                         f_scatter, f_retard, T_slices, n_pts, n_slices);
            boot_diff_scatter(:, b, :) = ds;
            boot_diff_retard(:,  b, :) = dr;
        end
    else
        for b = 1:B
            if opts.Verbose && mod(b, 25) == 0
                fprintf('[bootstrap_gam_difference]   resample %d / %d\n', b, B);
            end
            [ds, dr] = one_diff_resample(Ts, ns, Tc, nc, sp, cp, ...
                                         f_scatter, f_retard, T_slices, n_pts, n_slices);
            boot_diff_scatter(:, b, :) = ds;
            boot_diff_retard(:,  b, :) = dr;
        end
    end

    % Percentile bands per slice
    alpha = opts.CIAlpha;
    lo = 100 * alpha / 2;
    hi = 100 * (1 - alpha / 2);

    ci_diff_scatter = cell(n_slices, 1);
    ci_diff_retard  = cell(n_slices, 1);
    for sl = 1:n_slices
        ci_diff_scatter{sl} = prctile(boot_diff_scatter(:, :, sl), [lo, hi], 2);
        ci_diff_retard{sl}  = prctile(boot_diff_retard(:,  :, sl), [lo, hi], 2);
    end

    D.epv_grid        = epv_grid;
    D.ves_levels      = ves_levels;
    D.slice_quantiles = squant;
    D.diff_scatter    = diff_scatter;
    D.diff_retard     = diff_retard;
    D.ci_diff_scatter = ci_diff_scatter;
    D.ci_diff_retard  = ci_diff_retard;
    D.mask            = mask;
    D.nbootstrap      = B;
    D.ci_alpha        = alpha;

    if opts.Verbose
        fprintf('[bootstrap_gam_difference] Done (alpha = %.2f).\n', alpha);
    end
end


% =========================================================================
% One paired resample: draw severe and control independently, refit both,
% return per-slice difference predictions (n_pts x n_slices) for each metric.
% =========================================================================
function [ds, dr] = one_diff_resample(Ts, ns, Tc, nc, sp, cp, ...
                                      f_scatter, f_retard, T_slices, n_pts, n_slices)
    Tsb = Ts(randi(ns, ns, 1), :);
    Tcb = Tc(randi(nc, nc, 1), :);

    sev_s = fitrgam(Tsb, f_scatter, ...
        'NumTreesPerPredictor', sp.NumTrees, 'MaxNumSplitsPerPredictor', sp.MaxSplits, ...
        'InitialLearnRateForPredictors', sp.LearnRate);
    sev_r = fitrgam(Tsb, f_retard, ...
        'NumTreesPerPredictor', sp.NumTrees, 'MaxNumSplitsPerPredictor', sp.MaxSplits, ...
        'InitialLearnRateForPredictors', sp.LearnRate);
    ctl_s = fitrgam(Tcb, f_scatter, ...
        'NumTreesPerPredictor', cp.NumTrees, 'MaxNumSplitsPerPredictor', cp.MaxSplits, ...
        'InitialLearnRateForPredictors', cp.LearnRate);
    ctl_r = fitrgam(Tcb, f_retard, ...
        'NumTreesPerPredictor', cp.NumTrees, 'MaxNumSplitsPerPredictor', cp.MaxSplits, ...
        'InitialLearnRateForPredictors', cp.LearnRate);

    ds = zeros(n_pts, n_slices);
    dr = zeros(n_pts, n_slices);
    for sl = 1:n_slices
        ds(:, sl) = predict(sev_s, T_slices{sl}) - predict(ctl_s, T_slices{sl});
        dr(:, sl) = predict(sev_r, T_slices{sl}) - predict(ctl_r, T_slices{sl});
    end
end