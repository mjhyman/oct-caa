function results = compare_isosurfaces_severe_control(T_severe, T_control, varargin)
% COMPARE_ISOSURFACES_SEVERE_CONTROL  Compare severe vs. control GAM slice
%   curves (y = optical property, x = EPVS SWP, fanned at the 10th/50th/
%   90th percentile of vessel SWP) on a SHARED, aligned grid so the two
%   groups are actually comparable point-for-point.
%
%   Two alignment choices matter and are made explicit rather than
%   implicit here:
%
%   1. EPVS SWP axis. Severe and control will generally have different
%      epv_swp ranges. By default this function evaluates BOTH fitted
%      GAMs only over the INTERSECTION of the two groups' ranges
%      ('GridRange','intersection'), so no plotted point is ever a
%      silent extrapolation beyond what a group's own data supports.
%      'GridRange','union' extends both curves across the full combined
%      range instead; any region outside a group's own observed range is
%      genuine model extrapolation for that group and should be reported
%      as such if used (the density mask will flag it, but a reader
%      skimming the figure will not know that without a caption note).
%
%   2. Vessel SWP fanning levels. By default the 10th/50th/90th
%      percentile levels are computed from the POOLED (severe + control)
%      ves_swp distribution ('SharedPercentiles', true), so "p50 vessel
%      SWP" is the same actual ves_swp value in both curves. Setting
%      'SharedPercentiles', false computes percentiles independently per
%      group -- this compares each group's own low/median/high vessel
%      proximity, which is a different (and easily misread) comparison
%      since "p50" would then correspond to different absolute ves_swp
%      values in each group.
%
%   Each group's GAM complexity (MaxSplits/NumTrees) is tuned/fit
%   independently -- there's no reason the optimal complexity should
%   match across two datasets with different N and dispersion. Only the
%   prediction grid and (by default) the slice levels are shared.
%
% REQUIRED
%   T_severe, T_control : tables with 'scattering','retardance','ves_swp','epv_swp'
%
% NAME-VALUE PAIRS
%   'TuneParams'        Auto-tune each group's MaxSplits/NumTrees      (default: true)
%   'MaxSplitsGrid'     Grid searched per group if tuning               (default: [5 10 20 40])
%   'NumTreesGrid'      Grid searched per group if tuning               (default: [50 100 200])
%   'KFold'             CV folds used for tuning                        (default: 5)
%   'MaxSplits'         [severe, control] fixed values if not tuning    (default: [10 10])
%   'NumTrees'          [severe, control] fixed values if not tuning    (default: [100 100])
%   'LearnRate'         fitrgam learning rate (both groups)              (default: 0.1)
%   'MaxSample'         Max voxels per group used for fitting            (default: 1e5)
%   'NumGridPts'        Points along the shared EPVS SWP axis            (default: 100)
%   'GridRange'         'intersection' or 'union' of epv_swp ranges      (default: 'intersection')
%   'SliceQuantiles'    Vessel SWP percentiles for fanning                (default: [0.10 0.50 0.90])
%   'SharedPercentiles' Use pooled ves_swp for slice levels               (default: true)
%   'MinDensity'        Min voxel count per slice bin to trust a curve    (default: 5)
%   'BootstrapCI'       Compute bootstrap CIs on slice curves             (default: true)
%   'NBootstrap'        Bootstrap resamples                               (default: 100)
%   'CIAlpha'           CI alpha                                          (default: 0.05)
%   'SevereLabel'       Legend label for the severe group                 (default: 'Severe')
%   'ControlLabel'      Legend label for the control group                (default: 'Control')
%   'TitleStr'          Label for titles / filenames                     (default: 'Severe_vs_Control')
%   'dirout'            If non-empty, saves figures + tuning CSVs         (default: '')
%   'Verbose'           Print progress                                    (default: true)
%
% OUTPUT results struct:
%   .severe / .control    each: .gam_pair, .tuning, .ves_levels,
%                          .slice_scatter, .slice_retard, .density_mask,
%                          .ci_scatter, .ci_retard, .T_fit, .pct_masked
%   .epv_grid              shared EPVS SWP axis used for both groups
%   .grid_range_used        [lo hi] actually used (post intersection/union)
%   .slice_quantiles         quantile fractions used
%   .shared_percentiles_used logical, whether ves_levels were pooled
%   .fig                     comparison figure handle
%
% DEPENDENCIES
%   Requires gam_slice_toolkit.m and tune_swp_psoct_gam.m on the path.

    p = inputParser();
    p.FunctionName = 'compare_isosurfaces_severe_control';

    addRequired(p,  'T_severe',          @istable);
    addRequired(p,  'T_control',         @istable);
    addParameter(p, 'TuneParams',        true,             @islogical);
    addParameter(p, 'MaxSplitsGrid',     [5 10 20 40],     @isnumeric);
    addParameter(p, 'NumTreesGrid',      [50 100 200],     @isnumeric);
    addParameter(p, 'KFold',             5,                @(x) isnumeric(x) && x >= 2);
    addParameter(p, 'MaxSplits',         [10 10],          @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'NumTrees',          [100 100],        @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'LearnRate',         0.1,              @(x) isnumeric(x) && x > 0 && x <= 1);
    addParameter(p, 'MaxSample',         1e5,              @(x) isnumeric(x) && x > 0);
    addParameter(p, 'NumGridPts',        100,              @(x) isnumeric(x) && x > 1);
    addParameter(p, 'GridRange',         'intersection',   @(x) any(strcmpi(x, {'intersection','union'})));
    addParameter(p, 'SliceQuantiles',    [0.10 0.50 0.90], @isnumeric);
    addParameter(p, 'SharedPercentiles', true,             @islogical);
    addParameter(p, 'MinDensity',        5,                @(x) isnumeric(x) && x >= 0);
    addParameter(p, 'BootstrapCI',       true,             @islogical);
    addParameter(p, 'NBootstrap',        100,              @(x) isnumeric(x) && x >= 10);
    addParameter(p, 'CIAlpha',           0.05,             @(x) isnumeric(x) && x > 0 && x < 1);
    addParameter(p, 'SevereLabel',       'Severe');
    addParameter(p, 'ControlLabel',      'Control');
    addParameter(p, 'TitleStr',          'Severe_vs_Control');
    addParameter(p, 'dirout',            '');
    addParameter(p, 'Verbose',           true,             @islogical);

    parse(p, T_severe, T_control, varargin{:});
    opts = p.Results;

    toolkit = gam_slice_toolkit();

    % ---------------------------------------------------------------
    % Clean each group
    % ---------------------------------------------------------------
    Tw_severe  = toolkit.clean_and_subsample(T_severe,  opts.MaxSample, opts.Verbose);
    Tw_control = toolkit.clean_and_subsample(T_control, opts.MaxSample, opts.Verbose);

    % ---------------------------------------------------------------
    % Tune (or use fixed) complexity per group
    % ---------------------------------------------------------------
    if opts.TuneParams
        if opts.Verbose
            fprintf('[compare_isosurfaces_severe_control] Tuning severe GAM...\n');
        end
        tuning_severe = tune_swp_psoct_gam(Tw_severe, ...
            'MaxSplitsGrid', opts.MaxSplitsGrid, 'NumTreesGrid', opts.NumTreesGrid, ...
            'KFold', opts.KFold, 'MaxSample', opts.MaxSample, ...
            'TitleStr', [opts.TitleStr '_severe'], 'dirout', opts.dirout, 'Verbose', opts.Verbose);
        MaxSplits_severe = tuning_severe.recommended.MaxSplits;
        NumTrees_severe  = tuning_severe.recommended.NumTrees;

        if opts.Verbose
            fprintf('[compare_isosurfaces_severe_control] Tuning control GAM...\n');
        end
        tuning_control = tune_swp_psoct_gam(Tw_control, ...
            'MaxSplitsGrid', opts.MaxSplitsGrid, 'NumTreesGrid', opts.NumTreesGrid, ...
            'KFold', opts.KFold, 'MaxSample', opts.MaxSample, ...
            'TitleStr', [opts.TitleStr '_control'], 'dirout', opts.dirout, 'Verbose', opts.Verbose);
        MaxSplits_control = tuning_control.recommended.MaxSplits;
        NumTrees_control  = tuning_control.recommended.NumTrees;
    else
        tuning_severe = []; tuning_control = [];
        MaxSplits_severe  = opts.MaxSplits(1); NumTrees_severe  = opts.NumTrees(1);
        MaxSplits_control = opts.MaxSplits(2); NumTrees_control = opts.NumTrees(2);
    end

    % ---------------------------------------------------------------
    % Fit final GAMs
    % ---------------------------------------------------------------
    gam_severe  = toolkit.fit_gam_pair(Tw_severe,  MaxSplits_severe,  NumTrees_severe,  opts.LearnRate);
    gam_control = toolkit.fit_gam_pair(Tw_control, MaxSplits_control, NumTrees_control, opts.LearnRate);

    % ---------------------------------------------------------------
    % Shared EPVS SWP grid
    % ---------------------------------------------------------------
    range_severe  = [min(Tw_severe.epv_swp),  max(Tw_severe.epv_swp)];
    range_control = [min(Tw_control.epv_swp), max(Tw_control.epv_swp)];

    switch lower(opts.GridRange)
        case 'intersection'
            lo = max(range_severe(1), range_control(1));
            hi = min(range_severe(2), range_control(2));
            if lo >= hi
                error('compare_isosurfaces_severe_control:noOverlap', ...
                    ['Severe and control epv_swp ranges do not overlap ', ...
                     '(severe [%.4g %.4g] vs control [%.4g %.4g]). Cannot build a shared ', ...
                     'grid without extrapolation. Consider ''GridRange'',''union'' and report ', ...
                     'the extrapolated regions explicitly.'], ...
                    range_severe(1), range_severe(2), range_control(1), range_control(2));
            end
        case 'union'
            lo = min(range_severe(1), range_control(1));
            hi = max(range_severe(2), range_control(2));
    end
    epv_grid = linspace(lo, hi, opts.NumGridPts)';

    if opts.Verbose
        fprintf(['[compare_isosurfaces_severe_control] Shared EPVS SWP grid (%s): ', ...
                 '[%.4g, %.4g]  (severe range [%.4g %.4g], control range [%.4g %.4g])\n'], ...
                opts.GridRange, lo, hi, range_severe(1), range_severe(2), ...
                range_control(1), range_control(2));
    end

    % ---------------------------------------------------------------
    % Vessel SWP slice levels
    % ---------------------------------------------------------------
    if opts.SharedPercentiles
        pooled_ves        = [Tw_severe.ves_swp; Tw_control.ves_swp];
        ves_levels         = quantile(pooled_ves, opts.SliceQuantiles);
        ves_levels_severe  = ves_levels;
        ves_levels_control = ves_levels;
        if opts.Verbose
            fprintf('[compare_isosurfaces_severe_control] Pooled vessel SWP slice levels: %s\n', ...
                    mat2str(ves_levels, 4));
        end
    else
        ves_levels_severe  = quantile(Tw_severe.ves_swp,  opts.SliceQuantiles);
        ves_levels_control = quantile(Tw_control.ves_swp, opts.SliceQuantiles);
        if opts.Verbose
            fprintf(['[compare_isosurfaces_severe_control] WARNING: SharedPercentiles=false -- ', ...
                     'severe slice levels %s and control slice levels %s correspond to the ', ...
                     'same percentile label but DIFFERENT absolute vessel SWP values.\n'], ...
                    mat2str(ves_levels_severe, 4), mat2str(ves_levels_control, 4));
        end
    end

    % ---------------------------------------------------------------
    % Predict slices for each group on the shared epv_grid
    % ---------------------------------------------------------------
    [slice_scatter_severe,  slice_retard_severe]  = toolkit.predict_slices(gam_severe,  epv_grid, ves_levels_severe);
    [slice_scatter_control, slice_retard_control] = toolkit.predict_slices(gam_control, epv_grid, ves_levels_control);

    mask_severe  = toolkit.slice_density_mask(Tw_severe,  epv_grid, ves_levels_severe,  opts.MinDensity);
    mask_control = toolkit.slice_density_mask(Tw_control, epv_grid, ves_levels_control, opts.MinDensity);

    pct_masked_severe  = 100 * mean(mask_severe(:));
    pct_masked_control = 100 * mean(mask_control(:));
    if opts.Verbose
        fprintf(['[compare_isosurfaces_severe_control] Slice-point coverage below MinDensity=%d: ', ...
                 'severe %.1f%%, control %.1f%%\n'], opts.MinDensity, pct_masked_severe, pct_masked_control);
        if abs(pct_masked_severe - pct_masked_control) > 20
            fprintf(['[compare_isosurfaces_severe_control] WARNING: coverage is asymmetric between ', ...
                     'groups by more than 20 percentage points. Regions where one group''s curve is ', ...
                     'masked and the other''s is not cannot be compared -- that gap reflects data ', ...
                     'density, not necessarily biology.\n']);
        end
    end

    % ---------------------------------------------------------------
    % Bootstrap CIs per group (optional)
    % ---------------------------------------------------------------
    ci_scatter_severe = {}; ci_retard_severe = {};
    ci_scatter_control = {}; ci_retard_control = {};
    if opts.BootstrapCI
        [ci_scatter_severe, ci_retard_severe] = toolkit.bootstrap_slices(Tw_severe, ...
            MaxSplits_severe, NumTrees_severe, opts.LearnRate, epv_grid, ves_levels_severe, ...
            opts.NBootstrap, opts.CIAlpha, opts.Verbose);
        [ci_scatter_control, ci_retard_control] = toolkit.bootstrap_slices(Tw_control, ...
            MaxSplits_control, NumTrees_control, opts.LearnRate, epv_grid, ves_levels_control, ...
            opts.NBootstrap, opts.CIAlpha, opts.Verbose);
    end

    % ---------------------------------------------------------------
    % Plot comparison
    % ---------------------------------------------------------------
    fig = toolkit.plot_slices_compare(epv_grid, ...
        ves_levels_severe,  slice_scatter_severe,  slice_retard_severe,  mask_severe,  ci_scatter_severe,  ci_retard_severe,  opts.SevereLabel, ...
        ves_levels_control, slice_scatter_control, slice_retard_control, mask_control, ci_scatter_control, ci_retard_control, opts.ControlLabel, ...
        opts.SliceQuantiles, opts.TitleStr, opts.dirout);

    % ---------------------------------------------------------------
    % Package
    % ---------------------------------------------------------------
    results.severe.gam_pair      = gam_severe;
    results.severe.tuning        = tuning_severe;
    results.severe.ves_levels    = ves_levels_severe;
    results.severe.slice_scatter = slice_scatter_severe;
    results.severe.slice_retard  = slice_retard_severe;
    results.severe.density_mask  = mask_severe;
    results.severe.ci_scatter    = ci_scatter_severe;
    results.severe.ci_retard     = ci_retard_severe;
    results.severe.T_fit         = Tw_severe;
    results.severe.pct_masked    = pct_masked_severe;

    results.control.gam_pair      = gam_control;
    results.control.tuning        = tuning_control;
    results.control.ves_levels    = ves_levels_control;
    results.control.slice_scatter = slice_scatter_control;
    results.control.slice_retard  = slice_retard_control;
    results.control.density_mask  = mask_control;
    results.control.ci_scatter    = ci_scatter_control;
    results.control.ci_retard     = ci_retard_control;
    results.control.T_fit         = Tw_control;
    results.control.pct_masked    = pct_masked_control;

    results.epv_grid                 = epv_grid;
    results.grid_range_used          = [lo, hi];
    results.slice_quantiles          = opts.SliceQuantiles;
    results.shared_percentiles_used  = opts.SharedPercentiles;
    results.fig                      = fig;

    if opts.Verbose
        fprintf('[compare_isosurfaces_severe_control] Done: "%s".\n', opts.TitleStr);
    end
end
