function results = fit_isosurface_by_dataset(T, varargin)
% FIT_ISOSURFACE_BY_DATASET  Fit a tuned GAM to ONE dataset (one subject /
%   stage / region combination) and plot scattering & retardance slice
%   curves: y = optical property, x = EPVS SWP, fanned at the 10th/50th/
%   90th percentile of vessel SWP.
%
%   GAM complexity (MaxSplits, NumTrees) is auto-tuned per dataset via
%   tune_swp_psoct_gam (selected by cross-validated R^2) unless disabled.
%   Run once per subject/stage/region; loop externally across your set of
%   datasets to build a panel of these.
%
% REQUIRED
%   T : table with 'scattering', 'retardance', 'ves_swp', 'epv_swp'
%
% NAME-VALUE PAIRS
%   'TuneParams'      Auto-tune MaxSplits/NumTrees via CV R^2       (default: true)
%   'MaxSplitsGrid'   Grid searched if TuneParams=true               (default: [5 10 20 40])
%   'NumTreesGrid'    Grid searched if TuneParams=true               (default: [50 100 200])
%   'KFold'           CV folds used for tuning                       (default: 5)
%   'MaxSplits'       Fixed value if TuneParams=false                (default: 10)
%   'NumTrees'        Fixed value if TuneParams=false                (default: 100)
%   'LearnRate'       fitrgam learning rate                          (default: 0.1)
%   'MaxSample'       Max voxels used for fitting                    (default: 1e5)
%   'NumGridPts'      Points along the EPVS SWP axis                 (default: 100)
%   'SliceQuantiles'  Vessel SWP percentiles for fanning              (default: [0.10 0.50 0.90])
%   'MinDensity'      Min voxel count per slice bin to trust a curve  (default: 5)
%   'BootstrapCI'     Compute bootstrap CIs on slice curves           (default: true)
%   'NBootstrap'      Bootstrap resamples                             (default: 100)
%   'CIAlpha'         CI alpha                                        (default: 0.05)
%   'TitleStr'        Label for titles / filenames                    (default: '')
%   'dirout'          If non-empty, saves figures + tuning CSV        (default: '')
%   'Verbose'         Print progress                                   (default: true)
%
% OUTPUT results struct:
%   .gam_pair          fitted scattering/retardance GAMs + complexity used
%   .tuning            output of tune_swp_psoct_gam (or [] if TuneParams=false)
%   .epv_grid          EPVS SWP axis used for the slice curves
%   .ves_levels        vessel SWP values held fixed for each slice
%   .slice_quantiles   quantile fractions used
%   .slice_scatter     cell{n_slices,1} scattering ~ epv_swp curves
%   .slice_retard      cell{n_slices,1} retardance ~ epv_swp curves
%   .density_mask      n_pts x n_slices logical, true = under MinDensity
%   .ci_scatter        cell{n_slices,1} Nx2 bootstrap CI (or {} if disabled)
%   .ci_retard         cell{n_slices,1} Nx2 bootstrap CI (or {} if disabled)
%   .fig               figure handle
%   .T_fit             (possibly subsampled) table used for fitting
%
% NOTE ON COMPARABILITY
%   This function builds its EPVS SWP grid from THIS dataset's own
%   min/max, so its x-axis is not automatically aligned with any other
%   dataset's grid -- two calls to this function are not directly
%   comparable/overlayable. For a severe-vs-control comparison on a
%   shared, aligned grid, use compare_isosurfaces_severe_control instead.
%
% DEPENDENCIES
%   Requires gam_slice_toolkit.m and tune_swp_psoct_gam.m on the path.

    p = inputParser();
    p.FunctionName = 'fit_isosurface_by_dataset';

    addRequired(p,  'T',              @istable);
    addParameter(p, 'TuneParams',     true,           @islogical);
    addParameter(p, 'MaxSplitsGrid',  [5 10 20 40],   @isnumeric);
    addParameter(p, 'NumTreesGrid',   [50 100 200],   @isnumeric);
    addParameter(p, 'KFold',          5,              @(x) isnumeric(x) && x >= 2);
    addParameter(p, 'MaxSplits',      10,             @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NumTrees',       100,            @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'LearnRate',      0.1,            @(x) isnumeric(x) && x > 0 && x <= 1);
    addParameter(p, 'MaxSample',      1e5,            @(x) isnumeric(x) && x > 0);
    addParameter(p, 'NumGridPts',     100,            @(x) isnumeric(x) && x > 1);
    addParameter(p, 'SliceQuantiles', [0.10 0.50 0.90], @isnumeric);
    addParameter(p, 'MinDensity',     5,              @(x) isnumeric(x) && x >= 0);
    addParameter(p, 'BootstrapCI',    true,           @islogical);
    addParameter(p, 'NBootstrap',     100,            @(x) isnumeric(x) && x >= 10);
    addParameter(p, 'CIAlpha',        0.05,           @(x) isnumeric(x) && x > 0 && x < 1);
    addParameter(p, 'TitleStr',       '');
    addParameter(p, 'dirout',         '');
    addParameter(p, 'Verbose',        true,           @islogical);

    parse(p, T, varargin{:});
    opts = p.Results;

    toolkit = gam_slice_toolkit();

    T_work = toolkit.clean_and_subsample(T, opts.MaxSample, opts.Verbose);

    % ---------------------------------------------------------------
    % Tune (or use fixed) GAM complexity
    % ---------------------------------------------------------------
    if opts.TuneParams
        if opts.Verbose
            fprintf('[fit_isosurface_by_dataset] Tuning GAM complexity for "%s"...\n', opts.TitleStr);
        end
        tuning = tune_swp_psoct_gam(T_work, ...
            'MaxSplitsGrid', opts.MaxSplitsGrid, ...
            'NumTreesGrid',  opts.NumTreesGrid,  ...
            'KFold',         opts.KFold,         ...
            'MaxSample',     opts.MaxSample,     ...
            'TitleStr',      opts.TitleStr,      ...
            'dirout',        opts.dirout,        ...
            'Verbose',       opts.Verbose);
        MaxSplits = tuning.recommended.MaxSplits;
        NumTrees  = tuning.recommended.NumTrees;
    else
        tuning    = [];
        MaxSplits = opts.MaxSplits;
        NumTrees  = opts.NumTrees;
    end

    % ---------------------------------------------------------------
    % Fit final GAM pair at the chosen complexity
    % ---------------------------------------------------------------
    gam_pair = toolkit.fit_gam_pair(T_work, MaxSplits, NumTrees, opts.LearnRate);

    % ---------------------------------------------------------------
    % Grid + slice levels (this dataset's own range/percentiles)
    % ---------------------------------------------------------------
    epv_grid   = linspace(min(T_work.epv_swp), max(T_work.epv_swp), opts.NumGridPts)';
    ves_levels = quantile(T_work.ves_swp, opts.SliceQuantiles);

    [slice_scatter, slice_retard] = toolkit.predict_slices(gam_pair, epv_grid, ves_levels);
    density_mask = toolkit.slice_density_mask(T_work, epv_grid, ves_levels, opts.MinDensity);

    if opts.Verbose
        fprintf('[fit_isosurface_by_dataset] %.1f%% of slice points below MinDensity=%d.\n', ...
                100 * mean(density_mask(:)), opts.MinDensity);
    end

    % ---------------------------------------------------------------
    % Bootstrap CIs (optional)
    % ---------------------------------------------------------------
    ci_scatter = {};
    ci_retard  = {};
    if opts.BootstrapCI
        [ci_scatter, ci_retard] = toolkit.bootstrap_slices(T_work, MaxSplits, NumTrees, ...
            opts.LearnRate, epv_grid, ves_levels, opts.NBootstrap, opts.CIAlpha, opts.Verbose);
    end

    % ---------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------
    fig = toolkit.plot_slices_single(epv_grid, ves_levels, opts.SliceQuantiles, ...
        slice_scatter, slice_retard, density_mask, ci_scatter, ci_retard, ...
        opts.TitleStr, opts.dirout);

    % ---------------------------------------------------------------
    % Package
    % ---------------------------------------------------------------
    results.gam_pair        = gam_pair;
    results.tuning          = tuning;
    results.epv_grid        = epv_grid;
    results.ves_levels      = ves_levels;
    results.slice_quantiles = opts.SliceQuantiles;
    results.slice_scatter   = slice_scatter;
    results.slice_retard    = slice_retard;
    results.density_mask    = density_mask;
    results.ci_scatter      = ci_scatter;
    results.ci_retard       = ci_retard;
    results.fig             = fig;
    results.T_fit           = T_work;

    if opts.Verbose
        fprintf('[fit_isosurface_by_dataset] Done: "%s".\n', opts.TitleStr);
    end
end
