% TEST_GAM_SLICE_PIPELINE  Synthetic-data smoke test for:
%   tune_swp_psoct_gam.m, fit_isosurface_by_dataset.m,
%   compare_isosurfaces_severe_control.m (+ gam_slice_toolkit.m)
%
% Generates fake severe/control voxel tables with a known, deliberately
% simple structure (linear trend + an interaction bump + noise), then runs
% every code path: tuning sweep, single-dataset slice fitting, and the
% severe-vs-control aligned-grid comparison, with bootstrap CIs on.
%
% This checks that the pipeline RUNS end-to-end without error and that
% outputs have sane shapes -- it does not and cannot validate the real
% biological structure of your actual data.

close all;
clear;
rng(42);
Verbose = true;

fprintf('\n========================================\n');
fprintf('STEP 0: Generate synthetic data\n');
fprintf('========================================\n');

n_control = 3000;
n_severe  = 6000;   % severe: more tissue/voxels, wider proximity ranges

% Control: lower vessel/EPVS burden, narrower ranges
ves_swp_control = abs(0.5 + 0.3 * randn(n_control, 1));
epv_swp_control = abs(0.4 + 0.25 * randn(n_control, 1));

% Severe: higher burden, wider spread (overlapping but broader than control)
ves_swp_severe = abs(0.9 + 0.5 * randn(n_severe, 1));
epv_swp_severe = abs(0.8 + 0.45 * randn(n_severe, 1));

gen_optical = @(ves, epv) deal( ...
    1.0 + 0.8 * epv + 0.6 * ves .* epv + 0.05 * randn(size(ves)), ...   % scattering
    2.0 - 0.5 * epv + 0.3 * ves .* epv + 0.05 * randn(size(ves)));      % retardance

[scattering_control, retardance_control] = gen_optical(ves_swp_control, epv_swp_control);
[scattering_severe,  retardance_severe]  = gen_optical(ves_swp_severe,  epv_swp_severe);

T_control = table(ves_swp_control, epv_swp_control, scattering_control, retardance_control, ...
    'VariableNames', {'ves_swp', 'epv_swp', 'scattering', 'retardance'});
T_severe = table(ves_swp_severe, epv_swp_severe, scattering_severe, retardance_severe, ...
    'VariableNames', {'ves_swp', 'epv_swp', 'scattering', 'retardance'});

fprintf('  control: N=%d, ves_swp range [%.3f %.3f], epv_swp range [%.3f %.3f]\n', ...
        height(T_control), min(T_control.ves_swp), max(T_control.ves_swp), ...
        min(T_control.epv_swp), max(T_control.epv_swp));
fprintf('  severe:  N=%d, ves_swp range [%.3f %.3f], epv_swp range [%.3f %.3f]\n', ...
        height(T_severe), min(T_severe.ves_swp), max(T_severe.ves_swp), ...
        min(T_severe.epv_swp), max(T_severe.epv_swp));

dirout = '/projectnb/npbssmic/ns/CAA/gam_tests/';

pass = true;

fprintf('\n========================================\n');
fprintf('STEP 1: tune_swp_psoct_gam\n');
fprintf('========================================\n');
try
    tuning = tune_swp_psoct_gam(T_severe, ...
        'MaxSplitsGrid', [5 10], 'NumTreesGrid', [50 100], 'KFold', 3, ...
        'TitleStr', 'SmokeTest_Severe', 'dirout', dirout, 'Verbose', Verbose);

    assert(isfield(tuning, 'recommended'), 'tuning.recommended missing');
    assert(isfield(tuning.recommended, 'MaxSplits'), 'recommended.MaxSplits missing');
    assert(isfield(tuning.recommended, 'NumTrees'), 'recommended.NumTrees missing');
    assert(istable(tuning.complexity_table), 'complexity_table should be a table');
    assert(istable(tuning.density_table), 'density_table should be a table');
    assert(height(tuning.complexity_table) == 4, 'complexity_table should have 2x2=4 rows for the test grid');
    fprintf('  PASS: tune_swp_psoct_gam ran; recommended MaxSplits=%d, NumTrees=%d\n', ...
            tuning.recommended.MaxSplits, tuning.recommended.NumTrees);
catch e
    pass = false;
    fprintf('  FAIL: tune_swp_psoct_gam -- %s\n', e.message);
    for i = 1:numel(e.stack)
        fprintf('        at %s line %d\n', e.stack(i).name, e.stack(i).line);
    end
end

fprintf('\n========================================\n');
fprintf('STEP 2: fit_isosurface_by_dataset (severe only)\n');
fprintf('========================================\n');
try
    res_single = fit_isosurface_by_dataset(T_severe, ...
        'MaxSplitsGrid', [5 10], 'NumTreesGrid', [50 100], 'KFold', 3, ...
        'NumGridPts', 40, 'NBootstrap', 20, 'BootstrapCI', true, ...
        'TitleStr', 'SmokeTest_Severe_Single', 'dirout', dirout, 'Verbose', Verbose);

    assert(numel(res_single.epv_grid) == 40, 'epv_grid length wrong');
    assert(numel(res_single.ves_levels) == 3, 'ves_levels count wrong');
    assert(numel(res_single.slice_scatter) == 3, 'slice_scatter cell count wrong');
    assert(numel(res_single.slice_retard) == 3, 'slice_retard cell count wrong');
    for i = 1:3
        assert(numel(res_single.slice_scatter{i}) == 40, 'slice_scatter{%d} length wrong', i);
        assert(numel(res_single.slice_retard{i}) == 40, 'slice_retard{%d} length wrong', i);
    end
    assert(isequal(size(res_single.density_mask), [40, 3]), 'density_mask shape wrong');
    assert(numel(res_single.ci_scatter) == 3, 'ci_scatter cell count wrong');
    assert(isequal(size(res_single.ci_scatter{1}), [40, 2]), 'ci_scatter{1} shape wrong');
    fprintf('  PASS: fit_isosurface_by_dataset ran; all output shapes correct\n');
catch e
    pass = false;
    fprintf('  FAIL: fit_isosurface_by_dataset -- %s\n', e.message);
    for i = 1:numel(e.stack)
        fprintf('        at %s line %d\n', e.stack(i).name, e.stack(i).line);
    end
end

fprintf('\n========================================\n');
fprintf('STEP 3: fit_isosurface_by_dataset (control only)\n');
fprintf('========================================\n');
try
    res_single_control = fit_isosurface_by_dataset(T_control, ...
        'MaxSplitsGrid', [5 10], 'NumTreesGrid', [50 100], 'KFold', 3, ...
        'NumGridPts', 40, 'NBootstrap', 20, 'BootstrapCI', true, ...
        'TitleStr', 'SmokeTest_Control_Single', 'dirout', dirout, 'Verbose', Verbose);
    fprintf('  PASS: fit_isosurface_by_dataset (control) ran\n');
catch e
    pass = false;
    fprintf('  FAIL: fit_isosurface_by_dataset (control) -- %s\n', e.message);
end

fprintf('\n========================================\n');
fprintf('STEP 4: compare_isosurfaces_severe_control (intersection grid, pooled percentiles)\n');
fprintf('========================================\n');
try
    cmp = compare_isosurfaces_severe_control(T_severe, T_control, ...
        'MaxSplitsGrid', [5 10], 'NumTreesGrid', [50 100], 'KFold', 3, ...
        'NumGridPts', 40, 'NBootstrap', 20, 'BootstrapCI', true, ...
        'GridRange', 'intersection', 'SharedPercentiles', true, ...
        'TitleStr', 'SmokeTest_Compare', 'dirout', dirout, 'Verbose', Verbose);

    assert(numel(cmp.epv_grid) == 40, 'compare epv_grid length wrong');
    assert(isequal(cmp.severe.ves_levels, cmp.control.ves_levels), ...
           'SharedPercentiles=true should give identical ves_levels for both groups');
    assert(cmp.grid_range_used(1) >= min(T_severe.epv_swp), 'grid lo should not go below severe min under intersection');
    assert(cmp.grid_range_used(2) <= max(T_severe.epv_swp) + 1e-9, 'grid hi should not exceed severe max under intersection');
    assert(cmp.grid_range_used(1) >= min(T_control.epv_swp) - 1e-9, 'grid lo should not go below control min under intersection');
    assert(cmp.grid_range_used(2) <= max(T_control.epv_swp) + 1e-9, 'grid hi should not exceed control max under intersection');
    fprintf('  PASS: compare_isosurfaces_severe_control (intersection) ran; shared grid=[%.4f %.4f]\n', ...
            cmp.grid_range_used(1), cmp.grid_range_used(2));
    fprintf('        severe pct_masked=%.1f%%, control pct_masked=%.1f%%\n', ...
            cmp.severe.pct_masked, cmp.control.pct_masked);
catch e
    pass = false;
    fprintf('  FAIL: compare_isosurfaces_severe_control (intersection) -- %s\n', e.message);
    for i = 1:numel(e.stack)
        fprintf('        at %s line %d\n', e.stack(i).name, e.stack(i).line);
    end
end

fprintf('\n========================================\n');
fprintf('STEP 5: compare_isosurfaces_severe_control (union grid, independent percentiles)\n');
fprintf('========================================\n');
try
    cmp2 = compare_isosurfaces_severe_control(T_severe, T_control, ...
        'MaxSplitsGrid', [5 10], 'NumTreesGrid', [50 100], 'KFold', 3, ...
        'NumGridPts', 30, 'BootstrapCI', false, ...
        'GridRange', 'union', 'SharedPercentiles', false, ...
        'TitleStr', 'SmokeTest_Compare_Union', 'dirout', dirout, 'Verbose', Verbose);

    assert(~isequal(cmp2.severe.ves_levels, cmp2.control.ves_levels), ...
           'SharedPercentiles=false with different distributions should give different ves_levels');
    assert(cmp2.grid_range_used(1) <= min(T_severe.epv_swp) + 1e-9, 'union grid lo should reach down to at least the lower group min');
    assert(cmp2.grid_range_used(2) >= max(T_severe.epv_swp) - 1e-9 || cmp2.grid_range_used(2) >= max(T_control.epv_swp) - 1e-9, ...
           'union grid hi should reach up to at least one group''s max');
    fprintf('  PASS: compare_isosurfaces_severe_control (union, independent percentiles) ran\n');
catch e
    pass = false;
    fprintf('  FAIL: compare_isosurfaces_severe_control (union) -- %s\n', e.message);
end

fprintf('\n========================================\n');
fprintf('STEP 6: Edge case -- disjoint ranges should raise a clear error, not crash oddly\n');
fprintf('========================================\n');
try
    T_far = T_control;
    T_far.epv_swp = T_far.epv_swp + 100;  % shove control epv_swp far away from severe's range
    try
        compare_isosurfaces_severe_control(T_severe, T_far, ...
            'TuneParams', false, 'MaxSplits', [10 10], 'NumTrees', [100 100], ...
            'BootstrapCI', false, 'GridRange', 'intersection', ...
            'TitleStr', 'SmokeTest_NoOverlap', 'Verbose', false);
        pass = false;
        fprintf('  FAIL: expected an error for non-overlapping ranges but none was raised\n');
    catch e2
        if ~isempty(strfind(e2.message, 'do not overlap'))
            fprintf('  PASS: correctly raised a clear no-overlap error\n');
        else
            pass = false;
            fprintf('  FAIL: raised an error, but not the expected no-overlap message: %s\n', e2.message);
        end
    end
catch e
    pass = false;
    fprintf('  FAIL: edge-case test setup itself errored -- %s\n', e.message);
end

fprintf('\n========================================\n');
if pass
    fprintf('ALL SMOKE TESTS PASSED\n');
else
    fprintf('SOME SMOKE TESTS FAILED -- see FAIL lines above\n');
end
fprintf('========================================\n');
