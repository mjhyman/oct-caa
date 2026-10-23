%% Add percentage change from the WM median to an existing shells CSV
%{
Reads a CSV written by measure_parenchyma_shells.m, adds:
  mus_med, ret_med : subject/region white-matter median used
  mus_pct, ret_pct : (value - median) / median * 100
  mus_ok,  ret_ok  : raw value below the bound
and writes a new CSV with "_pct" appended (the input is not overwritten).

Requires add_percent_change.m in the same folder or on the path.
%}

clear; clc;

%% Settings
data_dir    = '/projectnb/npbssmic/ns/CAA/';
csv_in      = fullfile(data_dir, 'parenchyma_shells_2026-10-07.csv');   % update
median_file = fullfile(data_dir, 'median_white_matter_values_14-Nov-2025.mat');
mus_upper   = 25;   % raw mus bound
ret_upper   = 45;   % raw retardance bound

assert(isfile(csv_in),      'CSV not found: %s', csv_in);
assert(isfile(median_file), 'Median file not found: %s', median_file);

% Helper lives next to this script
addpath(fileparts(mfilename('fullpath')));

%% Read, convert, write
opts = detectImportOptions(csv_in);
opts = setvartype(opts, {'subject','severity','region','object_id','class'}, 'string');
T = readtable(csv_in, opts);
fprintf('Read %d rows from %s\n', height(T), csv_in);

M = load(median_file, 'parench_median');
T = add_percent_change(T, M.parench_median, mus_upper, ret_upper);

[p, f] = fileparts(csv_in);
csv_out = fullfile(p, [f '_pct.csv']);
writetable(T, csv_out);
fprintf('Wrote %s\n', csv_out);
fprintf('Rows failing bounds: mus %d, ret %d\n', nnz(~T.mus_ok), nnz(~T.ret_ok));