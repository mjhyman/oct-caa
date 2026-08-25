function [xy_out, se_out, info] = bin_swp(pair, N)
% bin_swp - Average y within disjoint x-windows, returning at most N points.
%
% INPUTS:
%   pair (n x 2) : column 1 = x, column 2 = y
%   N (integer)  : number of bins
%
% OUTPUTS:
%   xy_out (M x 2) : [bin_midpoint, y_bin_mean] for bins that hold data
%   se_out (M x 1) : standard error of y in each returned bin, NaN when the
%                    bin holds a single observation
%   info (struct)  : nDropped, nZeros, nOutliers, nEmptyBins, edges
%
% Fixes relative to the previous version:
%   1. y was filtered for zeros BEFORE the outlier mask was built, so the
%      mask (length numel(y_after_zeros)) was then applied to the unfiltered
%      x. MATLAB permits a logical index shorter than the array, so this ran
%      silently and mispaired every observation after the first removed zero.
%   2. y_bin_mean was shrunk by (y_bin_mean ~= 0) before valid_bins was
%      built, leaving it shorter than y_bin_count -> the size mismatch.
%   3. An empty bin and a bin whose mean is genuinely 0 were indistinguishable.
%      Emptiness is now decided by the count, which is what it means.
%   4. Single-observation bins reported SE = 0 (implying perfect precision)
%      instead of NaN.

narginchk(2, 2);
if size(pair, 2) ~= 2
    error('bin_swp:shape', 'pair must be [n x 2]; got [%d x %d].', size(pair,1), size(pair,2));
end

%% --- Cast to single to reduce memory, keep x and y locked together ---
x = single(pair(:,1));
y = single(pair(:,2));
clear pair

nStart = numel(x);

% Non-finite: isnan alone leaves Inf, and one Inf makes every bin edge Inf
keep = isfinite(x) & isfinite(y);
x = x(keep);
y = y(keep);
nDropped = nStart - numel(x);

%% --- Remove zeros and IQR outliers in ONE joint mask ---
% Bounds are computed on the non-zero y, matching the original intent.
nz = (y ~= 0);
if ~any(nz)
    error('bin_swp:allZero', 'All y values are zero.');
end

Q1  = prctile(y(nz), 25);
Q3  = prctile(y(nz), 75);
IQRv = Q3 - Q1;
lower_bound = Q1 - 1.5 * IQRv;
upper_bound = Q3 + 1.5 * IQRv;

inRange   = (y >= lower_bound) & (y <= upper_bound);
valid_pts = nz & inRange;          % single mask, applied to both vectors

nZeros    = nnz(~nz);
nOutliers = nnz(nz & ~inRange);

x = x(valid_pts);
y = y(valid_pts);

if isempty(x)
    error('bin_swp:empty', 'No points remain after zero and outlier removal.');
end

%% --- Bin edges ---
xmin = min(x);
xmax = max(x);
if xmax == xmin
    error('bin_swp:degenerate', 'All x values are identical (%g).', xmin);
end

edges = linspace(double(xmin), double(xmax), N + 1);
[~, ~, bin] = histcounts(double(x), edges);

% Drop anything outside the edges (bin == 0)
inBin = bin > 0;
bin   = bin(inBin);
y     = y(inBin);

%% --- Accumulate at an explicit [N 1] size ---
% Vectorised sums rather than accumarray(..., @mean) / (..., @std): the
% function-handle form groups values one bin at a time and is far slower on
% millions of points. y is centred first so the one-pass variance stays
% numerically stable in single precision.
bin = double(bin(:));               % accumarray wants double subscripts
yd  = double(y(:));

y_bin_count = accumarray(bin, 1, [N 1], @sum, 0);

yBar = mean(yd);
yc   = yd - yBar;
s1   = accumarray(bin, yc,    [N 1], @sum, 0);
s2   = accumarray(bin, yc.^2, [N 1], @sum, 0);

y_bin_mean = nan(N, 1);
y_bin_se   = nan(N, 1);

has1 = y_bin_count > 0;
has2 = y_bin_count > 1;

y_bin_mean(has1) = yBar + s1(has1) ./ y_bin_count(has1);

varY = (s2(has2) - (s1(has2).^2) ./ y_bin_count(has2)) ./ (y_bin_count(has2) - 1);
y_bin_se(has2) = sqrt(max(varY, 0)) ./ sqrt(y_bin_count(has2));

%% --- Apply one mask to midpoints, means, and SEs ---
% y_bin_mean is still [N 1] here, so this can no longer mismatch.
valid_bins    = has1 & ~isnan(y_bin_mean);
bin_midpoints = (edges(1:end-1) + edges(2:end)) / 2;

xy_out = [bin_midpoints(valid_bins)', y_bin_mean(valid_bins)];
se_out = y_bin_se(valid_bins);

info = struct('nDropped',   nDropped, ...
              'nZeros',     nZeros, ...
              'nOutliers',  nOutliers, ...
              'nEmptyBins', N - nnz(valid_bins), ...
              'edges',      edges);
end