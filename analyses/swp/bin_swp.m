function [xy_out, se_out] = bin_swp(pair, N)
% Efficient scatter plot: averages all points in disjoint x-windows,
% returning at most N points
% INPUTS:
%   pair (Nx2 matrix): [x; y]
%   N (integer): max number of plotted points (bins)
%   xlab (string): x-axis label
%   ylab (string): y-axis label
%   tit (string): title of figure
%   plt_flag (bool): boolean for plotting within function
% OUTPUT:
%   xy_out (Mx2): matrix of [x_bin_mean, y_bin_mean] pairs
    
% Cast to single to reduce memory
x = single(pair(:,1));
y = single(pair(:,2));
clear pair

% Removing outliers using IQR method
Q1 = prctile(y, 25); % First quartile
Q3 = prctile(y, 75); % Third quartile
IQR = Q3 - Q1;       % Interquartile range
lower_bound = Q1 - 1.5 * IQR; % Lower bound
upper_bound = Q3 + 1.5 * IQR; % Upper bound

% Remove outliers
valid_pts = (y >= lower_bound) & (y <= upper_bound);
x = x(valid_pts);
y = y(valid_pts);

% Retrieve maximum of x-axis data
xmin = min(x);
xmax = max(x);

% Protect against degenerate data
if xmax == xmin
    error('All x values are identical.');
end

% Calculate edges for binning
edges = linspace(xmin, xmax, N+1);

% Bin data
[~, ~, bin] = histcounts(x, edges);

% Change "bin" from double to smaller data type (if applicable)
bin_max = max(bin(:));
if bin_max <= intmax("uint8")
    bin = uint8(bin);
elseif bin_max <= intmax("uint16")
    bin = uint16(bin);
elseif bin_max <= intmax("uint32")
    bin = uint32(bin);
elseif bin_max <= intmax("uint64")
    bin = uint64(bin);
elseif bin_max <= realmax("single")
    bin = single(bin);
end

% Remove data outside bin range (bin==0)
valid_pts = bin > 0;
bin = bin(valid_pts);
x = x(valid_pts);
y = y(valid_pts);

% Compute mean y for each bin using accumarray
y_bin_mean = accumarray(bin(:), y(:), [N 1], @mean, single(0));

%%% Calculate the standard error for y values in each bin
% standard deviation
y_bin_std = accumarray(bin(:), y(:), [N 1], @std, single(0));
% number of samples
y_bin_count = accumarray(bin(:), 1, [N 1], @sum);
% standard error
y_bin_se = y_bin_std ./ sqrt(y_bin_count);

% Remove bins with no points (NaN means no data in that bin)
valid_bins = ~isnan(y_bin_mean);

% Calculate midpoints of the bins for the x-axis points
bin_midpoints = (edges(1:end-1) + edges(2:end)) / 2;

% Construct output with midpoints for x values
xy_out = [bin_midpoints(valid_bins)', y_bin_mean(valid_bins)];

% Standard Error output
se_out = y_bin_se(valid_bins);

end