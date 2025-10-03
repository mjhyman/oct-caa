function scatter_subset(pair, xlab, ylab, tit, th, d)
% Efficient scatter plot: averages all points in disjoint x-windows of width d above threshold th
%
% INPUTS:
%   pair (Nx2 matrix): [EPVS density, optical property]
%   xlab (string): x-axis label
%   ylab (string): y-axis label
%   tit (string): title of figure
%   th (double): minimum threshold for EPVS density
%   d  (double): window width

    if nargin < 6
        error('All six inputs are required.');
    end

    % Remove pairs below threshold
    x = pair(:,1);
    y = pair(:,2);
    keep = x >= th;
    x = x(keep);
    y = y(keep);

    if isempty(x)
        error('No data points remain after thresholding.');
    end

    % Define bin edges
    x_min = min(x);
    x_max = max(x);
    edges = x_min:d:x_max;
    if edges(end) < x_max
        edges = [edges x_max];
    end

    % Bin data using histcounts
    [~, ~, bin] = histcounts(x, edges);

    % Remove data outside bin range (bin==0)
    valid_pts = bin > 0;
    bin = bin(valid_pts);
    x = x(valid_pts);
    y = y(valid_pts);

    % Compute mean x and y for each bin using accumarray
    x_bin_mean = accumarray(bin, x, [], @mean);
    y_bin_mean = accumarray(bin, y, [], @mean);

    % Scatter plot
    figure;
    scatter(x_bin_mean, y_bin_mean, 60, 'b', 'filled');
    xlabel(xlab);
    ylabel(ylab);
    title(tit);
    set(gca, 'fontsize', 20);
end