function [x_cent, y_mean, y_stds, y_sems] = window_avg2(data, window_size)
%WINDOW_AVG Computes window-averaged statistics over x-axis in disjoint windows.
%
% INPUTS:
%   data        - Nx2 matrix: column 1 = x, column 2 = y
%   window_size - scalar, width of window on the x-axis
%
% OUTPUTS:
%   x_cent   - column vector of center x values for each window
%   y_mean   - column vector of mean y values within each window
%   y_stds   - column vector of y standard deviation within each window
%   y_sems   - column vector of standard error within each window
%
% Example:
%   [x_cent, y_mean, y_stds, y_sems] = window_avg(data, 0.5);

if nargin < 2
    error('Requires data and window_size as inputs.');
end

x = data(:,1);
y = data(:,2);

x_min = min(x);
x_max = max(x);

% Calculate window edges
edges = x_min : window_size : x_max;
if edges(end) < x_max
    edges = [edges, x_max]; % Ensure last edge includes all data
end
n_windows = length(edges) - 1;

% Pre-allocate
x_cent  = nan(n_windows, 1);
y_mean  = nan(n_windows, 1);
y_stds  = nan(n_windows, 1);
y_sems  = nan(n_windows, 1);

for i = 1:n_windows
    start_x = edges(i);
    end_x = edges(i+1);
    center = (start_x + end_x) / 2;

    if i < n_windows
        in_window = (x >= start_x) & (x < end_x);
    else
        in_window = (x >= start_x) & (x <= end_x); % Include right edge on last
    end

    y_window = y(in_window);

    if ~isempty(y_window)
        x_cent(i)  = center;
        y_mean(i)  = mean(y_window);
        y_stds(i)  = std(y_window);
        y_sems(i)  = y_stds(i) / sqrt(length(y_window));
    end
end

% Remove any nan windows (where no points were found)
valid = ~isnan(y_mean);
x_cent = x_cent(valid);
y_mean = y_mean(valid);
y_stds = y_stds(valid);
y_sems = y_sems(valid);

end