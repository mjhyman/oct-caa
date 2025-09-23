function [x_cent, y_mean, y_stds, y_sems] = window_avg(data, window_size)
% Window average over x-axis, averaging y-values within disjoint windows.
%
% Inputs:
%   data        - Nx2 matrix: column 1 = x, column 2 = y
%   window_size - width of the window on the x-axis
%
% Outputs:
%   x_cent   - center x value of each window
%   y_mean   - average y value within each window
%   y_stds   - standard deviation of y within each window
%   y_sems   - standard error of y within each window

x = data(:,1);
y = data(:,2);

x_min = min(x);
x_max = max(x);

x_cent = [];
y_mean = [];
y_stds = [];
y_sems = [];

for start_x = x_min : window_size : x_max
    end_x = start_x + window_size;
    center = (start_x + end_x) / 2;

    in_window = x >= start_x & x < end_x;
    y_window = y(in_window);

    if ~isempty(y_window)
        x_cent(end+1) = center;
        y_mean(end+1) = mean(y_window);
        y_stds(end+1)  = std(y_window);
        y_sems(end+1)  = std(y_window) / sqrt(length(y_window));
    end
end
end
