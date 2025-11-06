%% Omit outliers from repeated measures within each subject/region
function clean_data = omit_outlier(raw_data, th)
% Remove outliers, compute median
% INPUTS:
%   raw_data (vector): vector of values to be parsed
%   th (float): upper limit to apply threshold

% Remove NaNs first
raw_data = raw_data(~isnan(raw_data));

% Apply upper-limit threshold
raw_data = raw_data(raw_data <= th);

% Outlier removal using IQR
try
    Q1 = prctile(real(raw_data), 25);
    Q3 = prctile(real(raw_data), 75);
    IQR = Q3 - Q1;
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
catch
    error('\nFailed to apply the 1.5*IQR thresholding\n')
end

% Remove outliers from data
clean_data = raw_data(real(raw_data) >= lower_bound & ...
                      real(raw_data) <= upper_bound);
end