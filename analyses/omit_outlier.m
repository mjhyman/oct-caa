%% Omit outliers from repeated measures within each subject/region
function clean_data = omit_outlier(raw_data)
% Remove outliers, compute median

% Remove NaNs first
raw_data = raw_data(~isnan(raw_data));

% Outlier removal using IQR
try
    Q1 = prctile(real(raw_data), 25);
    Q3 = prctile(real(raw_data), 75);
    IQR = Q3 - Q1;
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
catch
    error('wrongo')
end

% Remove outliers from data
clean_data = raw_data(real(raw_data) >= lower_bound & ...
                      real(raw_data) <= upper_bound);
% Then compute mean on clean data
clean_data = median(real(clean_data));
end