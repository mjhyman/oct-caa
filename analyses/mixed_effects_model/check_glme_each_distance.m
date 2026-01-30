function check_glme_each_distance(lme, tbl, numericPredictors,...
                              nlThresh, lintest_dir, op, dist, reg)
% Generates diagnostic plots and calculates nonlinearity scores.
%
% Inputs:
%   lme      - LinearMixedModel object (from fitlme)
%   tbl      - Table used to fit the model
%   numericPredictors (cell array) - numeric predictors
%   nlThresh - Threshold for nonlinearity score (e.g., 0.05). Default: 0.05
%   lintest_dir (string) - output directory for storing plots
%   op (string) - optical property name (scattering or retardance)
%   dist (uint) - distance in microns
%   reg (string) - region (front or occip)

%% Create output directory for storing plots
% Concatenate output directory
ddir = fullfile(lintest_dir, string(op), string(reg), string(dist));
% Check if directory exists. If not, then create
if ~exist(ddir, 'dir')
    mkdir(ddir);
end
fprintf('\nTesting GLME for %s, %s\n', op, reg);

%% Extract model info
fittedVals = fitted(lme);
residVals = residuals(lme);
fe = fixedEffects(lme);
feNames = lme.Formula.FELinearFormula.PredictorNames;

%% Check normality of residuals with lilliefors test
% Perform Lilliefors test for normality on residuals
[h, pValue] = lillietest(residVals);
% Print to console, including optical property, region, distance
if h == 1
    fprintf('\tResiduals are not normally distributed (p = %.4f)\n', pValue);
else
    fprintf('\tResiduals are normally distributed (p = %.4f)\n', pValue);
end
% Save the result of Lilliefors test to ddir
save(fullfile(ddir, 'lilliefors_result.mat'), 'h', 'pValue');

%% Create QQ plot of residuals
figure;
qqplot(residVals);
title('QQ Plot of Residuals');
grid on;
% Save the QQ plot
saveas(gcf, fullfile(ddir, 'qq_plot_residuals.png'));

%% 1. Residuals vs. Fitted
figure;
scatter(fittedVals, residVals, 'filled');
xlabel('Fitted Values');
ylabel('Residuals');
title('Residuals vs Fitted Values');
yline(0, '--k');
grid on;
% Save the Residuals vs Fitted plot
saveas(gcf, fullfile(ddir, 'residuals_vs_fitted.png'));

%% 2. Residuals vs Each Numeric Predictor + Nonlinearity Score
fprintf('\nNonlinearity Scores:\n');
for i = 1:numel(numericPredictors)
    predictor = numericPredictors{i};
    x = tbl.(predictor);

    % Sort for smooth plotting
    [xSorted, sortIdx] = sort(x);
    residSorted = residVals(sortIdx);
    smoothed = smooth(xSorted, residSorted, 0.2, 'loess');

    % Compute nonlinearity score
    nlScore = mean(smoothed.^2);
    flag = nlScore > nlThresh;

    % Print score
    if flag
        fprintf('  %s: %.4f  **\n', predictor, nlScore);
    else
        fprintf('  %s: %.4f\n', predictor, nlScore);
    end

    % Plot
    figure;
    scatter(x, residVals, 'filled');
    hold on;
    plot(xSorted, smoothed, 'r-', 'LineWidth', 2);
    xlabel(predictor);
    ylabel('Residuals');

    if flag
        title(sprintf('Residuals vs %s (Nonlinearity Score: %.4f*)',...
                       predictor, nlScore));
    else
        title(sprintf('Residuals vs %s (Nonlinearity Score: %.4f)',...
                       predictor, nlScore));
    end

    yline(0, '--k');
    grid on;

    % Save figure to ddir
    saveas(gcf, fullfile(ddir, sprintf('residuals_vs_%s.png', predictor)));
end

%% 3. Component + Residual (Partial Residual) Plots
for i = 1:numel(feNames)
    predictor = feNames{i};

    if ismember(predictor, numericPredictors)
        x = tbl.(predictor);
        coef = fe(i); % fixed effect coefficient
        partialRes = coef * x + residVals;

        figure;
        scatter(x, partialRes, 'filled');
        xlabel(predictor);
        ylabel('Component + Residual');
        title(['Component + Residual Plot for ', predictor]);
        grid on;
        % Save figure to ddir
        saveas(gcf, fullfile(ddir,sprintf('partial_residuals_%s.png',...
                            predictor)));
    end
end

%% Plot Residuals    
figure;
plotResiduals(lme, 'fitted');
grid on;
title('Residuals')
% Save figure to ddir
saveas(gcf, fullfile(ddir, 'residuals_plot.png'));
end
