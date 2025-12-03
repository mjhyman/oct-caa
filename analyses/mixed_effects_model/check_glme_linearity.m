function check_glme_linearity(lme, tbl, nlThresh, tstr)
% checkLinearityGLME_withNLScore(lme, tbl, nlThresh)
% Generates diagnostic plots and calculates nonlinearity scores.
%
% Inputs:
%   lme      - LinearMixedModel object (from fitlme)
%   tbl      - Table used to fit the model
%   nlThresh - Threshold for nonlinearity score (e.g., 0.05). Default: 0.05

    if nargin < 3
        nlThresh = 0.05;
    end

    % Extract model info
    fittedVals = fitted(lme);
    residVals = residuals(lme);
    fe = fixedEffects(lme);
    feNames = lme.Formula.FELinearFormula.PredictorNames;

    % Identify numeric predictors
    predictorVars = tbl.Properties.VariableNames;
    numericPredictors = predictorVars(varfun(@isnumeric, tbl,...
                                      'OutputFormat', 'uniform'));

    %% 1. Residuals vs. Fitted
    figure;
    scatter(fittedVals, residVals, 'filled');
    xlabel('Fitted Values');
    ylabel('Residuals');
    title({'Residuals vs Fitted Values',tstr});
    yline(0, '--k');
    grid on;

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
            title({sprintf('Residuals vs %s (Nonlinearity Score: %.4f*)',...
                predictor, nlScore),tstr});
        else
            title({sprintf('Residuals vs %s (Nonlinearity Score: %.4f)',...
                predictor, nlScore),tstr});
        end

        yline(0, '--k');
        grid on;
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
            title({['Component + Residual Plot for ', predictor],tstr});
            grid on;
        end
    end

    %% Plot Residuals    
    figure;
    plotResiduals(lme, 'fitted');
    grid on;
    title('Residuals')
end
