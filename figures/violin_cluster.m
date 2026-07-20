function h = violin_cluster(dataCell, groupLabels, subGroupLabels, ...
                            hexColors, useLog, yLab, fontSize)
    % violin_cluster  Grouped violin plots with embedded box-and-whisker
    %                 and within-cluster Cliff's delta brackets.
    %
    % Expects exactly 2 sub-groups per main group (e.g. "EPVS" and "Ves").
    % Cliff's delta is computed for each main-group pair (EPVS vs Ves)
    % and displayed as a bracket above each cluster.
    %
    % Inputs:
    %   dataCell       - Cell array of data vectors (one per condition)
    %   groupLabels    - Cell array of main-group labels (one per cell)
    %                    e.g. {'Control','Control','Mild','Mild','Severe','Severe'}
    %   subGroupLabels - Cell array of sub-group labels (one per cell)
    %                    e.g. {'EPVS','Ves','EPVS','Ves','EPVS','Ves'}
    %   hexColors      - Cell array of hex color strings (one per sub-group)
    %   useLog         - Boolean for logarithmic y-axis
    %   yLab           - String for the Y-axis label
    %   fontSize       - Numeric font size (e.g., 12)

    % ------------------------------------------------------------------ %
    % 1. Unique groups / sub-groups
    % ------------------------------------------------------------------ %
    uniqueGroups = unique(string(groupLabels), 'stable');
    uniqueSubs   = unique(string(subGroupLabels), 'stable');

    numTotalCells = numel(dataCell);
    numSub        = numel(uniqueSubs);
    numMain       = numel(uniqueGroups);

    if numSub ~= 2
        error('violin_cluster: Cliff''s delta bracket requires exactly 2 sub-groups per cluster (found %d).', numSub);
    end

    % ------------------------------------------------------------------ %
    % 2. KDE + box statistics per cell
    % ------------------------------------------------------------------ %
    numKDEPoints = 512;
    kdeCurves    = cell(numTotalCells, 1);

    for i = 1:numTotalCells
        x = dataCell{i}(:);
        x = x(~isnan(x));

        if useLog
            x = log10(x(x > 0));
        end

        % Support must be strictly outside all data points.
        % Use actual min/max and pad by 0.1% of the range (or a small
        % absolute value if the range is zero).
        xMin  = min(x);
        xMax  = max(x);
        xPad  = max((xMax - xMin) * 0.001, 1e-6);

        [f, xi] = ksdensity(x, 'NumPoints', numKDEPoints, ...
                               'Support', [xMin - xPad, xMax + xPad]);

        f  = (f  / max(f));
        f = f(:)';
        xi = xi(:)';

        if useLog
            xi = 10.^xi;
        end

        % Box statistics on original (non-log) data
        rawX   = dataCell{i}(:);
        rawX   = rawX(~isnan(rawX));
        q1     = quantile(rawX, 0.25);
        med    = median(rawX);
        q3     = quantile(rawX, 0.75);
        iqrVal = q3 - q1;
        rawLo  = quantile(rawX, 0.001);
        rawHi  = quantile(rawX, 0.999);
        whisLo = max(rawLo, q1 - 1.5 * iqrVal);
        whisHi = min(rawHi, q3 + 1.5 * iqrVal);

        kdeCurves{i} = struct('xi', xi, 'f', f, ...
                              'q1', q1, 'med', med, 'q3', q3, ...
                              'whisLo', whisLo, 'whisHi', whisHi);
    end

    % ------------------------------------------------------------------ %
    % 3. Compute x-positions (needed for brackets before plotting)
    % ------------------------------------------------------------------ %
    groupGap = 1.0;
    subGap   = 0.3;

    xPositions = zeros(1, numTotalCells);
    for i = 1:numTotalCells
        mainIdx        = find(strcmp(uniqueGroups, string(groupLabels{i})));
        subIdx         = find(strcmp(uniqueSubs,   string(subGroupLabels{i})));
        xPositions(i)  = (mainIdx * groupGap) + (subIdx * subGap);
    end

    % ------------------------------------------------------------------ %
    % 4. Cliff's delta within each main group (sub1 vs sub2)
    % ------------------------------------------------------------------ %
    maxN     = 100000;
    cliffD   = zeros(1, numMain);
    cliffStr = cell(1, numMain);

    for m = 1:numMain
        % Find the two cells belonging to this main group
        idx1 = find(strcmp(string(groupLabels), uniqueGroups(m)) & ...
                    strcmp(string(subGroupLabels), uniqueSubs(1)));
        idx2 = find(strcmp(string(groupLabels), uniqueGroups(m)) & ...
                    strcmp(string(subGroupLabels), uniqueSubs(2)));

        if isempty(idx1) || isempty(idx2)
            warning('violin_cluster: could not find both sub-groups for group "%s". Skipping Cliff''s delta.', uniqueGroups(m));
            cliffStr{m} = '';
            continue;
        end

        x1 = dataCell{idx1}(:);  x1 = x1(~isnan(x1));
        x2 = dataCell{idx2}(:);  x2 = x2(~isnan(x2));

        % Subsample for efficiency
        x1 = x1(randperm(numel(x1), min(numel(x1), maxN)));
        x2 = x2(randperm(numel(x2), min(numel(x2), maxN)));

        n1 = numel(x1);
        n2 = numel(x2);

        rnks       = tiedrank([x1; x2]);
        U          = sum(rnks(1:n1)) - n1*(n1+1)/2;
        cliffD(m)  = (2*U / (n1*n2)) - 1;
        cliffStr{m} = sprintf('d=%.2f (%s)', cliffD(m), ...
                               interpretCliffsDelta(abs(cliffD(m))));
    end

    % ------------------------------------------------------------------ %
    % 5. Plotting
    % ------------------------------------------------------------------ %
    figure('Color', 'w', 'Units', 'inches', 'Position', [0 0 6 5]);
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [6 5], ...
             'PaperPosition', [0 0 6 5]);
    hold on;

    violinWidth = 0.35;
    boxW        = 0.08;

    subHandles = gobjects(1, numSub);

    % Track y-extent of whiskers for bracket placement
    yWhisMax = -Inf;

    for i = 1:numTotalCells
        subIdx  = find(strcmp(uniqueSubs,   string(subGroupLabels{i})));

        xPos = xPositions(i);
        c    = validatecolor(hexColors{subIdx});
        kc   = kdeCurves{i};

        yWhisMax = max(yWhisMax, kc.whisHi);

        % --- Violin patch ----------------------------------------------
        fScaled = kc.f * violinWidth;
        xPatch  = [xPos + fScaled,  fliplr(xPos - fScaled)];
        yPatch  = [kc.xi,            fliplr(kc.xi)];
        p = patch(xPatch, yPatch, c, ...
                  'FaceAlpha', 0.4, ...
                  'EdgeColor', c, ...
                  'LineWidth', 1);

        if ~isgraphics(subHandles(subIdx))
            subHandles(subIdx) = p;
        end

        % --- Whiskers --------------------------------------------------
        plot([xPos xPos], [kc.whisLo, kc.q1], '-', ...
             'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        plot([xPos xPos], [kc.q3, kc.whisHi], '-', ...
             'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);

        % --- IQR box ---------------------------------------------------
        boxH = max(kc.q3 - kc.q1, eps);
        rectangle('Position', [xPos - boxW, kc.q1, 2*boxW, boxH], ...
                  'FaceColor', [1 1 1], ...
                  'EdgeColor', [0.2 0.2 0.2], ...
                  'LineWidth', 1.5);

        % --- Median line -----------------------------------------------
        plot([xPos - boxW, xPos + boxW], [kc.med, kc.med], '-', ...
             'Color', [0.2 0.2 0.2], 'LineWidth', 2.5);
    end

    % ------------------------------------------------------------------ %
    % 6. Cliff's delta brackets (one per main group)
    % ------------------------------------------------------------------ %
    % Estimate a sensible y-range for padding
    allData  = vertcat(dataCell{:});
    allData  = allData(~isnan(allData));
    yRange   = quantile(allData, 0.999) - quantile(allData, 0.001);

    bracketPad = yRange * 0.05;
    bracketH   = yRange * 0.02;
    bracketY   = yWhisMax + bracketPad;

    for m = 1:numMain
        if isempty(cliffStr{m}), continue; end

        % x-positions of the two violins in this cluster
        idx1 = strcmp(string(groupLabels), uniqueGroups(m)) & ...
                    strcmp(string(subGroupLabels), uniqueSubs(1));
        idx2 = strcmp(string(groupLabels), uniqueGroups(m)) & ...
                    strcmp(string(subGroupLabels), uniqueSubs(2));

        x1 = xPositions(idx1);
        x2 = xPositions(idx2);

        % Horizontal bar
        plot([x1, x2], [bracketY, bracketY], '-k', 'LineWidth', 1.2);
        % Left tick
        plot([x1, x1], [bracketY - bracketH, bracketY], '-k', 'LineWidth', 1.2);
        % Right tick
        plot([x2, x2], [bracketY - bracketH, bracketY], '-k', 'LineWidth', 1.2);
        % Label
        text((x1 + x2) / 2, bracketY + bracketH * 0.5, cliffStr{m}, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment',   'bottom', ...
             'FontSize', 9, ...
             'FontName', 'Arial', ...
             'Color', [0.3 0.3 0.3]);
    end

    % ------------------------------------------------------------------ %
    % 7. Axis formatting
    % ------------------------------------------------------------------ %
    h = gca;

    tickPos = (1:numMain) * groupGap + (((1 + numSub) / 2) * subGap);

    set(h, 'XTick',      tickPos, ...
           'XTickLabel', uniqueGroups, ...
           'TickDir',    'out', ...
           'FontName',   'Arial', ...
           'FontSize',   fontSize);

    box off;

    if useLog
        set(h, 'YScale', 'log');
    end

    % Expand y-axis to accommodate brackets and labels
    yTop = bracketY + yRange * 0.1;
    ylim([quantile(allData, 0.001) - bracketPad, yTop]);

    ylabel(yLab, 'FontName', 'Arial', 'FontSize', fontSize);
    % legend(subHandles, uniqueSubs, 'Location', 'northwest', ...
    %        'FontName', 'Arial', 'FontSize', fontSize);

    grid on;
    hold off;

    % Print summary to command window
    fprintf('\nCliff''s delta (EPVS vs Vessel within each cluster):\n');
    for m = 1:numMain
        if ~isempty(cliffStr{m})
            fprintf('  %s: %s\n', uniqueGroups(m), cliffStr{m});
        end
    end
    fprintf('Magnitude: |d|<0.147 negligible, 0.147-0.33 small, 0.33-0.474 medium, >0.474 large\n\n');
end


% ------------------------------------------------------------------ %
% Local helper: Cliff's delta interpretation
% ------------------------------------------------------------------ %
function s = interpretCliffsDelta(d)
    if d < 0.147
        s = 'negligible';
    elseif d < 0.33
        s = 'small';
    elseif d < 0.474
        s = 'medium';
    else
        s = 'large';
    end
end