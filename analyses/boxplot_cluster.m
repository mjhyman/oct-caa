function h = boxplot_cluster(dataCell, groupLabels, subGroupLabels,...
                            hexColors, useLog, yLab, ylims, fontSize)
    % Inputs:
    %   useLog   - Boolean (true/false) for logarithmic y-axis
    %   yLab     - String for the Y-axis label
    %   fontSize - Numeric value for font size (e.g., 12)

    % 1. Determine unique groups and subgroups
    uniqueGroups = unique(string(groupLabels), 'stable');
    uniqueSubs = unique(string(subGroupLabels), 'stable');
    
    numTotalCells = numel(dataCell);
    numSub = numel(uniqueSubs);
    numMain = numel(uniqueGroups);
    
    % 2. Calculate Stats (5th, 25th, 50th, 75th, 95th percentiles)
    stats = zeros(5, numTotalCells);
    for i = 1:numTotalCells
        % prctile handles large vectors efficiently
        stats(:, i) = prctile(dataCell{i}, [5, 25, 50, 75, 95]);
    end

    % 3. Plotting logic
    figure('Color', 'w'); hold on;
    width = 0.2;      
    groupGap = 1.0;   
    subGap = 0.3;     
    
    subHandles = gobjects(1, numSub); 
    
    for i = 1:numTotalCells
        mainIdx = find(strcmp(uniqueGroups, groupLabels{i}));
        subIdx = find(strcmp(uniqueSubs, subGroupLabels{i}));
        
        xPos = (mainIdx * groupGap) + (subIdx * subGap);
        s = stats(:, i);
        c = validatecolor(hexColors{subIdx});
        
        % Draw Whisker
        line([xPos xPos], [s(1) s(5)], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
        
        % Draw Box (Patch)
        px = [xPos-width/2, xPos+width/2, xPos+width/2, xPos-width/2];
        py = [s(2), s(2), s(4), s(4)];
        p = patch(px, py, c, 'FaceAlpha', 1, 'EdgeColor', 'k');
        
        if isempty(subHandles(subIdx)) || ~isgraphics(subHandles(subIdx))
            subHandles(subIdx) = p;
        end
        
        % Draw Median
        line([xPos-width/2 xPos+width/2], [s(3) s(3)], 'Color', 'k', 'LineWidth', 2);
    end

    % 4. Formatting Axes and Font
    h = gca;
    tickPos = (1:numMain) * groupGap + (((1 + numSub) / 2) * subGap);
    
    set(h, 'XTick', tickPos, ...
           'XTickLabel', uniqueGroups, ...
           'TickDir', 'out', ...
           'FontName', 'Arial', ...
           'FontSize', fontSize);
    
    % Apply Logarithmic Scale if requested
    if useLog
        set(h, 'YScale', 'log');
    end
    
    % Labels and Legend
    ylabel(yLab, 'FontName', 'Arial', 'FontSize', fontSize);
    % legend(subHandles, uniqueSubs, 'Location', 'northwest',...
    %        'FontName', 'Arial', 'FontSize', fontSize);
    ylim(ylims);
    
    grid on;
end