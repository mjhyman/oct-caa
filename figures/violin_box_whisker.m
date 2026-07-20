function stats_table = violin_box_whisker(data, glabels, tstr, ylims, ystr, colors)
% violin_box_whisker - Violin plots with embedded boxplots and Cliff's delta
%
% Inputs:
%   data        - cell array of data vectors, one per group
%   glabels     - cell array of strings for x-axis labels
%   tstr        - string for plot title
%   ylims       - y-axis limits [min max]
%   ystr        - string for y-axis label
%   colors      - nGroups x 3 RGB matrix of colors
%
% Output:
%   stats_table - MATLAB table with pairwise Cliff's delta; export with:
%                 writetable(stats_table, 'results.xlsx', 'Sheet', 'name')
%
% Example:
%   data = {randn(1e6,1), randn(1e6,1)+1, randn(1e6,1)+0.5};
%   T = violin_box_whisker(data, {'Control','Mild','Severe'}, 'Test', [0 5], 'Value', lines(3));

if ~exist('tiedrank', 'file')
    error('violin_box_whisker requires the Statistics and Machine Learning Toolbox (tiedrank missing).');
end

%% --- Remove NaN and define nGroups first ---
data    = cellfun(@(x) x(~isnan(x(:))), data, 'UniformOutput', false);
nGroups = numel(data);

%% --- Input validation ---
if nGroups < 2
    error('violin_box_whisker requires at least 2 groups.');
end
if size(colors, 1) < nGroups
    error('colors must have at least %d rows, one per group.', nGroups);
end
for i = 1:nGroups
    if isempty(data{i})
        error('Group %d ("%s") is empty after NaN removal.', i, glabels{i});
    end
end
if iscell(tstr)
    tstr = tstr{1};
end

violinWidth = 0.35;
nKDE        = 512;
maxN        = 100000;  % cap for computational efficiency only

figure;
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 5]);
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [6, 5], 'PaperPosition', [0, 0, 6, 5]);
hold on;

%% -----------------------------------------------------------------------
%  1. DRAW VIOLINS
%% -----------------------------------------------------------------------
yMax = -Inf;
yMin =  Inf;

for i = 1:nGroups

    x   = data{i}(:);
    xLo = quantile(x, 0.001);
    xHi = quantile(x, 0.999);
    yMax = max(yMax, xHi);
    yMin = min(yMin, xLo);

    %% --- KDE ---
    [f, xi] = ksdensity(x, 'NumPoints', nKDE);
    f  = f(:)';
    xi = xi(:)';
    f  = f / max(f) * violinWidth;

    %% --- Box statistics ---
    q1      = quantile(x, 0.25);
    med     = median(x, 'omitnan');
    q3      = quantile(x, 0.75);
    iqr_val = q3 - q1;
    whislo  = max(xLo, q1 - 1.5 * iqr_val);
    whishi  = min(xHi, q3 + 1.5 * iqr_val);

    %% --- Violin patch ---
    xPatch = [i + f,  fliplr(i - f)];
    yPatch = [xi,     fliplr(xi)];
    patch(xPatch, yPatch, colors(i,:), ...
        'FaceAlpha', 0.4, 'EdgeColor', colors(i,:), 'LineWidth', 1);

    %% --- Whiskers ---
    plot([i i], [whislo q1],    '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    plot([i i], [q3    whishi], '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);

    %% --- IQR box ---
    boxW = 0.08;
    boxH = max(iqr_val, eps);
    rectangle('Position', [i - boxW, q1, 2*boxW, boxH], ...
        'FaceColor', [1 1 1], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5);

    %% --- Median line ---
    plot([i - boxW, i + boxW], [med med], '-', ...
        'Color', [0.2 0.2 0.2], 'LineWidth', 2.5);

end

%% -----------------------------------------------------------------------
%  2. COMPUTE CLIFF'S DELTA (subsampled for efficiency)
%% -----------------------------------------------------------------------
dataSub = cellfun(@(x) x(randperm(numel(x), min(numel(x), maxN))), ...
    data, 'UniformOutput', false);

pairs    = nchoosek(1:nGroups, 2);
nPairs   = size(pairs, 1);
cliffs_d = zeros(nPairs, 1);

for k = 1:nPairs
    i  = pairs(k, 1);
    j  = pairs(k, 2);
    x1 = dataSub{i}(:);
    x2 = dataSub{j}(:);
    n1 = numel(x1);
    n2 = numel(x2);
    allXY       = [x1; x2];
    rnks        = tiedrank(allXY);
    U           = sum(rnks(1:n1)) - n1*(n1+1)/2;
    cliffs_d(k) = (2*U / (n1*n2)) - 1;
end

%% -----------------------------------------------------------------------
%  3. BUILD OUTPUT TABLE
%% -----------------------------------------------------------------------
comp_col   = cell(nPairs, 1);
cliffs_col = cell(nPairs, 1);
interp_col = cell(nPairs, 1);
notes_col  = repmat({'positive d = group i tends to be larger than group j'}, nPairs, 1);

for k = 1:nPairs
    i = pairs(k, 1);
    j = pairs(k, 2);
    comp_col{k}   = sprintf('%s vs %s', glabels{i}, glabels{j});
    cliffs_col{k} = sprintf('%.3f', cliffs_d(k));
    interp_col{k} = interpreteCliffsDelta(abs(cliffs_d(k)));
end

stats_table = table(comp_col, cliffs_col, interp_col, notes_col, ...
    'VariableNames', {'Comparison', 'Cliffs_d', 'Interpretation', 'Notes'});

property_col = repmat({tstr}, height(stats_table), 1);
stats_table  = addvars(stats_table, property_col, 'Before', 'Comparison', ...
    'NewVariableNames', 'Property');

disp(stats_table(:, {'Property', 'Comparison', 'Cliffs_d', 'Interpretation'}));
fprintf('\nCliffs delta: |d|<0.147 negligible, 0.147-0.33 small, 0.33-0.474 medium, >0.474 large\n');
fprintf('NOTE: voxel-level pseudoreplication — interpret as descriptive, not inferential\n');

%% -----------------------------------------------------------------------
%  4. BRACKETS WITH CLIFF'S DELTA
%% -----------------------------------------------------------------------
yRange      = yMax - yMin;
bracketPad  = yRange * 0.05;
bracketStep = yRange * 0.12;  % increased to accommodate text labels
bracketH    = yRange * 0.02;
currentY    = yMax + bracketPad;

for k = 1:nPairs
    i    = pairs(k, 1);
    j    = pairs(k, 2);
    d    = cliffs_d(k);
    dstr = sprintf('d=%.2f (%s)', d, interpreteCliffsDelta(abs(d)));

    plot([i, j], [currentY, currentY],            '-k', 'LineWidth', 1.2);
    plot([i, i], [currentY - bracketH, currentY], '-k', 'LineWidth', 1.2);
    plot([j, j], [currentY - bracketH, currentY], '-k', 'LineWidth', 1.2);
    text((i+j)/2, currentY + bracketH * 0.5, dstr, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', 9, ...
        'Color',    [0.3 0.3 0.3]);

    currentY = currentY + bracketStep;
end

%% -----------------------------------------------------------------------
%  5. AXES FORMATTING
%% -----------------------------------------------------------------------
xlim([0.5, nGroups + 0.5]);
ylim([ylims(1), max(ylims(2), currentY + bracketPad)]);
xticks(1:nGroups);
xticklabels(glabels);
ylabel(ystr);
title(tstr);
box off;
set(gca, 'FontSize', 13, 'LineWidth', 1.2);

hold off;
end


%% -----------------------------------------------------------------------
%  LOCAL HELPER: Cliff's delta interpretation
%% -----------------------------------------------------------------------
function s = interpreteCliffsDelta(d)
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