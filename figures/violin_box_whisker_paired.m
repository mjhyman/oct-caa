function [stats_table, fig] = violin_box_whisker_paired(data, glabels, tstr, ylims, ystr, varargin)
% violin_box_whisker_paired - Paired violin plots (EPVS / vessel) per group,
% with embedded boxplots, optional log-space rendering, and two families of
% effect sizes:
%
%   BETWEEN-GROUP  group i vs group j, within a structure (independent
%                  samples)  -> Cliff's delta
%   WITHIN-GROUP   EPVS vs vessel, within a group. If the two vectors are
%                  voxel-aligned (SWP), this is a PAIRED contrast and uses
%                  the paired dominance statistic d_paired = P(E>V) - P(E<V)
%                  plus the matched-pairs rank-biserial r. If they are not
%                  aligned (e.g. object volumes, different N), it falls back
%                  to unpaired Cliff's delta and says so in the table.
%
% Required inputs:
%   data     - nGroups x 2 cell array. Column 1 = EPVS, column 2 = vessel.
%   glabels  - 1 x nGroups cell of group names, e.g. {'Control','Severe'}
%   tstr     - plot title (also written to the Property column)
%   ylims    - y-axis limits [min max], IN THE UNITS OF `data`. With
%              'LogScale' true these are log-transformed for you.
%   ystr     - y-axis label
%
% Name-value options:
%   'Colors'       2 x 3 RGB, row 1 = EPVS, row 2 = vessel.
%                  Default #1964B0 (EPVS), #DB5829 (vessel).
%   'StructLabels' 1 x 2 cell of structure names. Default {'EPVS','Vessel'}
%   'Paired'       logical scalar or 1 x nGroups. True means column 1 and
%                  column 2 of that group are matched row-for-row.
%                  Default: auto-detected from equal vector lengths.
%   'LogScale'     true to work in log10 space (default false). The data are
%                  log10-transformed BEFORE the KDE and box statistics, so
%                  density and quartiles are estimated in log space; the y
%                  ticks are then relabelled in original units. This is not
%                  the same as set(gca,'YScale','log'), which would leave the
%                  KDE estimated in linear space and drawn on a warped axis.
%                  Non-positive values are dropped, with a warning.
%   'DateStamp'    true (default) to print the run date in the figure corner.
%                  Set false for final manuscript figures.
%
% Outputs:
%   stats_table - one row per comparison. Date_generated records the run date
%                 and Scale records linear vs log10, so exported sheets stay
%                 traceable. Export with writetable(stats_table, ...)
%   fig         - handle to the figure, for dated saving by the caller.
%
% Example:
%   T = violin_box_whisker_paired(volCell, {'Control','Severe'}, ...
%           'Object volume - front', [1e2 1e7], 'Object volume (\mum^3)', ...
%           'LogScale', true, 'Paired', false);

if ~exist('tiedrank', 'file')
    error('violin_box_whisker_paired requires the Statistics and Machine Learning Toolbox (tiedrank missing).');
end

%% --- Option parsing ---
ip = inputParser;
ip.FunctionName = 'violin_box_whisker_paired';
addParameter(ip, 'Colors',       []);
addParameter(ip, 'StructLabels', {'EPVS','Vessel'});
addParameter(ip, 'Paired',       []);
addParameter(ip, 'LogScale',     false);
addParameter(ip, 'DateStamp',    true);
parse(ip, varargin{:});

colors   = ip.Results.Colors;
slabels  = ip.Results.StructLabels;
isPaired = ip.Results.Paired;
useLog   = logical(ip.Results.LogScale);
showDate = logical(ip.Results.DateStamp);

if isempty(colors)
    colors = [hex2rgb01('#1964B0');   % EPVS
              hex2rgb01('#DB5829')];  % vessel
end
if iscell(tstr), tstr = tstr{1}; end

%% --- Run date ---
runDate = datestr(now, 'yyyy-mm-dd');      %#ok<*TNOW1,*DATST>

%% --- Shape check ---
if ~iscell(data) || size(data, 2) ~= 2
    error('data must be an nGroups x 2 cell array: column 1 = EPVS, column 2 = vessel.');
end
nGroups = size(data, 1);
nStruct = 2;

if nGroups < 1
    error('data must contain at least one group.');
end
if numel(glabels) ~= nGroups
    error('glabels has %d entries but data has %d groups.', numel(glabels), nGroups);
end
if size(colors, 1) < nStruct
    error('Colors must have at least 2 rows (EPVS, vessel).');
end

% double + column (ksdensity-style maths and tiedrank both need double)
data = cellfun(@(x) double(x(:)), data, 'UniformOutput', false);

%% --- Pairing: auto-detect unless told ---
autoPaired = false(1, nGroups);
for i = 1:nGroups
    autoPaired(i) = numel(data{i,1}) == numel(data{i,2});
end
if isempty(isPaired)
    isPaired = autoPaired;
elseif isscalar(isPaired)
    isPaired = repmat(logical(isPaired), 1, nGroups);
else
    isPaired = logical(isPaired(:))';
end
if numel(isPaired) ~= nGroups
    error('Paired must be scalar or have one entry per group.');
end
for i = 1:nGroups
    if isPaired(i) && ~autoPaired(i)
        error(['Group "%s" was declared paired but EPVS (n=%d) and vessel (n=%d) ' ...
               'differ in length.'], glabels{i}, numel(data{i,1}), numel(data{i,2}));
    end
end

%% --- Validity masking, then optional log10 transform ---
% Non-finite always dropped; non-positive additionally dropped on log scale.
% For paired groups the mask is applied JOINTLY so rows stay aligned.
for i = 1:nGroups
    v1 = isfinite(data{i,1});
    v2 = isfinite(data{i,2});
    if useLog
        v1 = v1 & (data{i,1} > 0);
        v2 = v2 & (data{i,2} > 0);
    end

    if isPaired(i)
        keep  = v1 & v2;
        nDrop = numel(keep) - nnz(keep);
        if nDrop > 0
            warning('violin_box_whisker_paired:droppedRows', ...
                'Group "%s": dropped %d of %d matched rows (non-finite%s).', ...
                glabels{i}, nDrop, numel(keep), tern(useLog, ' or non-positive', ''));
        end
        data{i,1} = data{i,1}(keep);
        data{i,2} = data{i,2}(keep);
    else
        reportDrop(glabels{i}, slabels{1}, v1, useLog);
        reportDrop(glabels{i}, slabels{2}, v2, useLog);
        data{i,1} = data{i,1}(v1);
        data{i,2} = data{i,2}(v2);
    end

    for s = 1:nStruct
        if isempty(data{i,s})
            error('Group "%s", %s is empty after removing invalid values.', ...
                  glabels{i}, slabels{s});
        end
    end

    if useLog
        data{i,1} = log10(data{i,1});
        data{i,2} = log10(data{i,2});
    end
end

%% --- Axis limits follow the same transform ---
if useLog
    if any(ylims <= 0)
        error('With LogScale true, ylims must be positive (they are given in data units).');
    end
    ylims = log10(ylims);
end

%% --- Layout constants ---
offset      = 0.22;
violinWidth = 0.18;
boxW        = 0.045;
nKDE        = 512;      % output points on the violin outline
nBins       = 2048;     % histogram bins for the binned KDE (see fastKDE)
maxN        = 100000;   % subsample cap for rank statistics only

xpos = zeros(nGroups, nStruct);
for i = 1:nGroups
    xpos(i,1) = i - offset;
    xpos(i,2) = i + offset;
end

fig = figure('Name', sprintf('%s  (%s)', tstr, runDate), 'NumberTitle', 'off');
set(fig, 'Units', 'inches', 'Position', [0, 0, 7, 5.5]);
set(fig, 'PaperUnits', 'inches', 'PaperSize', [7, 5.5], 'PaperPosition', [0, 0, 7, 5.5]);
hold on;

%% -----------------------------------------------------------------------
%  1. DRAW VIOLINS  (in log space when useLog, so the KDE is a log-space KDE)
%% -----------------------------------------------------------------------
yMax    = -Inf;
yMin    =  Inf;
hViolin = gobjects(nStruct, 1);

for i = 1:nGroups
    for s = 1:nStruct

        x   = data{i,s};
        xc  = xpos(i,s);

        % --- All five order statistics from ONE sorted pass ---
        qq      = quantile(x, [0.001 0.25 0.5 0.75 0.999]);
        xLo     = qq(1);
        q1      = qq(2);
        med     = qq(3);
        q3      = qq(4);
        xHi     = qq(5);
        iqr_val = q3 - q1;
        whislo  = max(xLo, q1 - 1.5 * iqr_val);
        whishi  = min(xHi, q3 + 1.5 * iqr_val);
        yMax = max(yMax, xHi);
        yMin = min(yMin, xLo);

        % --- KDE: binned, O(n) in the data instead of O(n * nKDE) ---
        [f, xi] = fastKDE(x, nKDE, xLo, xHi, iqr_val, nBins);
        f  = f / max(f) * violinWidth;

        % --- Violin patch ---
        h = patch([xc + f, fliplr(xc - f)], [xi, fliplr(xi)], colors(s,:), ...
            'FaceAlpha', 0.4, 'EdgeColor', colors(s,:), 'LineWidth', 1);
        if i == 1
            hViolin(s) = h;
        end

        % --- Whiskers, IQR box, median ---
        plot([xc xc], [whislo q1],    '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        plot([xc xc], [q3    whishi], '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        rectangle('Position', [xc - boxW, q1, 2*boxW, max(iqr_val, eps)], ...
            'FaceColor', [1 1 1], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5);
        plot([xc - boxW, xc + boxW], [med med], '-', ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 2.5);
    end
end

%% -----------------------------------------------------------------------
%  2. EFFECT SIZES
%% -----------------------------------------------------------------------
% All statistics here are rank-based, so log10 (monotonic) leaves Value
% unchanged. Median_diff is the exception: on a log scale it is a median
% log10 RATIO, not a difference. Flagged in the Notes column.
rs = RandStream('mt19937ar', 'Seed', 0);

needSub = false(nGroups, nStruct);
if nGroups >= 2
    needSub(:) = true;                 % every cell feeds a between-group test
end
for i = 1:nGroups
    if ~isPaired(i)
        needSub(i,:) = true;           % unpaired within-group fallback
    end
end

dataSub = cell(nGroups, nStruct);
for i = 1:nGroups
    for s = 1:nStruct
        if needSub(i,s)
            dataSub{i,s} = drawSub(rs, data{i,s}, maxN);
        end
    end
end

% --- 2a. Between-group Cliff's delta, within each structure ---
if nGroups >= 2
    pairs = nchoosek(1:nGroups, 2);
else
    pairs = zeros(0, 2);
end
nPairs   = size(pairs, 1);
cliffs_d = zeros(nPairs, nStruct);

for s = 1:nStruct
    for k = 1:nPairs
        cliffs_d(k,s) = cliffsDeltaUnpaired(dataSub{pairs(k,1), s}, ...
                                            dataSub{pairs(k,2), s});
    end
end

% --- 2b. Within-group EPVS vs vessel ---
within_d    = zeros(nGroups, 1);
within_r    = nan(nGroups, 1);
within_med  = nan(nGroups, 1);
within_stat = cell(nGroups, 1);

for i = 1:nGroups
    if isPaired(i)
        diffs = data{i,1} - data{i,2};
        within_d(i)   = (nnz(diffs > 0) - nnz(diffs < 0)) / numel(diffs);
        within_med(i) = median(diffs);

        dsub = drawSub(rs, diffs, maxN);
        dsub = dsub(dsub ~= 0);
        if isempty(dsub)
            within_r(i) = 0;
        else
            rk   = tiedrank(abs(dsub));
            Wpos = sum(rk(dsub > 0));
            Wneg = sum(rk(dsub < 0));
            within_r(i) = (Wpos - Wneg) / (Wpos + Wneg);
        end
        within_stat{i} = 'd_paired (matched dominance)';
    else
        within_d(i)    = cliffsDeltaUnpaired(dataSub{i,1}, dataSub{i,2});
        within_stat{i} = 'Cliffs d (UNPAIRED - samples not matched)';
    end
end

%% -----------------------------------------------------------------------
%  3. BUILD OUTPUT TABLE
%% -----------------------------------------------------------------------
scaleStr = tern(useLog, 'log10', 'linear');
medNote  = tern(useLog, 'positive = EPVS larger; Median_diff is a median log10 ratio (10^d = fold change)', ...
                        'positive = EPVS larger; benchmarks borrowed from Cliffs d');

Contrast = {}; Structure = {}; Comparison = {}; Statistic = {}; Interp = {}; Notes = {};
Value = []; RankBiserial = []; MedianDiff = []; N1 = []; N2 = [];

for s = 1:nStruct
    for k = 1:nPairs
        i = pairs(k,1); j = pairs(k,2);
        Contrast{end+1,1}     = 'between-group';   %#ok<*AGROW>
        Structure{end+1,1}    = slabels{s};
        Comparison{end+1,1}   = sprintf('%s vs %s', glabels{i}, glabels{j});
        Statistic{end+1,1}    = 'Cliffs d (unpaired)';
        Value(end+1,1)        = cliffs_d(k,s);
        RankBiserial(end+1,1) = NaN;
        MedianDiff(end+1,1)   = NaN;
        Interp{end+1,1}       = interpreteCliffsDelta(abs(cliffs_d(k,s)));
        N1(end+1,1)           = numel(data{i,s});
        N2(end+1,1)           = numel(data{j,s});
        Notes{end+1,1}        = 'positive d = first group larger; rank-based, unchanged by log';
    end
end

for i = 1:nGroups
    Contrast{end+1,1}     = 'within-group';
    Structure{end+1,1}    = sprintf('%s vs %s', slabels{1}, slabels{2});
    Comparison{end+1,1}   = glabels{i};
    Statistic{end+1,1}    = within_stat{i};
    Value(end+1,1)        = within_d(i);
    RankBiserial(end+1,1) = within_r(i);
    MedianDiff(end+1,1)   = within_med(i);
    Interp{end+1,1}       = interpreteCliffsDelta(abs(within_d(i)));
    N1(end+1,1)           = numel(data{i,1});
    N2(end+1,1)           = numel(data{i,2});
    if isPaired(i)
        Notes{end+1,1} = medNote;
    else
        Notes{end+1,1} = 'samples not matched; descriptive only';
    end
end

for k = 1:nPairs
    i = pairs(k,1); j = pairs(k,2);
    if isPaired(i) && isPaired(j)
        Contrast{end+1,1}     = 'within-group change';
        Structure{end+1,1}    = sprintf('%s vs %s', slabels{1}, slabels{2});
        Comparison{end+1,1}   = sprintf('%s minus %s', glabels{i}, glabels{j});
        Statistic{end+1,1}    = 'delta of d_paired';
        Value(end+1,1)        = within_d(i) - within_d(j);
        RankBiserial(end+1,1) = NaN;
        MedianDiff(end+1,1)   = within_med(i) - within_med(j);
        Interp{end+1,1}       = '';
        N1(end+1,1)           = numel(data{i,1});
        N2(end+1,1)           = numel(data{j,1});
        Notes{end+1,1}        = 'how much the EPVS-vessel gap shifts with group';
    end
end

nRows       = numel(Value);
Date_gen    = repmat({runDate},  nRows, 1);
Property    = repmat({tstr},     nRows, 1);
Scale       = repmat({scaleStr}, nRows, 1);
stats_table = table(Date_gen, Property, Scale, Contrast, Structure, Comparison, ...
                    Statistic, Value, RankBiserial, MedianDiff, Interp, N1, N2, Notes, ...
    'VariableNames', {'Date_generated','Property','Scale','Contrast','Structure', ...
                      'Comparison','Statistic','Value','RankBiserial_r','Median_diff', ...
                      'Interpretation','N_1','N_2','Notes'});

fprintf('\n%s  |  generated %s  |  scale: %s\n', tstr, runDate, scaleStr);
disp(stats_table(:, {'Contrast','Structure','Comparison','Value','Interpretation'}));
fprintf('\nCliffs delta benchmarks: |d|<0.147 negligible, 0.147-0.33 small, 0.33-0.474 medium, >0.474 large\n');
fprintf('NOTE: voxel-level pseudoreplication - interpret as descriptive, not inferential\n');

%% -----------------------------------------------------------------------
%  4. BRACKETS
%% -----------------------------------------------------------------------
yRange      = yMax - yMin;
bracketPad  = yRange * 0.05;
bracketStep = yRange * 0.12;
bracketH    = yRange * 0.02;
currentY    = yMax + bracketPad;
withinColor = [0.35 0.35 0.35];

for i = 1:nGroups
    if isPaired(i)
        lbl = sprintf('d_{paired}=%.2f', within_d(i));
    else
        lbl = sprintf('d=%.2f (unpaired)', within_d(i));
    end
    drawBracket(xpos(i,1), xpos(i,2), currentY, bracketH, withinColor, lbl);
    currentY = currentY + bracketStep;
end

for s = 1:nStruct
    for k = 1:nPairs
        i = pairs(k,1); j = pairs(k,2);
        lbl = sprintf('%s: d=%.2f (%s)', slabels{s}, cliffs_d(k,s), ...
                      interpreteCliffsDelta(abs(cliffs_d(k,s))));
        drawBracket(xpos(i,s), xpos(j,s), currentY, bracketH, colors(s,:), lbl);
        currentY = currentY + bracketStep;
    end
end

%% -----------------------------------------------------------------------
%  5. AXES FORMATTING
%% -----------------------------------------------------------------------
ax = gca;
xlim([0.5, nGroups + 0.5]);
ylim([ylims(1), max(ylims(2), currentY + bracketPad)]);
xticks(1:nGroups);
xticklabels(glabels);
ylabel(ystr);
title(tstr);
legend(hViolin, slabels, 'Location', 'northwest', 'Box', 'off');
box off;
set(ax, 'FontSize', 13, 'LineWidth', 1.2);

% Relabel the log10 positions in original units (10^k, 2x10^k, 5x10^k).
% Ticks stop at the top of the data so the bracket band stays unlabelled.
if useLog
    yl = ylim;
    [tickPos, tickLab] = logTickLabels(yl(1), yl(2));
    if numel(tickPos) >= 2
        set(ax, 'YTick', tickPos, 'YTickLabel', tickLab, ...
                'TickLabelInterpreter', 'tex');
    end
end

if showDate
    text(1, -0.085, runDate, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 7, 'Color', [0.55 0.55 0.55], 'Interpreter', 'none');
end

hold off;
end


%% -----------------------------------------------------------------------
%  LOCAL HELPERS
%% -----------------------------------------------------------------------
function [tickPos, tickLab] = logTickLabels(loLog, hiLog)
% Tick positions in log10 space, labelled in original units. Uses 1-2-5 per
% decade for narrow ranges and decades only for wide ones.
dLo = floor(loLog);
dHi = ceil(hiLog);
if (dHi - dLo) >= 4
    mant = 1;
else
    mant = [1 2 5];
end

tickPos = [];
tickLab = {};
for e = dLo:dHi
    for m = mant
        v = log10(m) + e;
        if v >= loLog && v <= hiLog
            tickPos(end+1) = v;   %#ok<*AGROW>
            if m == 1
                tickLab{end+1} = sprintf('10^{%d}', e);
            else
                tickLab{end+1} = sprintf('%d\\times10^{%d}', m, e);
            end
        end
    end
end
end

function [f, xi] = fastKDE(x, nPoints, lo, hi, iqrv, nBins)
% Binned Gaussian KDE. Cost is O(n) for the histogram plus O(nBins) for the
% smoothing, versus a direct kernel sum's O(n * nPoints) - the single largest
% cost in this function at millions of voxels.
%
% The outline is clipped to [lo, hi] (the 0.1/99.9 percentiles), so very
% heavy tails no longer stretch the violin past the whiskers.
n = numel(x);

if iqrv > 0
    sigma = iqrv / 1.349;
else
    sigma = std(x);
end
if sigma <= 0 || ~isfinite(sigma)
    sigma = eps;
end
bw = 0.9 * sigma * n^(-1/5);

pad     = 3 * bw;
edges   = linspace(lo - pad, hi + pad, nBins + 1);
binW    = edges(2) - edges(1);
centers = edges(1:end-1) + binW/2;

counts = histcounts(x, edges);

bw = max(bw, 2 * binW);   % keep the kernel at least a couple of bins wide

kHalf = ceil(4 * bw / binW);
kx    = (-kHalf:kHalf) * binW;
kern  = exp(-0.5 * (kx / bw).^2);
kern  = kern / sum(kern);

dens = conv(counts, kern, 'same');

xi = linspace(lo, hi, nPoints);
f  = interp1(centers, dens, xi, 'linear', 0);
f  = max(f, 0);
if max(f) <= 0
    f = ones(size(f));
end
end

function y = drawSub(rs, x, k)
% Uniform subsample with replacement. For n >> k the difference from
% sampling without replacement is negligible, and this is O(k), not O(n).
n = numel(x);
if n <= k
    y = x;
else
    y = x(randi(rs, n, k, 1));
end
end

function d = cliffsDeltaUnpaired(x1, x2)
x1 = x1(:); x2 = x2(:);
n1 = numel(x1); n2 = numel(x2);
rnks = tiedrank([x1; x2]);
U    = sum(rnks(1:n1)) - n1*(n1+1)/2;
d    = (2*U / (n1*n2)) - 1;
end

function drawBracket(x1, x2, y, h, col, lbl)
plot([x1, x2], [y, y],     '-', 'Color', col, 'LineWidth', 1.2);
plot([x1, x1], [y - h, y], '-', 'Color', col, 'LineWidth', 1.2);
plot([x2, x2], [y - h, y], '-', 'Color', col, 'LineWidth', 1.2);
text((x1+x2)/2, y + h*0.5, lbl, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 9, 'Color', col);
end

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

function rgb = hex2rgb01(hexstr)
hexstr = strrep(hexstr, '#', '');
rgb = double([hex2dec(hexstr(1:2)), hex2dec(hexstr(3:4)), hex2dec(hexstr(5:6))]) / 255;
end

function reportDrop(glab, slab, valid, useLog)
nDrop = numel(valid) - nnz(valid);
if nDrop > 0
    warning('violin_box_whisker_paired:droppedValues', ...
        'Group "%s", %s: dropped %d of %d values (non-finite%s).', ...
        glab, slab, nDrop, numel(valid), tern(useLog, ' or non-positive', ''));
end
end

function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end