%% En-face views of each caa volume to check anatomical orientation
% For every subject x region: mean projections of mus along each of the
% three array dimensions, with the grey-matter mask outlined in red and the
% matrix origin (1,1) marked. Images are drawn exactly as stored (row 1 at
% the top, column 1 at the left), so you can read off which array edge is
% superior and which dimension is the sectioning (depth) axis.
%
% Output: one PNG per subject x region in a timestamped folder.

clearvars -except caa
close all;

out_base = '/projectnb/npbssmic/ns/CAA/figures/orientation_check/';
run_dir  = fullfile(out_base, char(datetime('now','Format','yyyyMMdd_HHmmss')));
mkdir(run_dir);

%% Load caa structs (skipped if already in workspace)
if ~exist('caa','var')
    mat_dir = '/projectnb/npbssmic/ns/CAA/';
    files = struct('caa6', "/caa6/caa6.mat", 'caa17', "/caa17/occip/caa17.mat", ...
                   'caa22',"/caa22/caa22.mat", 'caa25',"/caa25/caa25.mat", ...
                   'caa26',"/caa26/caa26.mat");
    fn = fieldnames(files);
    caa = struct();
    for s = 1:numel(fn)
        fprintf('Loading %s\n', upper(fn{s}));
        tmp = load(fullfile(mat_dir, files.(fn{s})));
        caa.(fn{s}) = tmp.(fn{s});
        clear tmp
    end
end

%% Plot
subs = fieldnames(caa);
for s = 1:numel(subs)
    sub = subs{s};
    regs = intersect(fieldnames(caa.(sub)), {'front','occip'}, 'stable');
    for r = 1:numel(regs)
        reg = regs{r};
        V = caa.(sub).(reg);
        if isfield(V, 'mus')
            I = single(V.mus); ilab = '\mu_s';
        else
            I = single(V.mask_wm); ilab = 'WM mask';
        end
        if isfield(V, 'mask_gm'), gm = logical(V.mask_gm); else, gm = false(size(I)); end
        fprintf('%s %s: size %s\n', upper(sub), reg, mat2str(size(I)));

        fh = figure('Color','w', 'Position',[50 50 1600 560]);
        tl = tiledlayout(fh, 1, 3, 'TileSpacing','compact', 'Padding','compact');
        for d = 1:3
            P = squeeze(mean(I, d, 'omitnan'));     % rows/cols = remaining dims, in order
            G = squeeze(any(gm, d));
            other = setdiff(1:3, d);

            ax = nexttile(tl);
            imagesc(ax, P);
            axis(ax, 'image'); colormap(ax, gray);
            lim = prctile(P(:), [1 99]);
            if lim(2) > lim(1), caxis(ax, lim); end
            hold(ax, 'on');
            if any(G(:)) && all(size(G) > 1)
                contour(ax, double(G), [0.5 0.5], 'r', 'LineWidth', 1);
            end
            plot(ax, 1, 1, 'o', 'MarkerFaceColor','y', 'MarkerEdgeColor','k', 'MarkerSize',9);
            text(ax, 1, 1, '  (1,1)', 'Color','y', 'FontWeight','bold', ...
                 'VerticalAlignment','top');

            ylabel(ax, sprintf('dim %d (row index, increases \\downarrow)', other(1)));
            xlabel(ax, sprintf('dim %d (column index, increases \\rightarrow)', other(2)));
            title(ax, sprintf('Projection along dim %d   [%d \\times %d]', ...
                  d, size(P,1), size(P,2)), 'FontWeight','normal');
        end
        title(tl, sprintf('%s %s  |  mean %s projections, grey matter outlined in red', ...
              upper(sub), reg, ilab), 'FontSize', 12);

        exportgraphics(fh, fullfile(run_dir, sprintf('orientation_%s_%s.png', sub, reg)), ...
                       'Resolution', 150);
        close(fh);
        clear V I gm
    end
end
fprintf('Saved to %s\n', run_dir);