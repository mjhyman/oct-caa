%% Histology stain vs EPVS-SWP, stratified by vessel-SWP percentile:
%% group-level (control / severe) alongside subject-level fits
% Overview:
%{
- loads the group-level struct  histo_gam      (GAM_histo_<tag>_<date>.mat)
  and the subject-level struct  histo_gam_subj (GAM_histo_subjects_<tag>_<date>.mat),
  both written by create_histo_swp_gam.m
- for each stain x region, draws ONE row of panels:
      [control (group)] [severe (group)] [CAA26] [CAA6] [CAA17] [CAA22] [CAA25]
  each showing fitted z-score vs EPVS-SWP at the 10th / 50th / 90th
  vessel-SWP percentile, same colours as the OCT slice figures:
      blue = 10th, purple = 50th, orange = 90th
- the group surfaces are evaluated on the SAME shared grid the subjects were
  evaluated on (stored in each subject entry), so group and subject curves are
  directly comparable; group support = union of its subjects' support masks
- writes a PNG/FIG/SVG per stain x region plus a CSV of every plotted slice
%}

%% Initialization
clc; close all;

%%% Directories (SCC)
ddir   = '/projectnb/npbssmic/ns/CAA/histology/gam';
figdir = '/projectnb/npbssmic/ns/CAA/figures/histology_gam/subject_vs_group';
if ~exist(figdir, 'dir'), mkdir(figdir); end

%%% Input files
group_file = 'GAM_histo_manuscript_18-Oct-2026.mat';
subj_file  = 'GAM_histo_subjects_manuscript_18-Oct-2026.mat';

%%% Subjects: display order, stage, and which group each belongs to
subj_order = {'CAA26','CAA6','CAA17','CAA22','CAA25'};
subj_stage = struct('CAA26','control','CAA6','control','CAA17','moderate', ...
                    'CAA22','severe','CAA25','severe');
group_members = struct('control', {{'CAA26','CAA6'}}, ...
                       'severe',  {{'CAA22','CAA25'}});

%%% Vessel-SWP percentiles and colours (blue = 10th, purple = 50th, orange = 90th)
ves_pcts     = [10 50 90];
slice_colors = [0.000 0.447 0.741;    % blue   - 10th
                0.494 0.184 0.556;    % purple - 50th
                0.850 0.325 0.098];   % orange - 90th

%%% Axis limits per stain (y) ; x is set from the region's shared grid
ylims_by_stain = struct('lhe',[-0.5 1], 'gfap',[-1 1], 'cd68',[-0.5 0.5]);
default_ylims  = [-1.5 1.5];

%%% Figure settings
fsize = 14;

%% Load
fprintf('Loading group struct:   %s\n', group_file)
G = load(fullfile(ddir, group_file));   G = G.histo_gam;
fprintf('Loading subject struct: %s\n', subj_file)
S = load(fullfile(ddir, subj_file));    S = S.histo_gam_subj;

region_map = struct('front','frontal', 'occip','occipital');
rcodes     = fieldnames(region_map);
stains     = fieldnames(S);

%% Plot
for ii = 1:numel(stains)
    st = stains{ii};
    if isfield(ylims_by_stain, st), ylims = ylims_by_stain.(st); else, ylims = default_ylims; end

    for rr = 1:numel(rcodes)
        rcode = rcodes{rr};
        reg   = region_map.(rcode);

        %%% Collect this stain x region's subject entries
        subj_entries = struct();
        for s = 1:numel(subj_order)
            key = sprintf('%s_%s', subj_order{s}, rcode);
            if isfield(S.(st), key), subj_entries.(subj_order{s}) = S.(st).(key); end
        end
        present = fieldnames(subj_entries);
        if isempty(present)
            fprintf('  %s %s: no subject-level fits, skipping\n', st, reg); continue
        end

        %%% Shared grid (identical across subjects in a region by construction)
        E0      = subj_entries.(present{1});
        epv_vec = E0.epv_vec(:)';  ves_vec = E0.ves_vec(:)';
        [EPV, VES] = meshgrid(epv_vec, ves_vec);
        Xg = [EPV(:), VES(:), EPV(:).*VES(:)];     % [epv, ves, interaction]
        xlims = [epv_vec(1) epv_vec(end)];

        %%% Group surfaces on the same grid, masked to union of member support
        grp_names = fieldnames(group_members);
        grp_entries = struct();
        for g = 1:numel(grp_names)
            gname = grp_names{g};
            gkey  = sprintf('%s_%s', gname, rcode);          % e.g. control_front
            if ~isfield(G, st) || ~isfield(G.(st), gkey)
                fprintf('  %s %s: no %s group fit\n', st, reg, gname); continue
            end
            mdl = get_gam_model(G.(st).(gkey));
            if isempty(mdl)
                warning('%s %s %s: no RegressionGAM found in group entry; skipping group panel', st, reg, gname);
                continue
            end
            F = reshape(predict(mdl, Xg), size(EPV));
            SUP = false(size(EPV));
            for m = 1:numel(group_members.(gname))
                sid = group_members.(gname){m};
                if isfield(subj_entries, sid), SUP = SUP | logical(subj_entries.(sid).support); end
            end
            if ~any(SUP(:)), SUP = true(size(EPV)); end     % no members present: show all
            F(~SUP) = NaN;
            grp_entries.(gname) = struct('fit', F, 'support', SUP);
        end

        %%% Figure: group panels first, then subjects
        panels = [cellfun(@(g) ['group:' g], fieldnames(grp_entries)', 'UniformOutput', false), ...
                  cellfun(@(s) ['subj:'  s], subj_order,                  'UniformOutput', false)];
        npan = numel(panels);

        figure('Position', [50 50 320*npan 420], 'Color', 'w');
        tl = tiledlayout(1, npan, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, sprintf('%s, %s: stain vs EPVS-SWP by vessel-SWP percentile', upper(st), reg), ...
              'FontSize', fsize + 2);

        rows = {};
        for p = 1:npan
            tok   = split(panels{p}, ':');  kind = tok{1}; name = tok{2};
            nexttile; hold on;

            switch kind
                case 'group'
                    F = grp_entries.(name).fit;  Flo = []; Fhi = [];
                    ttl = sprintf('%s (group)', name);  stage = name;
                case 'subj'
                    if ~isfield(subj_entries, name)
                        text(0.5, 0.5, sprintf('%s\nno %s volume', name, reg), 'Units','normalized', ...
                             'HorizontalAlignment','center', 'Color',[0.5 0.5 0.5], 'FontSize', fsize);
                        axis off; continue
                    end
                    E = subj_entries.(name);
                    F = E.fit;
                    if isfield(E,'do_bootstrap') && E.do_bootstrap, Flo = E.ci_lo; Fhi = E.ci_hi; else, Flo = []; Fhi = []; end
                    stage = subj_stage.(name);
                    ttl   = sprintf('%s (%s)', name, stage);
            end

            % Percentile rows chosen by POSITION among the in-support rows, so the
            % three slices are always distinct grid rows (value-based nearest-row
            % lookup can collapse two percentiles onto one row when support is
            % sparse, which hides the earlier-drawn trace under the later one).
            sup_rows = find(any(~isnan(F), 2));
            if isempty(sup_rows), axis off; continue; end
            pos = round(1 + (numel(sup_rows) - 1) * ves_pcts / 100);
            iy_all = sup_rows(pos);
            if numel(unique(iy_all)) < numel(iy_all)
                warning('%s %s %s: fewer than 3 distinct in-support vessel-SWP rows; slices overlap', st, reg, name);
            end

            % Draw the slice with the most in-support points first so a narrow
            % slice is never covered by a wide one.
            npts  = arrayfun(@(r) sum(~isnan(F(r,:))), iy_all);
            [~, draw_order] = sort(npts, 'descend');
            hl = gobjects(1, numel(ves_pcts));

            for q = draw_order(:)'
                iy = iy_all(q);
                y  = F(iy, :);
                m = ~isnan(y);
                if ~isempty(Flo)
                    lo = Flo(iy,:); hi = Fhi(iy,:);
                    mm = m & ~isnan(lo) & ~isnan(hi);
                    fill([epv_vec(mm) fliplr(epv_vec(mm))], [lo(mm) fliplr(hi(mm))], ...
                         slice_colors(q,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                end
                hl(q) = plot(epv_vec(m), y(m), '-', 'Color', slice_colors(q,:), 'LineWidth', 2.2, ...
                     'DisplayName', sprintf('ves-SWP %dth', ves_pcts(q)));

                for ix = find(m)
                    rows(end+1,:) = {kind, name, stage, st, reg, ves_pcts(q), ves_vec(iy), epv_vec(ix), y(ix)}; %#ok<AGROW>
                end
            end

            yline(0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
            xlim(xlims); ylim(ylims);
            title(ttl, 'FontSize', fsize);
            xlabel('EPVS-SWP');
            if p == 1
                ylabel(sprintf('%s (z-score)', upper(st)));
                legend(hl(isgraphics(hl)), 'Location', 'best');   % legend in 10/50/90 order
            end
            set(gca, 'FontSize', fsize, 'FontName', 'Arial');
        end

        outstem = fullfile(figdir, sprintf('histo_gam_%s_%s_subject_vs_group', st, reg));
        saveas(gcf, outstem, 'png');
        saveas(gcf, outstem, 'fig');
        exportgraphics(gcf, [outstem '.svg'], 'Resolution', 600);
        close(gcf);

        if ~isempty(rows)
            T = cell2table(rows, 'VariableNames', ...
                {'level','name','stage','stain','region','ves_pct','ves_swp','epv_swp','fit'});
            writetable(T, [outstem '.csv']);
        end
        fprintf('  Saved %s\n', outstem);
    end
end
fprintf('Done.\n')


%% ------------------------------------------------------------------------
%% Local function: pull the RegressionGAM out of a fit_swp_histo_gam result
function mdl = get_gam_model(stored)
mdl = [];
if isa(stored, 'RegressionGAM'), mdl = stored; return; end
if ~isstruct(stored), return; end
fn = fieldnames(stored);
for i = 1:numel(fn)
    v = stored.(fn{i});
    if isa(v, 'RegressionGAM'), mdl = v; return; end
    if isstruct(v)                                   % look one level deeper
        fn2 = fieldnames(v);
        for m = 1:numel(fn2)
            if isa(v.(fn2{m}), 'RegressionGAM'), mdl = v.(fn2{m}); return; end
        end
    end
end
end