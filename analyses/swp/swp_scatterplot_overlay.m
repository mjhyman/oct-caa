function swp_scatterplot_overlay(xy_cell, se_cell, labels, colors, xlab, ylab, tit,...
                                xlims, ylims, plt_dir, fname)
% Scatterplot function for optical property vs SWP
% INPUTS
%   xy_cell {cell}: each cell contains [N x 2] xy pairs of subject/group
%   labels {cell}: legend labels for each subject/group
%   colors {cell}: color specs for each subject/group (hex strings, e.g. '#3A7DC9')
%   xlab, ylab, tit (strings): axis labels and title
%   xlims, ylims (vector): axis limits
%   plt_dir, fname: figure save directory & filename

figure('Position', [533 6 1551 983]);
hold on;
for idx = 1:numel(xy_cell)
    % Convert hex to RGB
    c = sscanf(colors{idx}(2:end), '%2x%2x%2x', [1 3]) / 255;
    
    % Set NaN in SE to 0. This is artifact from taking the std of a single
    % data point, which results in #/0 -> NaN
    se_cell{idx}(isnan(se_cell{idx})) = 0;

    if ~isempty(xy_cell{idx})
        % Plot shaded SE region first so line draws on top
        x_fill = [xy_cell{idx}(:,1); flipud(xy_cell{idx}(:,1))];
        y_upper = xy_cell{idx}(:,2) + se_cell{idx}(:);
        y_lower = xy_cell{idx}(:,2) - se_cell{idx}(:);
        y_fill = [y_upper; flipud(y_lower)];
        fill(x_fill, y_fill, c, ...
            'FaceAlpha', 0.2, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % Plot main line on top
        plot(xy_cell{idx}(:,1), xy_cell{idx}(:,2),...
            'DisplayName', labels{idx},...
            'Color', c,'Linewidth',4);
    end
end

% Labels and figure properties
xlabel(xlab);
ylabel(ylab);
title(tit);
% legend('Location', 'best');
set(gca, 'fontsize', 50);
xlim(xlims); ylim(ylims);

% Save figure
grid on; box on;
ax = gca;
ax.GridAlpha = 1;
fname_pdf = strcat(fname,'.pdf');
fname_png = strcat(fname,'.png');
exportgraphics(ax,fullfile(plt_dir,fname_pdf),...
               'ContentType','vector','Resolution',600);
exportgraphics(ax,fullfile(plt_dir,fname_png),...
               'Resolution',600);
pause(1)
close;
end