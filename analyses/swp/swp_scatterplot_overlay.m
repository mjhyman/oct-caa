function swp_scatterplot_overlay(xy_cell, labels, colors, xlab, ylab, tit,...
                                xlims, ylims, plt_dir, fname)
% Scatterplot function for optical property vs SWP
% INPUTS
%   xy_cell {cell}: each cell contains [N x 2] xy pairs of subject/group
%   labels {cell}: legend labels for each subject/group
%   colors {cell}: color specs for each subject/group
%   xlab, ylab, tit (strings): axis labels and title
%   xlims, ylims (vector): axis limits
%   plt_dir, fname: figure save directory & filename

figure('Position', [956  -274   909   844]);
hold on;
for idx = 1:numel(xy_cell)
    if ~isempty(xy_cell{idx})
        scatter(xy_cell{idx}(:,1), xy_cell{idx}(:,2), ...
            50, colors{idx}, 'filled', 'DisplayName', labels{idx});
    end
end

xlabel(xlab);
ylabel(ylab);
title(tit);
set(gca, 'fontsize', 30);
lgd = legend('show');
lgd.FontSize = 14;
lgd.Location = 'southeast';
xlim(xlims); ylim(ylims);

% Save figure
fout = fullfile(plt_dir, fname);
saveas(gcf, fout, 'png')
pause(1)
close;

end