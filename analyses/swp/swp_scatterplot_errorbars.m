function swp_scatterplot_errorbars(xy, errors, xlab, ylab, tit,...
                                    xlims, ylims, plt_dir, fname)
% Scatterplot function for optical property vs SWP
% INPUTS
%   xy (matrix): xy pairs for optical property vs. SWP
%   err (vector): standard error at each point
%   xlab (string): x-axis label
%   ylab (string): y-axis label
%   tit (string): title string
%   xlims (vector): x-axis limits
%   ylims (vector): y-axis limits
%   plt_dir (string): output directory for saving figure
%   fname (string): output filename of figure

%%% Scatterplot of op vs. SWP
% Set figure size to half width and half height
figure('Position', [956  -274   909   844]);
% Scatterplot
scatter(xy(:,1), xy(:,2), 60, 'k', 'filled');
% errorbar(xy(:,1), xy(:,2), errors);
hold on;

%%% Error bar patch

% Create linear space
x = xy(:,1);
y = xy(:,2);
xq = linspace(min(x), max(x), size(x,1));
yq = interp1(x, y, xq, 'pchip');
errors_interp = interp1(x, errors, xq, 'pchip');
% Plot error bars
error_upper = yq + errors_interp;
error_lower = yq - errors_interp;
patch([xq fliplr(xq)], [error_upper fliplr(error_lower)],...
       'b','FaceAlpha', 1, 'EdgeColor', 'none');
%}

% Finalize plot
hold off;
xlabel(xlab);
ylabel(ylab);
title(tit);
set(gca, 'fontsize', 30);

% Set axis limits
xlim(xlims);
ylim(ylims);

% Save the figure
fout = fullfile(plt_dir, fname);
saveas(gcf,fout,'pdf')
saveas(gcf,fout,'fig')
pause(1)
close;

end