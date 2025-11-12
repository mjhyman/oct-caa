function swp_scatterplot(xy, xlab, ylab, tit, xlims, ylims, plt_dir, fname)
% Scatterplot function for optical property vs SWP
% INPUTS
%   xy (matrix): xy pairs
%   xlab (string): x-axis label
%   ylab (string): y-axis label
%   tit (string): title string
%   xlims (vector): x-axis limits
%   ylims (vector): y-axis limits
%   plt_dir (string): output directory for saving figure
%   fname (string): output filename of figure

% Set figure size to half width and half height
figure('Position', [956  -274   909   844]);

% Plot figure
scatter(xy(:,1), xy(:,2), 60, 'b', 'filled');
xlabel(xlab);
ylabel(ylab);
title(tit);
set(gca, 'fontsize', 30);

% Set axis limits
xlim(xlims);
ylim(ylims);

% Save the figure
fout = fullfile(plt_dir, fname);
saveas(gcf,fout,'png')
pause(1)
close;

end