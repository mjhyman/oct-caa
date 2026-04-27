function swp_scatterplot_fitted(xy, xlab, ylab, tit, xl, yl, plt_dir, fname)
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
figure('Position',  [533 6 1551 983]);

% Plot figure
scatter(xy(:,1), xy(:,2), 50, 'filled','MarkerFaceColor','#1964B0');
xlabel(xlab);
ylabel(ylab);
title(tit);
set(gca, 'fontsize', 30);

% Fit polynomial to xy and overlay
p = polyfit(xy(:,1), xy(:,2), 4); % Fit a second-degree polynomial
yfit = polyval(p, xy(:,1)); % Evaluate the polynomial at the x values
hold on; % Retain current plot
plot(xy(:,1), yfit,'LineWidth',4,'Color','k'); % Overlay the fitted curve

% Set axis limits
xlim(xl);
ylim(yl);

% Save the figure
ax = gca;
fname_pdf = strcat(fname,'.pdf');
fname_png = strcat(fname,'.png');
exportgraphics(ax,fullfile(plt_dir,fname_pdf),...
               'ContentType','vector','Resolution',600);
exportgraphics(ax,fullfile(plt_dir,fname_png),...
               'Resolution',600);
pause(1)
close;
end