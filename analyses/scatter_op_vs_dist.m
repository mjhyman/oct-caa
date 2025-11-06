%% Scatterplot function w/o error bars
function scatter_op_vs_dist(x, ves_op, ves_sem, epvs_op, epvs_sem, err_flag,...
                            xlab, ylab, yl, xt, yt, tit, dir_out, fname,...
                            psize)
% Scatter plot with error bars for optical property vs distance
% INPUTS:
%   - x (vector): x-axis positions
%   - ves_op (vector): optical propery for vessels
%   - ves_sem (vector): optical propery standard error for vessels
%   - epvs_op (vector): optical propery for EPVS
%   - epvs_sem (vector): optical propery standard error for EPVS
%   - err_flag (bool): true = errorbar plotting
%   - xlab (str): x-axis label
%   - ylab (str): y-axis label
%   - yl (vector): y-axis limits
%   - xt (array): x-axis ticks
%   - yt (array): y-axis ticks
%   - tit (str): figure title
%   - dir_out (str): output directory
%   - fname (str): output filename
%   - psize (int): point size

% Init figure
fig = figure('Position', [100, 100, 1000, 1000],'Resize', 'off');

% Scatterplot with error bars
if err_flag
    h1 = errorbar(x,ves_op,ves_sem,'k'); hold on;
    h2 = errorbar(x,epvs_op,epvs_sem,'r');
    set(h1,'MarkerSize',psize);
    set(h2,'MarkerSize',psize);
else
    scatter(x,ves_op,psize,'k','filled'); hold on;
    scatter(x,epvs_op,psize,'r','filled');
end

% y-axis limits, ticks, etc.
ylim(yl);
yticks(yt);
ytickformat('%.1f');
ylabel(ylab);

% x-axis label
xlabel(xlab);
xticks(xt);
xtickformat('%.0f');

% plot title
title(tit);

% Font and fontsize
fontname(fig, "Helvetica")
set(gca,'fontsize',40);

% Save output
fout = fullfile(dir_out,fname);
pause(1); saveas(gcf,fout); close;
end