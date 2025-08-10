%% Scatterplot function w/o error bars
function scatter_op_vs_dist(x, ves_op, ~, epvs_op, ~,...
                            xlab, ylab, yl, tit, dir_out, fname)
% Scatter plot with error bars for optical property vs distance
% INPUTS:
%   - x (vector): x-axis positions
%   - ves_op (vector): optical propery for vessels
%   - ves_sem (vector): optical propery standard error for vessels
%   - epvs_op (vector): optical propery for EPVS
%   - epvs_sem (vector): optical propery standard error for EPVS
%   - xlab (str): x-axis label
%   - ylab (str): y-axis label
%   - yl (vector): y-axis limits
%   - tit (str): figure title
%   - dir_out (str): output directory
%   - fname (str): output filename

    % Init figure
    figure('Position', [100, 100, 900, 900],'Resize', 'off');
    
    % Scatterplot with error bars
    scatter(x,ves_op,100,'k','filled');
    hold on;
    scatter(x,epvs_op,100,'r','filled');
    
    % Limits, Labels, and legend
    ylim(yl);
    xlabel(xlab); ylabel(ylab); title(tit);
    % legend({'Vessels','EPVS'});
    
    % Font and fontsize
    fontname("SansSerif")
    set(gca,'fontsize',24);
    
    % Save output
    fout = fullfile(dir_out,fname);
    pause(1); saveas(gcf,fout);
end