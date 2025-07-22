function scatter_histo_oct(hist,figdir)
%SCATTER_HISTO_OCT Scatter plots of retardance vs. histology
% This function is run after calling "corr_histo_oct"
%
%   INPUTS:
%       - hist (struct): histology structure
%       - xl (vector): x-axis limits
%       - yl (vector): y-axis limits
%       - figdir (str): output path

% Find all limits for x-axis and y-axis
x_epvs = [hist.histo_epvs];
y_epvs = [hist.ret_epvs];
x_ves = [hist.histo_ves];
y_ves= [hist.ret_epvs];
xmin = min([x_epvs(:);x_ves(:)]);
xmax = max([x_epvs(:);x_ves(:)]);
ymin = min([y_epvs(:);y_ves(:)]);
ymax = max([y_epvs(:);y_ves(:)]);
% Set the limits
xmin = xmin - 0.01;
xmax = xmax + 0.01;
ymin = ymin - 0.1;
ymax = ymax + 0.1;

% Iterate sections
for ii = 1:length(hist)
    % Extract vectors of data
    histo_epvs = hist(ii).histo_epvs;
    histo_ves = hist(ii).histo_ves;
    ret_epvs = hist(ii).ret_epvs;
    ret_ves = hist(ii).ret_ves;
    % Extract stats
    epvs_r = hist(ii).epvs_r;
    epvs_p = hist(ii).epvs_p;
    ves_r = hist(ii).ves_r;
    ves_p = hist(ii).ves_p;
    % Plot
    tit = hist(ii).baseName;
    fname = fullfile(figdir,strcat(tit,'_epvs_ves.jpg'));
    xl = [xmin,xmax];
    yl = [ymin, ymax];
    scat_plot(histo_epvs, ret_epvs, histo_ves, ret_ves, ...
              epvs_r,epvs_p,ves_r,ves_p,xl,yl,tit,fname)
end

function scat_plot(histo_epvs, ret_epvs, histo_ves, ret_ves, ...
                   epvs_r, epvs_p, ves_r, ves_p, xl,yl,tit,fname)
    % Overlay EPVS and vessels
    figure;
    set(gcf, 'Position',  [100, 100, 800, 800]);
    s1 = scatter(histo_epvs,ret_epvs,100,'filled','o','r'); hold on;
    s2 = scatter(histo_ves,ret_ves,100,'filled','o','b');
    set(gca,'FontSize',20);
    xlabel({'Normalized Relative Density of Gallyas Stain','(Unitless)'});
    ylabel('Retardance (\circ)');
    % xlim(xl); ylim(yl);
    % legend([s1,s2],{'EPVS','Vessels'})
    title(tit,'Interpreter','none')

    % Add text box of Spearman's rho statistics
    dim = [0.15 0.15 0.3 0.1];
    str1 = ['EPVS: p = ' num2str(epvs_p) ',  rho = ' num2str(epvs_r)];
    str2 = ['Vessel: p = ' num2str(ves_p) ',  rho = ' num2str(ves_r)];
    str = {str1,str2};
    a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    a.FontSize = 20;
    saveas(gcf,fname)
    pause(1);
    close;
end

end