function T = corr_region_distance(hist,xl,yl,figdir)
% Measure correlation b/w histology and OCT optical property
% This function is called after "corr_histo_oct," which measures the
% optical properties in both the histology and OCT. This function will
% combine the measurements across subjects and then report the rho and p
% for the combined values for each region
%
%   INPUTS:
%       - hist (struct): contains image, epvs mask, tissue border mask
%       - xl (vector): x-axis limits
%       - yl (vector): y-axis limits
%   OUTPUTS:
%       - hist (struct):
%           - combined (struct):
%           - front (struct):
%           - occip (struct):
%           - each contains the following:
%               - EPVS, vessel
%                   - distance, rho, p

%% Iterate over the subjects in the hist struct

% Struct to store mean values at distance
stats = struct();

% Identify # subjects with frontal and occipital
nfront = 0;
noccip = 0;
for ii = 1:numel(hist)
    name = hist(ii).baseName;
    if contains(name,'_frontal_')
        nfront = nfront + 1;
    else
        noccip = noccip + 1;
    end
end

% Initialize vectors for storing mean values
nmeas = numel(hist(ii).meas);

% Iterate over distances
for ii = 1:nmeas
    % Initialize vectors for frontal
    front_histo_epvs = [];
    front_histo_ves = [];
    front_ret_epvs = [];
    front_ret_ves = [];
    % Initialize vectors for occipital
    occip_histo_epvs = [];
    occip_histo_ves = [];
    occip_ret_epvs = [];
    occip_ret_ves = [];

    % iterate over samples
    for j=1:numel(hist)
        % Retrieve filename to identify region
        name = hist(j).baseName;
        %%% Retrieve local variables
        histo_epvs = hist(j).meas(ii).histo_epvs_mean;
        histo_ves = hist(j).meas(ii).histo_ves_mean;
        ret_epvs = hist(j).meas(ii).ret_epvs_mean;
        ret_ves = hist(j).meas(ii).ret_ves_mean;
        
        %%% Add to appropriate array based on region
        if contains(name,'_frontal_')
            front_histo_epvs(end+1) = histo_epvs;
            front_histo_ves(end+1) = histo_ves;
            front_ret_epvs(end+1) = ret_epvs;
            front_ret_ves(end+1) = ret_ves;
        else
            occip_histo_epvs(end+1) = histo_epvs;
            occip_histo_ves(end+1) = histo_ves;
            occip_ret_epvs(end+1) = ret_epvs;
            occip_ret_ves(end+1) = ret_ves;
        end
    end

    % Combine frontal and occipital
    stats(ii).histo_epvs = mean([front_histo_epvs, occip_histo_epvs]);
    stats(ii).ret_epvs = mean([front_ret_epvs, occip_ret_epvs]);
    stats(ii).histo_ves = mean([front_histo_epvs, occip_histo_epvs]);
    stats(ii).ret_ves = mean([front_ret_ves, occip_ret_ves]);
    % Frontal
    stats(ii).front_histo_epvs = mean(front_histo_epvs);
    stats(ii).front_ret_epvs = mean(front_ret_epvs);
    stats(ii).front_histo_ves = mean(front_histo_ves);
    stats(ii).front_ret_ves = mean(front_ret_ves);
    % Occipital
    stats(ii).occip_histo_epvs = mean(occip_histo_epvs);
    stats(ii).occip_ret_epvs = mean(occip_ret_epvs);
    stats(ii).occip_histo_ves = mean(occip_histo_ves);
    stats(ii).occip_ret_ves = mean(occip_ret_ves);
end

%% Measure correlation b/w EPVS & ves vs. OCT

%%% Combined frontal + occipital
% EPVS
histo_epvs = [stats.histo_epvs]';
ret_epvs = [stats.ret_epvs]';
% Vessels
histo_ves = [stats.histo_ves]';
ret_ves = [stats.ret_ves]';
% EPVS Spearman's
[comb_epvs_rho,comb_epvs_p] = corr(histo_epvs,ret_epvs,'Type','Spearman',...
                    'Rows','pairwise');
n_comb_epvs = length(histo_epvs);
% Vessel Spearman's
[comb_ves_rho,comb_ves_p] = corr(histo_ves,ret_ves,'Type','Spearman',...
                    'Rows','pairwise');
n_comb_ves = length(histo_ves);
% Scatter plots
tit = 'Occipital & Frontal';
fname = fullfile(figdir,'occip_front_epvs_ves.jpg');
scat_plot(histo_epvs, ret_epvs, histo_ves, ret_ves, xl, yl,tit,fname)

%%% Frontal
% EPVS
histo_epvs = [stats.front_histo_epvs]';
ret_epvs = [stats.front_ret_epvs]';
% Vessels
histo_ves = [stats.front_histo_ves]';
ret_ves = [stats.front_ret_ves]';
% EPVS Spearman's
[front_epvs_rho,front_epvs_p] = corr(histo_epvs,ret_epvs,'Type','Spearman',...
                    'Rows','pairwise');
n_front_epvs = length(histo_epvs);
% Vessel Spearman's
[front_ves_rho,front_ves_p] = corr(histo_ves,ret_ves,'Type','Spearman',...
                    'Rows','pairwise');
n_front_ves = length(histo_ves);
% Scatter plots
tit = 'Frontal';
fname = fullfile(figdir,'front_epvs_ves.jpg');
scat_plot(histo_epvs, ret_epvs, histo_ves, ret_ves, xl, yl,tit,fname)

%%% Occipital
% EPVS
histo_epvs = [stats.occip_histo_epvs]';
ret_epvs = [stats.occip_ret_epvs]';
% Vessels
histo_ves = [stats.occip_histo_ves]';
ret_ves = [stats.occip_ret_ves]';
% EPVS Spearman's
[occip_epvs_rho,occip_epvs_p] = corr(histo_epvs,ret_epvs,'Type','Spearman',...
                    'Rows','pairwise');
n_occip_epvs = length(histo_epvs);
% Vessel Spearman's
[occip_ves_rho,occip_ves_p] = corr(histo_ves,ret_ves,'Type','Spearman',...
                    'Rows','pairwise');
n_occip_ves = length(histo_ves);
% Scatter plots
tit = 'Occipital';
fname = fullfile(figdir,'occip_epvs_ves.jpg');
scat_plot(histo_epvs, ret_epvs, histo_ves, ret_ves, xl, yl,tit,fname)

%%% Create Table
rho = [comb_epvs_rho;comb_ves_rho;...
    front_epvs_rho;front_ves_rho;...
    occip_epvs_rho;occip_ves_rho];

p = [comb_epvs_p;comb_ves_p;...
    front_epvs_p;front_ves_p;...
    occip_epvs_p;occip_ves_p];

n = [n_comb_epvs; n_comb_ves;...
    n_front_epvs; n_front_ves;...
    n_occip_epvs; n_occip_ves];
T = table(p, rho, n, 'VariableNames', {'p','rho','n'},...
    'RowNames',{'Comb EPVS','Comb Ves','Front EPVS','Front Ves',...
                'Occip EPVS','Occip Ves'});


    function scat_plot(histo_epvs, ret_epvs, histo_ves, ret_ves, ...
                        xl, yl,tit,fname)
        % Overlay EPVS and vessels
        figure;
        set(gcf, 'Position',  [100, 100, 800, 800]);
        s1 = scatter(histo_epvs,ret_epvs,100,'filled','o','r'); hold on;
        s2 = scatter(histo_ves,ret_ves,100,'filled','o','b');
        xlim(xl); ylim(yl); set(gca,'FontSize',20);
        xlabel({'Normalized Relative Density of Gallyas Stain','(Unitless)'});
        ylabel('Retardance (\circ)');
        legend([s1,s2],{'EPVS','Vessels'})
        title(tit)
        saveas(gcf,fname)

    end

end