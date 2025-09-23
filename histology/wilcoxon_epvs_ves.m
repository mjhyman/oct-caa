function [stats] = wilcoxon_epvs_ves(histo,tail)
% Extract the average of the measurements surrounding the EPVS and control
% and perform Wilcoxon Ranked Sign Test. Measure from the edge of the EPVS
% or vessel to the first radius.
%   INPUTS
%       histo (struct): pathology struct. Each number entry corresponds to
%                       a different tissue section.
%       tail (str): the type of test (both, left, right)
%   OUTPUTS:
%       stats (cell): summary statistics for combine, front, occip,
%                   these include: p, h, w, z, median experimental, median
%                   control. each entry in the cell array is a table.

%% Extract EPVS and control measurements

%%% Determine number of frontal and occipital measurements
% Number of tissue sections for stain
nsec = numel(histo);
% Initialize counters for frontal and occipital
nfront = 0;
noccip = 0;
for ii = 1:nsec
    name = histo(ii).baseName;
    if contains(name,'_1_') || endsWith(name,'_1')
        nfront = nfront + 1;
    elseif contains(name,'_7_') || endsWith(name,'_7')
        noccip = noccip + 1;
    end
end
    
%%% Create vectors for storing experimental and control
% Vectors for combined frontal and occipital
epvs = [];
ctl = [];
% Frontal
epvs_f = [];
ctl_f = [];
% Occipital
epvs_o = [];
ctl_o = [];

%%% Iterate over number of radii measurements
for ii=1:nsec
    %%% Count number of radii measurements
    % Open all measurements from current tissue section
    data = histo(ii);
    % Determine whether frontal (_1_) or occipital (_7_)
    name = histo(ii).baseName;

    %%% Retrieve both frontal and occipital
    if contains(name,'_1_') || endsWith(name,'_1') || ...
            contains(name,'_7_') || endsWith(name,'_7')
        % Extract first sample surrounding EPVS
        epvs = [epvs, data.exp_mean];
        ctl = [ctl, data.ctl_mean];
    end

    %%% Separate by frontal and occipital
    if contains(name,'_1_') || endsWith(name,'_1')
        epvs_f = [epvs_f, data.exp_mean];
        ctl_f = [ctl_f, data.ctl_mean];
    elseif contains(name,'_7_') || endsWith(name,'_7')
        epvs_o = [epvs_o, data.exp_mean];
        ctl_o = [ctl_o, data.ctl_mean];
    end
end
%%% Perform Wilcoxon tests
% combined frontal and occipital
[comb] = wilco(epvs,ctl,tail);
% Measure frontal
[front] = wilco(epvs_f,ctl_f,tail);
% Measure occipital
[occip] = wilco(epvs_o,ctl_o,tail);

%%% Write to a table
T = table( ...
    [comb.p; front.p; occip.p], ...
    [comb.h; front.h; occip.h], ...
    [comb.w; front.w; occip.w], ...
    [comb.z; front.z; occip.z], ...
    [comb.n; front.n; occip.n], ...
    [comb.r; front.r; occip.r], ...
    [comb.exp_med; front.exp_med; occip.exp_med], ...
    [comb.ctl_med; front.ctl_med; occip.ctl_med], ...
    'VariableNames',{'p','h','w','z','n','effect_size','exp_med','ctl_med'},...
    'RowNames',{'comb','front','occip'});
stats = T;
end