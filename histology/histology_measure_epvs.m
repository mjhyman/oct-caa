function [outputArg1,outputArg2] = histology_measure_epvs(hist)
%HISTOLOGY_MEASURE_DONUT Measure controls + parenchyma around EPVS.
%   INPUTS:
%       - hist (struct): contains image, epvs mask, tissue border mask
%   OUTPUTS:
%       - hist (struct):
%           - exp (matrix): experimental measurements (around EPVS)
%           - ctl (matrix): control measurements (around EPVS)


%% Iterate over the subjects in the hist struct

nsub = length(fields(hist));

for ii = 1:nsub

%%% Measure parenchyma around EPVS


%%% Measure control 


end


end

