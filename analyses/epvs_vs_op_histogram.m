%% Matlab code to plot optical properties vs. N EPVS


%% Retrieve the # EPVS per volume
% Create a structure array with the desired fields
top = struct();
top.caa6 = caa6;
top.caa17 = caa17;
top.caa22 = caa22;
top.caa25 = caa25;
top.caa26 = caa26;

% List of structure names
sub_ids = fieldnames(top);

% Struct to store the number of EPVS
n_epvs = struct();

% Iterate over each structure in the list
for i = 1:length(sub_ids)
	% Retrieve top-level names of struct
    sub = sub_ids{i};
    tmp = top.(sub);
    
	% Retrieve regions of this sample
	regions = fieldnames(top.(sub));

    % Access the regions directly by field names
    for ii = 1:length(regions)
		% Retrieve epvs
        epvs = tmp.(regions{ii}).epvs;
        % Label connected components
        [~, n] = bwlabel(epvs);
		% Add # EPVS to struct
		n_epvs.(sub).(regions{ii}) = n;
    end
end

%% Create box/whisker plots


