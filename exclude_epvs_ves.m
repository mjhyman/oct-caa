function [seg_no_epvs] = exclude_epvs_ves(seg, epvs)
% EXCLUDE_EPVS returns the vessels devoid of EPVS
% Find overlap w/ EPVS, exclude these vessels.
%
% INPUTS:
%   seg (matrix): segmentation matrix
%   epvs (matrix): epvs matrix
% OUTPUTS:
%   seg_no_epvs (matrix): segmentation matrix devoid of vessels with EPVS

%%% Find pixels in vessels and EPVS
% discrete vessels
cc_seg = bwconncomp(seg, 26);
seg_lst = cc_seg.PixelIdxList;
% discrete EPVS
cc_epvs = bwconncomp(epvs, 26);
epvs_lst = cc_epvs.PixelIdxList;

%%% Initialize array to store indices of vessels to include
ves_include = true(length(seg_lst), 1);

%%% Create a linear index array for fast comparison
epvs_idx = unique(cell2mat(epvs_lst'));

%%% Find groups with overlap
for ii = 1:length(seg_lst)
    % Retrieve discrete vessel
    ves_tmp = seg_lst{ii};
    
    % Check for intersection using logical indexing
    overlap = any(ismember(ves_tmp, epvs_idx));
    
    % If any overlap is found, exclude this vessel
    if overlap
        ves_include(ii) = false;
    end
end

% Retrieve the vessel indices to keep
seg_idx = seg_lst(ves_include);

% Initialize a matrix of all zeros (the same size as the segmentation matrix)
seg_no_epvs = false(size(seg));

% Use logical indexing to assign the voxels to keep as one
for ii = 1:length(seg_idx)
    seg_no_epvs(seg_idx{ii}) = true;
end

end
