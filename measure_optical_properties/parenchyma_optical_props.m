function [mus_outer, ret_outer, ori_outer,...
        mus_inner, ret_inner, ori_inner,...
        parench_outer, parench_inner] = ...
    parenchyma_optical_props(seg, mus, ret, ori, se1, se2)
% Measure optical properties in parenchyma of all vessels
%   INPUTS:
%       seg (logical mat): vessels
%       mus (single mat): scattering coefficient matrix
%       ret (single mat): retardance matrix
%       ori (single mat): orientation matrix
%       se1 (handle): structuring element for inner dilation
%       se2 (handle): structuring element for outer dilation
%   OUTPUTS:
%       mus_outer (single array): average mus for parenchyma from XOR of
%                                   dilated vessels (outer)
%       ret_outer (single array): average ret for parenchyma from XOR of
%                                   dilated vessels (outer)
%       ori_std_outer (single array): std. dev. of orientation (outer)
%       mus_inner (single array): average mus for parenchyma from XOR of
%                                   original and dilated vessels (inner)
%       ret_inner (single array): average ret for parenchyma from XOR of
%                                   original and dilated vessels (inner)
%       ori_std_inner (single array): std. dev. of orientation (inner)
%   OUTLINE:
%       - Dilate vessels with se1
%           - Take XOR of of dilated and undilated vessels (inner)
%       - Dilate vessels with se2 (se2 > se1)
%           - xor of outer/inner = disjoint ring
%       - Measure scattering, retardance, orientation in xor output

% Keep 

% First case: XOR of dilated regions (seg_dil1 and seg_dil2)
seg_dil1 = imdilate(seg, se1);
seg_dil2 = imdilate(seg, se2);
parench_outer = xor(seg_dil1, seg_dil2);

% Second case: XOR of original and dilated regions (seg and seg_dil1)
parench_inner = xor(seg, seg_dil1);

% Measure properties for outer parenchyma ring
[mus_outer, ret_outer, ori_outer] =...
    measure_properties(parench_outer, mus, ret, ori);
% Measure properties for inner parenchyma ring
[mus_inner, ret_inner, ori_inner] =...
    measure_properties(parench_inner, mus, ret, ori);


% Helper function to measure mus, ret, and ori_std for a given parenchyma segmentation
function [mus, ret, ori_std] = measure_properties(parench, mus, ret, ori)
    % Create cell array of indices for each grouping of parenchyma
    cc = bwconncomp(parench, 26);
    idx = cc.PixelIdxList;

    % Measure mus & ret within parenchyma
    mus = rmmissing(cell2mat(cellfun(@(x) mean(x, 'omitnan'),...
                    cellfun(@(idx) mus(idx), idx, 'UniformOutput', false), ...
                    'UniformOutput', false)));
    
    ret = rmmissing(cell2mat(cellfun(@(x) mean(x, 'omitnan'),...
                    cellfun(@(idx) ret(idx), idx, 'UniformOutput', false), ...
                    'UniformOutput', false)));

    % Measure mean and std of orientation
    ori = deg2rad(ori);  % Convert orientation to radians and remove NaN
    ori_std = rmmissing(cell2mat(cellfun(@(x) circ_std(rmmissing(x)),...
        cellfun(@(idx) ori(idx), idx, 'UniformOutput', false), 'UniformOutput', false)));
end
end