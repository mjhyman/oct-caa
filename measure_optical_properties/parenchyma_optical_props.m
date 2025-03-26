function [mus_mat,ret_mat,ori_std] = ...
    parenchyma_optical_props(seg,mus,ret,ori,se)
% Measure optical properties in parenchyma of all vessels
%   INPUTS:
%       seg (logical mat): vessels
%       mus (single mat): scattering coefficient matrix
%       ret (single mat): retardance matrix
%       ori (single mat): orientation matrix
%       se (handle): structuring element for dilation
%   OUTPUTS:
%       mus_mat (single array): average mus for each vessel's parenchyma
%       ret_mat (single array): average ret for each vessel's parenchyma
%       orient_std (single array): standard deviation of orientation for
%                                   each parenchyma around vessel or epvs
%   OUTLINE:
%       - dilate vessels
%       - Take XOR of of dilated and undilated vessels
%       - Measure retardance within parenchyma
%       - Measure scattering within parenchyma

%%% Segment the parenchymal tissue of all vessels
% Dilate vessels
seg_dil = imdilate(seg,se);
% Take XOR of of dilated and undilated vessels
parench = xor(seg, seg_dil);

%%% Create cell array of indices for each grouping of parenchyma
% Retrieve all disjoint parenchyma groupings of voxels
cc = bwconncomp(parench,26);
% Cell array of lists of pixel indices for each disjoint EPVS group
idx = cc.PixelIdxList;

%%% Measure mus & ret within parenchyma
% Measure scattering for each discrete parenchyma
mus_mat = rmmissing(cell2mat(cellfun(@(x) mean(x,'omitnan'),...
                    cellfun(@(idx) mus(idx),idx,'UniformOutput',false),...
                    'UniformOutput',false)));
% Measure retardance for each discrete parenchyma
ret_mat = rmmissing(cell2mat(cellfun(@(x) mean(x,'omitnan'),...
                    cellfun(@(idx) ret(idx),idx,'UniformOutput',false),...
                    'UniformOutput',false)));

%%% Measure mean and std of orientation
% Convert orientation to radians and remove NaN
ori = deg2rad(ori);
% The "rmmissing" function removes the NaN from each cell array. This
% accounts for the NaN that were set in the caller script, which account
% for vessels, backgrounds, etc.
% Calculate circular standard deviation
ori_std = rmmissing(cell2mat(cellfun(@(x) circ_std(rmmissing(x)),...
                    cellfun(@(idx) ori(idx),idx,'UniformOutput',false),...
                    'UniformOutput',false)));
end