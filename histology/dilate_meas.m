function [pstain,pstain_mask] = dilate_meas(stain,mask,annot,se1,se2)
% Dilate the annotation, measure parenchyma
% INPUTS:
%   stain (2D image): stain image
%   mask (logical): binary mask of white matter
%   annot (logical): binary mask of annotations 
%   se1 (struct element): inner ring dilation structuring element
%   se2 (struct element): outter ring dilation structuring element
% OUTPUTS:
%   pstain (array): measurements within parenchyma

% Set the non-tissue pixels equal to NaN
stain(~mask) = NaN;
% Dilate epvs - inner ring
annot_inner = imdilate(annot,se1);
% Dilate epvs - outter ring
annot_outter = imdilate(annot,se2);
% Take XOR between dilated regions
pstain_mask = xor(annot_inner, annot_outter);
% Measure average within each discrete donut
cc = bwconncomp(pstain_mask);
pstain = cellfun(@(idx) mean(stain(idx),'omitnan'),...
                        cc.PixelIdxList);
end