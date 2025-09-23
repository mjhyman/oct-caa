function [epvs, ves] = dilate_meas_exclude_overlap( ...
                                stain,mask,epvs,ves,se1,se2)
% Dilate the annotation, measure parenchyma
% INPUTS:
%   stain (2D image): stain image
%   mask (logical): binary mask of white matter
%   epvs (logical): binary mask of epvs annotations
%   ves (logical): binary mask of vessel annotations
%   se1 (struct element): inner ring dilation structuring element
%   se2 (struct element): outter ring dilation structuring element
% OUTPUTS:
%   pstain (array): average stain value within parenchyma
%   pstain_mask (array): logical mask of measurements

%%% Set the non-tissue pixels equal to NaN
stain(~mask) = NaN;

%%% Dilate the EPVS + vessels
epvs_donut = dilate_xor(epvs, se1, se2);
ves_donut = dilate_xor(ves, se1, se2);

%%% Exclude overlap between donuts from both
overlap = and(epvs_donut, ves_donut);
epvs_donut = epvs_donut - overlap;
ves_donut = ves_donut - overlap;

%%% Measure average within each discrete donut
% EPVS
cc_epvs = bwconncomp(epvs_donut);
epvs = cellfun(@(idx) mean(stain(idx),'omitnan'),cc_epvs.PixelIdxList);
% Vessel
cc_ves = bwconncomp(ves_donut);
ves = cellfun(@(idx) mean(stain(idx),'omitnan'),cc_ves.PixelIdxList);

%% Function for performing dilation and exclusive or
    function [donut] = dilate_xor(annot, se1, se2)
        % Dilate epvs - inner ring
        annot_inner = imdilate(annot,se1);
        % Dilate epvs - outter ring
        annot_outter = imdilate(annot,se2);
        % Take XOR between dilated regions
        donut = xor(annot_inner, annot_outter);
    end

end