function [pmus, pret] = meas_optical_props(seg,mus,ret,thick,zvox,se)
% MEAS_OPTICAL_PROPS measure parenchymal optical properties
% The segmentation matrix is a complete 3D volume. The scattering
% coefficient matrix (mus) has the same number of depths as the
% segmentation matrix. The retardance matrix (ret) has a projection for
% each physical slice.
%
% Steps:
%   - z-axis projection of segmentation
%   - dilate projection
%   - take the difference b/w segmentation and dilation
%   - measure retardace just in the difference
%
% INPUTS:
%   seg (logical matrix): segmentation either EPVS or vessels
%   mus (single matrix): scattering matrix AIP (1 per slice)
%   ret (single matrix): retardance matrix AIP (1 per slice)
%   thick (uint): thickness of slice (microns)
%   zvox (uint): voxel dimensions (microns)
%   se (matrix): structuring element for dilation
% OUTPUTS:
%   pmus (single matrix): scattering coefficients in surrounding
%       parenchymal tissue
%   pret (single matrix): retardance in surrounding
%       parenchymal tissue
% 

%%% Identify z-axis range of segmentation
% z-axis maximum intensity projection of segmentation
zseg = max(seg,[],3);
% Matrix subscripts of segmentation
[~,~,iz] = ind2sub(size(seg),find(seg));
% Identify highest & lowest z-axis position
z_top = min(iz(:));
z_bot = max(iz(:));
% Calculate slice thickness in voxels
nvox = thick ./ zvox;
% Maximum slice
slice_max = size(ret,3);
% Calculate bottom slice
slice_bot = ceil(z_bot ./ nvox);
slice_bot = min(slice_bot,slice_max);
% Calculate top slice (add 1 to account for 1 indexing)
slice_top = floor(z_top ./ nvox) + 1;

%%% Scattering Coefficient
% Subset of scattering corresponding to segmentation
pmus = mus(:,:,slice_top:slice_bot);
% z-axis average intensity projection of scattering
zmus = mean(pmus,3);
% Dilate projection of segmentation
zseg_dil = imdilate(zseg,se);
% Difference b/w segmentation and dilation
tiss = xor(zseg, zseg_dil);
% Average scattering in parenchymal tissue
pmus = mean(zmus(tiss),"all");

%%% Retardance
% Recalculate z indices for retardance
ret_top = max(1,floor(slice_top / nvox));
ret_bot = min(size(seg,3),ceil(slice_bot / nvox));
% Subset of retardance corresponding to segmentation
pret = ret(:,:,ret_top:ret_bot);
% z-axis average intensity projection of retardance
zret = mean(pret,3);
% Difference b/w segmentation and dilation
tiss = xor(zseg, zseg_dil);
% Measure retardance in parenchymal tissue
pret = mean(zret(tiss),"all");
end