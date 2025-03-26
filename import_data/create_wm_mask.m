%%  Create WM mask from tissue mask and thresholded retardance 
function [wm, seg] = create_wm_mask(mask,ret,seg,th,fout,save_flag)
% INPUTS:
%   mask (matrix): tissue mask
%   ret (matrix): retardance matrix (full volume)
%   mus (matrix): scattering coefficient
%   seg (matrix): segmentation (not masked)
%   th (matrix): minimum voxel intensity threshold of retardance for
%               segmenting the white matter
%   fout (str): directory path and filename for output
%   res (str): resolution of mri
%   dtype (str): data type of mri
%   pflag (str): permute flag for mri
%   save_flag (logical): flag for saving output
% OUTPUTS:
%   wm (matrix): white matter mask
%   seg (matrix): segmentation thresholded with white matter

%%% Create WM mask from thresholded retardance
% Threshold retardance
ret(ret < th) = 0;
% convert to logical
ret = logical(ret);
% white = retardance AND tissue mask
wm = logical(ret .* mask);

%%% Apply WM mask to vessels
seg = logical(wm .* seg);

%%% Save WM mask to NIFTI
if save_flag
    %%% MRI properties
    pflag = 0;              % permute x,y dimensions flag
    res = [0.02,0.02,0.02]; % resolution in millimeters
    dtype = 'uchar';        % smallest data type in freeview
    % Set path
    mask_path = fullfile(fout);
    % Save MRI
    save_mri(wm,mask_path,res,dtype,pflag);
end

end