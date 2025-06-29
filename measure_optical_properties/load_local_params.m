%% Function to load local parameters for each subject/region
function [mus, ret, seg, seg_wm, ori, epvs] = load_local_params(sub, loc)
% Function to load the local parameters for a given subject and location.
% It dynamically loads the data based on subject and region.
% INPUTS:
%   sub (struct): structure for a specific subject
%   loc (string): region under analysis
% OUTPUTS:
%   mus (array): scattering coefficient
%   ret (array): retardance coefficient
%   seg (logical array): vascular segmentation
%   seg_wm (logical array): vascular segmentation in white matter
%   ori (array): orientation 
%   epvs (logical array): EPVS segmentation
    
% Load the common local parameters
mus = sub.(loc).mus;
ret = sub.(loc).ret_full;
seg = sub.(loc).seg;
ori = sub.(loc).orient;

% Create white matter segmentation of segmentation
wm = sub.(loc).mask_wm;
seg_wm = wm .* seg;

% Set epvs to empty if the subject does not have epvs data
if isfield(sub, loc) && isfield(sub.(loc), 'epvs')
    epvs = sub.(loc).epvs;
else
    epvs = [];
end
end