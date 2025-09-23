%% Update CAA6 occipital
% The z-axis appears flipped on each physical slice (5 depths each). The
% purpose of this script is to flip each depth.

% Directory to optical properties
dpath = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
            'optical_properties'];

%% Load all structures

% CAA 6 struct
fprintf('Loading CAA6\n')
caa6 = load(fullfile(dpath,'/caa6/caa6.mat'));
caa6 = caa6.caa6;
fprintf('Finished Loading CAA6\n')

%% Load updated white matter mask

%%% CAA 6 Occip
fprintf('Importing caa6 occip WM mask\n')
% import wm mask
mask_wm = fullfile(dpath,'caa6/occip/wm_mask_revised.nii');
mask_wm = MRIread(mask_wm,0,0);
% Keep just WM (wm = 1)
mask_wm = mask_wm.vol;
mask_wm = logical(mask_wm==1);
% Flip about the horizontal to align with scattering
mask_wm = flip(mask_wm,1);
% Add back to struct
caa6.occip.mask_wm = mask_wm;
fprintf('Finished updating caa6 occip WM mask\n')

%% Flip the Z-axis for each depth
% apply to white matter mask, mus, retardance

%%% mus
% Correct mus stack
mus = caa6.occip.mus;
mus = reverse_order(mus);
% Add back to struct
caa6.occip.mus = mus;
% Save output as nifti
fout = fullfile(dpath,'caa6/occip/mus_revised.nii');
save_mri(mus,fout,[0.02,0.02,0.02],'float',0);

%% Segmentation (yes)
% reverse order
seg = caa6.occip.seg;
seg = reverse_order(seg);
% Add to struct
caa6.occip.seg = seg;
% Save output as nifti
fout = fullfile(dpath,'caa6/occip/seg_revised.nii');
save_mri(seg,fout,[0.02,0.02,0.02],'uchar',0);

%%% Apply the white matter mask to segmentation
caa6.occip.seg_wm = seg .* mask_wm;

%% Save .MAT
% CAA 6
fprintf('Saving CAA6/n')
dout = fullfile(dpath,'/caa6/caa6.mat');
save(dout,'caa6','-v7.3');
fprintf('Finished updating CAA6')

%% Function to reverse the depths:

function vol = reverse_order(vol)
% Reverse the order of each physical slice
% INPUTS:
%   vol (3D): volume to be reversed
% OUTPUTS:
%   vol (3D): volume after reversal


%%% Depth 1-6 are continuous. The following depths are shifted every 5
p1 = vol(:,:,1:6);
p1 = p1(:,:,end:-1:1);
vol(:,:,1:6) = p1;
for ii = 7:5:500
    % take subset of vol stack
    try
        % Retrieve portion of stack
        stack = vol(:,:,ii:ii+4);
        % flip the stack
        stack = stack(:,:,end:-1:1);
        % Place back into original stack
        vol(:,:,ii:ii+4) = stack;
    catch
        % Retrieve portion of stack
        stack = vol(:,:,ii:end);
        % flip the stack
        stack = stack(:,:,end:-1:1);
        % Place back into original stack
        vol(:,:,ii:end) = stack;
    end
end
end