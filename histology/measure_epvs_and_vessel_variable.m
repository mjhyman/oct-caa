function hist = measure_epvs_and_vessel_variable(hist,...
                    radii, th, histmatch_flag, plot_flag)
%HISTOLOGY_MEASURE_DONUT Measure controls + parenchyma around EPVS.
%   INPUTS:
%       - hist (struct): contains image, epvs mask, tissue border mask
%       - mrad (uint): measurement radius (units = pixels)
%       - dpath (str): output path
%       - radii (array): radii of dilation (pixels)
%       - th (double): thickness of segmentation ring (pixels)
%       - histmatch_flag (bool): histogram matching
%       - plot_flag (bool): create figures for debugging
%   OUTPUTS:
%       - hist (struct):
%           - exp (matrix): experimental measurements (around EPVS)
%           - ctl (matrix): control measurements (uniform grid)


%% Histogram Matching (global normalization)
% Select first image as reference image for histogram matching. This will
% be used as the reference image. All other images will be normalized to
% this one.
ref_im = im2single(hist(1).image);
ref_mask = logical(hist(1).mask);
zmin = 0;
zmax = 1;

%% Iterate over the subjects in the hist struct

% Iterate subjects
for ii = 1:length(hist)
    %%% Retrieve local variables
    stain = im2single(hist(ii).image);
    epvs = logical(hist(ii).epvs);
    ves = logical(hist(ii).ves);
    mask = logical(hist(ii).mask);
    
    %%% Match histogram of current image to reference image
    if ~histmatch_flag
        stain_matched = stain;
    elseif ii ~= 1
        % Apply histogram match to pixels within mask
        stain_pixels = stain(mask);
        ref_pixels = ref_im(ref_mask);
        pixels_matched = imhistmatch(stain_pixels,ref_pixels);
        % Convert pixels back to image
        stain_matched = stain;
        stain_matched(mask) = pixels_matched;
    else
        stain_matched = ref_im;
    end
    
    %%% Perform z-score normalization
    masked_pixels = stain_matched(mask);
    mu = mean(masked_pixels(:), 'omitnan');
    sigma = std(masked_pixels(:), 'omitnan');
    z_stain = stain_matched;
    z_stain(mask) = (stain_matched(mask) - mu) / sigma;

    %%% Add the debugging figures to the struct
    hist(ii).stain = stain;
    hist(ii).stain_matched = stain_matched;
    hist(ii).z_stain = z_stain;

    %%% Find the z_stain limits
    if plot_flag
        zmin = min([zmin, min(z_stain(mask))]);
        zmax = max([zmax, max(z_stain(mask))]);
    end

    %%% Iterate over segmentation radii
    for j = 1:length(radii)
        fprintf('  --starting radius %d of %d\n', j, length(radii));
        %%% Create structuring elements for dilation
        % se1 is the inner dilation
        se1 = strel('disk',radii(j));
        % se2 = "th" voxels greater than > inner dilation
        se2 = strel('disk',radii(j)+th);
    
        %%% Dilate masks for EPVS, Vessel & measure parenchyma
        [stain_epvs,stain_ves] = dilate_meas_exclude_overlap(z_stain,mask,...
                                        epvs,ves,se1,se2);
        
        %%% Add experimental + control to structure
        hist(ii).rad(j).exp = stain_epvs;
        hist(ii).rad(j).ctl = stain_ves;
        hist(ii).rad(j).exp_mean = mean(stain_epvs,'omitnan');
        hist(ii).rad(j).ctl_mean = mean(stain_ves,'omitnan');        
    end
end
end