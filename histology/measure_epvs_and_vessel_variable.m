function hist = measure_epvs_and_vessel_variable(hist,...
                    radii, rad, pix)
%HISTOLOGY_MEASURE_DONUT Measure controls + parenchyma around EPVS.
%   INPUTS:
%       - hist (struct): contains image, epvs mask, tissue border mask
%       - mrad (uint): measurement radius (units = pixels)
%       - dpath (str): output path
%       - radii (array): vector of radii of dilation (pixels)
%       - rad (double): radius of segmentation ring (just the ring) (pixels)
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

%% Iterate over the subjects in the hist struct

% Iterate subjects
for ii = 1:length(hist)
    fprintf('Subject %d of %d\n', ii, length(hist));
    %%% Retrieve local variables
    stain = im2single(hist(ii).image);
    epvs = logical(hist(ii).epvs);
    ves = logical(hist(ii).ves);
    mask = logical(hist(ii).mask);
    
    %%% Match histogram of current image to reference image
    % First iteration, use this as reference image
    if ii == 1
        stain_matched = ref_im;
    else
        % Apply histogram match to pixels within mask
        stain_pixels = stain(mask);
        ref_pixels = ref_im(ref_mask);
        pixels_matched = imhistmatch(stain_pixels,ref_pixels);
        % Convert pixels back to image
        stain_matched = stain;
        stain_matched(mask) = pixels_matched;
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

    %%% Iterate over segmentation radii
    for j = 1:length(radii)
        fprintf('\t--radius %d of %d\n', j, length(radii));
        %%% Create structuring elements for dilation
        % se1 is the inner dilation
        se1 = strel('disk',radii(j));
        % se2 = "th" voxels greater than > inner dilation
        se2 = strel('disk',radii(j)+rad);
    
        %%% Dilate masks for EPVS, Vessel & measure parenchyma
        [stain_epvs,stain_ves] = dilate_meas_exclude_overlap(z_stain,mask,...
                                        epvs,ves,se1,se2);
        
        %%% Add experimental + control to structure
        % Create string name for segmentation ring size
        rad_name = strcat('rad',num2str(floor(pix.*rad)));
        % Create subfield name for specific radius
        rad_str = strcat('rad',num2str(floor(pix.*radii(j))));
        % Add to struct
        hist(ii).(rad_name).(rad_str).exp = stain_epvs;
        hist(ii).(rad_name).(rad_str).ctl = stain_ves;
        hist(ii).(rad_name).(rad_str).exp_mean = mean(stain_epvs,'omitnan');
        hist(ii).(rad_name).(rad_str).ctl_mean = mean(stain_ves,'omitnan');        
    end
end
end