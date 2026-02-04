function [hist] = histo_match_zscore(hist)
% HISTO_MATCH_ZSCORE: histogram match and z-score
%   INPUTS:
%       - hist (struct): Stain: image, epvs mask, tissue border mask
%   OUTPUTS:
%       - hist (struct): contains image, epvs mask, tissue border mask,
%                       stain_matched, z_stain
  
%% Histogram Matching (global normalization)
% Select first image as reference image for histogram matching. This will
% be used as the reference image. All other images will be normalized to
% this one.
ref_im = im2single(hist(1).image);
ref_mask = logical(hist(1).mask);

%% Iterate over each section within stain 

for ii = 1:length(hist)
    fprintf('Normalizing subject %d of %d\n', ii, length(hist));
    %%% Retrieve local variables
    stain = im2single(hist(ii).image);
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

    %%% Add back to struct
    hist(ii).stain_matched = stain_matched;
    hist(ii).z_stain = z_stain;
end
end