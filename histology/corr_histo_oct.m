function hist = corr_histo_oct(hist,radii,r)
% Measure correlation b/w histology and OCT optical property
% Create donut around EPVS or vessel, measure parenchyma and respective OCT
% value, measure correlation (Spearman's rho).
%
%   INPUTS:
%       - hist (struct): contains image, epvs mask, tissue border mask
%       - radii (array): inner radii of dilation
%       - r (uint): constant to dilate outer radius
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
    % iterate dilation inner radius
    for j=1:length(radii)
        %%% Create structuring elements for dilation
        se1 = strel('disk',radii(j));
        se2 = strel('disk',radii(j)+r);

        %%% Retrieve local variables
        stain = im2single(hist(ii).image);
        epvs = logical(hist(ii).epvs);
        ves = logical(hist(ii).ves);
        mask = logical(hist(ii).mask);
        ret = im2single(hist(ii).ret);        

        %%% Match histogram of current image to reference image
        if ii ~= 1
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
        
        %%% Dilate EPVS & vessel mask & measure parenchyma
        [stain_epvs,~] = dilate_meas(z_stain,mask,epvs,se1,se2);
        [stain_ves,~] = dilate_meas(z_stain,mask,ves,se1,se2);

        %%% Measure OCT retardance
        % Measure around EPVS
        ret_epvs = dilate_meas(ret,mask,epvs,se1,se2);
        % Measure around vessels
        ret_ves = dilate_meas(ret,mask,ves,se1,se2);
        
        %%% Add experimental + control to structure
        hist(ii).meas(j).histo_epvs = stain_epvs;
        hist(ii).meas(j).histo_ves = stain_ves;
        hist(ii).meas(j).histo_epvs_mean = mean(stain_epvs,'omitnan');
        hist(ii).meas(j).histo_ves_mean = mean(stain_ves,'omitnan');
        hist(ii).meas(j).ret_epvs = ret_epvs;
        hist(ii).meas(j).ret_ves = ret_ves;
        hist(ii).meas(j).ret_epvs_mean = mean(ret_epvs,'omitnan');
        hist(ii).meas(j).ret_ves_mean = mean(ret_ves,'omitnan');
    end
    %%% Measure correlation b/w EPVS & ves vs. OCT
    % Extract average pathology measurements at each radii
    histo_epvs = [hist(ii).meas.histo_epvs_mean];
    histo_ves = [hist(ii).meas.histo_ves_mean];
    ret_epvs = [hist(ii).meas.ret_epvs_mean];
    ret_ves = [hist(ii).meas.ret_ves_mean];
    % Add to struct
    hist(ii).histo_epvs = histo_epvs;
    hist(ii).histo_ves = histo_ves;
    hist(ii).ret_epvs = ret_epvs;
    hist(ii).ret_ves = ret_ves;
    
    % Compare EPVS
    [epvs_r,epvs_p] = corr(histo_epvs',ret_epvs','Type','Spearman',...
                        'Rows','pairwise');
    % Compare vessels
    [ves_r,ves_p] = corr(histo_ves',ret_ves','Type','Spearman',...
                        'Rows','pairwise');
    % Add to struct
    hist(ii).epvs_r = epvs_r;
    hist(ii).epvs_p = epvs_p;
    hist(ii).ves_r = ves_r;
    hist(ii).ves_p = ves_p;
end


end