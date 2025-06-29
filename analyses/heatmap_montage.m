function heatmap_montage(heatmap,segmentation,mask,fpath,fname,save_tif)
% INPUTS
% heatmap (matrix): heatmap matrix of metric
% segmentation (logical matrix): EPVS or vessel matrix
% mask (logical matrix): tissue mask (e.g. white matter, gray matter, etc.)
% fpath (string): output directory
% fname (string): output filename
% save_tif (logical): boolean for saving output figure

% Parameters
z_range = 1:size(heatmap, 3);  % Loop over all slices
tif_fname = fullfile(fpath,fname);

figure('Color', 'k'); colormap hot;

for z = z_range
    % Extract slices
    density_slice = heatmap(:,:,z);
    inclusion_slice = mask(:,:,z);
    epvs_slice = segmentation(:,:,z);

    % Normalize density for consistent color scaling
    norm_density = density_slice ./ max(heatmap(:));

    % Plot
    clf;
    imagesc(inclusion_slice); colormap gray; axis image off;
    hold on;
    h = imagesc(norm_density);
    colormap(gca, hot);
    set(h, 'AlphaData', norm_density > 0);

    % EPVS boundaries
    [B,~] = bwboundaries(epvs_slice);
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1);
    end

    title(sprintf('EPVS Density Overlay - Slice z = %d', z), ...
        'Color', 'w', 'FontSize', 14);
    colorbar;

    drawnow;

    % Save to multi-page TIF
    if save_tif
        frame = getframe(gcf);
        im = frame2im(frame);
        im = im2uint8(im);  % Ensure image is uint8 format

        if z == z_range(1)
            imwrite(im, tif_fname, 'tif', 'Compression', 'none', ...
                    'WriteMode', 'overwrite');
        else
            imwrite(im, tif_fname, 'tif', 'Compression', 'none', ...
                    'WriteMode', 'append');
        end
    end
end
end