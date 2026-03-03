%% Function to save .mat files and export TIFF
function save_epvs_heatmap(fbase, sub, loc, prefix, ...
                subsampled_volume, interpolated_volume, radius, p)
    
% Define full output path
fpath = fullfile(fbase, sub, loc);
if ~exist(fpath, 'dir')
    mkdir(fpath)
end

% Save subsampled volume
fout = fullfile(fpath, strcat(sub,'_',loc,'_',prefix,...
                '_radius_',num2str(radius),...
                '_exp_',num2str(p), ...
                '_subsample_heatmap.mat'));
save(fout, 'subsampled_volume', '-v7.3');

% Save interpolated volume
fout = fullfile(fpath, strcat(sub,'_',loc,'_',prefix,...
                '_radius_',num2str(radius),...
                '_exp_',num2str(p), ...
                '_interpolated_heatmap.mat'));
save(fout, 'interpolated_volume', '-v7.3');

% Export TIFF
fname = strcat(sub,'_',loc,'_',prefix,...
                '_radius_',num2str(radius),...
                '_exp_',num2str(p), ...
                '_interpolated_heatmap.tif');
try
    export_matrix_to_tiff(interpolated_volume, fpath, fname)
catch
    warning('Failed to export to TIF for %s %s', sub, loc);
end

end