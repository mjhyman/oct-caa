%% Measure the epvs density across each volume
subject_name = 'caa22';
region = 'front';

%% Settings
data_dir = '/projectnb/npbssmic/s/mhyman/CAA_data/matlab_structs/';
radius = 200;
fbase = '/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/';

%% Map subject names to file paths
subject_files = struct( ...
    'caa17', 'caa17.mat', ...
    'caa22', 'caa22.mat', ...
    'caa25', 'caa25.mat', ...
    'caa26', 'caa26.mat');

% Validate subject name
if ~isfield(subject_files, subject_name)
    error('Subject name "%s" not recognized.', subject_name);
end

% Check if heatmap already exists
fname = strcat(subject_name, '_', region, '_interpolated_heatmap.mat');
fname = fullfile(fbase, subject_name, region, fname);
if isfile(fname)
    fprintf('Heatmap already exists for %s %s. Skipping.\n', subject_name, region);
    return
end

% Load the subject data
file_path = fullfile(data_dir, subject_files.(subject_name));
fprintf('Loading subject %s\n', subject_name)
tmp = load(file_path);
subject_data = tmp.(subject_name);
fprintf('Finished loading %s\n', subject_name)

% Validate region
if ~isfield(subject_data, region)
    error('Region "%s" not found in subject "%s" data.', region, subject_name);
end

% Ensure EPVS field exists
if isfield(subject_data.(region), 'epvs')
    epvs = subject_data.(region).epvs;
    mask = subject_data.(region).mask_wm;
    fprintf('Starting %s %s \n', subject_name, region);
else
    fprintf('No EPVS data for %s %s. Skipping.\n', subject_name, region);
    return
end

% Calculate EPVS density
[subsampled_volume, interpolated_volume] = ...
    epvs_density(epvs, mask, radius);

% Save results
save_epvs_outputs(fbase, subject_name, region, ...
    subsampled_volume, interpolated_volume, radius);

% Save to Tif
fname = strcat(subject_name,'_',region,'_heatmap_','radius',radius,'.tif');
export_matrix_to_tiff(interpolated_volume, fbase, fname)

fprintf('Finished processing %s %s \n', subject_name, region);

%% Function to save .mat files and export TIFF
function save_epvs_outputs(fbase, sub, loc, ...
    subsampled_volume, interpolated_volume, radius)
    
% Define full output path
fpath = fullfile(fbase, sub, loc);
if ~exist(fpath, 'dir')
    mkdir(fpath)
end

% Save subsampled volume
fout = fullfile(fpath, strcat(sub, '_', loc, '_subsample_heatmap.mat'));
save(fout, 'subsampled_volume', '-v7.3');

% Save interpolated volume
fout = fullfile(fpath, strcat(sub, '_', loc, '_interpolated_heatmap.mat'));
save(fout, 'interpolated_volume', '-v7.3');

% Export TIFF
fname = strcat(sub, '_', loc, 'epvs_heatmap_radius', ...
    num2str(radius), '.tif');
try
    export_matrix_to_tiff(interpolated_volume, fpath, fname)
catch
    warning('Failed to export to TIF for %s %s', sub, loc);
end

end

%% Function to export matrix to TIFF
function export_matrix_to_tiff(A, fpath, filename)
    filename = fullfile(fpath, filename);

    if ndims(A) ~= 3
        error('Input must be a 3D matrix.');
    end

    t = Tiff(filename, 'w');

    tagstruct.ImageLength = size(A, 1);
    tagstruct.ImageWidth = size(A, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = size(A, 1);
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;

    for i = 1:size(A, 3)
        slice = uint16(65535 * mat2gray(A(:, :, i)));
        t.setTag(tagstruct);
        t.write(slice);
        if i < size(A, 3)
            t.writeDirectory();
        end
    end

    t.close();
end
