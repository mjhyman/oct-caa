%% Measure the epvs density across each volume
% Import each subjects .MAT struct
% Measure the EPVS density at each parenchymal voxel

%% Settings
% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/projectnb/npbssmic/s/mhyman/CAA_data/matlab_structs/';
% Search radius in voxels
radius = 200;
% Output directory
fbase = '/projectnb/npbssmic/s/mhyman/CAA_data/heatmaps/';

%% Start parallel pool
if isempty(gcp('nocreate'))
    parpool; % Start parallel pool with default settings
end

%% Define subject list and their file paths
subject_list = { ...
    struct('name', 'caa17', 'path', 'caa17.mat'), ...
    struct('name', 'caa22', 'path', 'caa22.mat'), ...
    struct('name', 'caa25', 'path', 'caa25.mat'), ...
    struct('name', 'caa26', 'path', 'caa26.mat') ...
};

%% Parallel Processing
parfor i = 1:length(subject_list)
    % Subject ID and file path
    sub = subject_list{i}.name;
    file_path = fullfile(data_dir, subject_list{i}.path);

    fprintf('Loading subject %s\n', sub)
    tmp = load(file_path);
    subject_data = tmp.(sub);
    fprintf('Finished loading %s\n', sub)

    % Iterate over regions
    regions = fieldnames(subject_data);
    for j = 1:length(regions)
        %%% Skip if heatmap already exists for subject/region
        % Retrieve tissue region
        loc = regions{j};
        % Create the filename
        fname = strcat(sub,'_',loc,'_interpolated_heatmap.mat');
        fname = fullfile(fbase,sub,loc,fname);
        % skip iteration if heatmap exists
        if isfile(fname)
            continue
        end
        
        %%% Ensure EPVS is within struct
        if isfield(subject_data.(loc), 'epvs')
            epvs = subject_data.(loc).epvs;
            mask = subject_data.(loc).mask_wm;
            fprintf('Starting %s %s \n', sub, loc);
        else
            continue
        end

        % Calculate EPVS density
        [subsampled_volume, interpolated_volume] = ...
            epvs_density(epvs, mask, radius);

        % Save outputs via function
        save_epvs_outputs(fbase, sub, loc, subsampled_volume,...
            interpolated_volume, radius);
    end
end

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