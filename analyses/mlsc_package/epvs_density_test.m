%% Test the function "epvs_density"

% generate a file to check if the script is actually running
fid = fopen(['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
            'tmp/run_success.txt'],'w');
fprintf(fid,'hello whatever, i dont care');
fclose(fid);

%% Settings
% Directory for loading seg, mus, ret, mask, epvs
data_dir = ['/autofs/cluster/octdata3/users/mjhyman/' ...
    'oct_caa_analyses/optical_properties'];
% Search radius in voxels
radius = 200;
% Output directory
fbase = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
    'optical_properties/heatmaps'];
% flag for importing data structs
flag_load_caa_structs = true;

%% Import Data

%%% Load the .MAT structs
if flag_load_caa_structs
    % CAA 17
    fprintf('Loading CAA17\n')
    caa17 = load(fullfile(data_dir,'/caa17/occip/caa17.mat'));
    caa17 = caa17.caa17;
    fprintf('Finished Loading CAA17\n')
    % CAA 22
    fprintf('Loading CAA22\n')
    caa22 = load(fullfile(data_dir,'/caa22/caa22.mat'));
    caa22 = caa22.caa22;
    fprintf('Finished Loading CAA22\n')
    % CAA 25
    fprintf('Loading CAA25\n')
    caa25 = load(fullfile(data_dir,'/caa25/caa25.mat'));
    caa25 = caa25.caa25;
    fprintf('Finished Loading CAA25\n')
    % CAA 26
    fprintf('Loading CAA26\n')
    caa26 = load(fullfile(data_dir,'/caa26/caa26.mat'));
    caa26 = caa26.caa26;
    fprintf('Finished Loading CAA26\n')
end

subjects = struct();
subjects.caa17 = caa17;
subjects.caa22 = caa22;
subjects.caa25 = caa25;
subjects.caa26 = caa26;

%% Iterate over subjects

% retrieve subject IDs
subs = fieldnames(subjects);

% Loop over each subject
for i = 1:length(subs)
    % subject ID
    sub = subs{i};  
    fprintf('STARTING sub %s\n',sub)
    % regions for this subject
    regions = fieldnames(subjects.(sub));
    % Loop over each region (e.g., 'front', 'occip')
    for j = 1:length(regions)
        % Region (front or occip)
        loc = regions{j};  
        fprintf('Starting %s %s \n', sub, loc);
        % Load EPVS + WM mask
        if isfield(subjects.(sub).(loc),'epvs')
            epvs = subjects.(sub).(loc).epvs;
            mask = subjects.(sub).(loc).mask_wm;        
        % skip for-loop iteration if EPVS is not present
        else
            continue
        end
        
        %%% Calculate EPVS density
        [subsampled_volume, interpolated_volume] =...
            epvs_density(epvs,mask,radius);
        
        %%% Save the outputs
        fprintf('Saving outputs\n')
        fpath = fullfile(fbase,sub,loc);
        % subsampled output filename
        fout = fullfile(fpath,strcat(sub,'_',loc,'_subsample_heatmap.mat'));
        save(fout,'subsampled_volume','-v7.3');
        % interpolated output filename
        fout = fullfile(fpath,strcat(sub,'_',loc,'_interpolated_heatmap.mat'));
        save(fout,'interpolated_volume','-v7.3');

        %%% Generate .TIF stack of interpolated volume
        fname = strcat(sub,'_',loc,'epvs_heatmap_radius',num2str(radius),'.tif');
        try
            export_matrix_to_tiff(interpolated_volume,fpath,fname)
        catch
            warning('Failed to export to Tif for %s %s',sub,loc);
        end
    end
end


%% Function to export matrix to TIFF
function export_matrix_to_tiff(A, fpath, filename)
% export 3D double matrix to TIF file
% INPUTS: 
%   A (double matrix): the heatmap file
%   fpath (double matrix): output file path
%   filename (double matrix): output filename
% OUTPUTS: saves TIF 

% Create the full output file
filename = fullfile(fpath,filename);

% Ensure input is a 3D matrix
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
tagstruct.SampleFormat = Tiff.SampleFormat.UInt; % Use UInt for uint16

for i = 1:size(A, 3)
    slice = uint16(65535 * mat2gray(A(:, :, i)));  % scale and convert
    t.setTag(tagstruct);
    t.write(slice);
    if i < size(A, 3)
        t.writeDirectory();
    end
end

t.close();
end
