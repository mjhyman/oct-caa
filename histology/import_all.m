function [stain] = import_all(stain_dir,stain_suffix,epvs_suffix,...
                                mask_suffix)
% INPUTS:
%   stain_dir (string): path to specific stain
%   stain_suffix (string): suffix of stain filename
%   epvs_suffix (string): suffix of EPVS annotations filename
%   mask_suffix (string): suffix of white matter mask (e.g. '_EPVS_Mask.tif')
%   
% OUTPUTS:
%   stain (struct): stain, epvs, and mask for each subject


% Set the folder containing your images
folder = stain_dir;
% List all EPVS mask files
fname = strcat('*',epvs_suffix);
mask_files = dir(fullfile(folder,fname));
% Initialize struct array
stain = struct();

for ii = 1:length(subdirs)
    fprintf('Importing subject %i out of %i\n',ii,length(subdirs))

    % Set the folder for subject/region
    folder = fullfile(hist_dir,subdirs{ii});
    % List all EPVS mask files
    mask_files = dir(fullfile(folder, '*_EPVS_Mask.tif'));

    % Extract base filename
    epvs_filename = mask_files(i).name;
    base_name = erase(epvs_filename,epvs_suffix);
    
    % Construct full file paths
    stain_path = fullfile(folder, [base_name stain_suffix]);
    epvs_path = fullfile(folder, epvs_filename);
    mask_path = fullfile(folder, [base_name mask_suffix]);
    
    % Load the images (assume both files exist)
    image = imread(stain_path);
    epvs = logical(imread(epvs_path));
    mask = logical(imread(mask_path));

    % Store in struct
    stain(i).baseName = base_name;
    stain(i).image = image;
    stain(i).epvs = epvs;
    stain(i).mask = mask;
end

end