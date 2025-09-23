function [stain] = import_stain(stain_dir,stain_suffix,epvs_suffix,...
                                ves_suffix, mask_suffix)
% INPUTS:
%   stain_dir (string): path to specific stain
%   stain_suffix (string): suffix of stain filename
%   epvs_suffix (string): suffix of EPVS annotations filename
%   ves_suffix (string): suffix of vessel annotations filename
%   mask_suffix (string): suffix of white matter mask (e.g. '_EPVS_Mask.tif')
%   
% OUTPUTS:
%   stain (struct): stain, epvs, and mask for each subject

% Set the suffix for identifying files to parse
base_suffix = ves_suffix;
% Set the folder containing your images
folder = stain_dir;
% List all EPVS mask files
fname = strcat('*',base_suffix);
mask_files = dir(fullfile(folder,fname));
% Initialize struct array
stain = struct();

for i = 1:length(mask_files)
    fprintf('Importing file %i out of %i\n',i,length(mask_files))
    % Extract base filename
    epvs_filename = mask_files(i).name;
    base_name = erase(epvs_filename,base_suffix);
    
    % Construct full file paths
    stain_path = fullfile(folder, [base_name, stain_suffix]);
    epvs_path = fullfile(folder, [base_name, epvs_suffix]);
    ves_path = fullfile(folder, [base_name, ves_suffix]);
    mask_path = fullfile(folder, [base_name, mask_suffix]);
    
    % Load the images (assume both files exist)
    image = imread(stain_path);
    epvs = logical(imread(epvs_path));
    ves = logical(imread(ves_path));
    mask = logical(imread(mask_path));

    %%% Invert the stain
    % Darker = more stain, but appears as lower number in Matlab
    % We want the opposite trend where higher number denotes more stain.
    image = 255 - image;

    % Store in struct
    stain(i).baseName = base_name;
    stain(i).image = image;
    stain(i).epvs = epvs;
    stain(i).ves = ves;
    stain(i).mask = mask;
end

end