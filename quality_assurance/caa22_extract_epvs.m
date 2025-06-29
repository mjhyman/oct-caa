%% Measure white matter parenchymal Optical Properties
%{
Purpose is to extract the EPVS from the segmentation, which contains some
of the EPVS.
%}

%% Prepare environment
clear; clc; close all;
% Add top-level directory + subdirectories
addpath(genpath(fullfile(pwd, '..')))
% Directory for loading seg, mus, ret, mask, epvs
data_dir = ['/autofs/cluster/octdata3/users/mjhyman/' ...
    'oct_caa_analyses/optical_properties'];
% WM mask directory (from Taylor)
wm_dir = ['/autofs/space/turtle_001/users/xz875/projects/' ...
          'Multi-resolution_Unets_Semi_OCT/prediction'];
% Voxel dimensions (microns) for all runs
res = [20,20,20]; % resolution in microns

%%% Flag for loading .MAT struct and creating WM masks
% flag for reloading the .MAT struct for each subject
flag_load_caa_structs = true;
% flag for importing updated wm masks
flag_import_wm_mask = true;

%% Load CAA22 .MAT struct
fprintf('Loading CAA22\n')
caa22 = load(fullfile(data_dir,'/caa22/caa22.mat'));
caa22 = caa22.caa22;
fprintf('Finished Loading CAA22\n')

%% Load the manually edited CAA22 front EPVS
fin = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
    'optical_properties/caa22/front/epvs_manual_latest.nii'];
epvs = MRIread(fin,0,0);
epvs = logical(epvs.vol);

%% Extract EPVS from segmentation
% Extract segmentation
seg = caa22.front.seg;
% Apply white matter mask
wm = caa22.front.mask_wm;
seg = seg .* wm;
% minimum number of voxels to keep
nmin = 200;
% maximum # voxels to keep
nmax = 25000;
% radius to dilate the epvs
r = 1;
fout = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
    'optical_properties/caa22/front/epvs_subsets/' ...
    'epvs_subset_200_25000_r1.nii'];

%%% Dilate - connect disjoint EPVS that are missing few voxels
% dilate segmentation to remove gaps
se = strel('sphere',r);
seg = imdilate(seg,se);

%%% Find connected components
% Define connectivity (6, 18, or 26 for 3D)
connectivity = 26; % Most common choice for volumetric data
% Find connected components
CC = bwconncomp(seg, connectivity);
% Measure number of voxels in each component
voxelCounts = cellfun(@numel, CC.PixelIdxList);
bin_edge = 0:1000:1e5;
figure; histogram(voxelCounts,bin_edge);

%%% Keep only the large components
largeComponentsIdx = find(voxelCounts >= nmin & voxelCounts <= nmax);
seg_filtered = ismember(labelmatrix(CC), largeComponentsIdx);
% recombine with manual annotations
epvs_combined = seg_filtered | epvs;

%%% Save the filtered to nifti
save_mri(epvs_combined,fout,[0.02,0.02,0.02],'uchar',0);

%% Filter EPVS+Seg
% Filter
seg_filtered = filter_pca(seg,10,26);
% Save output
fout = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
    'optical_properties/caa22/front/epvs_subsets/' ...
    'seg_pca_r10.nii'];
save_mri(seg_filtered,fout,[0.02,0.02,0.02],'uchar',0);

%% Extract EPVS based on radii
%{
% Define minimum width threshold (in voxels)
min_rad = 25; % e.g., 10 voxels
% Define output path
fout = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
    'optical_properties/caa22/front/epvs_subsets/' ...
    'epvs_min_rad_25.nii'];

% Define connectivity (6, 18, or 26 for 3D)
connectivity = 26;
% Create new blank matrix
seg_filtered = filter_minor_radius(seg,min_rad,connectivity);
% Visualize
figure; imagesc(seg_filtered(:,:,100));

%%% Recombine with manual EPVS and save
epvs_combined = seg_filtered | epvs;
save_mri(epvs_combined,fout,[0.02,0.02,0.02],'uchar',0);
%}

%% Function to filter vessels based on PCA

function seg_filtered = filter_pca(seg, r, connectivity)
% REMOVETHINCYLINDERS Removes thin cylindrical components with radius below r
%
%   BW_filtered = removeThinCylinders(BW, r, connectivity)
%
%   Inputs:
%       BW - 3D logical matrix
%       r - Radius threshold (components with radius < r will be removed)
%       connectivity - (Optional) Connectivity for connected components (default = 26)
%
%   Output:
%       BW_filtered - Filtered 3D logical matrix

    if nargin < 3
        connectivity = 26; % Default to full 3D connectivity
    end

    % Find connected components
    CC = bwconncomp(seg, connectivity);

    % Prepare label matrix
    L = labelmatrix(CC);

    % Initialize cylinder radius array
    cylinderRadii = zeros(CC.NumObjects, 1);

    % Loop through each component
    for i = 1:CC.NumObjects
        voxelIdx = CC.PixelIdxList{i};
        if numel(voxelIdx) < 3
            cylinderRadii(i) = 0; % Skip tiny components
            continue
        end

        % Convert indices to coordinates
        [x, y, z] = ind2sub(size(seg), voxelIdx);
        coords = [x, y, z];

        % Center the coordinates
        coords_centered = coords - mean(coords, 1);

        % Perform PCA to find main axis
        [coeff, ~, ~] = pca(coords_centered);

        % The first principal component is the cylinder's main axis
        mainAxis = coeff(:,1);

        % Project points onto the main axis
        projections = coords_centered * mainAxis;

        % Compute residuals (distance from main axis)
        projected_points = projections * mainAxis';
        residuals = coords_centered - projected_points;
        distances = sqrt(sum(residuals.^2, 2));

        % Estimate cylinder radius (mean distance to axis)
        cylinderRadii(i) = mean(distances);
    end

    % Retain components with radius >= r
    validComponentsIdx = find(cylinderRadii >= r);

    % Create filtered matrix
    seg_filtered = ismember(L, validComponentsIdx);
end



%% Functions to filter the segmentation by radii
function BW_filtered = filter_minor_radius(BW, nmin, connectivity)
% FILTERCOMPONENTSBYMINORRADIUS Retain connected components with minor radius >= nmin
%
%   BW_filtered = filterComponentsByMinorRadius(BW, nmin, connectivity)
%
%   Inputs:
%       BW - 3D logical matrix
%       nmin - Minimum minor axis radius (in voxels) to retain a component
%       connectivity - (Optional) Connectivity for connected components (default = 26)
%
%   Output:
%       BW_filtered - 3D logical matrix with filtered components

    if nargin < 3
        connectivity = 26; % Default to full 3D connectivity
    end

    % Find connected components
    CC = bwconncomp(BW, connectivity);

    % Prepare label matrix
    L = labelmatrix(CC);

    % Calculate minor radius for each component using arrayfun
    minorRadii = arrayfun(@(i) computeMinorRadius(CC.PixelIdxList{i}, size(BW)), 1:CC.NumObjects);

    % Find components that meet the minor radius threshold
    validComponentsIdx = find(minorRadii >= nmin);

    % Create filtered matrix
    BW_filtered = ismember(L, validComponentsIdx);

end

function minorRadius = computeMinorRadius(voxelIdx, matrixSize)
    if numel(voxelIdx) < 3
        minorRadius = 0;
        return
    end

    % Convert linear indices to subscripts
    [x, y, z] = ind2sub(matrixSize, voxelIdx);
    coords = [x y z];

    % Calculate covariance matrix
    covMatrix = cov(coords);

    % Calculate eigenvalues (principal axis variances)
    eigenvalues = eig(covMatrix);

    % The radii of the ellipsoid are sqrt of eigenvalues
    radii = sqrt(eigenvalues);

    % Minor radius is the smallest radius
    minorRadius = min(radii);
end
