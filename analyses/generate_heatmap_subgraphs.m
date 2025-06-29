%% Divide volume into grid and calculate metrics for each cube
% The purpose of this is to subdivide the entire tissue volume in cubes of
% dimensions [M, M, M] cubic millimeters and then generate the graph for
% each cube. The following step is to generate the heat map from these
% cubes.

%% Add top-level directory of code repository to path
clc; close all;
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
ncores = feature('numcores');
maxNumCompThreads(ncores);
% Flag for loading CAA structs (false if already in environment)
flag_load_caa_structs = false;
% Directory for loading seg, mus, ret, mask, epvs
data_dir = ['/autofs/cluster/octdata3/users/mjhyman/' ...
    'oct_caa_analyses/optical_properties'];

%%% Subvolume parameters
% Isotropic cube length (microns)
cube_side = 200;
% Size of each voxel (microns)
vox = [20, 20, 20];
% Compute number of voxels in x,y,z for each cube
n_x = floor(cube_side ./ vox(1));
n_y = floor(cube_side ./ vox(2));
n_z = floor(cube_side ./ vox(3));

%% Load the subject structs

% Load each subject
if flag_load_caa_structs
    % CAA 6
    fprintf('Loading CAA6\n')
    caa6 = load(fullfile(data_dir,'/caa6/caa6.mat'));
    caa6 = caa6.caa6;
    fprintf('Finished Loading CAA6\n')
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

%% Initialize directories, filenames, parameters

%%% Add all subjects to struct for ease
volumes = struct();
volumes.caa6 = caa6;
volumes.caa17 = caa17;
volumes.caa22 = caa22;
volumes.caa25 = caa25;
volumes.caa26 = caa26;
% retrieve subject IDs
subid = fieldnames(volumes);

%%% Struct for storing vascular heat map
heatmap = struct();

%% Create matrix of EPVS & Vascular heatmaps (volume fraction)

%%% Iterate over each subject
for ii = 1:length(subid)
    % retrieve local subject ID
    sub = subid{ii};
    % retrieve regions for this subject
    regions = fieldnames(volumes.(sub));
    % Iterate over tissue volume region(s)
    for j = 1:length(regions)
        % Retrieve local region
        reg = regions{j};
        fprintf('\nSTARTED HEATMAP -- %s %s\n',sub,reg)
        % Retrieve segmentation
        seg = logical(volumes.(sub).(reg).seg);
        % Retrieve tissue mask
        mask = logical(volumes.(sub).(reg).mask);
        % Retrieve optical properties (mus, ret, ori)
        mus = volumes.(sub).(reg).mus;
        ret = volumes.(sub).(reg).ret_full;
        ori = deg2rad(volumes.(sub).(reg).orient);
        % Set vessel voxels to NaN within optical properties
        mus(seg) = nan;
        ret(seg) = nan;
        ori(seg) = nan;
        % Calculate number of cubes in z dimension
        Nz = ceil(size(seg,3) ./ n_z);
        % Initialize vessel volume fraction matrix
        ves_vf_mat = zeros(size(seg,1), size(seg,2), Nz);
        % If EPVS exists, initialize heatmap matrix
        if isfield(volumes.(sub).(reg),'epvs')
            epvs = logical(volumes.(sub).(reg).epvs);
            epvs_vf_mat = zeros(size(seg,1), size(seg,2), Nz);
            % Set EPVS voxels to NaN within optical properties
            mus(epvs) = nan;
            ret(epvs) = nan;
            ori(epvs) = nan;
        end
        % Initialize optical properties matrices
        mus_mat = zeros(size(seg,1), size(seg,2), Nz);
        ret_mat = zeros(size(seg,1), size(seg,2), Nz);
        ori_mat = zeros(size(seg,1), size(seg,2), Nz);
        
        %% Iterate over the segmentation
        % Heatmap depth index - depth in matrix
        hm_z_idx = 0;
        % Iterate over the z-axis
        for z = 1:n_z:size(seg,3)
            % Iterate the heatmap depth index
            hm_z_idx = hm_z_idx + 1;
            % Iterate over rows
            for x = 1:n_x:size(seg,1)
                % Iterate over columns
                for y = 1:n_y:size(seg,2)
                    %% Crop matrices into cubes
                    % Initialize end indices for each axis
                    xf = x + n_x - 1;
                    yf = y + n_y - 1;
                    zf = z + n_z - 1;
                    % Take minimum of matrix dimensions and end indices
                    xf = min(xf, size(seg,1));
                    yf = min(yf, size(seg,2));
                    zf = min(zf, size(seg,3));
                    % Take cube from segmentation + mask
                    seg_cube = seg((x:xf), (y:yf), (z:zf));
                    mask_cube = mask((x:xf), (y:yf), (z:zf));
                    % Take cube from EPVS
                    if isfield(volumes.(sub).(reg),'epvs')
                        epvs_cube = epvs((x:xf), (y:yf), (z:zf));
                    end
                    
                    %% Mean of optical properties
                    % Take means of mus and retardance
                    mus_mean = mean(mus((x:xf),(y:yf),(z:zf)),'all','omitnan');
                    ret_mean = mean(ret((x:xf),(y:yf),(z:zf)),'all','omitnan');
                    % Calculate circular std dev. of orientation
                    ori_cube = ori((x:xf), (y:yf), (z:zf));
                    ori_std = circ_std(rmmissing(ori_cube(:)));
                    % Add to matrices
                    mus_mat((x:xf),(y:yf),hm_z_idx) = mus_mean;
                    ret_mat((x:xf),(y:yf),hm_z_idx) = ret_mean;
                    ori_mat((x:xf),(y:yf),hm_z_idx) = ori_std;
    
                    %% Calculate volume fraction if segmentation || EPVS
                    % Vasculature
                    if sum(seg_cube(:)) > 1
                        vf = sum(seg_cube(:)) ./ sum(mask_cube(:));
                        ves_vf_mat((x:xf), (y:yf), hm_z_idx) = vf;
                    end
                    % If EPVS is a field and within ROI
                    if isfield(volumes.(sub).(reg),'epvs') &&...
                        sum(epvs_cube(:)) > 1
                        vf = sum(epvs_cube(:)) ./ sum(mask_cube(:));
                        epvs_vf_mat((x:xf), (y:yf), hm_z_idx) = vf;
                    end
                end
            end
        end
       
        %% Add metrics to heatmap struct
        % Optical Properties
        heatmap.(sub).(reg).op.mus = mus_mat;
        heatmap.(sub).(reg).op.ret = ret_mat;
        heatmap.(sub).(reg).op.ori = ori_mat;
        % Volume fraction
        heatmap.(sub).(reg).ves.vf = ves_vf_mat;
        if isfield(volumes.(sub).(reg),'epvs')
            heatmap.(sub).(reg).epvs.vf = epvs_vf_mat;
        end
        
        %% Take maximum intensity projection of tissue mask  
        %%% Check to see if the mask matrix needs to be truncated
        if size(mask,1) > size(seg,1)
            mask = mask(1:size(seg,1),:,:);
        end
        if size(mask,2) > size(seg,2)
            mask = mask(:,1:size(seg,2),:);
        end
    
        % Initialize masks matrix
        masks = zeros(size(seg,1), size(seg,2), Nz);
    
        % Initialize zmin and zmax to index the tissue mask depths
        zmin = 1;
        zmax = zmin + n_z - 1;
        % Iterate over heatmap depths
        for z = 1:Nz
            % Take MIP
            masks(:,:,z) = max(mask(:,:,zmin:zmax),[],3);
            % Calculate z depth bounds for tissue mask
            zmin = zmin + n_z;
            zmax = min([(zmax + n_z), size(seg,3)]);
        end
        
        % Add mask to heatmap struct
        heatmap.(sub).(reg).mask = masks;
        fprintf('\nFINISHED HEATMAP -- %s %s\n',sub,reg)
    end
end
fprintf('Finished all subjects and regions\n')

% Save the heatmap struct
heat_out = append('heatmap_',num2str(cube_side),'.mat');
heat_out = fullfile(data_dir,'/heatmaps/', heat_out);
save(heat_out,'heatmap','-v7.3');

%% Generate heat maps - normalized across subjects
% Iterate over each metric, average across z axis, normalize
% the colorbar across all subjects for this metric.

% Load Heat map
% heat_out = append('heatmap_',num2str(cube_side),'.mat');
% heat_out = fullfile(data_dir,'/heatmaps/', heat_out);
% heatmap = load(heat_out);
% heatmap = heatmap.heatmap;
% subid = fieldnames(heatmap);


%%% Identify minimum / maximum of each metric
% Set the maximum to be the 95th percentile of each metric. This will
% scale the colorbar to account for outliers.
% ves_vf = vessel volume fraction, epvs_vf = EPVS volume fraction
% mus = scattering, ret = retardance, ori = orientation
ves_vf = [];
mus = [];
ret= [];
ori = [];

%%% Iterate through subjects and regions
for ii = 1:length(subid)
    % retrieve local subject ID
    sub = subid{ii};
    % retrieve regions for this subject
    regions = fieldnames(heatmap.(sub));
    % Iterate over tissue volume region(s)
    for j = 1:length(regions)
        % retrieve tissue mask
        mask = logical(heatmap.(subid{ii}).(regions{j}).mask);
        % retrieve heatmaps
        ves_vf_ = heatmap.(subid{ii}).(regions{j}).ves.vf;
        mus_ =  heatmap.(subid{ii}).(regions{j}).op.mus;
        ret_ =  heatmap.(subid{ii}).(regions{j}).op.ret;
        ori_ =  heatmap.(subid{ii}).(regions{j}).op.ori;
        % Take values that only lie within tissue mask
        ves_vf_ = ves_vf_(mask);
        mus_ = mus_(mask);
        ret_ = ret_(mask);
        ori_ = ori_(mask);
        % Add to struct
        ves_vf = [ves_vf; ves_vf_];
        mus = [mus; mus_];
        ret = [ret; ret_];
        ori = [ori; ori_];
    end
end



% Take real component of orientation
ori = real(ori);
% Minimum of each metric
epvs_vf_min = 0;
ves_vf_min = min(ves_vf);
mus_min = min(mus);
ret_min = min(ret);
ori_min = min(ori);
% Maximum = 95th percentile of each metric
epvs_vf_max = 5;
ves_vf_max = prctile(ves_vf,95);
mus_max = prctile(mus,95);
ret_max = prctile(ret,95);
ori_max = prctile(ori,95);

%%% Generate normalized heatmaps
% Variable for whether or not to invert the heatmap
flip_cbar = 0;
% Iterate subjects
for ii = 1:length(subid)
    % retrieve local subject ID
    sub = subid{ii};
    % retrieve regions for this subject
    regions = fieldnames(heatmap.(sub));
    % Iterate over tissue volume region(s)
    for j = 1:length(regions)    
        %%% Output filepath for figures
        roi_dir = strcat('ROI_',num2str(cube_side));
        heatmap_dir = fullfile(data_dir,'/heatmaps',sub,regions{j},roi_dir);
        if ~isfolder(heatmap_dir)
            mkdir(heatmap_dir);
        end

        %%% Import matrices
        masks = logical(heatmap.(subid{ii}).(regions{j}).mask);
        ves_vf = heatmap.(subid{ii}).(regions{j}).ves.vf;
        mus =  heatmap.(subid{ii}).(regions{j}).op.mus;
        ret =  heatmap.(subid{ii}).(regions{j}).op.ret;
        ori =  real(heatmap.(subid{ii}).(regions{j}).op.ori);
        % Take average across z dimension for metrics
        ves_vf = mean(ves_vf,3);
        mus = mean(mus,3);
        ret = mean(ret,3);
        ori = mean(ori,3);
        % Take maximum intensity projection of tissue mask
        masks = max(masks,[],3);
    
        %%% Plot first depth for each heatmap
        % Check if EPVS exists
        if isfield(heatmap.(subid{ii}).(regions{j}),'epvs')
            epvs_vf = heatmap.(subid{ii}).(regions{j}).epvs.vf;
            % Take average across z-dimension
            epvs_vf = mean(epvs_vf,3);
            % Create heatmap
            plot_save_heatmap([], epvs_vf, flip_cbar,...
                [epvs_vf_min, epvs_vf_max],...
                masks,'EPVS Volume Fraction','(unitless)',...
                heatmap_dir,'heatmap_epvs_vf')
        end
        % Vessel volume fraction
        plot_save_heatmap([], ves_vf, flip_cbar,...
            [ves_vf_min, ves_vf_max],...
            masks,'Vessel Volume Fraction','(unitless)',...
            heatmap_dir,'heatmap_ves_vf')
        % scattering coefficient
        plot_save_heatmap([], mus, flip_cbar, [mus_min, mus_max],...
            masks,'Scattering Coefficient','(cm^-^1)',...
            heatmap_dir,'heatmap_mus')
        % retardance
        plot_save_heatmap([], ret, flip_cbar, [ret_min, ret_max],...
            masks,'Retardance','(Degrees)',...
            heatmap_dir,'heatmap_ret')
        % Orientation
        plot_save_heatmap([], ori, flip_cbar, [ori_min, ori_max],...
            masks,'Circular Std. Dev. of Orientation','(Radians)',...
            heatmap_dir,'heatmap_ori')
    end
end
%}

%% Pathology heatmap: A-beta and p-tau
%{
% Initialize maximum/minimum of all heat maps
ab_max = 0;
ab_min = 1;
pt_max = 0;
pt_min = 1;
% Initialize struct to store pathology heatmap
path_heatmap = struct();
% Title above figures
tstr = {'A-beta','P-tau'};
% Types of stains
stains = {'ab','pt'};
% Filenames of each stain
stain_fname = struct();
stain_fname.AD_10382.ab = 'AD_10382_slice_8_Ab';
stain_fname.AD_10382.pt = 'AD_10382_slice_14_AT8';
stain_fname.AD_20832.ab = 'AD_20832_slice_8_Ab';
stain_fname.AD_20832.pt = 'AD_20832_slice_14_AT8';
stain_fname.AD_21354.ab = 'AD_21354_slice_8_Ab';
stain_fname.AD_21354.pt = 'AD_21354_slice_2_AT8';
stain_fname.AD_21424.ab = 'AD_21424_slice_14_Ab';
stain_fname.AD_21424.pt = 'AD_21424_slice_20_AT8';
stain_fname.CTE_6489.ab = 'CTE_6489_slice_14_Ab';
stain_fname.CTE_6489.pt = 'CTE_6489_slice_8_AT8';
stain_fname.CTE_6912.ab = 'CTE_6912_slice_8_Ab';
stain_fname.CTE_6912.pt = 'CTE_6912_slice_14_AT8';
stain_fname.CTE_7019.ab = 'CTE_7019_slice_8_Ab';
stain_fname.CTE_7019.pt = 'CTE_7019_slice_14_AT8';
stain_fname.CTE_7126.ab = 'CTE_7126_slice_14_Ab';
stain_fname.CTE_7126.pt = 'CTE_7126_slice_8_AT8';
stain_fname.NC_21499.ab = 'NC_21499_slice_14_Ab';
stain_fname.NC_21499.pt = 'NC_21499_slice_8_AT8';
stain_fname.NC_8095.ab = 'NC_8095_slice_8_Ab';
stain_fname.NC_8095.pt = 'NC_8095_slice_2_AT8';
% Minimum thresholds for segmenting plaques (after taking compliment)
stain_thresh = struct();
stain_thresh.AD_10382.ab = 0.2;
stain_thresh.AD_10382.pt = 0.2;
stain_thresh.AD_20832.ab = 0;
stain_thresh.AD_20832.pt = 0.10;
stain_thresh.AD_21354.ab = 0.23;
stain_thresh.AD_21354.pt = 0.11;
stain_thresh.AD_21424.ab = 0.15;
stain_thresh.AD_21424.pt = 0.15;
stain_thresh.CTE_6489.ab = 0.16;
stain_thresh.CTE_6489.pt = 0.13;
stain_thresh.CTE_6912.ab = 0.10;
stain_thresh.CTE_6912.pt = 0.13;
stain_thresh.CTE_7019.ab = 0.15;
stain_thresh.CTE_7019.pt = 0.15;
stain_thresh.CTE_7126.ab = 0.15;
stain_thresh.CTE_7126.pt = 0.15;
stain_thresh.NC_21499.ab = 0.10;
stain_thresh.NC_21499.pt = 0.10;
stain_thresh.NC_8095.ab = 0.10;
stain_thresh.NC_8095.pt = 0.10;
% Subject ID list containing staining
stain_subid = fieldnames(stain_fname);

for ii = 1:length(fields(stain_fname))
    for j = 1:length(stains)
        %%% A-Beta: load pathology, apply mask, calculate resolution
        fpath = fullfile(dpath, stain_subid{ii}, 'stain');
        % Filename of staining
        stain_name = stain_fname.(stain_subid{ii}).(stains{j});
        % Generate mask file name
        mask_name = append(stain_fname.(stain_subid{ii}).(stains{j}),'_mask');
        % Import the stain, mask, and calculate pxiel resolution
        [stain, mask, res] =...
            import_pathology(fpath, stain_name, mask_name);
        
        %%% generate heatmap from A-Beta staining
        th = stain_thresh.(stain_subid{ii}).(stains{j});
        [hm] = pathology_heatmap(res, cube_side, stain, th, mask);
        % Calculate maximum value of all heatmaps
        if strcmp(stains{j},'ab')
            ab_max = max(ab_max, max(hm(mask)));
            ab_min = min(ab_min, min(hm(mask)));
        else
            pt_max = max(pt_max, max(hm(mask)));
            pt_min = min(pt_min, min(hm(mask)));
        end

        %%% Add heatmap of A-beta and p-tau to struct
        path_heatmap.(stain_subid{ii}).(stains{j}).heatmap = hm;
        path_heatmap.(stain_subid{ii}).(stains{j}).mask = mask;
    end
end

%%% Plot and save the pathology heatmaps
for ii = 1:length(fields(stain_fname))
    for j = 1:length(stains)
        % Output filepath
        fpath = fullfile(dpath, stain_subid{ii}, 'stain');
        % Load the tissue mask
        mask = path_heatmap.(stain_subid{ii}).(stains{j}).mask;
        % Load the heatmap
        hm = path_heatmap.(stain_subid{ii}).(stains{j}).heatmap;
        % Figure Name
        stain_name = append(stain_fname.(stain_subid{ii}).(stains{j}),'_stain');
        figname = append(stain_name,'_heatmap_',num2str(cube_side));
        title_str = append(stain_subid{ii},' ', tstr{j});
        % Normalize each subvolume by largest dynamic range
        if strcmp(stains{j},'ab')
            hm = hm ./ ab_max;
            l = ab_min;
            u = ab_max;
        else
            hm = hm ./ pt_max;
            l = pt_min;
            u = pt_max;
        end
        % Plot heatmap
        plot_save_heatmap(1,hm,0,[0,1],mask,title_str,'(a.u.)',fpath,figname)
    end
end
%}

%% Plot and save the heat maps
function plot_save_heatmap(Ndepths, heatmaps, flip_cbar, colorbar_range,...
    masks, tstr, cbar_label, dpath, fname)
% PLOT_SAVE_HEATMAP use imagesc and set background = 0
% INPUT
%   Ndepths (int): number of depths in z dimension
%   heatmaps (double matrix): heatmaps of vascular metric
%   flip_cbar (logical): reverse the direction of the colorbar
%   colorbar_range (double array): [min, max]
%   masks (double): tissue mask (1=tissue, 0=other)
%   tstr (string): figure title
%   cbar_label (string): colorbar label
%   dpath (string): data directory path
%   fname (string): name of figure to save

%%% Set the number of depths to iterate for each heatmap
% If the number of depths is not specified, then set it equal to the number
% of z dimensions.
if isempty(Ndepths)
    Ndepths = size(heatmaps,3);
end
% Set fontsize for the heatmap figure
fontsize = 40;

%%% Iterate over frames in z dimension
for d = 1:Ndepths
    %%% Heatmap of d_th frame from the length density
    fh = figure();
    fh.WindowState = 'maximized';
    % If there are multiple heatmaps in the matrix
    if size(heatmaps,3) > 1
        heatmap = heatmaps(:,:,d);
    % Here it is just a single frame of a heatmap
    else
        heatmap = heatmaps;
    end
    % Initialize heatmap
    h = imagesc(heatmap);

    %%% Initialize colorbar
    % If the colorbar_range is passed in, then extract min & max
    if ~isempty(colorbar_range)
        cmap_min = colorbar_range(1);
        cmap_max = colorbar_range(2);
    % Otherwise, set limits from the current heatmap
    else
        % Min = Lowest value also greater than zero
        cmap_min = min(heatmap(heatmap(:)>0));
        % Max = Find the 95th percentile for upper limit
        cmap_max = prctile(heatmap(:), 95);
    end
    % Initialize the colormap limits
    cmap = jet(256);
    clim(gca, [cmap_min, cmap_max]);
    % Initialize colormap and colorbar
    if flip_cbar
        colormap(flipud(cmap));
    else
        colormap(cmap);
    end
    c = colorbar;

    %%% Apply tissue mask to the heatmap to remove background
    alpha_mask = double(masks(:,:,d));
    set(h, 'AlphaData', alpha_mask);

    %%% Configure figure parameters
    % Update title string with depth
    if size(heatmaps,3) > 1
        title_str = append(tstr, ' Depth ', num2str(d));
    else
        title_str = tstr;
    end
    title(title_str,'Interpreter','none');
    set(gca, 'FontSize', fontsize);
    % Label the colorbar    
    c.Label.String = cbar_label;
    % Offset colorbar label to the right of colorbar
    c.Label.Position = [10 (cmap_max - (cmap_max-cmap_min)/2)];
    c.Label.Rotation = 270;
    % Increase fontsize of colorbar
    c.FontSize = 40;
    % Remove x and y tick labels
    set(gca,'Yticklabel',[]);
    set(gca,'Xticklabel',[]);
    
    %%% Save figure as PNG
    % If there are multiple heatmaps in the matrix, save vascular heatmap
    if size(heatmaps,3) > 1
        fout = append(fname, '_', num2str(d));
    % Otherwise, save the pathology heatmap
    else
        fout = fname;
    end
    % If the colorbar is reversed then add suffix to filename
    if flip_cbar
        fout = append(fout, '_flip_cbar');
    end
    % Save figure
    fout = fullfile(dpath, fout);
    saveas(gca, fout,'png');
    pause(0.1)
    close;
end
end

%% Pathology heatmap
function [hm] = pathology_heatmap(pix, square_side, stain, th, mask)
%PATHOLOGY_HEATMAP generate heatmap with pathology
% Divide the pathology stain into isotropic squares, calculate the density
% of each square, add the density to a matrix.
%
% INPUTS:
%   pix (double array): pixel size (x,y) (microns)
%   square_side (uint): heatmap square dimension (microns)
%   stain (uint8 matrix): pathology stain
%   th (double): minimum threshold for segmenting AB or AT8 (after taking
%               complement of the staining)
%   mask (logical matrix): tissue mask for the pathology stain
% OUTPUTS:
%   hm (double matrix): masked heatmap of the pathology stain

%%% Verify that the stained image is scaled between [0,1]
% The inversion assumes that the image is scaled between [0,1], so this
% will return the incorrect value otherwise.
r = range(stain);
assert(0<=min(r), 'The stained image is not scaled between [0,1]');
assert(max(r)<=1, 'The stained image is not scaled between [0,1]');

%%% Translate cube size to pixel size
% Number of pixels in each dimension to create square
nx = floor(square_side ./ pix(1));
ny = floor(square_side ./ pix(2));

%%% Complement image and segment plaques (minimum threshold)
% Complement the image (you look nice, today)
stain = imcomplement(stain);
% Set pixels above background intensity equal to 1
stain(stain > th) = 1;
% Mask the stain to remove non-tissue pixels
stain(~mask) = 0;
% Set pixels < 1 equal to zero
stain(stain ~= 1) = 0;

%%% Generate heat map
% Initialize heatmap matrices for A-beta and p-tau
hm = zeros(size(stain,1), size(stain,2));
% Iterate over rows
for x = 1:nx:size(stain,1)
    % Iterate over columns
    for y = 1:ny:size(stain,2)
        %%% Crop segmentation into isotropic square
        % Initialize end indices for each axis
        xf = x + nx - 1;
        yf = y + ny - 1;
        % Take minimum of matrix dimensions and end indices
        xf = min(xf, size(stain,1));
        yf = min(yf, size(stain,2));
        % Take cube from segmentation
        path_square = stain((x:xf), (y:yf));

        %%% Calculate average plaque intensity in subvolume
        hm((x:xf), (y:yf)) = mean(path_square(:));
    end
end

% Mask the pathology heatmap with tissue mask
hm(~mask) = 0;
end

%% Read TIFF and extract image properties
function [stain, mask, res] = import_pathology(dpath, stain_name, mask_name)
%IMPORT_PATHOLOGY read the TIF file for the pathology (A-beta or p-tau)
% Retrieve TIFF metadata for the staining & compute the resolution of each
% pixel. Then, load the TIFF and the respective tissue mask. Rescale the
% staining to [0,1] to standardize across all images. Then, apply the
% tissue mask.
% INPUTS:
%   dpath (string): directory path to the staining TIFF file
%   stain_name (string): filename of staining TIFF file
%   mask_name (string): filename of tissue mask TIFF file
% OUTPUTS:
%   stain (double matrix): masked staining, scaled between [0,1]
%   mask (logical matrix): tissue mask for the pathology slide
%   bg (double matrix): background patch for normalizing the pathology.
%                       This is a region of the stain from the white matter
%                       that does not include any plaques.
%   res (double array): resolution of a pixel [x,y]

%%% Disable Tiff read warnings for unknown tags
id = 'imageio:tiffmexutils:libtiffWarning';
warning('off',id);

%%% Retrieve metadata and calculate resolution
% Read staining TIFF info
stain_file = append(stain_name, '.tif');
stain_file = fullfile(dpath, stain_file);
stain_metadata = imfinfo(stain_file);
% X and Y resolution (pixels / micron)
xres = stain_metadata.XResolution;
yres = stain_metadata.YResolution;
% Invert resolution (microns / pixel)
xres = 1./xres;
yres = 1./yres;
res = [xres, yres];

%%% Load TIFF file
% Read staining TIFF (RGB) & convert to grayscale
stain = Tiff(stain_file,'r');
stain = read(stain);
stain = rescale(stain, 0, 1);
stain = rgb2gray(stain);
% Load tissue mask
mask_fname = append(mask_name, '.tif');
mask_file = fullfile(dpath, mask_fname);
mask = Tiff(mask_file,'r');
mask = logical(read(mask));
end
