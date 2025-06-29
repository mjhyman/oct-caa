%% Divide volume into grid and calculate metrics for each cube
% The purpose of this is to subdivide the entire tissue volume in cubes of
% dimensions [M, M, M] cubic millimeters and then generate the graph for
% each cube. The following step is to generate the heat map from these
% cubes.

%% Add top-level directory of code repository to path
clear; clc; close all;
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

%% Load heat maps
heat_out = append('heatmap_',num2str(cube_side),'.mat');
heat_out = fullfile(data_dir,'/heatmaps/', heat_out);
heatmap = load(heat_out);
heatmap = heatmap.heatmap;
subid = fieldnames(heatmap);

%% Create 2D arrays of optical property vs. EPVS density
% 2xN Matrix for each optical property:
%   - top row = EPVS density from heatmap
%   - bottom row = optical property from heatmap
%   - column = pairwise observations from same tissue volume

%%% Initialize arrays to store pairs
epvs_mus = [];
epvs_ret = [];
epvs_ori = [];

%%% Iterate over each subject
for ii = 1:length(subid)
    % retrieve regions for this subject
    sub = subid{ii};
    regions = fieldnames(heatmap.(sub));
    % Iterate over tissue volume region(s)
    for j = 1:length(regions)
        % Skip if no EPVS is present
        reg = regions{j};
        if ~isfield(heatmap.(subid{ii}).(regions{j}),'epvs')
            fprintf('skipping %s %s\n',sub,reg)
            continue
        end
        fprintf('Measuring %s %s\n',sub,reg)
        % Load EPVS and optical properties
        epvs = heatmap.(subid{ii}).(regions{j}).epvs.vf;
        mus = heatmap.(subid{ii}).(regions{j}).op.mus;
        ret = heatmap.(subid{ii}).(regions{j}).op.ret;
        ori = heatmap.(subid{ii}).(regions{j}).op.ori;
        %% Iterate over heatmap patches
        % Iterate over the z-axis
        for z = 1:size(epvs,3)
            % Iterate over rows
            for x = 1:n_x:size(epvs,1)
                % Iterate over columns
                for y = 1:n_y:size(epvs,2)
                    %% Crop matrices into cubes
                    % Initialize end indices for each axis
                    xf = x + n_x - 1;
                    yf = y + n_y - 1;
                    % Take minimum of matrix dimensions and end indices
                    xf = min(xf, size(epvs,1));
                    yf = min(yf, size(epvs,2));
                    %%% Take cube from epvs and check if contains EPVS
                    epvs_cube = epvs((x:xf), (y:yf), z);
                    if ~all(epvs_cube(:) == 0)
                        % Take cube of remaining optical properties
                        mus_cube = mus((x:xf), (y:yf), z);
                        ret_cube = ret((x:xf), (y:yf), z);
                        ori_cube = ori((x:xf), (y:yf), z);
                        % Create pair of EPVS and optical property
                        mus_pair = [min(epvs_cube); min(mus_cube)];
                        ret_pair = [min(epvs_cube); min(ret_cube)];
                        ori_pair = [min(epvs_cube); min(ori_cube)];
                        % Add pair to matrix
                        epvs_mus = [epvs_mus,mus_pair];
                        epvs_ret = [epvs_ret,ret_pair];
                        epvs_ori = [epvs_ori,ori_pair];
                    end
                end
            end
        end
    end
end

%% Plot optical property vs. EPVS density

% scattering coefficient
x = epvs_mus(1,:); y = epvs_mus(2,:);
figure; scatter(x,y,20,'k','filled');
xlabel('EPVS Volume Fraction (percentage occupied by EPVS)')
ylabel('\mu_s (cm^-^1)')
title('\mu_s vs. EPVS Density');
set(gca,'fontsize',20)


% retardance
x = epvs_ret(1,:); y = epvs_ret(2,:);
figure; scatter(x,y,20,'k','filled');
xlabel('EPVS Volume Fraction (percentage occupied by EPVS)')
ylabel('retardance (degrees)')
title('Retardance vs. EPVS Density')
set(gca,'fontsize',20)

% orientation
x = epvs_ori(1,:); y = epvs_ori(2,:);
figure; scatter(x,y,20,'k','filled');
xlabel('EPVS Volume Fraction (percentage occupied by EPVS)')
ylabel('\sigma (radians)')
title('Circular \sigma of Orientation vs. EPVS Density')
set(gca,'fontsize',20)

%% Spearman's rho correaltion coefficient

% scattering coefficient
x = epvs_mus(1,:); y = epvs_mus(2,:);
mus_rho = corr(x',y','type','Spearman','rows','complete');

% retardance
x = epvs_ret(1,:); y = epvs_ret(2,:);
ret_rho = corr(x',y','type','Spearman','rows','complete');

% orientation
x = epvs_ori(1,:); y = epvs_ori(2,:);
ori_rho = corr(x',y','type','Spearman','rows','complete');
