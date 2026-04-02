%% Measure size-weighted proximity (SWP) of stains
% Purpose: measure SWP for each section and stain

%% Top-level settings
clc; close all; clear;

%%% Set maximum number of computational threads
maxNumCompThreads(28);

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/';
% Output directory for SWP data
swp_dir = fullfile(hdir, '/swp/');

%%% Flags for using histo-matched + z-score
% If false -> use just histogram matching
% If true -> use histogram matching + z-score 
zflag = true;
% Output file
t = datetime('today');
if zflag
    fname = strcat('histo_swp_hmatched_zscore_');
    swp_mat = fullfile(swp_dir,fname);
else
    fname = strcat('histo_swp_hmatched_');
    swp_mat = fullfile(swp_dir,fname);
end

%% Size weighted proximity (SWP) settings
%%% Radius search size
% Pixel size (um/pix)
pix = 1.3067;
% Radius (mm)
radius_mm = 4;
% Search radius in pixels
radius = radius_mm .* 1000 / pix;

%%% Other parameters
% Exponent p-value (constant)
p = 2;
% downsample factor for SWP
ds = 4;
% Parallel computing flag
parpool_flag = true;

%%% Output filename extension including radius_mm, ds
fname_ext = strcat('_voxelwise_radius',num2str(radius_mm),'mm_','ds',num2str(ds),...
                   '_',string(t),'.mat');

%% Import the SGE Task ID from the batch .SH script
% Load the SGE Task ID from the batch script
taskid = str2num(getenv('SGE_TASK_ID'));
fprintf('SGE Task ID: %i\n', taskid);
% If no task ID, then default to 3 for testing
if isempty(taskid)
    fprintf('Failed to retrieve SGE Task ID. Defaulting to task ID 3.\n');
    taskid = 1; % Default to GFAP if no task ID is found
end

%% Import stain according to SGE Task ID

if taskid == 1
    %%% LHE
    % Stain directory
    stain_dir = fullfile(hdir,'LHE/');
    % stain suffix
    stain_suffix = '_M_BackgroundMask.tif';
    % EPVS suffix
    epvs_suffix = '_EPVS_Mask.tif';
    % Vessel filename suffix
    ves_suffix = '_VESSEL_Mask.tif';
    % Mask suffix
    mask_suffix = '_mask.tif';
    % Import LHE stain
    fprintf('\nIMPORTING LHE\n')
    [lhe] = import_stain(stain_dir,stain_suffix,epvs_suffix, ...
                        ves_suffix,mask_suffix);
    % Histogram match and z-score normalization
    lhe = histo_match_zscore(lhe);
    % Measure LHE SWP
    [lhe] = iterate_sections(lhe, 'lhe', radius, p, ds);
    %%% Save output
    fout = strcat(swp_mat,'_LHE',fname_ext);
    save(fout,"lhe",'-v7.3')
elseif taskid == 2
    %%% CD68 Import Stains
    % Stain directory
    stain_dir = fullfile(hdir,'CD68/');
    % stain suffix
    stain_suffix = '_shrunk_deconv.tif';
    % EPVS suffix
    epvs_suffix = '_shrunk_EPVS.tif';
    % Vessel filename suffix
    ves_suffix = '_shrunk_VESSEL.tif';
    % Mask suffix
    mask_suffix = '_shrunk_mask.tif';
    % Import CD68 stain
    fprintf('\nIMPORTING CD68\n')
    [cd68] = import_stain(stain_dir,stain_suffix,epvs_suffix, ...
                        ves_suffix,mask_suffix);
    % Histogram match and z-score normalization
    cd68 = histo_match_zscore(cd68);
    % Measure CD68 SWP
    [cd68] = iterate_sections(cd68, 'cd68', radius, p, ds);
    %%% Save output
    fout = strcat(swp_mat,'_CD68',fname_ext);
    save(fout,"cd68",'-v7.3')
else
    %%% GFAP Import Stains
    % Stain directory
    stain_dir = fullfile(hdir,'GFAP/');
    % stain suffix
    stain_suffix = '_shrunk_deconv.tif';
    % EPVS suffix
    epvs_suffix = '_shrunk_EPVS.tif';
    % Vessel filename suffix
    ves_suffix = '_shrunk_VESSEL.tif';
    % Mask suffix
    mask_suffix = '_shrunk_mask.tif';
    % Import GFAP stain
    fprintf('\nIMPORTING GFAP\n')
    [gfap] = import_stain(stain_dir,stain_suffix,epvs_suffix, ...
                        ves_suffix,mask_suffix);
    % Histogram match and z-score normalization
    gfap = histo_match_zscore(gfap);
    % Measure gfap SWP
    [gfap] = iterate_sections(gfap, 'gfap', radius, p, ds);
    %%% Save output
    fout = strcat(swp_mat,'_GFAP',fname_ext);
    save(fout,"gfap",'-v7.3')
end

%% Function to iterate section

function [stain] = iterate_sections(stain,sname,radius,p,ds)
% ITERATE_SECTIONS measure SWP for each section
% INPUTS:
%   stain (struct): stain struct
%   sname (string): name of stain
%   radius (int): search radius in pixels
%   p (int): the exponent
%   ds (int): downsample factor (measure every "ds" pixel)
%   parpool_flag (bool): whether to use parallel processing (default false)
% OUTPUTS:

% Iterate over sections
for ii = 1:length(stain)
    % Retrieve epvs, vessels, mask
    epvs = stain(ii).epvs; 
    mask = stain(ii).mask;
    ves = stain(ii).ves;

    % Within mask, set vessels to 0 to exclude from SWP
    mask(ves) = 0;
    fprintf('\nMeasuring SWP for stain %s section %i/%i\n',...
            sname,ii,length(stain));
    
    %%% Calculate EPVS SWP (subsampled and interpolated)
    [~,~,epvs_swp] = planar_size_weighted_proximity(epvs, mask, ves, radius, p,'all');
    stain(ii).epvs_swp = epvs_swp;

    %%% Calculate Vessel SWP (subsampled and interpolated)
    % swap the epvs and vessel input arguments
    [~,~,ves_swp] = planar_size_weighted_proximity(ves, mask, epvs, radius, p,'all');
    stain(ii).ves_swp = ves_swp;

    %%% Add parameters to struct (radius, p, and ds)
    stain(ii).radius = radius;
    stain(ii).p = p;
    stain(ii).downsample = ds;
    % Print to console
    fprintf('\nFinished SWP for stain %s section %i/%i\n',...
            sname,ii,length(stain));
end
end