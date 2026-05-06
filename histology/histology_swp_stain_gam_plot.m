%% Plot Differences in the Histology GAM b/w severities
% Overview:
%{
- import main GAM struct (this generated a GAM for each stain and section)
- Compare each slice (i.e. severe - control, moderate - control, etc.)
%}

%% Initialization
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
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
ncores = feature('numcores');
maxNumCompThreads(ncores);

%% Top-level settings
%%% Directories (SCC)
% Input directory
ddir='/projectnb/npbssmic/ns/CAA/histology/gam';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/figures/histology_gam/';

%%% Import GAM
gam = load(fullfile(ddir,'GAM_histo_01-May-2026.mat'));
gam = gam.histo_gam;

%%% Output figure settings
% Date for file output
t = datetime('today');
% Font size for figures
fsize = 30;

%% Rename sections to severity
% CAA26 = control
% CAA6 = mild
% CAA17 = mild
% CAA22/25 = sever
% The second number: 1 = frontal, 7 = occipital

% 1. Get the top-level stains (cd68, gfap, lhe)
stains = fieldnames(gam);
gam_remapped = struct();

for i = 1:numel(stains)
    stainType = stains{i};
    % Pre-initialize stain field to avoid assignment errors
    gam_remapped.(stainType) = struct('frontal', struct(), 'occipital', struct());
    
    sections = fieldnames(gam.(stainType));
    
    for j = 1:numel(sections)
        oldName = sections{j};
        data = gam.(stainType).(oldName);
        
        % --- REGION IDENTIFICATION ---
        % Checks for 'front', '_1_', or ending in '_1'
        if contains(oldName, 'front') || contains(oldName, '_1_') || endsWith(oldName, '_1')
            region = 'frontal';
        elseif contains(oldName, 'occip') || contains(oldName, '_7_') || endsWith(oldName, '_7')
            region = 'occipital';
        else
            % Skip if it's a metadata field or something else
            continue; 
        end
        
        % --- SEVERITY MAPPING ---
        % Order matters here: Check specific subjects before 'combined'
        if contains(oldName, 'CAA26') || contains(oldName, 'CAA_26')
            sev = 'control';
        elseif contains(oldName, 'CAA6') || contains(oldName, 'CAA_6')
            sev = 'mild';
        elseif contains(oldName, 'CAA17') || contains(oldName, 'CAA_17')
            sev = 'moderate';
        elseif contains(oldName, 'CAA22') || contains(oldName, 'CAA_22')
            sev = 'severe_CAA22';
        elseif contains(oldName, 'CAA25') || contains(oldName, 'CAA_25')
            sev = 'severe_CAA25';
        elseif strcmpi(oldName, 'severe_front') || strcmpi(oldName, 'severe_occip')
            sev = 'severe_combined';
        else
            continue;
        end
        
        % --- ASSIGNMENT ---
        gam_remapped.(stainType).(region).(sev) = data;
    end
end

%% Compare the GAM slices for each stain + region separately

% Retrieve stain names
stains = fields(gam_remapped);
regions = fields(gam_remapped.cd68);
ylims = [-1,1];

% Iterate Stains
for ii = 1:numel(fields(gam_remapped))
    % Iterate Regions
    for j = 1:2
        % print to console
        fprintf('\n----%s %s----\n',stains{ii},regions{j})
        dirout = fullfile(figdir,stains{ii},regions{j});
        % Perform Slice Comparison
        histo_swp_compare_gam_slices(...
            gam_remapped.(stains{ii}).(regions{j}),...
                          stains{ii},...
                          regions{j},...
                          ylims, dirout);
    end
end
