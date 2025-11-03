%% Import matlab structs of donuts
%{
Purpose of this script:
- import structs
- wilcoxon signed rank test at each radii
%}

%% Top-level settings
clear; clc; close all;

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/histology/figures';

%%% Output directory to save output
stat_sheet = fullfile(hdir,'stats/histo_iterated_stats_23Oct2025.xlsx');

%% Import the stains
cd68 = load(fullfile(hdir,"cd68_rings_15-Oct-2025.mat"));
gfap = load(fullfile(hdir,"gfap_rings_15-Oct-2025.mat"));
lhe = load(fullfile(hdir,"lhe_rings_15-Oct-2025.mat"));

% Combine all three stains into a single struct for ease of access
stains = struct();
stains.cd68 = cd68.cd68;
stains.gfap = gfap.gfap;
stains.lhe = lhe.lhe;

%% Iterate over each stain and radius

% Number of stains for each radius size
nsmall = length(fields(stain(1).rad40));
nlarge = length(fields(stain(1).rad100));

% Initialize struct for storing stats
stats = struct();

% Iterate over stains
stain_name = fields(stains);

for ii = 1:length(stain_name)
    % Extract local stain
    stain = stains.(stain_name{ii});
    % Number of tissue sections in stain
    nsec = length(stain);

    % Iterate over small radii
    for j = 1:small
        stats.rad40
    end
    pause(0.1)
end


%% One-sided wilcoxon signed rank test
stain_names = fieldnames(stains);
radius_types = {'rad40', 'rad100'}; % Main radius structs
result_rows = {}; % For results table

%%% State which tail to use for test
% Hypothesize myelin rarefaction around EPVS
%   --> (med(EPVS) - med(control)) < 0
%   --> left tailed
tails.cd68 = 'right';
tails.gfap = 'right';
tails.lhe = 'left';

for s = 1:numel(stain_names)
    stain = stain_names{s};
    section_struct = stains.(stain); % Array of sections
    
    for rt = 1:numel(radius_types)
        radtype = radius_types{rt};
        
        % Get all subradius names from the first section
        subrads = fields(section_struct(1).(radtype));
        
        for sub = 1:numel(subrads)
            subrad = subrads{sub};
            exp_vec = [];
            ctl_vec = [];

            for sec = 1:numel(section_struct)
                % If subradius exists & has needed fields
                if isfield(section_struct(sec).(radtype), subrad)
                    rad_struct = section_struct(sec).(radtype).(subrad);
                    if isfield(rad_struct,'exp_mean') &&...
                            isfield(rad_struct,'ctl_mean')
                        exp_vec(end+1) = rad_struct.exp_mean;
                        ctl_vec(end+1) = rad_struct.ctl_mean;
                    end
                end
            end
            
            if ~isempty(exp_vec) && ~isempty(ctl_vec) &&...
                    numel(exp_vec)==numel(ctl_vec)
                [p, h, stats] = signrank(exp_vec, ctl_vec,...
                                        'method','approximate',...
                                        'tail',tails.(stain));
                radius_value = sscanf(subrad, 'rad%d');
                med_exp = median(exp_vec);
                med_ctl = median(ctl_vec);
                result_rows(end+1,:) = {stain, radtype, radius_value, ...
                                        med_exp, med_ctl, p,...
                                        stats.signedrank, stats.zval,...
                                        numel(exp_vec)};
            end
        end
    end
end

stat_tbl = cell2table(result_rows, ...
    'VariableNames', {'Stain', 'RadiusType', 'Radius', 'MedEPVS',...
                    'MedVessel', 'P_Value',...
                    'SignedRankStatistic','Z', 'NumPairs'});

writetable(stat_tbl, stat_sheet);