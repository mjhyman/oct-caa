%% Import matlab structs of donuts
%{
Purpose of this script:
- import structs
- wilcoxon signed rank test at each radii
%}

%% Top-level settings
clc; close all;
clearvars -except cd68 gfap lhe

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/histology/figures';

%%% Flags for using histo-matched + z-score
% If false -> use just histogram matching
% If true -> use histogram matching + z-score 
zflag = false;
if zflag
    % Output directory/filename to save output spreadsheet
    stat_sheet = fullfile(hdir,['stats/histo_iterated_stats_' ...
                                'hmatched_zscore_21Nov2025.xlsx']);
else
    % Output directory/filename to save output spreadsheet
    stat_sheet = fullfile(hdir,['stats/histo_iterated_stats_' ...
                                'hmatched_21Nov2025.xlsx']);
end

%% Import the stains
if ~exist('lhe','var') && ~exist('gfap','var') && ~exist('cd68','var')
    % Load stains from files
    fprintf('\nLoading CD68\n')
    cd68 = load(fullfile(hdir,"cd68_rings_21-Nov-2025.mat"));
    fprintf('\nLoading GFAP\n')
    gfap = load(fullfile(hdir,"gfap_rings_21-Nov-2025.mat"));
    fprintf('\nLoading LHE\n')
    lhe = load(fullfile(hdir,"lhe_rings_21-Nov-2025.mat"));
    % Combine all three stains into a single struct for ease of access
    stains = struct();
    stains.cd68 = cd68.cd68;
    stains.gfap = gfap.gfap;
    stains.lhe = lhe.lhe;
end

%% One-sided wilcoxon signed rank test

%%% Initialization
stain_names = fieldnames(stains);   % retrieve names of all stains
radius_types = {'rad40', 'rad100'}; % Define the radii sizes that were used
result_rows = {}; % initialize cell array for results table

%%% State which tail to use for test
% Hypothesize myelin rarefaction around EPVS
%   --> (med(EPVS) - med(control)) < 0
%   --> left tailed
tails.cd68 = 'right';
tails.gfap = 'right';
tails.lhe = 'left';
% Set all as both
% tails.cd68 = 'both';
% tails.gfap = 'both';
% tails.lhe = 'both';

for s = 1:numel(stain_names)
    stain = stain_names{s};
    section_struct = stains.(stain); % Array of sections
    %%% Iterate radii
    for rt = 1:numel(radius_types)
        radtype = radius_types{rt};
        % Get all subradius names from the first section
        subrads = fields(section_struct(1).(radtype));
        % Iterate subjects
        for sub = 1:numel(subrads)
            subrad = subrads{sub};
            exp_vec = [];
            ctl_vec = [];
            % Iterate sections
            for sec = 1:numel(section_struct)
                % If subradius exists & has needed fields
                if isfield(section_struct(sec).(radtype), subrad)
                    rad_struct = section_struct(sec).(radtype).(subrad);
                    % Add section means to exp_vec and ctl_vec
                    if zflag
                        exp_vec(end+1) = rad_struct.exp_mean;
                        ctl_vec(end+1) = rad_struct.ctl_mean;
                    else
                        exp_vec(end+1) = rad_struct.exp_hmatched_mean;
                        ctl_vec(end+1) = rad_struct.ctl_hmatched_mean;
                    end
                end
            end
            
            % Statistical test 
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