%% Export Histology Donut Measurements to Spreadsheet
%{
The histopathology measurements are stored in matlab structs. The purpose
of this is to export these to spreadsheets for each stain so that the
Bayesian statistical model can be used in R.

%}

%% Top-level settings
clear; clc; close all;
% Get the current folder
currentFolder = pwd;
% Move one directory up
parentFolder = fileparts(currentFolder);
% Add the parent folder and all its subfolders to the MATLAB search path
addpath(genpath(parentFolder));

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/';
% Directory to save output spreadsheet
sheet_dir = '/projectnb/npbssmic/ns/CAA/histology/stats/';
dt = datetime('now', 'Format', 'yyyy-MMM-dd');
fname = strcat('histogram_matched_donuts_',string(dt),'.xlsx');
output_file = fullfile(sheet_dir,fname);

%%% Load the stain and segmentation structs
lhe = load(fullfile(hdir,'lhe_rings_21-Nov-2025.mat')); lhe = lhe.lhe;
gfap = load(fullfile(hdir,'gfap_rings_21-Nov-2025.mat')); gfap = gfap.gfap;
cd68 = load(fullfile(hdir,'cd68_rings_21-Nov-2025.mat')); cd68 = cd68.cd68;


%% Export to Spreadsheet
% Make sheet for each stain, each sheet will have:
% Columnns: group(experimental/control),
%           subID,
%           stage (control = 0, mild=1, mod=2, severe=3),
%           region (front/occip),
%           distance, stain(z-score)

%% Export to Spreadsheet (Scalar Averages Edition)
clc;
stain_names = {'lhe','gfap','cd68'};

% Start fresh: delete existing file if it exists
if exist(output_file, 'file'); delete(output_file); end

% Combine stains for ease of access
stains_struct = struct('lhe', lhe, 'gfap', gfap, 'cd68', cd68);

for ii = 1:numel(stain_names)
    this_stain = stain_names{ii};
    sections = stains_struct.(this_stain);
    
    % Initialize cell array to collect data for this sheet
    % Columns: Groups, subjectID, Stage, Region, Distance, stain_value
    all_data = {}; 
    
    for j = 1:numel(sections)
        name = sections(j).baseName;
        
        % --- Subject & Stage Mapping ---
        if contains(name,'caa6','IgnoreCase',true) || contains(name,'caa_6','IgnoreCase',true)
            new_name = 'caa6'; stage = 1;
        elseif contains(name,'caa17','IgnoreCase',true) || contains(name,'caa_17','IgnoreCase',true)
            new_name = 'caa17'; stage = 2;
        elseif contains(name,'caa22','IgnoreCase',true) || contains(name,'caa_22','IgnoreCase',true)
            new_name = 'caa22'; stage = 3;
        elseif contains(name,'caa25','IgnoreCase',true) || contains(name,'caa_25','IgnoreCase',true)
            new_name = 'caa25'; stage = 3;
        else
            new_name = 'caa26'; stage = 0;
        end

        % --- Region Mapping ---
        tokens = regexp(name, 'CAA_?\d+_(\d+)', 'tokens');
        region = 'unknown'; 
        if ~isempty(tokens)
            regionNum = str2double(tokens{1}{1});
            if regionNum == 1, region = 'front';
            elseif regionNum == 7, region = 'occip'; end
        end

        % --- Process rad40 Substruct ---
        if isfield(sections(j), 'rad40')
            rad_fields = fieldnames(sections(j).rad40);
            
            for k = 1:numel(rad_fields)
                rad_name = rad_fields{k};
                distance = str2double(regexp(rad_name, '\d+', 'match', 'once'));
                
                % 1. Extract values and force to COLUMN vectors using (:)
                exp_vals = sections(j).rad40.(rad_name).exp_hmatched(:);
                ctl_vals = sections(j).rad40.(rad_name).ctl_hmatched(:);
                
                % 2. Process Experimental Data
                if ~isempty(exp_vals)
                    n_exp = numel(exp_vals);
                    % Duplicate metadata for every entry in the column vector
                    rows_exp = [repmat({"Experimental"}, n_exp, 1), ...
                                repmat({new_name}, n_exp, 1), ...
                                repmat({stage}, n_exp, 1), ...
                                repmat({region}, n_exp, 1), ...
                                repmat({distance}, n_exp, 1), ...
                                num2cell(exp_vals)]; % Value column
                    all_data = [all_data; rows_exp]; 
                end
                
                % 3. Process Control Data
                if ~isempty(ctl_vals)
                    n_ctl = numel(ctl_vals);
                    % Duplicate metadata for every entry in the column vector
                    rows_ctl = [repmat({"Control"}, n_ctl, 1), ...
                                repmat({new_name}, n_ctl, 1), ...
                                repmat({stage}, n_ctl, 1), ...
                                repmat({region}, n_ctl, 1), ...
                                repmat({distance}, n_ctl, 1), ...
                                num2cell(ctl_vals)]; % Value column
                    all_data = [all_data; rows_ctl]; 
                end
            end
        end
    end
    
    % --- Create Table and Write to Excel ---
    if ~isempty(all_data)
        T = cell2table(all_data, 'VariableNames', ...
            {'Groups', 'subjectID', 'Stage', 'Region', 'Distance', 'stain_value'});
        
        % Write the current stain as its own tab
        writetable(T, output_file, 'Sheet', this_stain);
        fprintf('Sheet "%s" complete. Rows: %d\n', this_stain, size(T,1));
    else
        fprintf('Warning: No valid data found for stain: %s\n', this_stain);
    end
end

fprintf('Export finished! File saved as: %s\n', output_file);