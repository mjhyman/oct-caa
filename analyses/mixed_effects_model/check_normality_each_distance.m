%% Check normality of percentage difference dataset
% Import the spreadsheet w/ the table of the table of organized data
% Plot each distance (separate vessel + EPVS)
% Use lillyfors test at each distance (separate vessel + EPVS)
% If all are normal, then we can use a GLME

%% Add top-level directory of code repository to path
clc; close all; clear;
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
% Directory containing seg, mus, ret, mask, epvs structs
data_dir = '/projectnb/npbssmic/ns/CAA/';
% Directory for storing output
lintest_dir = '/projectnb/npbssmic/ns/CAA/glme/linearity_tests/';

%% Import spreadsheet w/ table of organized percentage difference data
fprintf('\nImporting Data\n')
% Load the percentage difference data from the spreadsheet
scattering = readtable(fullfile(data_dir,'caa_all_radii_40um_donut_14-01-2026'),...
                        'Sheet','scattering');
% Load percentage difference retardance
retardance = readtable(fullfile(data_dir,'caa_all_radii_40um_donut_14-01-2026'),...
                        'Sheet','retardance');
%% Combine both tables into a structure for easier access
data.scattering = scattering;
data.retardance = retardance;

%%% Ensure correct categorical and double
f = fieldnames(data);
for ii = 1 : numel(f)
    % Convert groups, region, and subjecID to categorical
    data.(f{ii}).Groups = categorical(data.(f{ii}).Groups);
    data.(f{ii}).Region = categorical(data.(f{ii}).Region);
    data.(f{ii}).subjectID = categorical(data.(f{ii}).subjectID);
    % Convert Stage, distance, opticalproperty to doubles
    data.(f{ii}).Stage = double(data.(f{ii}).Stage);
    data.(f{ii}).distance = double(data.(f{ii}).distance);
    data.(f{ii}).OpticalProperty = double(data.(f{ii}).OpticalProperty);
    
    %%% Convert distance to normalized range
    % d_max = max(data.(f{ii}).distance);
    % d_min = min(data.(f{ii}).distance);
    % range_dist = d_max - d_min;
    % tmp = (data.(f{ii}).distance - d_min) / range_dist;
    % data.(f{ii}).distance = double(tmp);

    %%% Log-transform optical property
    data.(f{ii}).OpticalProperty = log(data.(f{ii}).OpticalProperty);
end

%% Test GLME at each distance
% Define optical properites
op = {'scattering','retardance'};
% Define distances
distances = 40 : 40 : 480;
% define regions
regions = {'front','occip'};
% structure for storing results of lilliefors
normal_results = struct();

%%% Define the GLME formula
% fixed effects = Groups, disease stage, region
% random intercept = subjectID
fml = 'OpticalProperty ~ Groups + Stage + (1 | subjectID)';
%%% Define the numeric predictors
numericPredictors = {'Stage'};
% struct for storing results of random and fixed effects
effects = struct();

%%% Iterate optical property (scattering or retardance)
for k = 1:2
    % Extract respective table
    tbl = data.(op{k});
    %%% Iterate distances
    for ii = 1:length(distances)
        %%% Iterate regions (front or occip)
        for j = 1:length(regions)
            %% Define GLME
            %%% Check the GLME Linearity
            % Create conditions for subsetting
            condition_distance = (tbl.distance == distances(ii));
            condition_region = (tbl.Region == regions(j));
        
            % Display the conditions for debugging
            disp(['\nSubsetting with Region: ' char(regions(j))]);
        
            % Combine conditions to create the subset
            tbl_subset = tbl(condition_distance & condition_region, :);

            %%% Fit the GLME model
            fprintf('\nFitting Model\n')
            glme = fitglme(tbl_subset, fml);
            
            %%% Test Linearity
            check_glme_each_distance(glme, tbl_subset,...
                                numericPredictors,...
                                0.05, lintest_dir,...
                                op{k}, distances(ii), regions{j});
    
            %%% Measure fixed effects
            [~, ~, fe_stats] = fixedEffects(glme);
            [~,~, re_stats] = randomEffects(glme);
            % Store to struct
            effects.op{k}.distances(ii).regions{j}.fe_stats = fe_stats;
            effects.op{k}.distances(ii).regions{j}.re_stats = re_stats;
            
            %% Normal distribution test
            %{
            % Extract data at distance and region
            optical_property = tbl{tbl.distance == distances(ii) &...
                    strcmp(tbl.Region, regions{j}),'OpticalProperty'};   
    
            % Histogram data for visualization
            figure;
            histogram(optical_property);
            title(['Histogram of Optical Property at ' regions{j} ...
                    ' for Distance ' num2str(distances(ii))]);
            xlabel('Optical Property');
            ylabel('Frequency');
            
            %%% Perform the Lilliefors test for normality
            [h, p] = lillietest(tmp);            
            % Define distance string for struct
            dist_str = strcat('dist',string(distances(ii)));
            normal_results.(op{k}).(regions{j}).(dist_str).h = h;
            normal_results.(op{k}).(regions{j}).(dist_str).p = p;
            %}
        end
    end
end
% Save struct to the data_dir
save(fullfile(data_dir,'log_transform_glme_stats.mat'),"effects",'-v7.3');


%% Iterate optical properties and regions and test GLME
% This will roll each distance into one model

% Define optical properites
op = {'scattering','retardance'};
% define regions
regions = {'front','occip'};
%%% Define the GLME formula
% fixed effects = Groups, disease stage, region
% random intercept = subjectID
fml = 'OpticalProperty ~ Groups + Stage + distance + (1 | subjectID)';
%%% Define the numeric predictors
numericPredictors = {'Stage','distance'};

%%% Iterate optical property (scattering or retardance)
for k = 1:2
    % Extract respective table
    tbl = data.(op{k});
    %%% Iterate regions (front or occip)
    for j = 1:length(regions)
        %%% Check the GLME Linearity
        % Create conditions for subsetting
        condition_region = (tbl.Region == regions(j));
    
        % Display the conditions for debugging
        disp(['\nSubsetting with Region: ' char(regions(j))]);
    
        % Combine conditions to create the subset
        tbl_subset = tbl(condition_region, :);          

        %%% Fit the GLME model
        fprintf('\nFitting Model\n')
        glme = fitglme(tbl_subset, fml);
        
        %%% Test Linearity
        % check_glme_linearity(glme, tbl_subset, numericPredictors,...
        %                     0.05, lintest_dir,...
        %                     op{k}, regions{j});

        %%% Measure fixed effects
        [~, ~, fe_stats] = fixedEffects(glme);
        [~,~, re_stats] = randomEffects(glme);
        pause(1)

    end
end

