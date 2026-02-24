%% Measure size-weighted proximity (SWP) of stains
% Purpose: plot stain vs. log(SWP) for each section and average
% Overview:
%{
- import each stain's SWP struct
- take log of each SWP
- plot sain vs. log(SWP) for each stain
- bin for each subject
- average for each subject
%}

%% Top-level settings
% clc; clear; close all;

%%% Directories (SCC)
% Input directory
hdir='/projectnb/npbssmic/ns/CAA/histology/swp';
% Directory to save output figures
figdir = '/projectnb/npbssmic/ns/CAA/histology/figures';

%%% Add top-level directory to path
addpath(genpath('/projectnb/npbssmic/s/mhyman/oct-caa'));

%%% Filenames
lhe_fname = 'histo_swp_hmatched_zscore__LHE_04-Feb-2026.mat';
gfap_fname = 'histo_swp_hmatched_zscore__GFAP__04-Feb-2026.mat';
cd68_fname = 'histo_swp_hmatched_zscore__CD68_04-Feb-2026.mat';

% number of bins along x-axis
nbins = 100;

% Font size
fsize = 30;

%% Import stains

% Import LHE
lhe = load(fullfile(hdir, lhe_fname));
lhe = lhe.lhe;
% Import GFAP
gfap = load(fullfile(hdir, gfap_fname));
gfap = gfap.gfap;
% Import CD68
cd68 = load(fullfile(hdir, cd68_fname));
cd68 = cd68.cd68;

% Add all stains to the struct "stain"
stains = struct();
stains.cd68 = cd68;
stains.gfap = gfap;
stains.lhe = lhe;

%% Take log(SWP) of each section
stain_names = fields(stains);

% Iterate stains
for ii = 1:numel(stain_names)
    % Get the current stain name
    stain = stains.(stain_names{ii});
    % print stain name to console
    fprintf('\nRunning stain %s',stain_names{ii})
    % Iterate sections within stain
    for j = 1:numel(stain)
        % print section # to console
        fprintf('\n\tSection %i / %i\n',j,numel(stain))
        % Take log of SWP
        log_swp = log(stain(j).swp);
        
        %%% Create matrix of stain vs. log(swp) within mask
        % Retrieve mask, vessels, EPVS
        mask = stain(j).mask;
        ves = stain(j).ves;
        epvs = stain(j).epvs;
        % Set the EPVS and vessels to 0
        mask(ves) = 0;
        mask(epvs) = 0;
        
        %%% Remove log_swp values outside of mask bounds
        % Set log_swp values < 1 to be equal to 1
        log_swp(log_swp<1) = 1;
        % Set log(swp) values outside mask to 0
        log_swp(~mask) = 0;
        % Retrieve indices of mask that are parenchyma tissue
        tissueIndices = find(mask);
        log_swp = log_swp(tissueIndices);

        %%% Remove the log_swp values naer 1. Many subjects have 1-2 points
        % that are outliers and cause erroneous results when binning since
        % the next nearest point may be at log(swp) == 6. The binning
        % algorithm will just copy the value at log(swp)==1 up to the next
        % data point.
        % Set threshold for outliers
        swp_min = 1.1;
        % Find indices of log_swp below threshold
        outlierIndices = find(log_swp < swp_min);
        % Remove log_swp below threshold
        log_swp(outlierIndices) = [];
        
        %%% Extract the stain pixels that correspond to log_swp
        % Find the stain z-scores corresponding to the masked log_swp
        stain_y = stain(j).z_stain(tissueIndices);
        % Remove stain pixels with outlier indices
        stain_y(outlierIndices) = [];
        
        %%% Create array to store the heatmap pairs
        % Add heatpair back to struct
        heatpair = [log_swp, stain_y];
        stains.(stain_names{ii})(j).heatpair = heatpair;

        %%% Bin the heatpair
        binned_heatpair = bin_swp(heatpair, nbins);
        stains.(stain_names{ii})(j).binned_heatpair = binned_heatpair;
            
        %% Debugging figures
        % x-axis limits
        xmin = 1;
        xmax = 16;
        xlims = [xmin,xmax];
        % Name of tissue section
        sname = stains.(stain_names{ii})(j).baseName;
        figure;

        %%% Plot stain vs. log(swp)
        % Downsample heatpair to 1000 points for scatterplot
        ds = floor(size(heatpair,1) ./ 1e3);
        heatpair_ds = heatpair(1:ds:end,:);
        
        % Plot stain vs. log(swp) for the current section
        subplot(2,1,1);
        scatter(heatpair_ds(:,1), heatpair_ds(:,2), 'filled');
        xlim(xlims);
        xlabel('log(SWP)');
        ylabel('Stain Intensity (z-score)');
        title(sprintf('%s - %s',stain_names{ii},sname),'Interpreter','none');

        %%% Plot binned stain vs. log(swp)
        subplot(2,1,2);
        scatter(binned_heatpair(:,1),binned_heatpair(:,2),'filled');
        xlim(xlims);
        xlabel('log(SWP)');
        ylabel('Stain Intensity (z-score)');
        title(sprintf('Binned %s - %s', stain_names{ii}, sname),'Interpreter','none');
        fname = strcat(stain_names{ii},'_',sname,'_heatpair_binned_subplots');
        pause(0.5);
        saveas(gcf, fullfile(figdir, fname),'jpg');
        close;

    end
end

%% Average region (frontal and occipital separate)

% Create struct for storing averages
stain_avg = struct();
% Create struct for storing statistics + model of fit
stats = struct();

for ii = 1:numel(stain_names)
    %%% Names
    % Get the current stain name
    stain = stains.(stain_names{ii});
    % Retrieve "baseName" field from section across all sections
    section_names = {stain.baseName}; 
    
    %%% Separate by section (frontal, occipital)
    % Find the names of the sections for frontal (1) and occipital (7)
    front_names = section_names(contains(section_names, '_1_') |...
        endsWith(section_names, '_1'));
    occip_names = section_names(contains(section_names, '_7_') |...
        endsWith(section_names, '_7'));
    % Separate front/occip and combine across subjects
    front = {stain(contains(section_names, front_names)).heatpair};
    occip = {stain(contains(section_names, occip_names)).heatpair};
    % Concatenate each subject's cell into single matrix
    front = vertcat(front{:});
    occip = vertcat(occip{:});

    % Separate x-axis and y-axis for each region
    front_x = reshape(front(:,1:2:end),[],1);
    front_y = reshape(front(:,2:2:end),[],1);
    occip_x = reshape(occip(:,1:2:end),[],1);
    occip_y = reshape(occip(:,2:2:end),[],1);
    
    %%% Sort combined matrix
    % Sort x-axis and apply same sorting indices to y-axis
    [front_x, sortIdxFront] = sort(front_x);
    [occip_x, sortIdxOccip] = sort(occip_x);
    front_y = front_y(sortIdxFront);
    occip_y = occip_y(sortIdxOccip);   

    %%% Bin the matrices
    % Bin the combined subjects along x and y
    combined_xy = [front_x, front_y; occip_x, occip_y];
    binned_combined = bin_swp(combined_xy, nbins);
    % Bin the front and occipital data separately
    binned_front = bin_swp(front, nbins);
    binned_occip = bin_swp(occip, nbins);
    
    % Store binned data
    stain_avg.(stain_names{ii}).binned_combined = binned_combined;
    stain_avg.(stain_names{ii}).binned_front = binned_front;
    stain_avg.(stain_names{ii}).binned_occip = binned_occip;
end

%% Plot average (all subjects combined)
% Labels and such
xlab = 'log(SWP)';
ylab = 'Stain Intensity (z-score)';
xlims = [1,16];
ylims = [-0.5, 0.25];
% Initialize statistics

for ii = 1:numel(stain_names)
    %%% Retrieve
    binned_combined = stain_avg.(stain_names{ii}).binned_combined;
    binned_front = stain_avg.(stain_names{ii}).binned_front;
    binned_occip = stain_avg.(stain_names{ii}).binned_occip;
    
    %%% Fit curve to data (combined, frontal, occip)
    [f_comb, stat_comb] = fit(binned_combined(:,1),binned_combined(:,2),'poly9');
    [f_front, stat_front] = fit(binned_front(:,1),binned_front(:,2),'poly9');
    [f_occip, stat_occip] = fit(binned_occip(:,1),binned_occip(:,2),'poly9');
    % Store to struct
    stats.(stain_names{ii}).binned_combined.fit = f_comb;
    stats.(stain_names{ii}).binned_combined.stats = stat_comb;
    stats.(stain_names{ii}).binned_front.fit = f_front;
    stats.(stain_names{ii}).binned_front.stats = stat_front;
    stats.(stain_names{ii}).binned_occip.fit = f_occip;
    stats.(stain_names{ii}).binned_occip.stats = stat_occip;

    %%% Plots
    % Combined
    figure;
    scatter(binned_combined(:,1),binned_combined(:,2),[],'k','filled');
    hold on; p = plot(f_comb);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title([stain_names{ii},' combined']);
    xlabel(xlab); ylabel(ylab); xlim(xlims); ylim(ylims);
    set(gca,'fontsize',fsize)
    pause(1)
    fname = strcat(stain_names{ii},'_front_occip_all_subs_combined');
    saveas(gcf, fullfile(figdir, fname),'jpg'); close

    % Frontal
    figure; scatter(binned_front(:,1),binned_front(:,2),[],'k','filled');
    hold on; p = plot(f_front);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title([stain_names{ii},' frontal']);
    xlabel(xlab); ylabel(ylab); xlim(xlims); ylim(ylims);   
    set(gca,'fontsize',fsize)
    pause(1)
    fname = strcat(stain_names{ii},'_front_all_subs_combined');
    saveas(gcf, fullfile(figdir, fname),'jpg'); close

    % Occipital
    figure; scatter(binned_occip(:,1),binned_occip(:,2),[],'k','filled');
    hold on; p = plot(f_occip);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title([stain_names{ii},' occipital']);
    xlabel(xlab); ylabel(ylab); xlim(xlims); ylim(ylims);
    set(gca,'fontsize',fsize)
    pause(1)
    fname = strcat(stain_names{ii},'_occip_all_subs_combined');
    saveas(gcf, fullfile(figdir, fname),'jpg'); close
end

%% Average region by severity (keep frontal and occipital separate)
%%% Separate by severity
% Control (CAA26)
% Mild (CAA6)
% Moderate (CAA17 occipital)
% Severe (CAA22 + CAA25)

% Clear stage names
clear ctl_front_names ctl_occip_names mild_front_names mild_occip_names
clear mod_occip_names sev_front_names sev_occip_names

% Stage cell array
stages = {'ctl','mild','mod','sev'};

% Create struct for storing averages
stage_avg = struct();

for ii = 1:numel(stain_names)
    %% Assign stain name to respective region and severity
    % Get current stain and section base names
    stain = stains.(stain_names{ii});
    section_names = {stain.baseName};
    
    % Helper: boolean mask for section index (1 = frontal, 7 = occipital)
    is_front = contains(section_names, '_1_') | endsWith(section_names, '_1');
    is_occip = contains(section_names, '_7_') | endsWith(section_names, '_7');
    
    front_names = section_names(is_front);
    occip_names = section_names(is_occip);
    
    % Helper function: match any of a set of patterns (case-insensitive)
    matchAny = @(names, patterns) ...
        any(cell2mat(cellfun(@(p) contains(names, p, 'IgnoreCase', true), patterns, 'UniformOutput', false)), 2);
    
    % Define disease-stage patterns
    patterns.ctl  = {'CAA26', 'CAA_26'};
    patterns.mild = {'CAA6', 'CAA_6'};
    patterns.mod  = {'CAA17', 'CAA_17'};
    patterns.sev  = {'CAA22', 'CAA_22', 'CAA25', 'CAA_25'};
    
    % Initialize struct
    section_stage_names = struct();
    
    % Populate with programmatic loop to avoid repetition
    stages = fieldnames(patterns);
    for s = 1:numel(stages)
        stage = stages{s};
        patt = patterns.(stage);
        % match against front and occip lists
        front_mask = matchAny(front_names(:), patt);
        occip_mask = matchAny(occip_names(:), patt);
        if any(front_mask)
            section_stage_names.(stage).front = front_names(front_mask);
        end
        if any(occip_mask)
            section_stage_names.(stage).occip = occip_names(occip_mask);
        end
    end
    
    %% Iterate Severity + bin combined, frontal, occipital
    % Cell array of title string for disease stage
    stage_names = {'control','mild','moderate','severe'};

    for j = 1:numel(stages)
        %%% Retrieve the frontal section names for this disease stage
        % Frontal
        try
            front_names = section_stage_names.(stages{j}).front;
        catch
            front_names = [];
        end
        % Occipital
        try
            occip_names = section_stage_names.(stages{j}).occip;
        catch
            occip_names = [];
        end
        
        % Check if both name arrays are empty (stain unavailable)
        if isempty(front_names) && isempty(occip_names)
            continue;  % Skip if both are empty
        end
    
        % Initialize front and occip arrays
        front = [];
        occip = [];
    
        % Check if frontal heatpairs exist
        if ~isempty(front_names)
            front_cells = stain(contains(section_names, front_names));
            if any(~cellfun(@isempty, {front_cells.heatpair}))
                front = {front_cells(~cellfun(@isempty, {front_cells.heatpair})).heatpair};
                front = vertcat(front{:});
            end
        end
    
        % Check if occipital heatpairs exist
        if ~isempty(occip_names)
            occip_cells = stain(contains(section_names, occip_names));
            if any(~cellfun(@isempty, {occip_cells.heatpair}))
                occip = {occip_cells(~cellfun(@isempty, {occip_cells.heatpair})).heatpair};
                occip = vertcat(occip{:});
            end
        end
        
        %%% Compute front, occip, combined (as applicable)
        % Initialize binned data containers
        binned_combined = [];
        binned_front = [];
        binned_occip = [];
        % Check and process frontal data
        if ~isempty(front)
            % Reshape
            front_x = reshape(front(:,1:2:end), [], 1);
            front_y = reshape(front(:,2:2:end), [], 1);
            % Sort frontal data
            [front_x, sortIdxFront] = sort(front_x);
            front_y = front_y(sortIdxFront);
            % Bin the frontal data
            binned_front = bin_swp(front, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_front = binned_front;
            % Fit line to data
            [f_front, stat_front] = fit(binned_front(:,1),binned_front(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).front.fit = f_front;
            stats.(stain_names{ii}).(stages{j}).front.stats = stat_front;
        end
        % Check and process occipital data
        if ~isempty(occip)
            % Reshape
            occip_x = reshape(occip(:,1:2:end), [], 1);
            occip_y = reshape(occip(:,2:2:end), [], 1);
            % Sort occipital data
            [occip_x, sortIdxOccip] = sort(occip_x);
            occip_y = occip_y(sortIdxOccip);
            % Bin the occipital data
            binned_occip = bin_swp(occip, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_occip = binned_occip;
            % Fit line to data
            [f_occip, stat_occip] = fit(binned_occip(:,1),binned_occip(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).occip.fit = f_occip;
            stats.(stain_names{ii}).(stages{j}).occip.stats = stat_occip;
        end
        % Compute combined data only if at least one region has data
        if ~isempty(front) && ~isempty(occip)
            % Create combined matrix and bin
            combined_xy = [front_x, front_y; occip_x, occip_y];
            binned_combined = bin_swp(combined_xy, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_combined = binned_combined;
            % Stats of fit line
            [f_comb, stat_comb] = fit(binned_combined(:,1),binned_combined(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).comb.fit = f_comb;
            stats.(stain_names{ii}).(stages{j}).comb.stats = stat_comb;
        end
    end
end

%% Plot the average of each disease stage

xlims = [1,16];
ylims = [-1.5, 1];

% Iterate over stains
for ii = 1:numel(stain_names)
    % Iterate severities
    for j = 1:numel(stages)
        % Frontal
        if isfield(stage_avg.(stain_names{ii}).(stages{j}), 'binned_front')
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_front;
            f = stats.(stain_names{ii}).(stages{j}).front.fit;
            % Plot for frontal
            figure; scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel(xlab); ylabel(ylab); xlim(xlims); ylim(ylims);   
            title([stain_names{ii},' ',stage_names{j}, ' frontal']);
            set(gca, 'fontsize',fsize)
            fname = strcat(stain_names{ii},'_',stages{j},'_frontal');
            pause(1)
            saveas(gcf, fullfile(figdir, fname), 'jpg');
        end

        % Occipital
        if isfield(stage_avg.(stain_names{ii}).(stages{j}), 'binned_occip')
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_occip;
            f = stats.(stain_names{ii}).(stages{j}).occip.fit;
            % Plot for frontal
            figure; scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel(xlab); ylabel(ylab); xlim(xlims); ylim(ylims);   
            title([stain_names{ii},' ',stage_names{j}, ' occip']);
            set(gca, 'fontsize',fsize)
            fname = strcat(stain_names{ii},'_',stages{j},'_occip');
            pause(1)
            saveas(gcf, fullfile(figdir, fname), 'jpg');
        end
        
        % Combined
        if isfield(stage_avg.(stain_names{ii}).(stages{j}), 'binned_combined')
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_combined;
            f = stats.(stain_names{ii}).(stages{j}).comb.fit;
            % Plot for frontal
            figure; scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel(xlab); ylabel(ylab); xlim(xlims); ylim(ylims);   
            title([stain_names{ii},' ',stage_names{j}, ' combined']);
            set(gca, 'fontsize',fsize)
            fname = strcat(stain_names{ii},'_',stages{j},'_combined');
            pause(1)
            saveas(gcf, fullfile(figdir, fname), 'jpg');
        end
    end
end
