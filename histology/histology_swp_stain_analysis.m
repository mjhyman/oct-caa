%% Measure size-weighted proximity (SWP) of stains
% Purpose: plot stain vs. log(SWP) for each section and average
% Overview:
%{
- import each stain's SWP struct
- plot stain vs. SWP (EPVS and vessel) for each stain
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
lhe_fname = 'histo_swp_hmatched_zscore__LHE_voxelwise_radius4mm_ds4_19-Mar-2026.mat';
gfap_fname = 'histo_swp_hmatched_zscore__GFAP_voxelwise_radius4mm_ds4_20-Mar-2026.mat';
cd68_fname = 'histo_swp_hmatched_zscore__CD68_voxelwise_radius4mm_ds4_20-Mar-2026.mat';

%%% Output figure settings
% x-axis binning number of points
nbins = 100;
% Date for file output
t = datetime('today');
% Font size for figures
fsize = 30;

%% Import stains

% % Import LHE
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

%% Create SWP/stain pair of each section
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
        % Extract vessel SWP and EPVS SWP
        ves_swp = stain(j).ves_swp;
        epvs_swp = stain(j).epvs_swp;        
        
        %%% Create matrix of stain vs. log(swp) within mask
        % Retrieve mask, vessels, EPVS
        mask = stain(j).mask;
        ves = stain(j).ves;
        epvs = stain(j).epvs;
        % Set the EPVS and vessels to 0
        mask(ves) = 0;
        mask(epvs) = 0;
        
        %%% Retrieve SWP values within parenchyma mask
        ves_swp = ves_swp(mask);
        epvs_swp = epvs_swp(mask);
               
        %%% Extract the stain pixels that correspond to log_swp
        % Find the stain z-scores corresponding to the masked log_swp
        stain_y = stain(j).z_stain(mask);
        
        %%% Create array to store the heatmap pairs
        % Add heatpair back to struct
        ves_swp_stain = [ves_swp, stain_y];
        epvs_swp_stain = [epvs_swp, stain_y];
        stains.(stain_names{ii})(j).ves_swp_stain = ves_swp_stain;
        stains.(stain_names{ii})(j).epvs_swp_stain = epvs_swp_stain;
            
        %% Debugging figures
        % x-axis limits
        xmin = 0;
        xmax = 25;
        xlims = [xmin,xmax];
        % Name of tissue section
        sname = stains.(stain_names{ii})(j).baseName;
        figure;

        %%% Plot stain vs. EPVS SWP and Vessel SWP
        % Bin to 100 points for scatterplot
        ves_swp_stain = bin_swp(ves_swp_stain, nbins);
        epvs_swp_stain = bin_swp(epvs_swp_stain, nbins);
        
        % Plot stain vs. Vessel SWP for the current section
        subplot(2,1,1);
        scatter(ves_swp_stain(:,1), ves_swp_stain(:,2), 'filled');
        xlim(xlims);
        xlabel('Vessel SWP');
        ylabel('Stain Intensity (z-score)');
        title(sprintf('%s vs. Vessel SWP - %s',stain_names{ii},sname),...
                        'Interpreter','none');
        % Plot stain vs. EPVS SWP for the current section
        subplot(2,1,2);
        scatter(epvs_swp_stain(:,1), epvs_swp_stain(:,2), 'filled');
        xlim(xlims);
        xlabel('EPVS SWP');
        ylabel('Stain Intensity (z-score)');
        title(sprintf('%s vs. EPVS SWP - %s',stain_names{ii},sname),...
                        'Interpreter','none');
        % Save settings
        fname = strcat(stain_names{ii},'_',sname,'_ves_epvs_stain_',string(t));
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
% Struct for iterating between ves_swp_stain & epvs_swp_stain
swp_cell = {'ves_swp_stain','epvs_swp_stain'};

for ii = 1:numel(stain_names)
    sprintf('Averaging subjects for each region for stain %i / %i',ii, numel(stain_names))
    % Iterate vessel or EPVS SWP
    for j = 1:2
        %%% Stain Names
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
        front = {stain(contains(section_names, front_names)).(swp_cell{j})};
        occip = {stain(contains(section_names, occip_names)).(swp_cell{j})};
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
        stain_avg.(stain_names{ii}).(swp_cell{j}).binned_combined = binned_combined;
        stain_avg.(stain_names{ii}).(swp_cell{j}).binned_front = binned_front;
        stain_avg.(stain_names{ii}).(swp_cell{j}).binned_occip = binned_occip;
    end
end

%% Plot average region (frontal and occipital separate)
% Labels and such
ylab = {'Stain Intensity','(z-score)'};
xlims = [0,30];
ylims = [-0.5,0.5];
screen = get(0,'ScreenSize'); % [left bottom width height]
fig_size = [1049, 6, 1034, 983];

for ii = 1:numel(stain_names)
    % Iterate the Vessel / EPVS SWPs
    for j = 1:2
        %%% Retrieve combined stains
        binned_combined = stain_avg.(stain_names{ii}).(swp_cell{j}).binned_combined;
        binned_front = stain_avg.(stain_names{ii}).(swp_cell{j}).binned_front;
        binned_occip = stain_avg.(stain_names{ii}).(swp_cell{j}).binned_occip;
        
        %%% Fit curve to data (combined, frontal, occip)
        [f_comb, stat_comb] = fit(binned_combined(:,1),binned_combined(:,2),'poly9');
        [f_front, stat_front] = fit(binned_front(:,1),binned_front(:,2),'poly9');
        [f_occip, stat_occip] = fit(binned_occip(:,1),binned_occip(:,2),'poly9');
        % Store to struct
        stats.(stain_names{ii}).(swp_cell{j}).binned_combined.fit = f_comb;
        stats.(stain_names{ii}).(swp_cell{j}).binned_combined.stats = stat_comb;
        stats.(stain_names{ii}).(swp_cell{j}).binned_front.fit = f_front;
        stats.(stain_names{ii}).(swp_cell{j}).binned_front.stats = stat_front;
        stats.(stain_names{ii}).(swp_cell{j}).binned_occip.fit = f_occip;
        stats.(stain_names{ii}).(swp_cell{j}).binned_occip.stats = stat_occip;
    end

    %%% Retrieve local data
    % Vessel
    ves_comb = stain_avg.(stain_names{ii}).(swp_cell{1}).binned_combined;
    ves_front = stain_avg.(stain_names{ii}).(swp_cell{1}).binned_front;
    ves_occip = stain_avg.(stain_names{ii}).(swp_cell{1}).binned_occip;
    ves_comb_fit = stats.(stain_names{ii}).(swp_cell{1}).binned_combined.fit;
    ves_front_fit = stats.(stain_names{ii}).(swp_cell{1}).binned_front.fit;
    ves_occip_fit = stats.(stain_names{ii}).(swp_cell{1}).binned_occip.fit;
    % EPVS
    epvs_comb = stain_avg.(stain_names{ii}).(swp_cell{2}).binned_combined;
    epvs_front = stain_avg.(stain_names{ii}).(swp_cell{2}).binned_front;
    epvs_occip = stain_avg.(stain_names{ii}).(swp_cell{2}).binned_occip;
    epvs_comb_fit = stats.(stain_names{ii}).(swp_cell{2}).binned_combined.fit;
    epvs_front_fit = stats.(stain_names{ii}).(swp_cell{2}).binned_front.fit;
    epvs_occip_fit = stats.(stain_names{ii}).(swp_cell{2}).binned_occip.fit;

    %%% Plot Combined
    fig = figure; set(fig,'Position',fig_size);
    % Vessel
    subplot(2,1,1)
    scatter(ves_comb(:,1),ves_comb(:,2),[],'k','filled');
    hold on; p = plot(ves_comb_fit);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title('Stain vs. Vessel SWP')
    xlabel('Vessel SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
    set(gca,'fontsize',fsize)
    % EPVS
    subplot(2,1,2)
    scatter(epvs_comb(:,1),epvs_comb(:,2),[],'k','filled');
    hold on; p = plot(epvs_comb_fit);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title('Stain vs. EPVS SWP')
    xlabel('EPVS SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
    % Super title
    sgtitle([stain_names{ii},' combined']);
    set(gca,'fontsize',fsize)
    pause(1);
    fname = strcat(stain_names{ii},'_combined_all_subs_combined_',string(t));
    saveas(gcf, fullfile(figdir, fname),'jpg'); close

    %%% Plot Front
    fig = figure; set(fig,'Position',fig_size);
    % Vessel
    subplot(2,1,1)
    scatter(ves_front(:,1),ves_front(:,2),[],'k','filled');
    hold on; p = plot(ves_front_fit);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title('Stain vs. Vessel SWP')
    xlabel('Vessel SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
    set(gca,'fontsize',fsize)
    % EPVS
    subplot(2,1,2)
    scatter(epvs_front(:,1),epvs_front(:,2),[],'k','filled');
    hold on; p = plot(epvs_front_fit);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title('Stain vs. EPVS SWP')
    xlabel('EPVS SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
    set(gca,'fontsize',fsize)
    % Super title
    sgtitle([stain_names{ii},' front']);
    set(gca,'fontsize',fsize)
    pause(1);
    fname = strcat(stain_names{ii},'_front_all_subs_combined_',string(t));
    saveas(gcf, fullfile(figdir, fname),'jpg'); close

    %%% Plot Occipital
    fig = figure; set(fig,'Position',fig_size);
    % Vessel
    subplot(2,1,1)
    scatter(ves_occip(:,1),ves_occip(:,2),[],'k','filled');
    hold on; p = plot(ves_occip_fit);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title('Stain vs. Vessel SWP')
    xlabel('Vessel SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
    set(gca,'fontsize',fsize)
    % EPVS
    subplot(2,1,2)
    scatter(epvs_occip(:,1),epvs_occip(:,2),[],'k','filled');
    hold on; p = plot(epvs_occip_fit);
    p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
    title('Stain vs. EPVS SWP')
    xlabel('EPVS SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
    set(gca,'fontsize',fsize)
    % Super title
    sgtitle([stain_names{ii},' occip']);
    set(gca,'fontsize',fsize)
    pause(1);
    fname = strcat(stain_names{ii},'_occip_all_subs_combined_',string(t));
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
        front_ves = [];
        occip_ves = [];
        front_epvs = [];
        occip_epvs = [];
    
        % Check if frontal heatpairs exist
        if ~isempty(front_names)
            front_cells = stain(contains(section_names, front_names));
            % Vessel SWP stain
            if any(~cellfun(@isempty, {front_cells.ves_swp_stain}))
                front_ves = {front_cells(~cellfun(@isempty, {front_cells.ves_swp_stain})).ves_swp_stain};
                front_ves = vertcat(front_ves{:});
            end
            % EPVS SWP stain
            if any(~cellfun(@isempty, {front_cells.epvs_swp_stain}))
                front_epvs = {front_cells(~cellfun(@isempty, {front_cells.epvs_swp_stain})).epvs_swp_stain};
                front_epvs = vertcat(front_epvs{:});
            end
        end
    
        % Check if occipital heatpairs exist
        if ~isempty(occip_names)
            occip_cells = stain(contains(section_names, occip_names));
            % Vessel SWP stain
            if any(~cellfun(@isempty, {occip_cells.ves_swp_stain}))
                occip_ves = {occip_cells(~cellfun(@isempty, {occip_cells.ves_swp_stain})).ves_swp_stain};
                occip_ves = vertcat(occip_ves{:});
            end
            % EPVS SWP stain
            if any(~cellfun(@isempty, {occip_cells.epvs_swp_stain}))
                occip_epvs = {occip_cells(~cellfun(@isempty, {occip_cells.epvs_swp_stain})).epvs_swp_stain};
                occip_epvs = vertcat(occip_epvs{:});
            end
        end
        
        %%% Compute front, occip, combined (as applicable)
        % Initialize binned data containers
        binned_front_ves = [];
        binned_occip_ves = [];
        binned_combined_ves = [];
        binned_front_epvs = [];
        binned_occip_epvs = [];
        binned_combined_epvs = [];
        % Check and process frontal data
        if ~isempty(front_ves)
            %%% Vessel
            % Reshape
            front_ves_x = reshape(front_ves(:,1:2:end), [], 1);
            front_ves_y = reshape(front_ves(:,2:2:end), [], 1);
            % Sort frontal data
            [front_ves_x, sortIdxFront] = sort(front_ves_x);
            front_ves_y = front_ves_y(sortIdxFront);
            % Bin the frontal data
            binned_front_ves = bin_swp(front_ves, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_front_ves = binned_front_ves;
            % Fit line to data
            [f_front_ves, stat_front_ves] = fit(binned_front_ves(:,1),binned_front_ves(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).front.ves.fit = f_front_ves;
            stats.(stain_names{ii}).(stages{j}).front.ves.stats = stat_front_ves;
            %%% EPVS
            % Reshape
            front_epvs_x = reshape(front_epvs(:,1:2:end), [], 1);
            front_epvs_y = reshape(front_epvs(:,2:2:end), [], 1);
            % Sort frontal data
            [front_epvs_x, sortIdxFront] = sort(front_epvs_x);
            front_epvs_y = front_epvs_y(sortIdxFront);
            % Bin the frontal data
            binned_front_epvs = bin_swp(front_epvs, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_front_epvs = binned_front_epvs;
            % Fit line to data
            [f_front_epvs, stat_front_epvs] = fit(binned_front_epvs(:,1),binned_front_epvs(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).front.epvs.fit = f_front_epvs;
            stats.(stain_names{ii}).(stages{j}).front.epvs.stats = stat_front_epvs;
        end
        % Check and process occipital data
        if ~isempty(occip_ves)
            %%% Vessel
            % Reshape
            occip_ves_x = reshape(occip_ves(:,1:2:end), [], 1);
            occip_ves_y = reshape(occip_ves(:,2:2:end), [], 1);
            % Sort frontal data
            [occip_ves_x, sortIdxoccip] = sort(occip_ves_x);
            occip_ves_y = occip_ves_y(sortIdxoccip);
            % Bin the occipal data
            binned_occip_ves = bin_swp(occip_ves, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_occip_ves = binned_occip_ves;
            % Fit line to data
            [f_occip_ves, stat_occip_ves] = fit(binned_occip_ves(:,1),binned_occip_ves(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).occip.ves.fit = f_occip_ves;
            stats.(stain_names{ii}).(stages{j}).occip.ves.stats = stat_occip_ves;
            %%% EPVS
            % Reshape
            occip_epvs_x = reshape(occip_epvs(:,1:2:end), [], 1);
            occip_epvs_y = reshape(occip_epvs(:,2:2:end), [], 1);
            % Sort occipal data
            [occip_epvs_x, sortIdxoccip] = sort(occip_epvs_x);
            occip_epvs_y = occip_epvs_y(sortIdxoccip);
            % Bin the occipal data
            binned_occip_epvs = bin_swp(occip_epvs, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_occip_epvs = binned_occip_epvs;
            % Fit line to data
            [f_occip_epvs, stat_occip_epvs] = fit(binned_occip_epvs(:,1),binned_occip_epvs(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).occip.epvs.fit = f_occip_epvs;
            stats.(stain_names{ii}).(stages{j}).occip.epvs.stats = stat_occip_epvs;
        end
        % Compute combined data only if at least one region has data
        if ~isempty(front_ves) && ~isempty(occip_ves)
            %%% Vessel
            % Create combined matrix and bin
            combined_xy = [front_ves_x, front_ves_y; occip_ves_x, occip_ves_y];
            binned_combined = bin_swp(combined_xy, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_combined_ves = binned_combined;
            % Stats of fit line
            [f_comb_ves, stat_comb_ves] = fit(binned_combined(:,1),binned_combined(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).comb.ves.fit = f_comb_ves;
            stats.(stain_names{ii}).(stages{j}).comb.ves.stats = stat_comb_ves;
            %%% EPVS
            % Create combined matrix and bin
            combined_xy = [front_epvs_x, front_epvs_y; occip_epvs_x, occip_epvs_y];
            binned_combined = bin_swp(combined_xy, nbins);
            stage_avg.(stain_names{ii}).(stages{j}).binned_combined_epvs = binned_combined;
            % Stats of fit line
            [f_comb_epvs, stat_comb_epvs] = fit(binned_combined(:,1),binned_combined(:,2),'poly9');
            stats.(stain_names{ii}).(stages{j}).comb.epvs.fit = f_comb_epvs;
            stats.(stain_names{ii}).(stages{j}).comb.epvs.stats = stat_comb_epvs;
        end
    end
end

%% Plot the average of each disease stage
% x- and y-axis limits
xlims = [0,25];
ylims = [-0.5, 0.5];
stage_avg.(stain_names{ii}).(stages{j}).binned_front_ves
% Iterate over stains
for ii = 1:numel(stain_names)
    % Iterate severities
    for j = 1:numel(stages)
        %%% Frontal
        if isfield(stage_avg.(stain_names{ii}).(stages{j}), 'binned_front_ves')
            fig = figure; set(fig,'Position',fig_size);
            %%% Vessel
            subplot(2,1,1)
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_front_ves;
            f = stats.(stain_names{ii}).(stages{j}).front.ves.fit;
            % Plot for frontal
            scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel('Vessel SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
            set(gca, 'fontsize',fsize)
            %%% EPVS
            subplot(2,1,2)
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_front_epvs;
            f = stats.(stain_names{ii}).(stages{j}).front.epvs.fit;
            % Plot for frontal
            scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel('EPVS SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
            set(gca, 'fontsize',fsize)
            %%% Save
            sgtitle([stain_names{ii},' ',stage_names{j}, ' frontal']);
            fname = strcat(stain_names{ii},'_',stages{j},'_frontal_',string(t));
            pause(1)
            saveas(gcf, fullfile(figdir, fname), 'jpg');
        end

        %%% Occipital
        if isfield(stage_avg.(stain_names{ii}).(stages{j}), 'binned_occip_ves')
            fig = figure; set(fig,'Position',fig_size);
            %%% Vessel
            subplot(2,1,1)
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_occip_ves;
            f = stats.(stain_names{ii}).(stages{j}).occip.ves.fit;
            % Plot for frontal
            scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel('Vessel SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
            set(gca, 'fontsize',fsize)
            %%% EPVS
            subplot(2,1,2)
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_occip_epvs;
            f = stats.(stain_names{ii}).(stages{j}).occip.epvs.fit;
            % Plot for frontal
            scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel('EPVS SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
            set(gca, 'fontsize',fsize)
            %%% Save
            sgtitle([stain_names{ii},' ',stage_names{j}, ' occipital']);
            fname = strcat(stain_names{ii},'_',stages{j},'_occip_',string(t));
            pause(1)
            saveas(gcf, fullfile(figdir, fname), 'jpg');
        end
        
        %%% Combined
        if isfield(stage_avg.(stain_names{ii}).(stages{j}), 'binned_combined_ves')
            fig = figure; set(fig,'Position',fig_size);
            %%% Vessel
            subplot(2,1,1)
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_combined_ves;
            f = stats.(stain_names{ii}).(stages{j}).comb.ves.fit;
            % Plot for frontal
            scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel('Vessel SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
            set(gca, 'fontsize',fsize)
            %%% EPVS
            subplot(2,1,2)
            tmp = stage_avg.(stain_names{ii}).(stages{j}).binned_combined_epvs;
            f = stats.(stain_names{ii}).(stages{j}).comb.epvs.fit;
            % Plot for frontal
            scatter(tmp(:,1),tmp(:,2),[],'k','filled');
            hold on; p = plot(f);
            p.LineWidth = 5; p.Color = 'r'; hold off; legend off;
            % Figure labels etc.
            xlabel('EPVS SWP'); ylabel(ylab); xlim(xlims); %ylim(ylims);
            set(gca, 'fontsize',fsize)
            %%% Save
            sgtitle([stain_names{ii},' ',stage_names{j}, ' combined']);
            fname = strcat(stain_names{ii},'_',stages{j},'_comb_',string(t));
            pause(1)
            saveas(gcf, fullfile(figdir, fname), 'jpg');
        end
    end
end
