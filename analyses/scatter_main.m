%% Scatterplot function for both within and across subjects
function scatter_main(subjects, fprops, offsets)
% Parse the parench struct, create arrays, scatterplot
% INPUTS
%   subjects (cell array): subjects to iterate over. set this to a single
%                           subject string to plot just one subject
%   fprops (struct): contains limits, x-axis, x-axis ticks, output folder,
%               radii to iterate, flag for plotting error bars
%   offsets (vector): median mus and retardance for each subject and
%                   region. This will be subtracted from each parenchyma
%                   measurement, which will make the measurements be the
%                   relative offset w.r.t. to the baseline for that tissue

%% Extract local variables from struct
% y-axis limits
mus_yl = fprops.mus_yl;
mus_yt = fprops.mus_yt;
ret_yl = fprops.ret_yl;
ret_yt = fprops.ret_yt;
ori_yl = fprops.ori_yl;
ori_yt = fprops.ori_yt;
% x-axis and tick marks
x = fprops.x;
xt = fprops.xt;
% Scatterplot point size
psize = fprops.psize;
% output folder & substring for filename
scat_out = fprops.scat_out;
subdir = fprops.subdir;
substr = fprops.substr;
% parenchyma structure stuff
parench = fprops.parench;
radii = fprops.radii;
% flag for including error bars
err_flag = fprops.err_flag;
% Upper limit for mus and retardance
mus_max = fprops.mus_max;
ret_max = fprops.ret_max;
rm_outlier = fprops.rm_outlier;


%% Initialize vectors: singel value for each distance
%%% Combined array for pairs of optical proprety vs. distance
% Mean
comb_mean_ves_mus = zeros(length(radii),1);
comb_mean_ves_ret = zeros(length(radii),1);
comb_mean_ves_ori = zeros(length(radii),1);
comb_mean_epvs_mus = zeros(length(radii),1);
comb_mean_epvs_ret = zeros(length(radii),1);
comb_mean_epvs_ori = zeros(length(radii),1);
% standard error of measurement (SEM)
comb_sem_ves_mus = zeros(length(radii),1);
comb_sem_ves_ret = zeros(length(radii),1);
comb_sem_ves_ori = zeros(length(radii),1);
comb_sem_epvs_mus = zeros(length(radii),1);
comb_sem_epvs_ret = zeros(length(radii),1);
comb_sem_epvs_ori = zeros(length(radii),1);

%%% Combined array for pairs of optical proprety vs. distance
% Mean
front_mean_ves_mus = zeros(length(radii),1);
front_mean_ves_ret = zeros(length(radii),1);
front_mean_ves_ori = zeros(length(radii),1);
front_mean_epvs_mus = zeros(length(radii),1);
front_mean_epvs_ret = zeros(length(radii),1);
front_mean_epvs_ori = zeros(length(radii),1);
% standard error of measurement (SEM)
front_sem_ves_mus = zeros(length(radii),1);
front_sem_ves_ret = zeros(length(radii),1);
front_sem_ves_ori = zeros(length(radii),1);
front_sem_epvs_mus = zeros(length(radii),1);
front_sem_epvs_ret = zeros(length(radii),1);
front_sem_epvs_ori = zeros(length(radii),1);

%%% Combined array for pairs of optical proprety vs. distance
% Mean
occip_mean_ves_mus = zeros(length(radii),1);
occip_mean_ves_ret = zeros(length(radii),1);
occip_mean_ves_ori = zeros(length(radii),1);
occip_mean_epvs_mus = zeros(length(radii),1);
occip_mean_epvs_ret = zeros(length(radii),1);
occip_mean_epvs_ori = zeros(length(radii),1);
% standard error of measurement (SEM)
occip_sem_ves_mus = zeros(length(radii),1);
occip_sem_ves_ret = zeros(length(radii),1);
occip_sem_ves_ori = zeros(length(radii),1);
occip_sem_epvs_mus = zeros(length(radii),1);
occip_sem_epvs_ret = zeros(length(radii),1);
occip_sem_epvs_ori = zeros(length(radii),1);

%% Iterate radii, subjects, optical property
% iterate distance from edge of epvs/vessel
for ii = 1:length(radii)
    % local radii name
    rname = string(radii(ii));
    %%% initialize arrays for combined
    comb_ves_mus = [];
    comb_ves_ret = [];
    comb_ves_ori = [];
    comb_epvs_mus = [];
    comb_epvs_ret = [];
    comb_epvs_ori = [];
    %%% initialize arrays for frontal
    front_ves_mus = [];
    front_ves_ret = [];
    front_ves_ori = [];
    front_epvs_mus = [];
    front_epvs_ret = [];
    front_epvs_ori = [];
    %%% initialize arrays for occiptial
    occip_ves_mus = [];
    occip_ves_ret = [];
    occip_ves_ori = [];
    occip_epvs_mus = [];
    occip_epvs_ret = [];
    occip_epvs_ori = [];

    %% iterate subjects
    for j = 1:length(subjects)
        % local subject ID
        sub = subjects{j};
        % retrieve regions for this subject
        regions = fieldnames(parench.(sub));
        
        %% iterate regions
        for k = 1:length(regions)
            % region for this iteration
            reg = regions{k};
            % Retrieve vessel measurements (mus, ret, ori)
            ves = parench.(sub).(reg).(rname).outter.ves;

            %%% Apply upper limit threshold & 1.5*IQR outlier removal
            if rm_outlier
                mus = omit_outlier(ves.pmus, mus_max);
                ret = omit_outlier(ves.pret, ret_max);
            else
                mus = ves.pmus;
                ret = ves.pret;
            end
            ori = ves.pori;

            %%% Subtract median - standardize measurements
            if ~isempty(offsets)
                % Subtract the median optical property to standardize
                med_mus = offsets.(sub).(reg).med.mus;
                med_ret = offsets.(sub).(reg).med.ret;
                mus = mus - med_mus;
                ret = ret - med_ret;
            end

            %%% Add data to combined vectors
            comb_ves_mus = [comb_ves_mus, mus];
            comb_ves_ret = [comb_ves_ret, ret];
            comb_ves_ori = [comb_ves_ori, ori];
            %%% Add to frontal or occipital
            if strcmp(reg,'front')
                front_ves_mus = [front_ves_mus, mus];
                front_ves_ret = [front_ves_ret, ret];
                front_ves_ori = [front_ves_ori, ori];
            elseif strcmp(reg,'occip')
                occip_ves_mus = [occip_ves_mus, mus];
                occip_ves_ret = [occip_ves_ret, ret];
                occip_ves_ori = [occip_ves_ori, ori];
            end

            %%% retrieve EPVS measurements (if exist)
            if isfield(parench.(sub).(reg).(rname).outter, 'epvs')
                % Create local
                epvs = parench.(sub).(reg).(rname).outter.epvs;
                
                %%% Apply upper limit threshold & 1.5*IQR outlier removal
                if rm_outlier
                    mus = omit_outlier(epvs.pmus, mus_max);
                    ret = omit_outlier(epvs.pret, ret_max);
                else
                    mus = epvs.pmus;
                    ret = epvs.pret;
                end
                ori = epvs.pori;

                %%% Subtract median - standardize measurements
                if ~isempty(offsets)
                    % Subtract the median optical property to standardize
                    med_mus = offsets.(sub).(reg).med.mus;
                    med_ret = offsets.(sub).(reg).med.ret;
                    mus = mus - med_mus;
                    ret = ret - med_ret;
                end

                %%% Add to combined
                comb_epvs_mus = [comb_epvs_mus, mus];
                comb_epvs_ret = [comb_epvs_ret, ret];
                comb_epvs_ori = [comb_epvs_ori, ori];
                %%% Add to Frontl or occipital
                if strcmp(reg,'front')
                    front_epvs_mus = [front_epvs_mus, mus];
                    front_epvs_ret = [front_epvs_ret, ret];
                    front_epvs_ori = [front_epvs_ori, ori];
                elseif strcmp(reg,'occip')
                    occip_epvs_mus = [occip_epvs_mus, mus];
                    occip_epvs_ret = [occip_epvs_ret, ret];
                    occip_epvs_ori = [occip_epvs_ori, ori];
                end
            end
        end
    end

    %% Take average across subjects at distance
    
    %%% Combined - take average across subjects
    % Add average optical property to main array
    comb_mean_ves_mus(ii) = mean((comb_ves_mus));
    comb_mean_ves_ret(ii) = mean((comb_ves_ret));
    comb_mean_ves_ori(ii) = real(mean(comb_ves_ori));
    comb_mean_epvs_mus(ii) = mean((comb_epvs_mus));
    comb_mean_epvs_ret(ii) = mean((comb_epvs_ret));
    comb_mean_epvs_ori(ii) = real(mean(comb_epvs_ori));

    %%% Front - take average across subjects
    % Add average optical property to main array
    front_mean_ves_mus(ii) = mean((front_ves_mus));
    front_mean_ves_ret(ii) = mean((front_ves_ret));
    front_mean_ves_ori(ii) = real(mean(front_ves_ori));
    front_mean_epvs_mus(ii) = mean((front_epvs_mus));
    front_mean_epvs_ret(ii) = mean((front_epvs_ret));
    front_mean_epvs_ori(ii) = real(mean(front_epvs_ori));

    %%% Occip - take average across subjects
    % Add average optical property to main array
    occip_mean_ves_mus(ii) = mean((occip_ves_mus));
    occip_mean_ves_ret(ii) = mean((occip_ves_ret));
    occip_mean_ves_ori(ii) = real(mean(occip_ves_ori));
    occip_mean_epvs_mus(ii) = mean((occip_epvs_mus));
    occip_mean_epvs_ret(ii) = mean((occip_epvs_ret));
    occip_mean_epvs_ori(ii) = real(mean(occip_epvs_ori));
end

%% Scatterplots of each region (and combined) on the same y-axis limits

%%% COMBINED scatterplots (EPVS and ves)
% Mus
xlab = 'Distance (\mum)';
ylab = '\mus (cm^-^1)';
tit = 'Combined: \mus vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('COMBINED_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_mus, comb_sem_ves_mus,...
                comb_mean_epvs_mus, comb_sem_epvs_mus,err_flag,...
                xlab, ylab, mus_yl, xt, mus_yt, tit, dir_out, fname, psize)
% Retardance
xlab = 'Distance (\mum)';
ylab = 'Retardance (degrees)';
tit = 'Combined: Retardance vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('COMBINED_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_ret, comb_sem_ves_ret,...
                comb_mean_epvs_ret, comb_sem_epvs_ret,err_flag,...
                xlab, ylab, ret_yl, xt, ret_yt, tit, dir_out, fname, psize)
% Orientation
% xlab = 'Distance (\mum)';
% ylab = '\sigma_{orientation} (radians)';
% tit = 'Combined: CSDO vs. Distance';
% dir_out = fullfile(scat_out,subdir);
% fname = strcat('COMBINED_epvs_ves_ori_vs_distance',substr,'.png');
% scatter_op_vs_dist(x, comb_mean_ves_ori, comb_sem_ves_ori,...
%                 comb_mean_epvs_ori, comb_sem_epvs_ori,err_flag,...
%                 xlab, ylab, ori_yl, xt, ori_yt, tit, dir_out, fname, psize)

%%% FRONTAL scatterplots (EPVS and ves)
% Mus
xlab = 'Distance (\mum)';
ylab = '\mus (cm^-^1)';
tit = 'Frontal: \mus vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('FRONT_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_mus, front_sem_ves_mus,...
                front_mean_epvs_mus, front_sem_epvs_mus,err_flag,...
                xlab, ylab, mus_yl, xt, mus_yt,tit, dir_out, fname, psize)
% Retardance
xlab = 'Distance (\mum)';
ylab = 'Retardance (degrees)';
tit = 'Frontal: Retardance vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('FRONT_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_ret, front_sem_ves_ret,...
                front_mean_epvs_ret, front_sem_epvs_ret,err_flag,...
                xlab, ylab, ret_yl,xt, ret_yt,tit, dir_out, fname, psize)
% Orientation
% xlab = 'Distance (\mum)';
% ylab = '\sigma_{orientation} (radians)';
% tit = 'Frontal: CSDO vs. Distance';
% dir_out = fullfile(scat_out,subdir);
% fname = strcat('FRONT_epvs_ves_ori_vs_distance',substr,'.png');
% scatter_op_vs_dist(x, front_mean_ves_ori, front_sem_ves_ori,...
%                 front_mean_epvs_ori, front_sem_epvs_ori,err_flag,...
%                 xlab, ylab, ori_yl,xt, ori_yt,tit, dir_out, fname, psize)

%%% OCCIPITAL scatterplots (EPVS and ves)
% Mus
xlab = 'Distance (\mum)';
ylab = '\mus (cm^-^1)';
tit = 'Occipital: \mus vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('OCCIP_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_mus, occip_sem_ves_mus,...
                occip_mean_epvs_mus, occip_sem_epvs_mus,err_flag,...
                xlab, ylab, mus_yl, xt, mus_yt, tit, dir_out, fname, psize)
% Retardance
xlab = 'Distance (\mum)';
ylab = 'Retardance (degrees)';
tit = 'Occipital: Retardance vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('OCCIP_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_ret, occip_sem_ves_ret,...
                occip_mean_epvs_ret, occip_sem_epvs_ret,err_flag,...
                xlab, ylab, ret_yl,xt, ret_yt,tit,dir_out, fname, psize)
% Orientation
% xlab = 'Distance (\mum)';
% ylab = '\sigma_{orientation} (radians)';
% tit = 'Occipital: CSDO vs. Distance';
% dir_out = fullfile(scat_out,subdir);
% fname = strcat('OCCIP_epvs_ves_ori_vs_distance',substr,'.png');
% scatter_op_vs_dist(x, occip_mean_ves_ori, occip_sem_ves_ori,...
%                 occip_mean_epvs_ori, occip_sem_epvs_ori,err_flag,...
%                 xlab, ylab, ori_yl,xt, ori_yt,tit,dir_out, fname, psize)

%% Print min and max to console

% MUS min/max across vessels and EPVS
mus_min = min([front_mean_ves_mus;occip_mean_ves_mus;...
               front_mean_epvs_mus;occip_mean_epvs_mus;]);
mus_max = max([front_mean_ves_mus;occip_mean_ves_mus;...
               front_mean_epvs_mus;occip_mean_epvs_mus;]);

% RET min/max across vessels and EPVS
ret_min = min([front_mean_ves_ret;occip_mean_ves_ret;...
               front_mean_epvs_ret;occip_mean_epvs_ret;]);
ret_max = max([front_mean_ves_ret;occip_mean_ves_ret;...
               front_mean_epvs_ret;occip_mean_epvs_ret;]);

% Print to console
fprintf('\nMUS: min = %f, max = %f',mus_min,mus_max)
fprintf('\nRET: min = %f, max = %f',ret_min,ret_max)
end

%% Scatterplot function w/o error bars
function scatter_op_vs_dist(x, ves_op, ves_sem, epvs_op, epvs_sem, err_flag,...
                            xlab, ylab, yl, xt, yt, tit, dir_out, fname,...
                            psize)
% Scatter plot with error bars for optical property vs distance
% INPUTS:
%   - x (vector): x-axis positions
%   - ves_op (vector): optical propery for vessels
%   - ves_sem (vector): optical propery standard error for vessels
%   - epvs_op (vector): optical propery for EPVS
%   - epvs_sem (vector): optical propery standard error for EPVS
%   - err_flag (bool): true = errorbar plotting
%   - xlab (str): x-axis label
%   - ylab (str): y-axis label
%   - yl (vector): y-axis limits
%   - xt (array): x-axis ticks
%   - yt (array): y-axis ticks
%   - tit (str): figure title
%   - dir_out (str): output directory
%   - fname (str): output filename
%   - psize (int): point size

% Init figure
fig = figure('Position', [100, 100, 1000, 1000],'Resize', 'off');

% Scatterplot with error bars
if err_flag
    h1 = errorbar(x,ves_op,ves_sem,'k'); hold on;
    h2 = errorbar(x,epvs_op,epvs_sem,'r');
    set(h1,'MarkerSize',psize);
    set(h2,'MarkerSize',psize);
else
    scatter(x,ves_op,psize,'k','filled'); hold on;
    scatter(x,epvs_op,psize,'r','filled');
end

% y-axis limits, ticks, etc.
ylim(yl);
yticks(yt);
ytickformat('%.1f');
ylabel(ylab);

% x-axis label
xlabel(xlab);
xticks(xt);
xtickformat('%.0f');

% plot title
title(tit);

% Font and fontsize
fontname(fig, "Helvetica")
set(gca,'fontsize',40);

% Save output
fout = fullfile(dir_out,fname);
saveas(gcf,fout);
pause(0.5)
close;
end