%% Scatterplot function for both within and across subjects
function scatter_main(subjects, fprops)
% Parse the parench struct, create arrays, scatterplot
% INPUTS
%   subjects (cell array): subjects to iterate over. set this to a single
%                           subject string to plot just one subject
%   fprops (struct): contains limits, x-axis, x-axis ticks, output folder,
%               radii to iterate, flag for plotting error bars

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

%%
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

%%% Iterate distance, subjects, optical property
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
            % Apply upper limit threshold & 1.5*IQR outlier removal
            if rm_outlier
                mus = omit_outlier(ves.pmus, mus_max);
                ret = omit_outlier(ves.pret, ret_max);
            else
                mus = ves.pmus;
                ret = ves.pret;
            end
            ori = ves.pori;
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
                % Apply upper limit threshold & 1.5*IQR outlier removal
                if rm_outlier
                    mus = omit_outlier(epvs.pmus, mus_max);
                    ret = omit_outlier(epvs.pret, ret_max);
                else
                    mus = epvs.pmus;
                    ret = epvs.pret;
                end
                ori = epvs.pori;
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
xlab = 'Distance (\mum)';
ylab = '\sigma_{orientation} (radians)';
tit = 'Combined: CSDO vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('COMBINED_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_ori, comb_sem_ves_ori,...
                comb_mean_epvs_ori, comb_sem_epvs_ori,err_flag,...
                xlab, ylab, ori_yl, xt, ori_yt, tit, dir_out, fname, psize)

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
xlab = 'Distance (\mum)';
ylab = '\sigma_{orientation} (radians)';
tit = 'Frontal: CSDO vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('FRONT_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_ori, front_sem_ves_ori,...
                front_mean_epvs_ori, front_sem_epvs_ori,err_flag,...
                xlab, ylab, ori_yl,xt, ori_yt,tit, dir_out, fname, psize)

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
xlab = 'Distance (\mum)';
ylab = '\sigma_{orientation} (radians)';
tit = 'Occipital: CSDO vs. Distance';
dir_out = fullfile(scat_out,subdir);
fname = strcat('OCCIP_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_ori, occip_sem_ves_ori,...
                occip_mean_epvs_ori, occip_sem_epvs_ori,err_flag,...
                xlab, ylab, ori_yl,xt, ori_yt,tit,dir_out, fname, psize)

end