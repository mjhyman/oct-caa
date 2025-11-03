function subject_scatter_op_vs_dist(x, parench, sub,...
                                    mus_yl, ret_yl, ori_yl,...
                                    mus_yt, ret_yt, ori_yt,...
                                    xt, dir_out)

%% Scatterplot of optical property vs. distance (disjoint ring)
% y-axis = average optical property  for epvs or vessel
% x-axis = distance of ring from edge of epvs or vessel
% repeat this for all 3 properties (mus, ret, ori)
% This will only use the "outer" ring, which is disjoint from the annotated
% vessel/epvs

err_flag = false;


%%% Initialization
% Retrieve radii within struct
radii = fieldnames(parench.(sub).front);

%%% Combined array for pairs of optical proprety vs. distance
comb_mean_ves_mus = zeros(length(radii),1);
comb_mean_ves_ret = zeros(length(radii),1);
comb_mean_ves_ori = zeros(length(radii),1);
comb_mean_epvs_mus = zeros(length(radii),1);
comb_mean_epvs_ret = zeros(length(radii),1);
comb_mean_epvs_ori = zeros(length(radii),1);

%%% Combined array for pairs of optical proprety vs. distance
front_mean_ves_mus = zeros(length(radii),1);
front_mean_ves_ret = zeros(length(radii),1);
front_mean_ves_ori = zeros(length(radii),1);
front_mean_epvs_mus = zeros(length(radii),1);
front_mean_epvs_ret = zeros(length(radii),1);
front_mean_epvs_ori = zeros(length(radii),1);

%%% Combined array for pairs of optical proprety vs. distance
occip_mean_ves_mus = zeros(length(radii),1);
occip_mean_ves_ret = zeros(length(radii),1);
occip_mean_ves_ori = zeros(length(radii),1);
occip_mean_epvs_mus = zeros(length(radii),1);
occip_mean_epvs_ret = zeros(length(radii),1);
occip_mean_epvs_ori = zeros(length(radii),1);

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

    %%% Iterate over the regions
    regs = fields(parench.(sub));
    for k = 1:length(regs)
        % Retrieve vessel measurements (mus, ret, ori)
        reg = regs{k};
        ves = parench.(sub).(reg).(rname).outter.ves;
        
        %%% Add data to combined vectors
        % Take average of mus, ret, ori
        comb_ves_mus = [comb_ves_mus, omit_outlier(ves.pmus)];
        comb_ves_ret = [comb_ves_ret, omit_outlier(ves.pret)];
        comb_ves_ori = [comb_ves_ori, omit_outlier(ves.pori)];
        %%% Add to frontal or occipital
        if strcmp(reg,'front')
            front_ves_mus = [front_ves_mus, omit_outlier(ves.pmus)];
            front_ves_ret = [front_ves_ret, omit_outlier(ves.pret)];
            front_ves_ori = [front_ves_ori, omit_outlier(ves.pori)];
        elseif strcmp(reg,'occip')
            occip_ves_mus = [occip_ves_mus, omit_outlier(ves.pmus)];
            occip_ves_ret = [occip_ves_ret, omit_outlier(ves.pret)];
            occip_ves_ori = [occip_ves_ori, omit_outlier(ves.pori)];
        end
    
        %%% retrieve EPVS measurements (if exist)
        if isfield(parench.(sub).(reg).(rname).outter, 'epvs')
            % Create local
            epvs = parench.(sub).(reg).(rname).outter.epvs;
            %%% Add to combined
            comb_epvs_mus = [comb_epvs_mus, omit_outlier(epvs.pmus)];
            comb_epvs_ret = [comb_epvs_ret, omit_outlier(epvs.pret)];
            comb_epvs_ori = [comb_epvs_ori, omit_outlier(epvs.pori)];
            %%% Add to Frontl or occipital
            if strcmp(reg,'front')
                front_epvs_mus = [front_epvs_mus, omit_outlier(epvs.pmus)];
                front_epvs_ret = [front_epvs_ret, omit_outlier(epvs.pret)];
                front_epvs_ori = [front_epvs_ori, omit_outlier(epvs.pori)];
            elseif strcmp(reg,'occip')
                occip_epvs_mus = [occip_epvs_mus, omit_outlier(epvs.pmus)];
                occip_epvs_ret = [occip_epvs_ret, omit_outlier(epvs.pret)];
                occip_epvs_ori = [occip_epvs_ori, omit_outlier(epvs.pori)];
            end
        end
    end

    %% Take average across subjects at distance
    %%% Combined - take average across subjects
    % Add average optical property to main array
    comb_mean_ves_mus(ii) = mean(comb_ves_mus,'omitnan');
    comb_mean_ves_ret(ii) = mean(comb_ves_ret,'omitnan');
    comb_mean_ves_ori(ii) = real(mean(comb_ves_ori,'omitnan'));
    comb_mean_epvs_mus(ii) = mean(comb_epvs_mus,'omitnan');
    comb_mean_epvs_ret(ii) = mean(comb_epvs_ret,'omitnan');
    comb_mean_epvs_ori(ii) = real(mean(comb_epvs_ori,'omitnan'));

    %%% Front - take average across subjects
    % Add average optical property to main array
    front_mean_ves_mus(ii) = mean(front_ves_mus,'omitnan');
    front_mean_ves_ret(ii) = mean(front_ves_ret,'omitnan');
    front_mean_ves_ori(ii) = real(mean(front_ves_ori,'omitnan'));
    front_mean_epvs_mus(ii) = mean(front_epvs_mus,'omitnan');
    front_mean_epvs_ret(ii) = mean(front_epvs_ret,'omitnan');
    front_mean_epvs_ori(ii) = real(mean(front_epvs_ori,'omitnan'));

    %%% Occip - take average across subjects
    % Add average optical property to main array
    occip_mean_ves_mus(ii) = mean(occip_ves_mus,'omitnan');
    occip_mean_ves_ret(ii) = mean(occip_ves_ret,'omitnan');
    occip_mean_ves_ori(ii) = real(mean(occip_ves_ori,'omitnan'));
    occip_mean_epvs_mus(ii) = mean(occip_epvs_mus,'omitnan');
    occip_mean_epvs_ret(ii) = mean(occip_epvs_ret,'omitnan');
    occip_mean_epvs_ori(ii) = real(mean(occip_epvs_ori,'omitnan'));
end

%% Scatterplots
%%% COMBINED scatterplots (EPVS and ves)
substr = strcat('_',sub);
% Mus
xlab = 'Distance (\mum)';
ylab = '\mus (cm^-^1)';
tit = 'Combined: \mus vs. Distance';
fname = strcat('COMBINED_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_mus, [], comb_mean_epvs_mus, [], ...
                    err_flag,...
                    xlab, ylab, mus_yl, xt, mus_yt,...
                    tit, dir_out, fname)
% Retardance
ylab = 'Retardance (degrees)';
tit = 'Combined: Retardance vs. Distance';
fname = strcat('COMBINED_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_ret, [],...
                    comb_mean_epvs_ret, [],...
                    err_flag,...
                    xlab, ylab, ret_yl, xt, ret_yt,...
                    tit, dir_out, fname, err_flag)
% Orientation
ylab = '\sigma_{orientation} (radians)';
tit = 'Combined: CSDO vs. Distance';
fname = strcat('COMBINED_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, comb_mean_ves_ori, [],...
                    comb_mean_epvs_ori, [],...
                    err_flag,...
                    xlab, ylab, ori_yl, xt, ori_yt,...
                    tit, dir_out, fname, err_flag)

%%% FRONTAL scatterplots (EPVS and ves)
% Mus
ylab = '\mus (cm^-^1)';
tit = 'Frontal: \mus vs. Distance';
fname = strcat('FRONT_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_mus, [],...
                    front_mean_epvs_mus, [],...
                    err_flag,...
                    xlab, ylab, mus_yl, xt, mus_yt,...
                    tit, dir_out, fname, err_flag)
% Retardance
ylab = 'Retardance (degrees)';
tit = 'Frontal: Retardance vs. Distance';
fname = strcat('FRONT_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_ret, [],...
                    front_mean_epvs_ret, [],...
                    err_flag,...
                    xlab, ylab, ret_yl, xt, ret_yt,...
                    tit, dir_out, fname, err_flag)
% Orientation
ylab = '\sigma_{orientation} (radians)';
tit = 'Frontal: CSDO vs. Distance';
fname = strcat('FRONT_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, front_mean_ves_ori, [],...
                    front_mean_epvs_ori, [],...
                    err_flag,...
                    xlab, ylab, ori_yl, xt, ori_yt,...
                    tit, dir_out, fname, err_flag)

%%% OCCIPITAL scatterplots (EPVS and ves)
% Mus
ylab = '\mus (cm^-^1)';
tit = 'Occipital: \mus vs. Distance';
fname = strcat('OCCIP_epvs_ves_mus_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_mus, [],...
                    occip_mean_epvs_mus, [],...
                    err_flag,...
                    xlab, ylab, mus_yl, xt, mus_yt,...
                    tit, dir_out, fname, err_flag)
% Retardance
xlab = 'Distance (\mum)';
ylab = 'Retardance (degrees)';
tit = 'Occipital: Retardance vs. Distance';
fname = strcat('OCCIP_epvs_ves_ret_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_ret, [],...
                    occip_mean_epvs_ret, [],...
                    err_flag,...
                    xlab, ylab, ret_yl, xt, ret_yt,...
                    tit, dir_out, fname, err_flag)
% Orientation
ylab = '\sigma_{orientation} (radians)';
tit = 'Occipital: CSDO vs. Distance';
fname = strcat('OCCIP_epvs_ves_ori_vs_distance',substr,'.png');
scatter_op_vs_dist(x, occip_mean_ves_ori, [],...
                    occip_mean_epvs_ori, [],...
                    err_flag,...
                    xlab, ylab, ori_yl, xt, ori_yt,...
                    tit, dir_out, fname, err_flag)

end