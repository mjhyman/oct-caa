%% LMM for each test
function [stats, p, exp_mus, exp_ret, exp_ori,...
            ctl_mus, ctl_ret, ctl_ori,...
            tbl_mus, tbl_ret, tbl_sori] =...
        lmm_test_diff_subid(test_idx, parench, n_ctl, n_exp)
% LMM_TEST Create linear mixed effects model to compare optical props
% This test will measure across all subjects and compare:
%   - vessels w/ EPVS vs. vessels w/o EPVS 
% INPUTS:
%   test_idx (int): index for the respective statistical test
%   parench (struct): parenchyma optical properties struct
%   n_ves (int): number of vessel measurements
%   n_epvs (int): number of EPVS measurements
% OUTPUTS:
%   stats (struct): contains the FixedEffects output tables for each
%                   optical property
%   data_exp (struct): experimental data arrays for each optical property
%   data_cntrl (struct): control data arrays for each optical property
%   p (struct): p-values for: mus, retardance, mean orientation, std dev of
%                       orienation


% Arrays for control (ctl) and experimental (exp)
ctl_mus = zeros(n_ctl,1);
ctl_ret = zeros(n_ctl,1);
ctl_ori = zeros(n_ctl,1);
ctl_reg = cell(n_ctl,1);
exp_mus = zeros(n_exp,1);
exp_ret = zeros(n_exp,1);
exp_ori = zeros(n_exp,1);
exp_reg = cell(n_exp,1);
% Index to track location in arrays of data (ctl_* and exp_*)
cidx = 1;
eidx = 1;
% Subject ID vectors
ctl_subid = zeros(n_ctl,1);
exp_subid = zeros(n_exp,1);
% Counter for tissue ID
tiss_idx = 1;

% Iterate over subjects in parench
subs = fields(parench);
for ii = 1:length(subs)
    % retrieve subject
    sub = subs{ii};
    % retrieve regions for subject
    regs = fields(parench.(sub));
    % iterate over regions
    for j = 1:length(regs)
        % retrieve region
        reg = regs{j};
        % retrieve optical properties (vessels and epvs)
        segmentations = fields(parench.(sub).(reg).rad2.inner);
        %%% Test1: include all subjects
        if test_idx == 1
            % Iterate over segmentations
            for k=1:length(segmentations)
                % retrieve segmentation
                seg = segmentations{k};
                [mus,ret,sori,n] = retrieve_op(parench,sub,reg,seg);
                % Add to vessel or epvs
                if strcmp(seg,'ves')
                    [ctl_mus,ctl_ret,ctl_ori,ctl_subid,ctl_reg,cidx] =...
                        copy_to_vector(mus, ret, sori,...
                                        ctl_mus, ctl_ret,...
                                        ctl_ori, ctl_subid,...
                                        ctl_reg, cidx, n, tiss_idx,reg);
                    % Iterate tissue sample counter
                    tiss_idx = tiss_idx + 1;
                elseif strcmp(seg,'epvs')
                    [exp_mus,exp_ret,exp_ori,exp_subid,exp_reg,eidx] =...
                        copy_to_vector(mus, ret, sori,...
                                        exp_mus, exp_ret,...
                                        exp_ori, exp_subid,...
                                        exp_reg, eidx, n, tiss_idx,reg);
                    % Iterate tissue sample counter
                    tiss_idx = tiss_idx + 1;
                end
            end
        %%% Test2: only examine subjects w/ EPVS
        elseif test_idx == 2
            if isfield(parench.(sub).(reg),'epvs')
                % Iterate over segmentations
                for k=1:length(segmentations)
                    % retrieve segmentation
                    seg = segmentations{k};
                    [mus,ret,sori,n] = retrieve_op(parench,sub,reg,seg);
                    % Add to vessel or epvs
                    if strcmp(seg,'ves')
                        [ctl_mus,ctl_ret,ctl_ori,ctl_subid,ctl_reg,cidx] =...
                            copy_to_vector(mus, ret, sori,...
                                        ctl_mus, ctl_ret,...
                                        ctl_ori, ctl_subid,...
                                        ctl_reg, cidx, n, tiss_idx,reg);
                        % Iterate tissue sample counter
                        tiss_idx = tiss_idx + 1;
                    elseif strcmp(seg,'epvs')
                        [exp_mus,exp_ret,exp_ori,exp_subid,exp_reg,eidx] =...
                        copy_to_vector(mus, ret, sori,...
                                        exp_mus, exp_ret,...
                                        exp_ori, exp_subid,...
                                        exp_reg, eidx, n, tiss_idx,reg);
                        % Iterate tissue sample counter
                        tiss_idx = tiss_idx + 1;
                    end
                end
            end
        %%% Test3: All subjects. only measure vessels w/o EPVS
        elseif test_idx == 3
            % Retrieve parenchyma scattering and retardance
            seg = 'ves';
            [mus,ret,sori,n] = retrieve_op(parench,sub,reg,seg);
            % Identify cases w/ EPVS and label as experimental
            if isfield(parench.(sub).(reg),'epvs')
                [exp_mus,exp_ret,exp_ori,exp_subid,exp_reg,eidx] =...
                        copy_to_vector(mus, ret, sori,...
                                        exp_mus, exp_ret,...
                                        exp_ori, exp_subid,...
                                        exp_reg, eidx, n, tiss_idx,reg);
                % Iterate tissue sample counter
                tiss_idx = tiss_idx + 1;
            else
                [ctl_mus,ctl_ret,ctl_ori,ctl_subid,ctl_reg,cidx] =...
                        copy_to_vector(mus, ret, sori,...
                                        ctl_mus, ctl_ret,...
                                        ctl_ori, ctl_subid,...
                                        ctl_reg, cidx, n, tiss_idx,reg);
                % Iterate tissue sample counter
                tiss_idx = tiss_idx + 1;
            end
        end
    end
end

%% Create a table for fitting the LME model
% Column 1 = group label
% Column 2 = vascular metric value
% Create the group labels
g_exp = repmat({'experimental'},[length(exp_mus),1]);
g_cnrtl = repmat({'control'},[length(ctl_mus),1]);
g_exp_cnrtl = vertcat(g_exp, g_cnrtl);
% Scattering coefficient array
mus_exp_cnrtl = [exp_mus; ctl_mus];
% Retardance array
ret_exp_cnrtl = [exp_ret; ctl_ret];
% Circular standard deviation of orientation array
ori_exp_cnrtl = [exp_ori; ctl_ori];
% Combine the subject IDs into column vector
subids_exp_cnrtl = vertcat(exp_subid, ctl_subid);
% Combine the brain region indices into column vector
reg_exp_cntrl = vertcat(exp_reg,ctl_reg);
% Table (group labels, subjectID, vascular metric values)
tbl_mus = table(g_exp_cnrtl,reg_exp_cntrl,subids_exp_cnrtl,mus_exp_cnrtl,...
      'VariableNames',{'Groups','Region','subID','OpticalProperty'});
tbl_ret = table(g_exp_cnrtl,reg_exp_cntrl,subids_exp_cnrtl,ret_exp_cnrtl,...
      'VariableNames',{'Groups','Region','subID','OpticalProperty'});
tbl_sori = table(g_exp_cnrtl,reg_exp_cntrl,subids_exp_cnrtl,ori_exp_cnrtl,...
      'VariableNames',{'Groups','Region','subID','OpticalProperty'});
% Take real part of complex numbers
tbl_sori.OpticalProperty = real(tbl_sori.OpticalProperty);

%%% Fit linear mixed-effects model (LME) (table)
% Define the model:
%   response = optical property array
%   random effect (intercept/subject): subID of tissue volume
%   fixed effect: Groups (experimental or control)
fml = 'OpticalProperty ~ Groups + (1 | subID)';
% Fit the model for scattering
lme_mus = fitlme(tbl_mus,fml);
% Fit the model for retardance
lme_ret = fitlme(tbl_ret,fml);
% Fit the model for circular std. dev. of orientation
lme_sori = fitlme(tbl_sori,fml);

%%% Estimates of fixed effects 
[~,~,stats_mus] = fixedEffects(lme_mus);
[~,~,stats_ret] = fixedEffects(lme_ret);
[~,~,stats_sori] = fixedEffects(lme_sori);
% Add stats tables to struct
stats = struct();
stats.mus = stats_mus;
stats.ret = stats_ret;
stats.sori = stats_sori;
% Store the p-value for the stats tests
p = struct();
p.mus = stats_mus{2,6};
p.ret = stats_ret{2,6};
p.sori = stats_sori{2,6};
end