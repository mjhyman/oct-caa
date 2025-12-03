%% LMM for each test
function [stats, p, exp_mus, exp_ret, exp_ori,...
            ctl_mus, ctl_ret, ctl_ori,...
            tbl_mus, tbl_ret, tbl_sori] =...
        lmm_linearity(check_lin, flag_subs, parench, rad, n_ctl, n_exp)
% LMM_TEST Create linear mixed effects model to compare optical props
% This test will measure across all subjects and compare:
%   - vessels w/ EPVS vs. vessels w/o EPVS 
% INPUTS:
%   check_lin (logical): flag for checking linearity assumptions
%   flag_subs (int): flag for how to increment the subject ID assignment
%       0 -> iterate subject ID for each tissue volume and EPVS/vessel.
%            The vessels and EPVS within the same tissue volume will have
%            different subject IDs.
%       1 -> iterate subject ID for each tissue volume. Both the vessel and
%            EPVS within the same volume will have the same subID
%   parench (struct): parenchyma optical properties struct
%   rad (string): string indicating EPVS ring radius for structure
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
        segmentations = fields(parench.(sub).(reg).(rad).outter);
        % Iterate over segmentations
        for k=1:length(segmentations)
            % retrieve segmentation
            seg = segmentations{k};
            [mus,ret,ori,n] = retrieve_op(parench,sub,reg,rad,seg);
            % Add to vessel or epvs
            if strcmp(seg,'ves')
                [ctl_mus,ctl_ret,ctl_ori,ctl_subid,ctl_reg,cidx] =...
                    copy_to_vector(mus, ret, ori,...
                                    ctl_mus, ctl_ret,...
                                    ctl_ori, ctl_subid,...
                                    ctl_reg, cidx, n, tiss_idx, reg);
                % Iterate tissue sample counter
                if ~flag_subs
                    tiss_idx = tiss_idx + 1;
                end
            elseif strcmp(seg,'epvs')
                [exp_mus,exp_ret,exp_ori,exp_subid,exp_reg,eidx] =...
                    copy_to_vector(mus, ret, ori,...
                                    exp_mus, exp_ret,...
                                    exp_ori, exp_subid,...
                                    exp_reg, eidx, n, tiss_idx, reg);
                if ~flag_subs
                    tiss_idx = tiss_idx + 1;
                end
            end
        end
        % Iterate tissue sample counter. Vessel and EPVS will have the same
        % subject index
        if flag_subs
            tiss_idx = tiss_idx + 1;
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

%%% Add the log transformation of the optical property
% mus
log_mus = log(tbl_mus.OpticalProperty);
log_mus(log_mus == -inf) = 0;
% retardance
log_ret = log(tbl_ret.OpticalProperty);
log_ret(log_ret == -inf) = 0;
% orientation
log_ori = log(tbl_sori.OpticalProperty);
log_ori(log_ori== -inf) = 0;
% Add to table
tbl_mus.logOpticalProperty = log_mus;
tbl_ret.logOpticalProperty = log_ret;
tbl_sori.logOpticalProperty = log_ori;

%% Fit general linear mixed-effects model (LME) (table)
% Define the model:
%   response = optical property array
%   random effect (intercept/subject): subID of tissue volume
%   fixed effect: Groups (experimental or control)
fprintf('Fitting LME for each optical property')
fml = 'OpticalProperty ~ Groups*Region + (1 | subID)';
% Define as log transformation
fml = 'logOpticalProperty ~ Groups*Region + (1 | subID)';
% Fit the model for scattering
lme_mus = fitglme(tbl_mus,fml);
% Fit the model for retardance
lme_ret = fitglme(tbl_ret,fml);
% Fit the model for circular std. dev. of orientation
lme_sori = fitglme(tbl_sori,fml);

%%% Check linearity assumptions of GLME
if check_lin
    fprintf('Checking GLME linearity assumptions\n')
    % Mus
    tstr = 'Scattering Coefficient';
    check_glme_linearity(lme_mus, tbl_mus, 0.05, tstr);
    % retardance
    tstr = 'Retardance';
    check_glme_linearity(lme_ret, tbl_ret, 0.05, tstr);
    % orientation
    % tstr = 'Orientation';
    % check_glme_linearity(lme_sori, tbl_sori, 0.05, tstr);
end

%%% Estimates of fixed effects 
fprintf('measuring fixed effects')
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