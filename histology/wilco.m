function [stats] = wilco(exp,ctl,tail)
% WILCO: perform Wilcoxon test
%   INPUTS
%       exp (vector): experimental values (EPVS)
%       ctl (vector): control values (vessel or parenchyma)
%       tail (str): the type of test (both, left, right)
%   OUTPUTS
%       w (float): wilcoxon signed rank test statistic
%       z (float): z-statistics
%       p (float): p-value
%       h (int): result of hypothesis test

% The approximate method is typically used for large samples (>15). It is
% used here just to compute the Z-score
[~,~,stats_signrank] = signrank(exp,ctl,'method','approximate',...
    'tail',tail);
w = stats_signrank.signedrank;
z = stats_signrank.zval;
% The exact method is used to compute the exact p-value since there are few
% number of samples
[p,h,~] = signrank(exp,ctl,'tail',tail);

%%% Add to struct
stats.p = p;
stats.h = h;
stats.w = w;
stats.z = z;
stats.exp_med = median(exp);
stats.ctl_med = median(ctl);
% Calculate effect size
d = exp - ctl;
n = sum(d ~= 0);
stats.r = z ./ sqrt(n);
stats.n = n;
end