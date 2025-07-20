%% Statistical Tests for Histology
%{
Each tissue has multiple samples from the parenchyma surrounding the EPVS
(experimental) or at a uniform sampling of the parenchyma (control). The
collaborator recommended including random sampling across the parenchyma,
rather than just around the non-pathological vessels. The average
experimental and average control are measured for each tissue. This was
performed for the following:
- CD68
- GFAP
- LHE

To compare the mean of the EPVS to the mean of the tissue (exp. vs.
control):

** Wilcoxon Signed Rank Test is used:
    Non-parametric statistical test
    Null Hypothesis: (X-Y) comes from distribution w/ median = 0
    Assumptions: X,Y from continuous distribution symmetric about median
    [p,h,stats] = signrank(x,y)

The Gallyas stains were co-registered to the respective OCT retardance
depth in the OCT volumes. The EPVS and vessel annotations are also
co-registered between the stain and the retardance. The concentric circle
(donut) method is used to measure the mean Gallyas and mean retardance
surrounding the EPVS and vessels. To measure the correlation b/w the
Gallyas stain and the retardance, the

** Spearman's Rho Correlation Coefficient
    Returns pairwise correlation coefficient b/w pair of inputs
    [rho, pval] = corr(X,Y,'Type','Spearman','pairwise')


%}

%% Scattering - CD68

% Extract exp + control from all subjects & perform Wilcoxon
[p,h,w,z] = wilcoxon_epvs_ctl(cd68);
% Print the stats
fprintf('CD68: p = %f, h = %d\n',p,h)
fprintf('CD68: exp_med = %f, ctl_med = %f\n',median(epvs),median(ves))
fprintf('CD68: W = %f\n',w);
fprintf('CD68: Z-score = %f\n',z);

%% Scattering - GFAP
% Extract exp + control from all subjects & perform Wilcoxon
[p,h,w,z] = wilcoxon_epvs_ctl(gfap);
% Print the stats
fprintf('GFAP: p = %f, h = %d\n',p,h)
fprintf('GFAP: exp_med = %f, ctl_med = %f\n',median(epvs),median(ves))
fprintf('CD68: W = %f\n',w);
fprintf('CD68: Z-score = %f\n',z);

%% Retardance - LHE (myelin rarefaction surrounding EPVS vs. tissue)
% Extract exp + control from all subjects & perform Wilcoxon
[p,h,w,z] = wilcoxon_epvs_ctl(lhe);
% Print the stats
fprintf('LHE: p = %f, h = %d\n',p,h)
fprintf('LHE: exp_med = %f, ctl_med = %f\n',median(epvs),median(ves))
fprintf('CD68: W = %f\n',w);
fprintf('CD68: Z-score = %f\n',z);

%% Function to extract exp. + control from struct
function [p,h,w,z] = wilcoxon_epvs_ctl(histo)
% Extract the average of the measurements surrounding the EPVS and control
%   INPUTS
%       histo (struct): pathology struct. Each number entry corresponds to
%                       a different tissue section.
%   OUTPUTS:
%       w (float): wilcoxon signed rank test statistic
%       z (float): z-statistics
%       p (float): p-value
%       h (int): result of hypothesis test

% Create vectors for storing experimental and control
nsec = length(histo);
epvs = zeros(nsec,1);
ves = zeros(nsec,1);
% Iterate over all tissue sections
for ii = 1:nsec
    % Open all measurements from current tissue section
    data = histo(ii).meas;
    % Extract first sample surrounding EPVS
    epvs(ii) = data(1).exp_mean;
    ves(ii) = data(1).ctl_mean;
end

%%% Wilcoxon Signed-Rank test
% The approximate method is typically used for large samples (>15). It is
% used here just to compute the Z-score
[~,~,stats] = signrank(epvs,ves,'method','approximate');
w = stats.signedrank;
z = stats.zval;
% The exact method is used to compute the exact p-value since there are few
% number of samples
[p,h,~] = signrank(epvs,ves);

end