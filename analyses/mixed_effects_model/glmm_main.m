%% Test linear mixed effects model
%{
The purpose of this script is to analyze the optical properties in the
parenchymal tissue surrounding the EPVS and the vasculature. This covers
just the following cases:
- CAA6 (frontal + occip)
- CAA 17 occipital
- CAA 22 (frontal + occip)
- CAA 25 (frontal + occip)
- CAA 26 (frontal + occip)

Outline:
- IMPORT struct containing optical properties (mus + retardance)
- Linearity assumptions for LMM
- LME model
%}

%% Prepare environment
clc; close all;
addpath(genpath('/projectnb/npbssmic/s/mhyman/oct-caa'));

% Directory for loading seg, mus, ret, mask, epvs
data_dir = '/projectnb/npbssmic/ns/CAA/';
fig_out = '/projectnb/npbssmic/ns/CAA/figures/';

%%% 40 um donut
% load the parenchymal optical properties
load(fullfile(data_dir, ...
    'parenchyma_optical_properties_40um_thick_03Nov2025.mat'));
% String indicating EPVS ring radius (in voxels) to access structure
rad = 'rad2';

%% Iterate distances, bootstrap, test for linearity



%% General Linear Mixed Effects Models

%%% Perform LMM Test 1 (same subID for exp/cont within volume)
% Flag for incrementing the subject index just for the tissue volume
flag_subs = true;
[stats1, p1, exp_mus1, exp_ret1, exp_sori1,...
    ctl_mus1, ctl_ret1, ctl_sori1,...
    tbl_mus, tbl_ret, tbl_sori] =...
        lmm_linearity(check_lin, flag_subs, parench, rad, n_ctl, n_exp);
% Export the tables to CSV
fout = fullfile(data_dir,'lmm_test_40um_thick_09Oct2025.xlsx');
writetable(tbl_mus,fout,'Sheet','scattering');
writetable(tbl_ret,fout,'Sheet','retardance');
writetable(tbl_sori,fout,'Sheet','orientation');
