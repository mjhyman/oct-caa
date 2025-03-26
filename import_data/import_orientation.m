%% Import the orientation volumes for the CAA dataset.
%{
For each subject:
- import struct
- import orientation
- align orientation to scattering
- save struct to:
    /autofs/cluster/octdata2/users/mjhyman/...
        oct_caa_analyses/optical_properties/[subject]
%}

% TODO: have only saved CAA6
% need to save all the others

%% Prepare environment
clc; close all;
% Add top-level directory
d = pwd;
addpath(fullfile(pwd));
addpath('/autofs/cluster/octdata3/users/mjhyman/oct-caa/freesurfer');
% Output directory for optical properties
fname = 'optical_properties.mat';
fout_base = ['/autofs/cluster/octdata3/users/mjhyman/oct_caa_analyses/' ...
    'optical_properties'];
% Flag for loading the .MAT files
load_mats = false;

%% Import the .MAT files for each subject
if load_mats
    % CAA 6
    fprintf('Loading CAA6\n');
    load(fullfile(fout_base,'caa6/caa6.mat'));
    fprintf('Finished loading CAA 6\n');
    % CAA 17
    fprintf('Loading CAA17\n');
    load(fullfile(fout_base,'caa17/occip/caa17.mat'));
    fprintf('Finished loading CAA17 \n');
    % CAA 22
    fprintf('Loading CAA22\n');
    load(fullfile(fout_base,'caa22/caa22.mat'));
    fprintf('Finished loading CAA17 \n');
    % CAA 25
    fprintf('Loading CAA25\n');
    load(fullfile(fout_base,'caa25/caa25.mat'));
    fprintf('Finished loading CAA17 \n');
    % CAA 26
    fprintf('Loading CAA26\n');
    load(fullfile(fout_base,'caa26/caa26.mat'));
    fprintf('Finished loading CAA17 \n');
end
%% Parameters for scattering and retardance
% MRIread parameters
% header only flag
h_flag = 0;
% permute x,y dimensions flag
p_flag = 0;

%% CAA 6 Front:
% import orientation
orient = fullfile(fout_base,'caa6/front/orientation_reg2mus.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);
% Add orientation to .MAT
caa6.front.orient = orient;

%% CAA 6 Occip:
% import orientation
orient = fullfile(fout_base,'caa6/occip/orientation_full_stack.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);

% Resize orientation
mus = caa6.occip.mus;
orient = imresize3(orient,[size(mus,1),size(mus,2),size(mus,3)]);

% Add to data struct and save
caa6.occip.orient = orient;
save(fullfile(fout_base,'/caa6/caa6.mat'),'caa6','-v7.3');

%% CAA 17 occip
% Import orientation
orient = fullfile(fout_base,'caa17/occip/orientation_full_stack.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);

% Resize orientation
mus = caa17.occip.mus;
orient = imresize3(orient,[size(mus,1),size(mus,2),size(mus,3)]);

% Debugging
ret = caa17.occip.ret_full;
mask = caa17.occip.mask;
figure; imagesc(mus(:,:,1)); title('mus - 1');
figure; imagesc(orient(:,:,1)); title('orient - 1');
figure; imagesc(ret(:,:,1)); title('ret - 1');
figure; imagesc(mask(:,:,1)); title('mask - 1');

% Add to data struct and save
caa17.occip.orient = orient;
save(fullfile(fout_base,'/caa17/occip/caa17.mat'),'caa17','-v7.3')

%% CAA 22 front
% import orientation
orient = fullfile(fout_base,'caa22/front/orientation_full_stack.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);

% Resize orientation
mus = caa22.front.mus;
orient = imresize3(orient,[size(mus,1),size(mus,2),size(mus,3)]);

% Debugging
ret = caa22.front.ret_full;
mask = caa22.front.mask;
figure; imagesc(mus(:,:,1)); title('mus - 1'); clim([0,50])
figure; imagesc(orient(:,:,1)); title('orient - 1');
figure; imagesc(ret(:,:,1)); title('ret - 1');
figure; imagesc(mask(:,:,1)); title('mask - 1');

% Add to data struct and save
caa22.front.orient = orient;

%% CAA 22 occip:
% Import orientation
orient = fullfile(fout_base,'caa22/occip/orientation_full_stack.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);

% Add to data struct and save
mus = caa22.occip.mus;
orient = imresize3(orient,[size(mus,1),size(mus,2),size(mus,3)]);

% Debugging
ret = caa22.occip.ret_full;
mask = caa22.occip.mask;
figure; imagesc(mus(:,:,1)); title('mus - 1'); %clim([0,10])
figure; imagesc(orient(:,:,1)); title('orient - 1');
figure; imagesc(ret(:,:,1)); title('ret - 1');
figure; imagesc(mask(:,:,1)); title('mask - 1');

% Add to data struct
caa22.occip.orient = orient;

% Save struct
save(fullfile(fout_base,'/caa22/caa22.mat'),'caa22','-v7.3');

%% CAA 25 front:
% Import orientation
orient = fullfile(fout_base,'caa25/front/orientation_full_stack.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);

% Resize orient
mus = caa25.front.mus;
orient = imresize3(orient,[size(mus,1),size(mus,2),size(mus,3)]);

% Debugging
ret = caa25.front.ret_full;
mask = caa25.front.mask;
figure; imagesc(mus(:,:,1)); title('mus - 1'); clim([0,10])
figure; imagesc(orient(:,:,1)); title('orient - 1');
figure; imagesc(ret(:,:,1)); title('ret - 1');
figure; imagesc(mask(:,:,1)); title('mask - 1');

% Add to data struct and save
caa25.front.orient = orient;

%%% CAA 25 occip
% Import orientation
orient = fullfile(fout_base,'caa25/occip/orientation_full_stack.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);

% Resize orientation
mus = caa25.occip.mus;
orient = imresize3(orient,[size(mus,1),size(mus,2),size(mus,3)]);

% Debugging
ret = caa25.occip.ret_full;
mask = caa25.occip.mask;
figure; imagesc(mus(:,:,1)); title('mus - 1'); %clim([0,10])
figure; imagesc(orient(:,:,1)); title('orient - 1');
figure; imagesc(ret(:,:,1)); title('ret - 1');
figure; imagesc(mask(:,:,1)); title('mask - 1');

% Add to data struct
caa25.occip.orient = orient;

% Save struct
save(fullfile(fout_base,'/caa25/caa25.mat'),'caa25','-v7.3');

%%% CAA 26 front:
% Import orientation
orient = fullfile(fout_base,'caa26/front/orientation_full_stack.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);

% Resize orientation
mus = caa26.front.mus;
orient = imresize3(orient,[size(mus,1),size(mus,2),size(mus,3)]);

% Debugging
mus = caa26.front.mus;
ret = caa26.front.ret_full;
mask = caa26.front.mask;
figure; imagesc(mus(:,:,1)); title('mus - 1'); clim([0,10])
figure; imagesc(orient(:,:,1)); title('orient - 1');
figure; imagesc(ret(:,:,1)); title('ret - 1');
figure; imagesc(mask(:,:,1)); title('mask - 1');

% Add to data struct and save
caa26.front.orient = orient;

%%% CAA 26 occip:
% Import orientation
orient = fullfile(fout_base,'caa26/occip/orientation_full_stack.nii');
orient = MRIread(orient,h_flag,p_flag);
orient = single(orient.vol);

% Resize orientation
mus = caa26.occip.mus;
orient = imresize3(orient,[size(mus,1),size(mus,2),size(mus,3)]);

% Debugging
ret = caa26.occip.ret_full;
mask = caa26.occip.mask;
figure; imagesc(mus(:,:,1)); title('mus - 1'); %clim([0,10])
figure; imagesc(orient(:,:,1)); title('orient - 1');
figure; imagesc(ret(:,:,1)); title('ret - 1');
figure; imagesc(mask(:,:,1)); title('mask - 1');

% Add to data struct and save
caa26.occip.orient = orient;

% Save struct
save(fullfile(fout_base,'/caa26/caa26.mat'),'caa26','-v7.3');