%% Retrieve the optical properties
function [mus,ret,sori,n] = retrieve_op(parench,sub,reg,rad,seg)
% RETRIEVE_OP Retrieve mus, retardance, circular Std.Dev. orientation
% INPTUS:
%       parench (struct): structure containing the optical properties of
%                       parenchyma surrounding the EPVS or vessel
%       sub (string): subject ID
%       reg (string): region of brain (occipital vs. frontal)
%       rad (string): outer radius of segmentation surrounding vessel/EPVS
%       seg (string): segmentation ('EPVS' or 'ves')
% OUTPUTS:
%       mus (vector): scattering
%       ret (vector): retardance
%       sori (vector): circular std. dev. of orientation
%       n (vector): number of observations
op = parench.(sub).(reg).(rad).outter.(seg);
mus = rmmissing(op.pmus);
ret = rmmissing(op.pret);
sori = rmmissing(op.pori);
% count number of vessel measurements for subject
n = length(mus);
end