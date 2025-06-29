%% Add optical properties data to arrays
function [array_mus,array_ret,array_ori,array_subid,array_reg,idx] =...
         copy_to_vector(mus, ret, sori,...
                        array_mus, array_ret, array_ori, array_subid,...
                        array_reg, idx, n, tiss_idx, reg)
% COPY_TO_VECTOR move optical properties to respective arrays.
% This function will add the optical data the respective array. 
% INPUTS:
%   array_mus (vector): array of mus data
%   array_ret (vector): array of retardance data
%   array_ori (vector): array of circular std dev of orientation data
%   array_subid (vector): array of subject ID for each observation
%   array_reg (vector): brain region array
%   idx (int): index of location in data vector
%   n (int): number of observations for this tissue sample
%   tiss_idx (int): tissue sample sample index
%   reg (string): brain region (front or occip)
% OUTPUTS:
%   same definition as inputs

% Add vessel measurements to array
array_mus(idx:idx+n-1) = mus;
array_ret(idx:idx+n-1) = ret;
array_ori(idx:idx+n-1) = sori;
% Update subject ID array
array_subid(idx:idx+n-1) = ones(n,1) .* tiss_idx;
% Update region array. 1 = front. 2 = occip
array_reg(idx:idx+n-1,1) = {reg};
% Iterate counter
idx = idx + n;
end