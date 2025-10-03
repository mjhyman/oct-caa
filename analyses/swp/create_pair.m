%% Function to create pairs of epvs density and optical prop

function [mus_pair, ret_pair] = create_pair(epvs,mus,ret)

% Step 1: Find non-zero indices in epvs
[rows, cols, slices] = ind2sub(size(epvs), find(epvs > 0));

% Step 2: Extract corresponding values from scattering and retardance
mus_values = mus(sub2ind(size(mus), rows, cols, slices));
ret_values = ret(sub2ind(size(ret), rows, cols, slices));

% Step 3: Create pairs (non-zero epvs values with scattering or retardance)
mus_pair = [epvs(sub2ind(size(epvs), rows, cols, slices)), mus_values];
ret_pair = [epvs(sub2ind(size(epvs), rows, cols, slices)), ret_values];

end