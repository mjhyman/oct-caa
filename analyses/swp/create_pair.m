%% Function to create pairs of SWP metric and optical prop

function [mus_pair, ret_pair] = create_pair(swp,mus,ret)
    % Step 1: Find non-zero indices in swp
    non_zero_indices = find(swp > 0);

    % Step 2: Extract corresponding values from scattering and retardance
    mus_values = mus(non_zero_indices);
    ret_values = ret(non_zero_indices);

    % Step 3: Create pairs (non-zero swp values with scattering or retardance)
    swp_values = swp(non_zero_indices);

    % Combine into pairs
    mus_pair = [swp_values, mus_values];
    ret_pair = [swp_values, ret_values];
end