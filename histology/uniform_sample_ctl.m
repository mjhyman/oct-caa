function ctl = uniform_sample_ctl(stain, mask, epvs, sl, dl)
    % INPUTS:
    % stain: 2D single matrix (histology image)
    % mask: logical matrix (1 = tissue, 0 = background)
    % epvs: logical matrix (1 = EPVS, 0 = not EPVS)
    % sl: side length of square patch
    % dl: distance (stride) between patch centers
    
    % OUTPUT:
    % control_samples: cell array of sampled patches
    
    % Matrix size
    [rows, cols] = size(stain);
    
    % Combined exclusion mask
    valid_region = mask & ~epvs;  % 1 = acceptable sampling, 0 = exclude

    % Generate sampling grid
    [X, Y] = meshgrid(dl:dl:cols, dl:dl:rows);
    
    % Flatten grid points
    X = X(:);
    Y = Y(:);

    % Initialize control measurement output vector
    ny = floor(size(stain,1) / sl);
    nx = floor(size(stain,2) / sl);
    ctl = zeros(nx * ny,1);

    % Half window size
    hw = floor(sl / 2);

    count = 1;

    for idx = 1:length(X)
        cx = X(idx);
        cy = Y(idx);

        % Define patch bounds
        x_min = max(cx - hw, 1);
        x_max = min(cx + hw, cols);
        y_min = max(cy - hw, 1);
        y_max = min(cy + hw, rows);

        % Check if patch is fully within bounds
        if (x_max - x_min + 1 == sl) && (y_max - y_min + 1 == sl)
            % Measure control
            patch = stain(y_min:y_max, x_min:x_max);
            
            % Retain pixels within valid tissue (exclude EPVS/background)
            patch_mask = valid_region(y_min:y_max, x_min:x_max);
            patch(~patch_mask) = NaN;

            % Take average within patch (exclude NaN)
            ctl(count) = mean(patch(:),'omitnan');
            count = count + 1;
        end
    end
    % Only keep the set control measurements
    ctl = ctl(1:count-1);
    fprintf('Total control patches sampled: %d\n', count - 1);
end
