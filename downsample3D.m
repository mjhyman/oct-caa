function mat_ds = downsample3D(mat, ds)
    % downsample3D: Downsamples a 3D matrix by a given factor in all dimensions.
    % 
    % Inputs:
    %   mat - The original 3D matrix to be downsampled.
    %   ds - A scalar or a 3-element vector specifying the factor
    %        by which to downsample in each of the 3 dimensions.
    %        If a scalar is provided, it applies the same factor to all dimensions.
    %
    % Output:
    %   downsampledMatrix - The downsampled 3D matrix.

    % Ensure downsampleFactor is a 3-element vector
    if numel(ds) == 1
        ds = [ds, ds, ds];
    end

    % Get the size of the input matrix
    [szX, szY, szZ] = size(mat);
    
    % Calculate the new size for the downsampled matrix
    newSzX = floor(szX / ds(1));
    newSzY = floor(szY / ds(2));
    newSzZ = floor(szZ / ds(3));

    % Reshape the input matrix into 4D with dimensions corresponding to blocks
    reshapedMatrix = reshape(mat(1:newSzX*ds(1), 1:newSzY*ds(2), 1:newSzZ*ds(3)), ...
                             ds(1), [], ds(2), [], ds(3), []);

    % Compute the mean of each block along the appropriate dimensions
    mat_ds = mean(mean(mean(reshapedMatrix, 1), 3), 5);
end
