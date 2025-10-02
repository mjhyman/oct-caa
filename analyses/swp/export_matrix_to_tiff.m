%% Function to export matrix to TIFF
function export_matrix_to_tiff(A, fpath, filename)
    filename = fullfile(fpath, filename);

    if ndims(A) ~= 3
        error('Input must be a 3D matrix.');
    end

    t = Tiff(filename, 'w');

    tagstruct.ImageLength = size(A, 1);
    tagstruct.ImageWidth = size(A, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = size(A, 1);
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;

    for i = 1:size(A, 3)
        slice = uint16(65535 * mat2gray(A(:, :, i)));
        t.setTag(tagstruct);
        t.write(slice);
        if i < size(A, 3)
            t.writeDirectory();
        end
    end

    t.close();
end