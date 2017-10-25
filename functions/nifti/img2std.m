function std = img2std(img, hdr)

% convert 1-based indices into image matrix to standard space coordinates
% using sform matrix from NIFTI header

if hdr.sform_code > 0
    img = [img(:) - 1; 1];
    m = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
    std = m * img;
    std = std(1:3);
else
    error('img2std only works for NIFTI images with sform_code > 0');
end
