function img = std2img(std, hdr)

% convert standard space coordinates to 1-based indices into image matrix
% using sform matrix from NIFTI header

if hdr.sform_code > 0
    std = [std(:); 1];
    m = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
    img = inv(m) * std;
    img = img(1:3) + 1;
else
    error('std2img only works for NIFTI images with sform_code > 0');
end
