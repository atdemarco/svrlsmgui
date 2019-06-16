function write_nifti(hdr, img, fname, raw)

% WRITE_NIFTI  Write a NIFTI image
%
% WRITE_NIFTI(HDR, IMG, FNAME) write the NIFTI image specified in HDR and
% IMG into the file FNAME. It checks to make sure that the dimensions
% of the header and image match. This should work for 4D images.


if nargin < 4
    raw = 0;
end

% check whether header and image dimensions match
hdrdim = hdr.dim(2:4);
imgdim = size(img);
if length(imgdim) >= 3
    imgdim = imgdim(1:3);
elseif length(imgdim) == 2
    imgdim = [imgdim 1];
else
    error('Problem with dimensionality of image');
end
if any(hdrdim ~= imgdim)
    error('Dimensionality mismatch between header and image while writing %s', fname);
end

[pathstr, basename, ext] = fileparts(fname);
if strcmp(ext, '.nii')
  if ~strncmp(hdr.magic, 'n+1', 3)
    error('.nii extension but magic field not n+1');
  end
  if hdr.vox_offset < 352
    error('.nii extension but vox_offset < 352');
  end
  if exist(fname, 'file')
    delete(fname);
  end
  write_nifti_hdr(hdr, fname);
  fid = fopen(fname, 'a');
  status = fseek(fid, hdr.vox_offset, 'bof');
  if status == -1
    fseek(fid, 0, 'eof');
    already = ftell(fid);
    fwrite(fid, zeros(hdr.vox_offset - already, 1));
  end
  
elseif strcmp(ext, '.img')
  if strncmp(hdr.magic, 'n+1', 3)
    error('.img extension but magic field is n+1');
  end
  if hdr.vox_offset > 0
    error('.img extension but vox_offset > 0');
  end
  write_nifti_hdr(hdr, fullfile(pathstr, [basename '.hdr']));
  fid = fopen(fname, 'w');
else
  error('Invalid extension.');
end

switch hdr.datatype
    case 2 % DT_UNSIGNED_CHAR
        precision = 'uchar';
    case 4 % DT_SIGNED_SHORT
        precision = 'int16';
    case 8 % DT_SIGNED_INT
        precision = 'int32';
    case 16 % DT_FLOAT
        precision = 'float32';
    case 64 % DT_DOUBLE
        precision = 'double';
    case 512 % DT_UINT16
        precision = 'uint16'; % added for itksnap files 10/11/17 AD
    otherwise
        error('Unsupported data type; datatype = %d\n', h.datatype);
end

if hdr.scl_slope ~= 0 && ~raw
    img = img ./ hdr.scl_slope - hdr.scl_inter;
end

fwrite(fid, img(:), precision);
fclose(fid);
