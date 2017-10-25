function analyze_to_nifti(fname)

% ANALYZE_TO_NIFTI(FNAME)
%
% Converts the ANALYZE image FNAME to .nii (one file) format. If no FNAME
% is provided, converts every .img file in the current directory.
%
% The header information is set up based on what SPM5 does if you reslice.
% The origin is set to the center.


if nargin == 0
  d = dir('*.img');
  if length(d) >= 1
    for i = 1:length(d)
      fprintf('%s\n', d(i).name);
      do(d(i).name);
    end
  end
else
  do(fname);
end


function do(fname)

[hdr, img] = read_nifti(fname);

nhdr = make_nifti_hdr(hdr.datatype, hdr.dim(2:4), abs(hdr.pixdim(2:4)));

x = (nhdr.dim(2) / 2 - 0.5) * nhdr.pixdim(2);
y = (nhdr.dim(3) / 2 - 0.5) * nhdr.pixdim(3);
z = (nhdr.dim(4) / 2 - 0.5) * nhdr.pixdim(4);

nhdr.pixdim(1) = -1;
nhdr.vox_offset = 352;
nhdr.qform_code = 2;
nhdr.sform_code = 2;
nhdr.quatern_b = 0;
nhdr.quatern_c = 1;
nhdr.quatern_d = 0;
nhdr.qoffset_x = x;
nhdr.qoffset_y = -y;
nhdr.qoffset_z = -z;
nhdr.srow_x = [-nhdr.pixdim(2) 0 0 x];
nhdr.srow_y = [0 nhdr.pixdim(3) 0 -y];
nhdr.srow_z = [0 0 nhdr.pixdim(4) -z];
nhdr.magic = 'n+1 ';
nhdr.magic(4) = 0;

fname((end - 2):end) = 'nii';

write_nifti(nhdr, img, fname);
