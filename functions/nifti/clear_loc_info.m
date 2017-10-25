function h = clear_loc_info(h)

% H = CLEAR_LOC_INFO(H)
%
% Strips location and orientation info from NIFTI headers. If run on an
% ANALYZE file, it would zero the origin field.


h.pixdim(1) = 0;
h.qform_code = 0;
h.sform_code = 0;
h.quatern_b = 0;
h.quatern_c = 0;
h.quatern_d = 0;
h.qoffset_x = 0;
h.qoffset_y = 0;
h.qoffset_z = 0;
h.srow_x = zeros(1, 4);
h.srow_y = zeros(1, 4);
h.srow_z = zeros(1, 4);
h.originator = zeros(1, 5);
h.origin = 0;
