function write_nifti_hdr(h, fname)

% WRITE_NIFTI_HDR  Write NIFTI header
%
% WRITE_NIFTI_HDR(H, FNAME) writes the header structure H into the file
% FNAME. The validity of H is not checked.

[pathstr, basename, ext] = fileparts(fname);
if strcmp(ext, '.img')
  fname = fullfile(pathstr, [basename '.hdr']);
end

if strcmp(ext, '.nii') && exist(fname, 'file')
  fid = fopen(fname, 'r+');
else
  fid = fopen(fname, 'w');
end

fwrite(fid, h.sizeof_hdr(1), 'int32');
fwrite(fid, padchar(h.data_type, 10), 'char');
fwrite(fid, padchar(h.db_name, 18), 'char');
fwrite(fid, h.extents(1), 'int32');
fwrite(fid, h.session_error(1), 'int16');
fwrite(fid, h.regular(1), 'char');
fwrite(fid, h.dim_info(1), 'char');

fwrite(fid, h.dim(1:8), 'int16');
fwrite(fid, h.intent_p1(1), 'float32');
fwrite(fid, h.intent_p2(1), 'float32');
fwrite(fid, h.intent_p3(1), 'float32');
fwrite(fid, h.intent_code(1), 'int16');
fwrite(fid, h.datatype(1), 'int16');
fwrite(fid, h.bitpix(1), 'int16');
fwrite(fid, h.slice_start(1), 'int16');
fwrite(fid, h.pixdim(1:8), 'float32');
fwrite(fid, h.vox_offset(1), 'float32');
fwrite(fid, h.scl_slope(1), 'float32');
fwrite(fid, h.scl_inter(1), 'float32');
fwrite(fid, h.slice_end(1), 'int16');
fwrite(fid, h.slice_code(1), 'uchar');
fwrite(fid, h.xyzt_units(1), 'uchar');
fwrite(fid, h.cal_max(1), 'float32');
fwrite(fid, h.cal_min(1), 'float32');
fwrite(fid, h.slice_duration(1), 'float32');
fwrite(fid, h.toffset(1), 'float32');
fwrite(fid, h.glmax(1), 'int32');
fwrite(fid, h.glmin(1), 'int32');

fwrite(fid, padchar(h.descrip, 80), 'char');
fwrite(fid, padchar(h.aux_file, 24), 'char');
fwrite(fid, h.qform_code(1), 'int16');
fwrite(fid, h.sform_code(1), 'int16');
fwrite(fid, h.quatern_b(1), 'float32');
fwrite(fid, h.quatern_c(1), 'float32');
fwrite(fid, h.quatern_d(1), 'float32');
fwrite(fid, h.qoffset_x(1), 'float32');
fwrite(fid, h.qoffset_y(1), 'float32');
fwrite(fid, h.qoffset_z(1), 'float32');
fwrite(fid, h.srow_x(1:4), 'float32');
fwrite(fid, h.srow_y(1:4), 'float32');
fwrite(fid, h.srow_z(1:4), 'float32');
fwrite(fid, padchar(h.intent_name, 16), 'char');
fwrite(fid, padchar(h.magic, 4), 'char');

if h.qform_code == 0 && h.sform_code == 0 && isfield(h, 'originator') && ...
        all(h.originator(1:3) >= 1) && all(h.originator(1:3) <= h.dim(2:4))
    warning('Writing ANALYZE originator over the top of NIFTI orientation fields');
    fseek(fid, 253, 'bof');
    fwrite(fid, h.originator, 'int16')'; % over the top of qform_code and the following few fields
end

fclose(fid);


function outstr = padchar(instr, n)

if length(instr) <= n
    outstr = char(zeros(1, n));
    outstr(1:length(instr)) = instr;
else
    outstr = instr(1:n);
end
