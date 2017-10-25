function h = read_nifti_hdr(fname)

% READ_NIFTI_HDR  Read NIFTI header file
%
% H = READ_NIFTI_HDR(FNAME) reads the header from the image file FNAME.
%
% It doesn't matter whether you specify the .hdr, .img or .nii file. If no
% extension is provided, then a .nii file is sought, then a .hdr file.
%
% This function can open files with byte swapping which differs from the
% local machine. Whether or not byte swapping differs, the field
% MACHINEFORMAT is set to 'ieee-le' or 'ieee-be' as appropriate. This is
% used by READ_NIFTI in reading the actual image data.
%
% Another nonstandard field written into the header is DATATYPESTR, which
% is a string describing the data type. It is only written for some common
% data types, otherwise there will be no such field.
%
% The ANALYZE originator is read into the ORIGINATOR field. If NIFTI
% location/orientation information is present, then this ORIGINATOR field
% will be garbage.
%
% Another nonstandard field is ORIGIN. If all location/orientation
% information suggests that this is a simple SPM file in RPI with no
% transformations, then ORIGIN is calculated from the NIFTI
% location/orientation information. Otherwise, if the ANALYZE originator
% seems to be valid, then ORIGIN is set based on ORIGINATOR. Otherwise,
% ORIGIN is set to the empty matrix.


didunzip = false;
if strcmpi(fname(end-2:end),'.gz') % it's a gzipped file
    didunzip = true;
    origfname = fname;
    fname = fname(1:end-3);
    unix_ec(['gunzip -c ' origfname]); % -c keeps the original file - we'll delete the unzipped version at the end
end
    

[pathstr, name, ext] = fileparts(fname);
switch ext
    case '.hdr'
        hdrfname = fname;
    case '.nii'
        hdrfname = fname;
    case '.img'
        hdrfname = fullfile(pathstr, [name '.hdr']);
    case ''
        if exist(fullfile(pathstr, [name '.nii'])) == 2
            hdrfname = fullfile(pathstr, [name '.nii']);
        elseif exist(fullfile(pathstr, [name '.hdr'])) == 2
            hdrfname = fullfile(pathstr, [name '.hdr']);
        else
            error('No file with valid extension found');
        end
    otherwise
        error('Unrecognized extension');
end
fid = fopen(hdrfname, 'r');

h.sizeof_hdr = fread(fid, 1, 'int32');
if h.sizeof_hdr ~= 348
    % try opening other-endian
    [fn, perm, machineformat] = fopen(fid);
    fclose(fid);
    switch machineformat
        case 'ieee-le'
            fid = fopen(hdrfname, 'r', 'ieee-be');
        case 'ieee-be'
            fid = fopen(hdrfname, 'r', 'ieee-le');
        otherwise
            error('Unknown machine format');
    end
    h.sizeof_hdr = fread(fid, 1, 'int32');
    if h.sizeof_hdr ~= 348 % still didn't work
        error('File %s does not appear to be ANALYZE or NIFTI', hdrfname);
    end
end
[fn, perm, h.machineformat] = fopen(fid);
% machineformat is not strictly part of the header, but it's used for
% opening the corresponding .img file

h.data_type = fread(fid, 10, '*char');
h.db_name = fread(fid, 18, '*char');
h.extents = fread(fid, 1, 'int32');
h.session_error = fread(fid, 1, 'int16');
h.regular = fread(fid, 1, '*char');
h.dim_info = fread(fid, 1, '*char');

h.dim = fread(fid, 8, 'int16')';
h.intent_p1 = fread(fid, 1, 'float32');
h.intent_p2 = fread(fid, 1, 'float32');
h.intent_p3 = fread(fid, 1, 'float32');
h.intent_code = fread(fid, 1, 'int16');
h.datatype = fread(fid, 1, 'int16');
% h.datatypestr isn't actually used for anything besides looking at
switch h.datatype
    case 2
        h.datatypestr = 'DT_UNSIGNED_CHAR';
    case 4
        h.datatypestr = 'DT_SIGNED_SHORT';
    case 8
        h.datatypestr = 'DT_SIGNED_INT';
    case 16
        h.datatypestr = 'DT_FLOAT';
end
h.bitpix = fread(fid, 1, 'int16')';
h.slice_start = fread(fid, 1, 'int16')';
h.pixdim = fread(fid, 8, 'float32')';
h.vox_offset = fread(fid, 1, 'float32');
h.scl_slope = fread(fid, 1, 'float32');
h.scl_inter = fread(fid, 1, 'float32');
h.slice_end = fread(fid, 1, 'int16');
h.slice_code = fread(fid, 1, 'char');
h.xyzt_units = fread(fid, 1, 'char');
h.cal_max = fread(fid, 1, 'float32');
h.cal_min = fread(fid, 1, 'float32');
h.slice_duration = fread(fid, 1, 'float32');
h.toffset = fread(fid, 1, 'float32');
h.glmax = fread(fid, 1, 'int32');
h.glmin = fread(fid, 1, 'int32');

h.descrip = fread(fid, 80, '*char');
h.aux_file = fread(fid, 24, '*char');
h.qform_code = fread(fid, 1, 'int16');
h.sform_code = fread(fid, 1, 'int16');
h.quatern_b = fread(fid, 1, 'float32');
h.quatern_c = fread(fid, 1, 'float32');
h.quatern_d = fread(fid, 1, 'float32');
h.qoffset_x = fread(fid, 1, 'float32');
h.qoffset_y = fread(fid, 1, 'float32');
h.qoffset_z = fread(fid, 1, 'float32');
h.srow_x = fread(fid, 4, 'float32')';
h.srow_y = fread(fid, 4, 'float32')';
h.srow_z = fread(fid, 4, 'float32')';
h.intent_name = fread(fid, 16, '*char');
h.magic = fread(fid, 4, '*char');

fseek(fid, 253, 'bof');
h.originator = fread(fid, 5, 'int16')'; % over the top of qform_code and the following few fields

% if the location and orientation information is standard for a normalized
% SPM image, then calculate the ANALYZE-style origin
if h.pixdim(1) == -1 && h.qform_code > 0 && h.sform_code > 0 && h.quatern_b == 0 && ...
        h.quatern_c == 1 && h.quatern_d == 0 && h.qoffset_x == h.srow_x(4) && ...
        h.qoffset_y == h.srow_y(4) && h.qoffset_z == h.srow_z(4) && ...
        all(h.srow_x(1:3) == [-h.pixdim(2) 0 0]) && all(h.srow_y(1:3) == [0 h.pixdim(3) 0]) && ...
        all(h.srow_z(1:3) == [0 0 h.pixdim(4)])
    h.origin = 1 - [h.qoffset_x / h.srow_x(1), ...
        h.qoffset_y / h.srow_y(2), h.qoffset_z / h.srow_z(3)];
% or if the originator field seems to be valid, use that instead
elseif all(h.originator(1:3) >= 1) && all(h.originator(1:3) <= h.dim(2:4))
    h.origin = h.originator(1:3);
else % otherwise we don't know what the origin is
    h.origin = [];
end

fclose(fid);

if didunzip
    delete(fname)
end