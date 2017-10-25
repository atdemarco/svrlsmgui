function h = make_nifti_hdr(varargin)

% MAKE_NIFTI_HEADER  Create a new NIFTI header structure
%
% H = MAKE_NIFTI_HDR(DATATYPE, DIM, PIXDIM, [ORIGIN])
% H = MAKE_NIFTI_HDR(DATATYPE, FNAME)
%
% Creates a new NIFTI header structure. The geometry can be supplied either
% explicitly with DIM, PIXDIM and (optionally) ORIGIN, or it can be copied
% from an existing image specified in FNAME.
%
% ORIGIN is ignored at the moment.


% process arguments
switch nargin
    case 2
        borrowgeom = true;
        datatype = varargin{1};
        fname = varargin{2};
        b = read_nifti_hdr(fname);
    case {3, 4}
        borrowgeom = false;
        datatype = varargin{1};
        dim = varargin{2};
        pixdim = varargin{3};
        if nargin == 4
            origin = varargin{4};
        else
            origin = [];
        end
    otherwise
        error('Usage error; type help make_nifti_hdr');
end

% create header
h.sizeof_hdr = 348;
h.data_type = '';
h.db_name = '';
h.extents = 0;
h.session_error = 0;
h.regular = 0;
h.dim_info = 0;

if borrowgeom
    h.dim = b.dim;
else
    if length(dim) > 7
        error('Too many dimensions');
    end
    h.dim = zeros(1, 8);
    h.dim(1) = length(dim);
    h.dim(2:(length(dim) + 1)) = dim;
end

h.intent_p1 = 0;
h.intent_p2 = 0;
h.intent_p3 = 0;
h.intent_code = 0;

if ischar(datatype)
    switch datatype
        case {'DT_UNSIGNED_CHAR', 'uchar'}
            h.datatype = 2;
            h.bitpix = 8;
        case {'DT_SIGNED_SHORT', 'int16'}
            h.datatype = 4;
            h.bitpix = 16;
        case {'DT_SIGNED_INT', 'int32'}
            h.datatype = 8;
            h.bitpix = 32;
        case {'DT_FLOAT', 'float32'}
            h.datatype = 16;
            h.bitpix = 32;
        otherwise
            error('Unsupported datatype');
    end
elseif length(datatype) == 1
    h.datatype = datatype;
    switch datatype
        case 2 % DT_UNSIGNED_CHAR
            h.bitpix = 8;
        case 4 % DT_SIGNED_SHORT
            h.bitpix = 16;
        case 8 % DT_SIGNED_INT
            h.bitpix = 32;
        case 16 % DT_FLOAT
            h.bitpix = 32;
        otherwise
            error('Bitpix unknown for this datatype');
    end
elseif length(datatype) == 2
    h.datatype = datatype(1);
    h.bitpix = datatype(2);
    warning('Inserting datatype and bitpix without checking validity');
end

if borrowgeom
    h.pixdim = b.pixdim;
else
    if length(pixdim) > 7
        error('Too many dimensions');
    end
    h.pixdim = zeros(1, 8);
    h.pixdim(1) = length(pixdim);
    h.pixdim(2:(length(pixdim) + 1)) = pixdim;
end

if borrowgeom
    h.slice_start = b.slice_start;
    h.slice_end = b.slice_end;
    h.slice_code = b.slice_code;
    h.xyzt_units = b.xyzt_units;
    h.slice_duration = b.slice_duration;
    h.toffset = b.toffset;
else
    h.slice_start = 0;
    h.slice_end = 0;
    h.slice_code = 0;
    h.xyzt_units = 0; % would be easy to set this to mm and sec
    h.slice_duration = 0;
    h.toffset = 0;
end
    
h.vox_offset = 0; % should be set when writing .nii
h.scl_slope = 1;
h.scl_inter = 0;
h.cal_max = 0;
h.cal_min = 0;
h.glmax = 0;
h.glmin = 0;

h.descrip = '';
h.aux_file = '';
h.intent_name = '';

if borrowgeom
    h.qform_code = b.qform_code;
    h.sform_code = b.sform_code;
    h.quatern_b = b.quatern_b;
    h.quatern_c = b.quatern_c;
    h.quatern_d = b.quatern_d;
    h.qoffset_x = b.qoffset_x;
    h.qoffset_y = b.qoffset_y;
    h.qoffset_z = b.qoffset_z;
    h.srow_x = b.srow_x;
    h.srow_y = b.srow_y;
    h.srow_z = b.srow_z;
    h.magic = b.magic;
else
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
    h.magic = '    '; % pretend to be ANALYZE
end
