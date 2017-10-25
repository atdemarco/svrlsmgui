function merge_images(varargin)

if nargin < 2
    error('Not enough input arguments');
end

remove_nans = true;

infnames = {};
count = 0;

for i = 1:(nargin - 1)
    if any(varargin{i} == '*')
        d = dir(varargin{i});
        for j = 1:length(d)
            if d(j).name(1) ~= '.'
                count = count + 1;
                infnames{count} = d(j).name;
            end
        end
    else
        count = count + 1;
        infnames{count} = varargin{i};
    end
end

maxdatatype = 0;
for i = 1:count
    hdr = read_nifti_hdr(infnames{i});
    
    if hdr.datatype > maxdatatype
        maxdatatype = hdr.datatype;
    end
end

switch maxdatatype
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
    otherwise
        error(sprintf('Unsupported data type; datatype = %d\n', hdr.datatype));
end

[firsthdr, firstimg] = read_nifti(infnames{1});
if remove_nans
  firstimg(isnan(firstimg)) = 0;
end
firsthdr.dim(1) = 4;
firsthdr.dim(5) = count;
firsthdr.datatype = maxdatatype;

write_nifti(firsthdr, firstimg, varargin{nargin});

fid = fopen(varargin{nargin}, 'a');
for i = 2:count
    [hdr, img] = read_nifti(infnames{i});
    if remove_nans
      img(isnan(img)) = 0;
    end
    
    if firsthdr.scl_slope ~= 0
        img = img ./ firsthdr.scl_slope - firsthdr.scl_inter;
    end

    fwrite(fid, img(:), precision);
end
fclose(fid);
