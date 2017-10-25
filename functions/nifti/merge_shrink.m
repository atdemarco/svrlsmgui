function merge_images(varargin)

if nargin < 2
    error('Not enough input arguments');
end

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

firsthdr = read_nifti_hdr(infnames{1});
firsthdr.dim(1) = 4;
firsthdr.dim(5) = count;
firsthdr.datatype = 4; % int16

write_nifti_hdr(firsthdr, varargin{nargin});

fid = fopen(varargin{nargin}, 'w');
for i = 1:count
    [hdr, img] = read_nifti(infnames{i});
    img(isnan(img)) = 0;
    img = round(img * 100);
    img = min(img, 32767);
    img = max(img, -32768);
    fwrite(fid, img(:), 'int16');
end
fclose(fid);
