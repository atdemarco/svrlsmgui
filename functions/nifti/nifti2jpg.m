function nifti2jpg(fname, skip, flip)

if nargin < 2
    skip = 1;
end
if nargin < 3
    flip = [0 0 1; 1 0 1]; % works for NIC T1s after dicom2analyze
end

[hdr, img] = read_nifti(fname);

vals = img(:);
vals = sort(vals(1:10:end));
anatmin = vals(round(length(vals) * 10 / 100));
anatmax = vals(round(length(vals) * 95 / 100));

for o = 1:3
    for i = 1:skip:hdr.dim(o + 1)
        if o == 1
            slice = squeeze(img(i, :, :));
        elseif o == 2
            slice = squeeze(img(:, i, :));
        else
            slice = squeeze(img(:, :, i));
        end
        if flip(1, o)
            slice = slice';
        end
        if flip(2, o)
            slice = flipud(slice);
        end
        slice = scaleimg(slice, anatmin, anatmax, 0.15, 1);
        imwrite(slice, sprintf('%d_%03d.jpg', o, i), 'JPEG', 'Quality', 50);
    end
end

function sd = scaleimg(img, imin, imax, omin, omax)

% scales an image such that values between imin and imax now lie between
% omin and omax. values above and below are clipped.

sd = omin + (img - imin) / (imax - imin) * (omax - omin);
sd(sd < omin) = omin;
sd(sd > omax) = omax;
