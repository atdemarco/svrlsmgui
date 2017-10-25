function [cimg, ctab, peaks] = cluster(fname, thresh, saveoutputs, minvol, nn_mode)

if nargin < 2, error('Not enough input arguments.'); end
if nargin < 3, saveoutputs = false; end
if nargin < 4, minvol = 0; end
if nargin < 5, nn_mode = 3; end

[hdr, img] = read_nifti(fname);

[fpath,filename,~]=fileparts(fname);
basename = fullfile(fpath,filename);

% when searching for deactivations, invert both the threshold and the image
if thresh < 0
  thresh = -thresh;
  img = -img;
  deactivations = true;
else
  deactivations = false;
end

% mask the image to zero voxels below threshold
img(isnan(img)) = 0;
mask = img >= thresh;
cimg = img .* mask;

% fill in clusters with negative integers
nextfillval = -1;
pos = find(cimg > 0); % find the first voxel in the first cluster
while ~isempty(pos)
  [x, y, z] = ind2sub(size(cimg), pos(1)); % get the coords of the first voxel
  cimg = fillcluster(cimg, x, y, z, nextfillval, nn_mode); % fill the cluster
  nextfillval = nextfillval - 1; % advance the fill value
  pos = find(cimg > 0); % find the first voxel in the next cluster
end

% characterize clusters
nclusters = -nextfillval - 1;
ctab = [];
peaks = [];
for c = 1:nclusters
  i = find(cimg == -c); % find all voxels in the cluster
  
  % volume
  nvoxels = length(i);
  volume = nvoxels * prod(abs(hdr.pixdim(2:4)));
  
  if volume >= minvol
    [x, y, z] = ind2sub(size(cimg), i); % get coords
    intensity = img(i); % get intensities
    [maxint, maxintidx] = max(intensity);

    % center of mass
    cmx = sum(x .* intensity) / sum(intensity);
    cmy = sum(y .* intensity) / sum(intensity);
    cmz = sum(z .* intensity) / sum(intensity);

    % coords of interest
    minextent = img2std([max(x), min(y), min(z)], hdr);
    maxextent = img2std([min(x), max(y), max(z)], hdr);
    peak = img2std([x(maxintidx), y(maxintidx), z(maxintidx)], hdr);
    cm = round(img2std([cmx, cmy, cmz], hdr) .* 10) / 10;

    % put everything together
    ctab = [ctab; c nvoxels volume maxint minextent' maxextent' peak' cm' x(maxintidx) y(maxintidx) z(maxintidx)];
  end
end

% if looking for deactivations then re-invert intensity values
if ~isempty(ctab)
  if deactivations
    ctab(:, 4) = -ctab(:, 4);
  end
end

% renumber clusters according to volume
if ~isempty(ctab)
  ctab = flipud(sortrows(ctab, 2));
  for c = 1:size(ctab, 1)
    cimg(cimg == -ctab(c, 1)) = c;
  end
  cimg(cimg < 0) = 0;
  
  ctab(:, 1) = (1:size(ctab, 1))'; % update first column of indices
  
  peaks = ctab(:, (end - 2):end);
  ctab = ctab(:, 1:(end - 3));
else
  cimg(cimg < 0) = 0;
end

% save cluster image if requested
if saveoutputs
  hdr2 = hdr;
  hdr2.scl_slope = 1;
  write_nifti(hdr2, cimg, [basename '_clustidx.nii']);
  outfid = fopen([basename '_clusttab.txt'], 'w');
end

% print the cluster table
headerline = 'idx\tnvox\tvol\tmaxint\tminx\tminy\tminz\tmaxx\tmaxy\tmaxz\tpeakx\tpeaky\tpeakz\tcmx\tcmy\tcmz\n';
if saveoutputs, fprintf(outfid, headerline); end
if ~isempty(ctab)
  for i = 1:size(ctab, 1)
    ctabline = sprintf('%d\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\n', ctab(i, :));
    if saveoutputs, fprintf(outfid, ctabline); end
  end
end

if saveoutputs, fclose(outfid); end


function img = fillcluster(img, x, y, z, fillval, nn_mode)
% fill the cluster with fillval starting from x, y, z
% note that even corner-to-corner connectivity is acceptable

[sx, sy, sz] = size(img);
img(x, y, z) = fillval;
stack = [x y z];

while ~isempty(stack)
  p = stack(1, :);
  
  for x = p(1) - 1 : p(1) + 1
    for y = p(2) - 1 : p(2) + 1
      for z = p(3) - 1 : p(3) + 1
        if x >= 1 && y >= 1 && z >= 1 && x <= sx && y <= sy && z <= sz
          if img(x, y, z) > 0
            if nn_mode == 3 || sqrt((x - p(1)) ^ 2 + (y - p(2)) ^ 2 + (z - p(3)) ^ 2) <= 1
              img(x, y, z) = fillval;
              stack = [stack; x y z];
            end
          end
        end
      end
    end
  end
  
  if size(stack, 1) == 1
    stack = [];
  else
    stack = stack(2:end, :);
  end
end
