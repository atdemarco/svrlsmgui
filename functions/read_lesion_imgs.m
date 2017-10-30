function [variables] = read_lesion_imgs(parameters, variables)
%% Read the lesion data and get the lesion volume of each subject
variables.lesion_vol = zeros(size(variables.one_score));
h = waitbar(0,'Reading lesioned images...');
for ni= 1 : numel(variables.SubjectID)
    waitbar(ni / length(variables.SubjectID)) % show progress.
    fname = fullfile(parameters.lesion_img_folder, [variables.SubjectID{ni}, '.nii']);
    if ~exist(fname,'file')
        error('Cannot find lesion image file: %s\n', fname);
        %continue;
    end
    vo = spm_vol(fname);
    tmp = spm_read_vols(vo);
    tmp(isnan(tmp(:))) = 0; % Denan the image.
    tmp = tmp > 0;  % Binarize
    Ldat(:,:,:,ni) = uint8(tmp);
    variables.lesion_vol(ni,1) = sum(tmp(:));
end

close(h) % close the waitbar...

variables.vo = vo;
variables.vo.name = 'NULL.nii';
variables.vo.dt = [64,0];

%% get a mask based on overlapping map and given mask image
mask_map = sum(Ldat, 4);

variables.l_idx = find(mask_map >= 1); % index of voxels with lesion on at least 1 subject
variables.m_idx = find(mask_map >= parameters.lesion_thresh);  % index of voxels that will be excluded from result

% write out mask_map to check that these numbers make sense...
variables.vo.fname = fullfile(variables.output_folder.base,['All lesion overlap n=' num2str(size(Ldat,4)) '.nii']);
spm_write_vol(variables.vo, mask_map);

% write out a mask of voxels that exceed the minimum lesion threshold requested
variables.vo.fname = fullfile(variables.output_folder.base,['Minimum lesion mask n=' num2str(parameters.lesion_thresh) '.nii']);
min_lesion_mask = mask_map >= parameters.lesion_thresh;
spm_write_vol(variables.vo, min_lesion_mask);

% get the voxel index within the generated mask
Mdat = reshape(Ldat, length(mask_map(:)), variables.SubNum).';
variables.lesion_dat = double(Mdat(:,variables.l_idx));

% Check the thresholded lesion data, remove subjects who have no survival voxel within the mask.
Ldat_tmp = double(Mdat(:,variables.m_idx));
sub_idx = find(sum(Ldat_tmp, 2) == 0);
variables.exclude_idx = sub_idx;
variables.SubNum = variables.SubNum - length(sub_idx);
if ~isempty(sub_idx)
    %fprintf('<strong> Warning: The following subjects have no surviving voxel after thresholding, \nand so will be excluded in the following analysis:\n</strong>')
    for ni = 1 : length(sub_idx)
        %fprintf('<strong>%s </strong>\n', variables.SubjectID{sub_idx(ni)})
        variables.lesion_dat(ni,:) = [];
        variables.one_score(ni) = [];
        variables.lesion_vol(ni) = [];
        variables.excluded_SubjectID{ni} = variables.SubjectID{sub_idx(ni)};
        if ~isempty(variables.covariates)
            variables.covariates(ni,:) = []; % remove whole rows.
        end
       variables.scorefiledata(ni,:) = [];
       variables.SubjectID(sub_idx(ni)) = [];
    end
    %fprintf('\n')
end

variables.neg_idx = [];
variables.pos_idx = [];

%% Write lesion volumes to text file
fid = fopen(fullfile(variables.output_folder.base,'Lesion Volumes.txt'),'w');
fprintf(fid,'Image Name\tLesion Volume (vox)');
for s = 1 : numel(variables.SubjectID) % one row at a time...
    fprintf(fid,'\n%s\t%d',variables.SubjectID{s},variables.lesion_vol(s));
end
fclose(fid);

end