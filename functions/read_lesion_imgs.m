function [variables] = read_lesion_imgs(parameters, variables)
%% Read the lesion data and get the lesion volume of each subject
variables.lesion_vol = zeros(size(variables.one_score));
for ni= 1 : numel(variables.SubjectID)
    fname = [variables.SubjectID{ni}, '.nii'];
    fullfname = fullfile(parameters.lesion_img_folder, fname);
    svrlsm_waitbar(parameters.waitbar,ni / length(variables.SubjectID),sprintf('Reading lesion file %s...',fname));
    % in the future, make the line above not interpret underscores (_) in file names as subscript character...
    if ~exist(fullfname,'file') % this should not happen since we've checked for missing files already...
        error('Cannot find lesion image file: %s\n', fullfname); 
    end
    vo = spm_vol(fullfname); % The true voxel intensities of the jth image are given by: val*V.pinfo(1,j) + V.pinfo(2,j)
    tmp = spm_read_vols(vo);
    tmp(isnan(tmp)) = 0; % Denan the image.
    
    if parameters.imagedata.do_binarize
        tmp = tmp > 0;  % Binarize if desired...
    end
    
    Ldat(:,:,:,ni) = uint8(tmp);
    variables.lesion_vol(ni,1) = sum(tmp(:));
    check_for_interrupt(parameters)
end
svrlsm_waitbar(parameters.waitbar,0,'');

variables.vo = vo;
variables.vo.name = 'NULL.nii';
%variables.vo.dt = [4,0]; 
variables.vo.dt = [64,0];

%% get a mask based on overlapping map and given mask image

mask_map = sum(Ldat, 4);
variables.l_idx = find(mask_map >= 1); % index of voxels with lesion on at least 1 subject
variables.m_idx = find(mask_map >= parameters.lesion_thresh);  % index of voxels that will be included from result

% verbose=false;
% if verbose
%     disp(['Report from read_lesion_imgs.m - numel l_idx = ' num2str(numel(variables.l_idx)) '... numel m_idx = ' num2str(numel(variables.m_idx))])
% end

% % Write out mask_map to check that these numbers make sense...
% variables.vo.fname = fullfile(variables.output_folder.base,['All lesion overlap n=' num2str(size(Ldat,4)) '.nii']);
% svrlsmgui_write_vol(variables.vo, mask_map);
% variables.files_created.all_lesion_overlap = variables.vo.fname;

% Write out a mask of voxels that exceed the minimum lesion threshold requested
variables.vo.fname = fullfile(variables.output_folder.base,['Minimum lesion mask n=' num2str(parameters.lesion_thresh) '.nii']);
min_lesion_mask = mask_map >= parameters.lesion_thresh;
svrlsmgui_write_vol(variables.vo, min_lesion_mask);
variables.files_created.min_lesion_mask = variables.vo.fname;

% Get the voxel index within the generated mask
Mdat = reshape(Ldat, length(mask_map(:)), variables.SubNum).';
variables.lesion_dat = double(Mdat(:,variables.l_idx));

% Check the thresholded lesion data, remove subjects who have no survival voxel within the mask.
Ldat_tmp = double(Mdat(:,variables.m_idx)); % only voxel indices above our requested threshold
exclude_sub_bool = sum(Ldat_tmp, 2) == 0;
sub_idx = find(exclude_sub_bool); % these indices have no voxels in our output mask, so remove them from the analysis
variables.exclude_idx = sub_idx;
variables.SubNum = variables.SubNum - length(sub_idx);
variables.excluded.novoxels={}; % even if nobody will be excluded.
if ~isempty(sub_idx)
    for ni = 1 : length(sub_idx)
        cursub_index = sub_idx(ni); 
        % I believe there was a bug here in the original svrlsm code where 
        % ni was used as an index instead of cursub_index aka sub_idx(ni)
        % for removing some of the variables.(asdf) fields below:
        variables.lesion_dat(cursub_index,:) = []; % < here 
        variables.one_score(cursub_index) = []; % < here 
        variables.lesion_vol(cursub_index) = []; % < and here 
        variables.excluded.novoxels{ni} = variables.SubjectID{cursub_index}; % < ok in original code
        variables.scorefiledata(cursub_index,:) = []; % < (not in original svrlsm code)
        variables.SubjectID(cursub_index) = [];  % < ok in original code
%         if ~isempty(variables.covariates)  % < (not in original svrlsm code)
%             variables.covariates(cursub_index,:) = []; % remove whole rows;  % < (not in original svrlsm code)
%         end
        if ~isempty(variables.covariatestable)  % < as of 6/7/18 we use covariables table to support nominal variables
            variables.covariatestable(cursub_index,:) = []; 
        end
    end
end

% Recompute map after subject removal and write out mask_map to check that these numbers make sense...
% This was moved below the subject removal on 2/24/18
mask_map = sum(Ldat(:,:,:,~exclude_sub_bool), 4); % non excluded subjects.
variables.vo.fname = fullfile(variables.output_folder.base,['All lesion overlap n=' num2str(numel(variables.SubjectID)) '.nii']);
svrlsmgui_write_vol(variables.vo, mask_map);
variables.files_created.all_lesion_overlap = variables.vo.fname;

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