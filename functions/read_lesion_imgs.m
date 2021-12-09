function variables = read_lesion_imgs(parameters, variables)
    if parameters.use_analysis_mask % Then analysis will only be conducted within voxels falling within the specified mask.
        svrlsm_waitbar(parameters.waitbar,0,'Reading mask file...');
        if ~exist(parameters.analysis_mask_file,'file'), error('Cannot find analysis mask file: %s\n', parameters.analysis_mask_file); end
        vo = spm_vol(parameters.analysis_mask_file); % The true voxel intensities of the jth image are given by: val*V.pinfo(1,j) + V.pinfo(2,j)
        tmp = spm_read_vols(vo);
        tmp(isnan(tmp)) = 0; % Denan 
    	tmp = tmp > 0;  % Binarize

        if parameters.imagedata.do_resample, [tmp,vo] = makeResizedImgIfDesired(parameters.imagedata.resample_to,tmp,vo); end % tmp is the original volume, vo is its header
        variables.AnalysisMask = tmp; % store.
        vo.fname = fullfile(variables.output_folder.base,'Analysis mask.nii');
        variables.files_created.analysis_mask = vo.fname;
        svrlsmgui_write_vol(vo,variables.AnalysisMask);
    end

    %% Ok now if an atlas was specified, then prepare that.
    if parameters.use_atlas_parcellation
        if ~exist(parameters.analysis_parcellation_file,'file'), error('Cannot find analysis parcellation file: %s\n', analysis_parcellation_file); end
        vo = spm_vol(parameters.analysis_parcellation_file);
        tmp = spm_read_vols(vo);
        if parameters.imagedata.do_resample, [tmp,vo] = makeResizedImgIfDesired(parameters.imagedata.resample_to,tmp,vo); end % tmp is the original volume, vo is its header
        variables.AnalysisParcellation = tmp; % store.
        vo.fname = fullfile(variables.output_folder.base,'Analysis Parcellation.nii');
        variables.files_created.analysis_parcellation = vo.fname;
        svrlsmgui_write_vol(vo,variables.AnalysisParcellation);
        variables.parcellation_vals = unique(variables.AnalysisParcellation);
        variables.parcellation_vals(variables.parcellation_vals==0) = []; % remove roi "0"
    end

    %% Read the lesion data and get the lesion volume of each subject
    variables.lesion_vol = zeros(size(variables.one_score));
    for ni = 1 : numel(variables.SubjectID)
        fname = [variables.SubjectID{ni}, '.nii'];
        fullfname = fullfile(parameters.lesion_img_folder, fname);
        svrlsm_waitbar(parameters.waitbar,ni / length(variables.SubjectID),sprintf('Reading lesion file %s...',strrep(fname,'_',''))); % replace underscores to skip subscript display
        if ~exist(fullfname,'file'), error('Cannot find lesion image file: %s\n', fullfname); end % this should not happen since we've checked for missing files already...
        vo = spm_vol(fullfname); % The true voxel intensities of the jth image are given by: val*V.pinfo(1,j) + V.pinfo(2,j)
        tmp = spm_read_vols(vo);
        tmp(isnan(tmp)) = 0; % Denan the image.

        if parameters.imagedata.do_binarize, tmp = tmp > 0; end % Binarize if desired...
        if parameters.imagedata.do_resample, [tmp,vo] = makeResizedImgIfDesired(parameters.imagedata.resample_to,tmp,vo); end % Does user want to resample the images prior to analysis? - implemented June 16 2019 - ad
        variables.lesion_vol(ni,1) = sum(tmp(:)); % Calculate lesion volume before any masking is conducted!
        if parameters.use_analysis_mask
            if ~all(size(variables.AnalysisMask) == size(tmp))
                error(['Image resolution of at least one lesion file (' fname ', ' num2str(size(tmp)) ') does not match the mask resolution (' num2str(size(variables.AnalysisMask)) '). Halting.'])
            end
            tmp(~variables.AnalysisMask) = 0; 
        end % zero out outside the mask

        if ~parameters.use_atlas_parcellation  % then business as usual - voxelwise
            Ldat(:,:,:,ni) = single(tmp); %#ok<AGROW>
        else
            if ~all(size(variables.AnalysisParcellation) == size(tmp)), error(['Image resolution of at least one lesion file (' fname ', ' num2str(size(tmp)) ') does not match the parcellation atlas resolution (' num2str(size(variables.AnalysisParcellation)) '). Halting.']), end
            % calculate the percent damage within each atlas ROI....
            Ldat(:,ni) = single(arrayfun(@(x) nnz(tmp&(variables.AnalysisParcellation==x)) / nnz((variables.AnalysisParcellation==x)),variables.parcellation_vals)); %#ok<AGROW>
        end
        check_for_interrupt(parameters)
    end
    
    svrlsm_waitbar(parameters.waitbar,0,'');

    variables.vo = vo;
    variables.vo.name = 'NULL.nii';
    variables.vo.dt = [64,0];

    %% Get a mask based on overlapping map and given mask image
    % June 25, 2021 updated this binarize behavior
    if parameters.imagedata.do_binarize  % we use ndims() to support the atlas stuff, which is 2D and not 4D as in voxelwise case
        mask_map = sum(Ldat, ndims(Ldat));%4); % classic behavior
    else
        mask_map = sum(Ldat~=0, ndims(Ldat)); %4); % now give a credit of "1" (i.e. true) for values at any voxel that isn't 0. this will prevent naughty masking behavior on continously valued (i.e. non-binarized) data - bug reported by Dr. Lorca puls - thanks Diego!
    end

    variables.l_idx = find(mask_map >= 1); % index of voxels with lesion on at least 1 subject
    variables.m_idx = find(mask_map >= parameters.lesion_thresh);  % index of voxels that will be included from result

    % Write out a mask of voxels that exceed the minimum lesion threshold requested
    variables.vo.fname = fullfile(variables.output_folder.base,['Minimum lesion mask n=' num2str(parameters.lesion_thresh) '.nii']);
    min_lesion_mask = mask_map >= parameters.lesion_thresh;
    svrlsmgui_write_vol(variables.vo, min_lesion_mask);
    variables.files_created.min_lesion_mask = variables.vo.fname;

    %% Get the voxel index within the generated mask - voxelwise
    Mdat = reshape(Ldat, length(mask_map(:)), variables.SubNum).';
    variables.lesion_dat = double(Mdat(:,variables.l_idx));

    % Check the thresholded lesion data, remove subjects who have no survival voxel within the mask.
    Ldat_tmp = double(Mdat(:,variables.m_idx)); % only voxel indices above our requested threshold
    exclude_sub_bool = sum(Ldat_tmp, 2) == 0;
    sub_idx = find(exclude_sub_bool); % these indices have no voxels in our output mask, so remove them from the analysis
    variables.exclude_idx = sub_idx;
    variables.SubNum = variables.SubNum - length(sub_idx);
    variables.excluded.novoxels={}; % even if nobody will be excluded.

    % This section updated on 8/8/18
    if ~isempty(sub_idx) % Then there are some subjects to remove from the analysis
        variables.lesion_dat(sub_idx,:) = [];
        variables.one_score(sub_idx) = [];
        variables.lesion_vol(sub_idx) = [];
        variables.excluded.novoxels(1:numel(sub_idx)) = variables.SubjectID(sub_idx);
        variables.scorefiledata(sub_idx,:) = [];
        variables.SubjectID(sub_idx) = [];
        if ~isempty(variables.covariatestable), variables.covariatestable(sub_idx,:) = []; end % < as of 6/7/18 we use covariables table to support nominal variables
    end

    %% Recompute map after subject removal and write out mask_map to check that these numbers make sense... this was moved below the subject removal on 2/24/18
    mask_map = sum(Ldat(:,:,:,~exclude_sub_bool), ndims(Ldat));%4); % non excluded subjects.
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

function [tmp,vo] = makeResizedImgIfDesired(resample_to,tmp,vo) % tmp is the original volume, vo is its header
    new_vox_size = resample_to .* ones(1,3); % make the new mm into a 1x3... 
    bb = spm_get_bbox(vo);
    vo.mat = spm_matrix([bb(1,:) 0 0 0 new_vox_size])*spm_matrix([-1 -1 -1]);
    vo.dim = ceil(vo.mat \ [bb(2,:) 1]' - 0.1)';
    vo.dim = vo.dim(1:3);
    new_image_dims = vo.dim;
    tmp = imresize3(single(tmp),new_image_dims,'nearest'); % don't interpolate any new values