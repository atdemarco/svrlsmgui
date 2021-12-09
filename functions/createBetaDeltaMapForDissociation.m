function handles = createBetaDeltaMapForDissociation(handles)
    handles = UpdateProgress(handles,'Dissocation analysis - computing conjunction and disjunction of real and null beta maps...',1);

    dissociation = handles.dissociation; % for convenience.
    %% Where we'll put the results of the dissocation ...
    dissociation.disjunction_svrb_file = fullfile(dissociation.output_folder.base,'svrb_disjunction_allperms.bin');
    dissociation.conjunction_svrb_file = fullfile(dissociation.output_folder.base,'svrb_conjunction_allperms.bin');
    
    %% What function will we use to derive the joint null distributions (and our real data) for conjunction and disjunction
    dissociation.disjunction_null_distribution_fcn = @(x,y) x-y; % simple difference
    dissociation.conjunction_null_distribution_fcn = @(x,y) max(x,y); % the negative tail now contains our conjunction results.
    %dissociation.conjunction_null_distribution_fcn = @(x,y) sign(x.*y) .* sqrt(x.^2 .* y.^2); %

    %% Now create a svrb delta map form the two big base "main effect" file beta map files
    for f = 1 : 2 % Grab a memmap handle to each
        permfile = dir(fullfile(fileparts(dissociation.maineffects{f}.parmsfile),'Voxwise*','*','pmu_beta_maps*.bin')); % This is for regular cluster correction
        if numel(permfile) == 0, permfile = dir(fullfile(fileparts(dissociation.maineffects{f}.parmsfile),'*','pmu_beta_maps*.bin')); end % This is for cfwer correction
        cur_bigfname = fullfile(permfile.folder,permfile.name);
        all_perm_data{f} = memmapfile(cur_bigfname,'Format','single'); %#ok<AGROW>
    end

    %% Make sure we only solve for feature/voxel locations for which we have data in both main effects analyses.
    idx_1_n = numel(dissociation.maineffects{1}.m_idx);
    idx_2_n = numel(dissociation.maineffects{2}.m_idx);
    [m_idx,ia,ib] = intersect(dissociation.maineffects{1}.m_idx,dissociation.maineffects{2}.m_idx);
    dissociation.m_idx = m_idx; % these are out .m_idx for the dissocation
    
    dissociation.orig_betas_A = dissociation.maineffects{1}.ori_beta_vals(ia); % need to use ia so we only get indices overlapping between the 2 analyses
    dissociation.orig_betas_B = dissociation.maineffects{2}.ori_beta_vals(ib); % need to use ib so we only get indices overlapping between the 2 analyses
    
    %% Compute the real data conjunction and disjunction map 
    dissociation.orig_betas_disjunction = feval(dissociation.disjunction_null_distribution_fcn,dissociation.orig_betas_A,dissociation.orig_betas_B); % apply the function we use to calculate the joint distribution...
    dissociation.orig_betas_conjunction = feval(dissociation.conjunction_null_distribution_fcn,dissociation.orig_betas_A,dissociation.orig_betas_B);
    
    idx1_frame_inds = 1:idx_1_n; % we'll offset these by permutation count and reuse these indices.
    idx2_frame_inds = 1:idx_2_n; % we'll offset these by permutation count and reuse these indices.

    if ~exist(dissociation.output_folder.base,'dir'), mkdir(dissociation.output_folder.base); end

    %% Create and save the new permutation data
    disjunct_fileID = fopen(dissociation.disjunction_svrb_file,'w');
    conjunct_fileID = fopen(dissociation.conjunction_svrb_file,'w');
    for p = 1 : dissociation.nperms
        thisperm_offset_1 = (p-1)*idx_1_n; % so for the first permutation our offset is 0.
        thisperm_offset_2 =  (p-1)*idx_2_n; % so for the first permutation our offset is 0.

        result1_thisperm_frame = all_perm_data{1}.Data(thisperm_offset_1 + idx1_frame_inds);
        result2_thisperm_frame = all_perm_data{2}.Data(thisperm_offset_2 + idx2_frame_inds);
        result1_thisperm_relvox = result1_thisperm_frame(ia); % extract overlapping voxel index data
        result2_thisperm_relvox = result2_thisperm_frame(ib); % extract overlapping voxels index data
        
        %% Calculate conjunction and disjunction map for each beta permutation
        disjunction_thisperm_relvox = feval(dissociation.disjunction_null_distribution_fcn,result1_thisperm_relvox,result2_thisperm_relvox); 
        fwrite(disjunct_fileID, disjunction_thisperm_relvox,'single'); % append to the big file

        conjunction_thisperm_relvox = feval(dissociation.conjunction_null_distribution_fcn,result1_thisperm_relvox,result2_thisperm_relvox); 
        fwrite(conjunct_fileID, conjunction_thisperm_relvox,'single'); % append to the big file

    end
    fclose(disjunct_fileID); % close big file
    fclose(conjunct_fileID); % close big file
    
    %% ok, now save some initial volumetric data to show our real data results, and our coverage, etc.
    saveInitialDissociationVolumes(dissociation)
    
    handles.dissociation = dissociation; % for convenience.
    
function saveInitialDissociationVolumes(dissociation)
    %% Save the real delta results before we return
    tmp = dissociation.vo; % this is the volume size info from the first analysis, which matches the second...

    %% Base A information - make sure we have indices all correct
    tmp.fname = fullfile(dissociation.output_folder.base,'A_valid_svrbs.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.m_idx) = dissociation.orig_betas_A;
    svrlsmgui_write_vol(tmp, outimg);

    tmp.fname = fullfile(dissociation.output_folder.base,'A_all_svrbs.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.maineffects{1}.m_idx) = dissociation.maineffects{1}.ori_beta_vals; % all original analysis indices... - values shuld overlap with the intersect version
    svrlsmgui_write_vol(tmp, outimg);

    %% Base B information - make sure we have indices all correct
    tmp.fname = fullfile(dissociation.output_folder.base,'B_valid_svrbs.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.m_idx) = dissociation.orig_betas_B;
    svrlsmgui_write_vol(tmp, outimg);

    tmp.fname = fullfile(dissociation.output_folder.base,'B_all_svrbs.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.maineffects{2}.m_idx) = dissociation.maineffects{2}.ori_beta_vals; % all original analysis indices... - values shuld overlap with the intersect version
    svrlsmgui_write_vol(tmp, outimg);

    %% Write out the lost-coverage volume (places where both analyses don't have valid indices)
    lostinds = setxor(dissociation.maineffects{1}.m_idx,dissociation.maineffects{2}.m_idx);% indices that are only found in one map result (so not included in dissociation analysis)
    tmp.fname = fullfile(dissociation.output_folder.base,'lost_coverage.nii');
    outimg = zeros(tmp.dim);
    if any(lostinds), outimg(lostinds) = 1; end % note these lost indices...
    svrlsmgui_write_vol(tmp, outimg);

    %% Now the real data disjunction map
    tmp.fname = fullfile(dissociation.output_folder.base,'A_disjunct_B_svrb.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.m_idx) = dissociation.orig_betas_disjunction;
    svrlsmgui_write_vol(tmp, outimg);

    %% Now the real data conjunction map
    tmp.fname = fullfile(dissociation.output_folder.base,'A_conjunct_B_svrb.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.m_idx) = dissociation.orig_betas_conjunction;
    svrlsmgui_write_vol(tmp, outimg);