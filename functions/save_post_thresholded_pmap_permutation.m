function save_post_thresholded_pmap_permutation(parameters,variables,cur_perm_image,include_mask,col)
    %% Perform the whole-brain voxelwise thresholding for this permutation
    cur_perm_image(~include_mask) = 0; % zero out voxels that didn't achieve our critical value...
    
    %% Write out inverted and non-inverted p images for this permutation (post-thresholding)
    labels = {'','_inv'};
    zerosmap = zeros(variables.vo.dim(1:3)); % make a zeros template....
    for L = labels
        if strcmp(L{1},'_inv')
            cur_perm_image = 1 - cur_perm_image; % invert (1-val) these indices
        end
        zerosmap(variables.m_idx) = cur_perm_image;
        variables.vo.fname = fullfile(variables.output_folder.cfwer,['Thresholded_pmap_perm' num2str(col) '_of_' num2str(parameters.PermNumVoxelwise) L{1} '.nii']);
        svrlsmgui_write_vol(variables.vo, zerosmap);
    end