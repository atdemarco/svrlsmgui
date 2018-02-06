function save_post_voxelwise_thresholded_permutation(parameters,variables,options,thresholded)
    switch parameters.tails
        case options.hypodirection{1} % 'one_positive'
            variables.vo.fname = fullfile(variables.output_folder.clusterwise,['pos_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
            spm_write_vol(variables.vo, thresholded.pos_threshed);
        case options.hypodirection{2} %'one_negative'
            variables.vo.fname = fullfile(variables.output_folder.clusterwise,['neg_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
            spm_write_vol(variables.vo, thresholded.neg_threshed);
        case options.hypodirection{3} %'two'
            variables.vo.fname = fullfile(variables.output_folder.clusterwise,['twotail_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
            spm_write_vol(variables.vo, thresholded.twotail_threshed);
    end
