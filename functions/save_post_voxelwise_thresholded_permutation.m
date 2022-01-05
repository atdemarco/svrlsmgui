function save_post_voxelwise_thresholded_permutation(parameters,variables,thresholded,f)
    switch parameters.tailshort
        case 'pos' % high scores bad
            variables.vo.fname = fullfile(variables.output_folder.clusterwise,['pos_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
            svrlsmgui_write_vol(variables.vo, thresholded.pos_threshed);
        case 'neg' % high scores good 
            variables.vo.fname = fullfile(variables.output_folder.clusterwise,['neg_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
            svrlsmgui_write_vol(variables.vo, thresholded.neg_threshed);
        case 'two' % two tail
            warning('temporarily disabled')
            variables.vo.fname = fullfile(variables.output_folder.clusterwise,['twotail_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
            svrlsmgui_write_vol(variables.vo, thresholded.twotail_threshed);
    end
