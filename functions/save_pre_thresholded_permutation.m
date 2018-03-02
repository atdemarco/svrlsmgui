function save_pre_thresholded_permutation(variables,parameters,templatevol)
    variables.vo.fname = fullfile(variables.output_folder.clusterwise,['UNthreshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
    svrlsmgui_write_vol(variables.vo, templatevol);