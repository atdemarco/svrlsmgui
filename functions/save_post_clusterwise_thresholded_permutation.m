function save_post_clusterwise_thresholded_permutation(variables,parameters,thresholded,largest_cur_cluster_size,permtype)
    out_map = remove_scatter_clusters(thresholded.testvol_thresholded, largest_cur_cluster_size-1);
    variables.vo.fname = fullfile(variables.output_folder.clusterwise,[permtype '_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumClusterwise) '_largest_cluster.nii']);
    svrlsmgui_write_vol(variables.vo, out_map);
