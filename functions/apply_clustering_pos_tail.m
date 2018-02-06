function [ctab,T,hypothdirection,survivingbetavals,clusterthresh] = apply_clustering_pos_tail(parameters,variables,real_beta_map_vol,sorted_clusters)
    hypothdirection = 'one-tailed, positive';
    pos_thresh_index = median([1 round((1-parameters.clusterwise_p) * parameters.PermNumClusterwise) parameters.PermNumClusterwise]); % row 9500 in 10000 permutations.

    % Read in the beta cutoff map (for positive 1 tailed)
    null_beta_cutoff_map_pos = spm_read_vols(spm_vol(fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (positive tail).nii')));
    one_tailed_beta_out_vol = real_beta_map_vol; % a copy to modify...

    % Do the voxelwise thresholding (for positive 1 tailed)
    one_tailed_beta_out_vol(real_beta_map_vol < null_beta_cutoff_map_pos) = 0; % zero out subthreshold values.
    survivingbetavals = numel(find(one_tailed_beta_out_vol(:))); % number of voxel beta values that survive permutation at the requested p value.

    % Write voxelwise thresholded (for positive 1 tailed)
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map.nii');
    spm_write_vol(variables.vo, one_tailed_beta_out_vol);

    % Cluster the resulting data.
    thresh=0;
    [cimg, ctab, peaks] = cluster(variables.vo.fname, thresh,1); %, saveoutputs);
    if ~isempty(ctab) % Calculate clusterwise P value for each cluster.
        clustertablefile = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map_clusttab.txt');
        T=readtable(clustertablefile);
        T.clusterP = nan(size(T,1),1);
        for r = 1 : size(T,1)
            cur_nvox = T.nvox(r);
            T.clusterP(r) = sum(cur_nvox < sorted_clusters)/numel(sorted_clusters);
        end
        U = [T(:,1) T(:,end) T(:,2:end-1)]; % Reorder columns
        delete(clustertablefile)
        writetable(U,fullfile(variables.output_folder.clusterwise,'Table of clusters.txt')); 
    else
        T=[];
    end

    % Perform clusterwise thresholding (for positive 1 tailed)
    clusterthresh = sorted_clusters(pos_thresh_index)-1;
    out_map = remove_scatter_clusters(one_tailed_beta_out_vol, clusterthresh);

    % Write clusterwise thresholded (for positive 1 tailed)
    variables.vo.fname = fullfile(variables.output_folder.clusterwise,'Thresholded by cluster size.nii');
    spm_write_vol(variables.vo, out_map);
