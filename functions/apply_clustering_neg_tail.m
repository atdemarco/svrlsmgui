function [ctab,T,hypothdirection,survivingbetavals,clusterthresh] = apply_clustering_neg_tail(parameters,variables,real_beta_map_vol,sorted_clusters)
    hypothdirection = 'one-tailed, negative';
    neg_thresh_index = median([1 round((1-parameters.clusterwise_p) * parameters.PermNumClusterwise) parameters.PermNumClusterwise]); % row 9500 in 10000 permutations.

    % Read in beta cutoff map (for negative 1 tailed)
    null_beta_cutoff_map_neg = spm_read_vols(spm_vol(fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (negative tail).nii')));

    % Do voxelwise thresholding (for negative 1 tail)
    one_tailed_beta_out_vol = real_beta_map_vol; % a copy to modify...
    one_tailed_beta_out_vol(real_beta_map_vol > null_beta_cutoff_map_neg) = 0; % zero out subthreshold values.
    survivingbetavals = numel(find(one_tailed_beta_out_vol(:))); % number of voxel beta values that survive permutation at the requested p value.

    % Write voxelwise thresholded (for negative 1 tailed) - the absolute values was added in v0.1 at Peter's request for display easy
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map (abs).nii');  % since v0.1 we abs() here
    spm_write_vol(variables.vo, abs(one_tailed_beta_out_vol)); % since v0.1 we abs() here

    % Cluster the resulting data.
    %thresh=-.01; % this is negative so we get "deactivations"
    thresh=.001; % now it's positive because we've abs'ed the beta map 
    [cimg, ctab, peaks] = cluster(variables.vo.fname, thresh,1); %, saveoutputs);

    if ~isempty(ctab)
        % Calculate clusterwise P value for each cluster.
        clustertablefile = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map (abs)_clusttab.txt'); % nb since v0.1 we abs() here
        T=readtable(clustertablefile);
        T.clusterP = nan(size(T,1),1);
        for r = 1 : size(T,1)
            cur_nvox = T.nvox(r);
            T.clusterP(r) = sum(cur_nvox < sorted_clusters)/numel(sorted_clusters);
        end
        U = [T(:,1) T(:,end) T(:,2:end-1)]; % Reorder columns

        % Write out reduced table with cluster p values...
        delete(clustertablefile)
        writetable(U,fullfile(variables.output_folder.clusterwise,'Table of clusters.txt')); % overwrite
    else
        T=[];
    end

    % Perform clusterwise thresholding (for negative 1 tailed)
    clusterthresh = sorted_clusters(neg_thresh_index)-1;
    out_map = remove_scatter_clusters(one_tailed_beta_out_vol, clusterthresh);

    % Write clusterwise thresholded (for negative 1 tailed)
    variables.vo.fname = fullfile(variables.output_folder.clusterwise,'Thresholded by cluster size.nii');
    spm_write_vol(variables.vo, out_map);