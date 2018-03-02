function variables = apply_clustering_neg_tail(parameters,variables,real_beta_map_vol,sorted_clusters)

    neg_thresh_index = median([1 round((1-parameters.clusterwise_p) * parameters.PermNumClusterwise) parameters.PermNumClusterwise]); % row 9500 in 10000 permutations.

    % Read in beta cutoff map (for negative 1 tailed)
    null_beta_cutoff_map_neg = spm_read_vols(spm_vol(variables.files_created.betamask));
    one_tailed_beta_out_vol = real_beta_map_vol; % a copy to modify...

    % Do voxelwise thresholding (for negative 1 tail)
    one_tailed_beta_out_vol(real_beta_map_vol > null_beta_cutoff_map_neg) = 0; % zero out subthreshold values.
    variables.clusterresults.survivingbetavals = numel(find(one_tailed_beta_out_vol(:))); % number of voxel beta values that survive permutation at the requested p value.

    % Write voxelwise thresholded (for negative 1 tailed) - the absolute values was added in v0.1 at Peter's request for display easy
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map (abs).nii');  % since v0.1 we abs() here
    svrlsmgui_write_vol(variables.vo, abs(one_tailed_beta_out_vol)); % since v0.1 we abs() here
    variables.files_created.voxelwisethresholdedbetas = variables.vo.fname;

    % Cluster the resulting data.
    outputPvals = true;
    if outputPvals
        volume_to_cluster = spm_read_vols(spm_vol(variables.files_created.thresholded_pmap)); % read in the voxelwise thresholded p map....
        % Thresholded P values (inv)_clusttab.txt
        [fpath,base,ext] = fileparts(fname);
        clustertablefname = fullfile(variables.output_folder.voxelwise,[base '_clusttab.txt']); % 'Voxelwise thresholded pval map_clusttab.txt');
        clusteredoutfname = fullfile(variables.output_folder.clusterwise,'Thresholded by cluster size (p values).nii');
%     else % original beta values.
%         variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map.nii'); % this is redundant with an earlier line...
%         volume_to_cluster = one_tailed_beta_out_vol;
%         clustertablefname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map_clusttab.txt'); 
%         clusteredoutfname = fullfile(variables.output_folder.clusterwise,'Thresholded by cluster size (betas).nii');
    end
    
    thresh=0; % .001; % now it's positive because we've abs'ed the beta map 
    [variables.clusterresults.cimg, variables.clusterresults.ctab, variables.clusterresults.peaks] = cluster(variables.vo.fname, thresh,1); %, saveoutputs);
    if ~isempty(variables.clusterresults.ctab) % Calculate clusterwise P value for each cluster.
        clustertablefile = clustertablefname; % nb since v0.1 we abs()'ed here
        T=readtable(clustertablefile);
        T.clusterP = nan(size(T,1),1);
        for r = 1 : size(T,1)
            cur_nvox = T.nvox(r);
            T.clusterP(r) = sum(cur_nvox < sorted_clusters)/numel(sorted_clusters);
        end
        U = [T(:,1) T(:,end) T(:,2:end-1)]; % Reorder columns
        delete(clustertablefile) % Write out reduced table with cluster p values...
        writetable(U,fullfile(variables.output_folder.clusterwise,'Table of clusters.txt')); % overwrite
    else
        T=[];
    end

    % Perform clusterwise thresholding (for negative 1 tailed)
    clusterthresh = sorted_clusters(neg_thresh_index)-1;
    out_map = remove_scatter_clusters(volume_to_cluster, clusterthresh); % one_tailed_beta_out_vol

    % Write clusterwise thresholded (for negative 1 tailed)
    variables.vo.fname = clusteredoutfname; % fullfile(variables.output_folder.clusterwise,'Thresholded by cluster size (betas).nii');
    svrlsmgui_write_vol(variables.vo, out_map);
