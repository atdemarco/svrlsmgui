function variables = apply_clustering_pos_tail(parameters,variables,real_beta_map_vol,sorted_clusters)

    pos_thresh_index = median([1 round((1-parameters.clusterwise_p) * parameters.PermNumClusterwise) parameters.PermNumClusterwise]); % row 9500 in 10000 permutations.

    % Read in the beta cutoff map (for positive 1 tailed)
    null_beta_cutoff_map_pos = spm_read_vols(spm_vol(variables.files_created.betamask));
    one_tailed_beta_out_vol = real_beta_map_vol; % a copy to modify...

    % Do the voxelwise thresholding (for positive 1 tailed)
    one_tailed_beta_out_vol(real_beta_map_vol < null_beta_cutoff_map_pos) = 0; % zero out subthreshold values.
    variables.clusterresults.survivingbetavals = numel(find(one_tailed_beta_out_vol(:))); % number of voxel beta values that survive permutation at the requested p value.

    % Write voxelwise thresholded (for positive 1 tailed)
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map.nii');
    svrlsmgui_write_vol(variables.vo, one_tailed_beta_out_vol);
    variables.files_created.voxelwisethresholdedbetas = variables.vo.fname;

    % Cluster the resulting data.
    outputPvals = true;
    if outputPvals
        fname = variables.files_created.thresholded_pmap; % what we'll cluster.
        voltype ='p values';
     else % original beta values.
        fname = variables.files_created.voxelwisethresholdedbetas; % what we'll cluster.
        voltype ='betas';
    end
    
    [fpath,base,ext] = fileparts(fname);
    volume_to_cluster = spm_read_vols(spm_vol(fname));
    clustertablefname = fullfile(variables.output_folder.voxelwise,[base '_clusttab.txt']);
    clusteredoutfname = fullfile(variables.output_folder.clusterwise,['Thresholded by cluster size (' voltype ').nii']);
    
    thresh=0;
    [variables.clusterresults.cimg, variables.clusterresults.ctab, variables.clusterresults.peaks] = cluster(variables.vo.fname, thresh,1);
    if ~isempty(variables.clusterresults.ctab) % Calculate clusterwise P value for each cluster.
        clustertablefile = clustertablefname;
        T=readtable(clustertablefile);
        T.clusterP = nan(size(T,1),1);
        for r = 1 : size(T,1)
            cur_nvox = T.nvox(r);
            T.clusterP(r) = sum(cur_nvox < sorted_clusters)/numel(sorted_clusters);
        end
        U = [T(:,1) T(:,end) T(:,2:end-1)]; % Reorder columns
        delete(clustertablefile)
        clustertable = fullfile(variables.output_folder.clusterwise,'Table of clusters.txt');
        writetable(U,clustertable);
        variables.files_created.clustertable = clustertable;
    else
        T=[];
        variables.files_created.clustertable = [];
    end

    variables.clusterresults.T = T;
    
    % Perform clusterwise thresholding (for positive 1 tailed)
    variables.clusterresults.clusterthresh = sorted_clusters(pos_thresh_index)-1;
    out_map = remove_scatter_clusters(volume_to_cluster, variables.clusterresults.clusterthresh);

    % Write clusterwise thresholded (for positive 1 tailed)
    variables.vo.fname = clusteredoutfname;
    svrlsmgui_write_vol(variables.vo, out_map);
    variables.files_created.clusteredfile = clusteredoutfname;