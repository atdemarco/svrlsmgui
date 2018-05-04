function variables = apply_clustering_one_tail(parameters,variables,real_beta_map_vol,sorted_clusters,tail)
    % This is the same for positive and negative tails.
    thresh_index = median([1 round((1-parameters.clusterwise_p) * parameters.PermNumClusterwise) parameters.PermNumClusterwise]); % row 9500 in 10000 permutations.

    % Read in the beta cutoff map (for positive 1-tailed)
    null_beta_cutoff_map = spm_read_vols(spm_vol(variables.files_created.betamask));
    one_tailed_beta_out_vol = real_beta_map_vol; % a copy to modify...

    % Do the voxelwise thresholding 
    switch tail
        case 'pos'
            one_tailed_beta_out_vol(real_beta_map_vol < null_beta_cutoff_map) = 0; % zero out subthreshold values.
        case 'neg'
            one_tailed_beta_out_vol(real_beta_map_vol > null_beta_cutoff_map) = 0; % zero out subthreshold values.
    end
    
    % Calculate how many volumes survive thresholding for display back to user.
    variables.clusterresults.survivingbetavals = numel(find(one_tailed_beta_out_vol(:))); % number of voxel beta values that survive permutation at the requested p value.

    % Write voxelwise thresholded volume - always abs() this for display purposes.
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map (abs).nii');
    svrlsmgui_write_vol(variables.vo, abs(one_tailed_beta_out_vol)); % note that negative betas will come out positive
    variables.files_created.voxelwisethresholdedbetas = variables.vo.fname;

    %% Cluster the pvalue file
    fname = variables.files_created.thresholded_pmap_inv; % we read in the inverted version (now it always exists) so extreme values are actually "peak" maxima ... we will re-vert them in our table...
    volume_to_cluster = spm_read_vols(spm_vol(fname));

    thresh=0;
    [cimg, variables.clusterresults.ctab, variables.clusterresults.peaks,saved] = cluster(fname,thresh,1);
    delete(saved.image); % since this saves in the wrong location and we re-generate this output anyway...
    if ~isempty(variables.clusterresults.ctab) % Calculate clusterwise P value for each cluster.
        T=readtable(saved.table); % read in the table we saved
        T.clusterP = nan(size(T,1),1); % reserve space
        possible_cluster_pvals = (1:numel(sorted_clusters))./numel(sorted_clusters);
        for r = 1 : size(T,1) % get cluster p value, but limit how extreme it can be (i.e. not p = 0)
            num_clusters_gt_cur_clust = max(1,sum(T.nvox(r) < sorted_clusters)); % if 0, then choose the smallest resolved p value (index 1)
            T.clusterP(r) = possible_cluster_pvals(num_clusters_gt_cur_clust); % old --> T.clusterP(r) = sum(cur_nvox < sorted_clusters)/numel(sorted_clusters); % < this can give a p of 0
        end
        U = [T(:,1) T(:,end) T(:,2:end-1)]; % Reorder columns
        U.maxintinv = U.maxint; % inverted max p values.
        U.maxint = 1 - U.maxint; % makes P values back to the normal direction where small is extreme.

        delete(saved.table) % so we can write an updated version.
        clustertable = fullfile(variables.output_folder.clusterwise,'Table of clusters.txt');
        writetable(U,clustertable);
        variables.files_created.clustertable = clustertable;
    else
        T=[];
        variables.files_created.clustertable = [];
    end

    variables.clusterresults.T = T;
    
    if isempty(T)
        warning('No clusters at all (no suprathreshold voxels).')
        T.clusterP = []; % so we don't get errors for not being able to access the field.
    elseif sum(T.clusterP <= parameters.clusterwise_p) == 0
        warning('Clusters, not none significant.')
    end
        
    % Perform clusterwise thresholding on our p_inv file we read in... use this for masking
    variables.clusterresults.clusterthresh = sorted_clusters(thresh_index)-1;
    %out_map = remove_scatter_clusters(volume_to_cluster, variables.clusterresults.clusterthresh);

    %% Write output
    
    %voxels_in_any_clusters = cimg>0;
    %voxels_in_significant_clusters = out_map > 0;
    
    %% Write files containing cluster indices
    
    % All cluster indices, regardless of significance
    variables.files_created.all_cluster_indices = fullfile(variables.output_folder.clusterwise,'All clust indices.nii');
    variables.vo.fname = variables.files_created.all_cluster_indices;
    %variables.vo.private.dat.scl_scope = 1;
    svrlsmgui_write_vol(variables.vo, cimg);
    
    % Only significant cluster indices
    variables.files_created.significant_cluster_indices = fullfile(variables.output_folder.clusterwise,'Significant clust indices.nii');
    variables.vo.fname = variables.files_created.significant_cluster_indices;
    cimg_only_significant = cimg;

    last_signif_clust_num = sum(T.clusterP <= parameters.clusterwise_p); %nb: sum is ok since they're sorted ascending in this column
    cimg_only_significant(cimg_only_significant > last_signif_clust_num) = 0; % zero out clusters whose indices are larger than the last significant cluster index
    svrlsmgui_write_vol(variables.vo, cimg_only_significant);

    %% Write files containing cluster-level p-values
    
    % All clusters, clusterwise thresholded with cluster p-values - not inverted
    variables.files_created.all_cluster_cluster_pvals = fullfile(variables.output_folder.clusterwise,'All clusts cluster pvals.nii');
    tmp = zeros(size(cimg));
    
    for c = 1 : numel(T.clusterP) % i.e. all cluster indices...
        tmp(cimg==c) = T.clusterP(c); 
    end
    
    variables.vo.fname = variables.files_created.all_cluster_cluster_pvals;
    svrlsmgui_write_vol(variables.vo, tmp);
    
    % All clusters, clusterwise thresholded with cluster p-values - inverted
    variables.files_created.all_cluster_cluster_pvals_inv = fullfile(variables.output_folder.clusterwise,'All clusts cluster pvals (inv).nii');
    tmp = zeros(size(cimg));
    for c = 1 : numel(T.clusterP) % i.e. all cluster indices...
        tmp(cimg==c) = 1-T.clusterP(c); % invert the p-value
    end
    variables.vo.fname = variables.files_created.all_cluster_cluster_pvals_inv;
    svrlsmgui_write_vol(variables.vo, tmp);
    
    % Just significant clusters with cluster p values - not inv
    variables.files_created.all_cluster_cluster_pvals = fullfile(variables.output_folder.clusterwise,'Signif clusts cluster pvals.nii');
    tmp = zeros(size(cimg));
    
    for c = 1 : last_signif_clust_num % only to the last significant cluster.
        tmp(cimg==c) = T.clusterP(c); % Non-inverted p-values
    end
    variables.vo.fname = variables.files_created.all_cluster_cluster_pvals;
    svrlsmgui_write_vol(variables.vo, tmp);
    
    % Just significant clusters with cluster p values - inv
    variables.files_created.all_cluster_cluster_pvals_inv = fullfile(variables.output_folder.clusterwise,'Signif clusts cluster pvals (inv).nii');
    tmp = zeros(size(cimg));
    for c = 1 : last_signif_clust_num % only to the last significant cluster.
        tmp(cimg==c) = 1 - T.clusterP(c); % Invert the p-values
    end
    variables.vo.fname = variables.files_created.all_cluster_cluster_pvals_inv;
    svrlsmgui_write_vol(variables.vo, tmp);
    
    %% Write files containing voxel-level p-values
%     % Clusterwise thresholded with voxel p-values - not inverted
%     variables.files_created.all_cluster_voxelwise_pvals = fullfile(variables.output_folder.clusterwise,'All clusters voxelwise p values.nii');
%     
%     % Clusterwise thresholded with voxel p-values - inverted
%     variables.files_created.all_cluster_voxelwise_pvals_inv = fullfile(variables.output_folder.clusterwise,'All clusters voxelwise p values (inv).nii');

    % Just significant clusters with voxelwise p values - inv
    cimg_only_significant_mask = cimg_only_significant>0; % inclusive mask
    variables.files_created.all_cluster_voxelwise_pvals_inv = fullfile(variables.output_folder.clusterwise,'Signif clusts vox pvals (inv).nii');
    p_inv = volume_to_cluster; % use our file we already read in, which is inverted already.
    p_inv(~cimg_only_significant_mask) = 0;    
    variables.vo.fname = variables.files_created.all_cluster_voxelwise_pvals_inv;    
    svrlsmgui_write_vol(variables.vo, p_inv);
    
    % Just significant clusters with voxelwise p values - not inv
    variables.files_created.all_cluster_voxelwise_pvals = fullfile(variables.output_folder.clusterwise,'Signif clusts vox pvals.nii');
    variables.vo.fname = variables.files_created.all_cluster_voxelwise_pvals;
    p_notinv = volume_to_cluster; % use our file we already read in, which is inverted already.
    p_notinv(~cimg_only_significant_mask) = 0;    
    p_notinv(cimg_only_significant_mask) = 1-p_notinv(cimg_only_significant_mask); % this should work.
    svrlsmgui_write_vol(variables.vo, p_notinv);