function variables = do_cfwer_clustering(handles,parameters,variables,all_perm_data,thresholded)
    % This function runs the clustering procedure on the output of a cfwer-thresholded volume
    % The code is based on a modified apply_clustering_one_tail()

    % Read in permutation-generated cluster sizes...
    clustervals = load(variables.files_created.largest_clusters);
    sorted_clusters = sort(clustervals.all_max_cluster_sizes);

    % Calculate how many voxels survive thresholding for display back to user.
    variables.clusterresults.survivingbetavals = -1; % update me...    % number of voxel beta values that survive permutation at the requested p value.

    %% Cluster the inv p-value file
    fname = variables.files_created.thresholded_pmap_inv; % we read in the inverted version (now it always exists) so extreme values are actually "peak" maxima ... we will re-vert them in our table...

    thresh=0;
    [cimg, variables.clusterresults.ctab, variables.clusterresults.peaks,saved] = cluster(fname,thresh,1);
    delete(saved.image); % since this saves in the wrong location and we re-generate this output anyway...
    if ~isempty(variables.clusterresults.ctab) % Calculate clusterwise P value for each cluster.
        T=readtable(saved.table);
        T.clusterP = nan(size(T,1),1);
        for r = 1 : size(T,1)
            cur_nvox = T.nvox(r);
            T.clusterP(r) = sum(cur_nvox < sorted_clusters)/numel(sorted_clusters);
        end
        U = [T(:,1) T(:,end) T(:,2:end-1)]; % Reorder columns
        U.maxintinv = U.maxint; % inverted max p values.
        U.maxint = 1 - U.maxint; % makes P values back to the normal direction where small is extreme.

        delete(saved.table) % so we can write an updated version.
        clustertable = fullfile(variables.output_folder.cfwer,'Table of clusters.txt');
        writetable(U,clustertable);
        variables.files_created.clustertable = clustertable;
    else
        T=[];
        variables.files_created.clustertable = [];
    end

    variables.clusterresults.T = T;

    if isempty(T)
        warning('No clusters at all.')
    elseif sum(T.clusterP <= parameters.clusterwise_p) == 0
        warning('No significant clusters.')
    end

    % Perform clusterwise thresholding on our p_inv file we read in... use this for masking
    variables.clusterresults.clusterthresh = -1; % update me -- this is undefined for cfwer. sorted_clusters(thresh_index)-1;

    %% Write output

    %% Write files containing cluster indices 

    % All cluster indices, regardless of significance
    variables.files_created.all_cluster_indices = fullfile(variables.output_folder.cfwer,'All clusts indices.nii');
    variables.vo.fname = variables.files_created.all_cluster_indices;
    svrlsmgui_write_vol(variables.vo, cimg);

    %% Write files containing cluster-level p-values

    % All clusters, clusterwise thresholded with cluster p-values - not inverted
    variables.files_created.all_cluster_cluster_pvals = fullfile(variables.output_folder.cfwer,'All clusts clust pvals.nii');
    tmp = zeros(size(cimg));
    if ~isempty(T)
        for c = 1 : numel(T.clusterP) % i.e. all cluster indices...
            tmp(cimg==c) = T.clusterP(c); 
        end
    end
    variables.vo.fname = variables.files_created.all_cluster_cluster_pvals;
    svrlsmgui_write_vol(variables.vo, tmp);

    % All clusters, clusterwise thresholded with cluster p-values - inverted
    variables.files_created.all_cluster_cluster_pvals_inv = fullfile(variables.output_folder.cfwer,'All clusts clust pvals (inv).nii');
    tmp = zeros(size(cimg));
    if ~isempty(T)
        for c = 1 : numel(T.clusterP) % i.e. all cluster indices...
            tmp(cimg==c) = 1-T.clusterP(c); % invert the p-value
        end
    end
    variables.vo.fname = variables.files_created.all_cluster_cluster_pvals_inv;
    svrlsmgui_write_vol(variables.vo, tmp);