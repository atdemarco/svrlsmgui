function variables = do_cfwer_clustering(handles,parameters,variables,all_perm_data,thresholded)
% This function runs the clustering procedure on the output of a cfwer-thresholded volume
% The code is based on a modified apply_clustering_one_tail()

% Do the voxelwise thresholding 
switch parameters.tailshort
    case {'pos','neg'}
%        one_tailed_beta_out_vol(real_beta_map_vol < null_beta_cutoff_map) = 0; % zero out subthreshold values.
    case 'two'
        error('not finished for two-tails')
%        one_tailed_beta_out_vol(real_beta_map_vol > null_beta_cutoff_map) = 0; % zero out subthreshold values.
end

% Read in permutation-generated cluster sizes...
clustervals = load(variables.files_created.largest_clusters);
sorted_clusters = sort(clustervals.all_max_cluster_sizes);
    
% Calculate how many voxels survive thresholding for display back to user.
variables.clusterresults.survivingbetavals = -1; % update me...    % number of voxel beta values that survive permutation at the requested p value.

%% Cluster the inv p-value file
fname = variables.files_created.thresholded_pmap_inv; % we read in the inverted version (now it always exists) so extreme values are actually "peak" maxima ... we will re-vert them in our table...
%volume_to_cluster = spm_read_vols(spm_vol(fname));

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
%out_map = remove_scatter_clusters(volume_to_cluster, variables.clusterresults.clusterthresh);

%% Write output

%% Write files containing cluster indices 

% All cluster indices, regardless of significance
variables.files_created.all_cluster_indices = fullfile(variables.output_folder.cfwer,'All clusters indices.nii');
variables.vo.fname = variables.files_created.all_cluster_indices;
svrlsmgui_write_vol(variables.vo, cimg);

% % Only significant cluster indices
% variables.files_created.significant_cluster_indices = fullfile(variables.output_folder.cfwer,'Significant clusters indices.nii');
% variables.vo.fname = variables.files_created.significant_cluster_indices;
% cimg_only_significant = cimg;
% last_signif_clust_num = sum(T.clusterP <= parameters.clusterwise_p); %nb: sum is ok since they're sorted ascending in this column
% cimg_only_significant(cimg_only_significant > last_signif_clust_num) = 0; % zero out clusters whose indices are larger than the last significant cluster index
% svrlsmgui_write_vol(variables.vo, cimg_only_significant);

%% Write files containing cluster-level p-values

% All clusters, clusterwise thresholded with cluster p-values - not inverted
variables.files_created.all_cluster_cluster_pvals = fullfile(variables.output_folder.cfwer,'All clusters cluster p values.nii');
tmp = zeros(size(cimg));
for c = 1 : numel(T.clusterP) % i.e. all cluster indices...
    tmp(cimg==c) = T.clusterP(c); 
end
variables.vo.fname = variables.files_created.all_cluster_cluster_pvals;
svrlsmgui_write_vol(variables.vo, tmp);

% All clusters, clusterwise thresholded with cluster p-values - inverted
variables.files_created.all_cluster_cluster_pvals_inv = fullfile(variables.output_folder.cfwer,'All clusters cluster p values (inv).nii');
tmp = zeros(size(cimg));
for c = 1 : numel(T.clusterP) % i.e. all cluster indices...
    tmp(cimg==c) = 1-T.clusterP(c); % invert the p-value
end
variables.vo.fname = variables.files_created.all_cluster_cluster_pvals_inv;
svrlsmgui_write_vol(variables.vo, tmp);

% % Just significant clusters with cluster p values - not inv
% variables.files_created.all_cluster_cluster_pvals = fullfile(variables.output_folder.cfwer,'Significant clusters cluster p values.nii');
% tmp = zeros(size(cimg));
% for c = 1 : last_signif_clust_num % only to the last significant cluster.
%     tmp(cimg==c) = T.clusterP(c); % Non-inverted p-values
% end
% variables.vo.fname = variables.files_created.all_cluster_cluster_pvals;
% svrlsmgui_write_vol(variables.vo, tmp);

% % Just significant clusters with cluster p values - inv
% variables.files_created.all_cluster_cluster_pvals_inv = fullfile(variables.output_folder.cfwer,'Significant clusters cluster p values (inv).nii');
% tmp = zeros(size(cimg));
% for c = 1 : last_signif_clust_num % only to the last significant cluster.
%     tmp(cimg==c) = 1 - T.clusterP(c); % Invert the p-values
% end
% variables.vo.fname = variables.files_created.all_cluster_cluster_pvals_inv;
% svrlsmgui_write_vol(variables.vo, tmp);

%% Write files containing voxel-level p-values

% % Just significant clusters with voxelwise p values - inv
% cimg_only_significant_mask = cimg_only_significant>0; % inclusive mask
% variables.files_created.all_cluster_voxelwise_pvals_inv = fullfile(variables.output_folder.clusterwise,'Significant cluster voxelwise p values (inv).nii');
% p_inv = volume_to_cluster; % use our file we already read in, which is inverted already.
% p_inv(~cimg_only_significant_mask) = 0;    
% variables.vo.fname = variables.files_created.all_cluster_voxelwise_pvals_inv;    
% svrlsmgui_write_vol(variables.vo, p_inv);

% % Just significant clusters with voxelwise p values - not inv
% variables.files_created.all_cluster_voxelwise_pvals = fullfile(variables.output_folder.clusterwise,'Significant cluster voxelwise p values.nii');
% variables.vo.fname = variables.files_created.all_cluster_voxelwise_pvals;
% p_notinv = volume_to_cluster; % use our file we already read in, which is inverted already.
% p_notinv(~cimg_only_significant_mask) = 0;    
% p_notinv(cimg_only_significant_mask) = 1-p_notinv(cimg_only_significant_mask); % this should work.
% svrlsmgui_write_vol(variables.vo, p_notinv);