function [ctab,T,hypothdirection,survivingbetavals,clusterthresh] = apply_clustering_two_tailed(parameters,variables,real_beta_map_vol,sorted_clusters)
        hypothdirection = 'two-tailed';
        two_tailed_thresh_index = median([1 round((1-(parameters.clusterwise_p/2)) * parameters.PermNumClusterwise) parameters.PermNumClusterwise]); % row 9750 in 10000 permutations.

        % Read in the beta cutoff map (for two tailed upper and lower tail)
        null_beta_cutoff_map_twotail_upper = spm_read_vols(spm_vol(fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (two tail, upper).nii')));
        null_beta_cutoff_map_twotail_lower = spm_read_vols(spm_vol(fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (two tail, lower).nii')));

        % Do voxelwise thresholding (for two tailed both tails)
        two_tailed_beta_out_vol = real_beta_map_vol; % a copy to modify...
        two_tailed_beta_out_vol((real_beta_map_vol > 0) & (real_beta_map_vol < null_beta_cutoff_map_twotail_upper)) = 0; % zero out subthreshold values.
        two_tailed_beta_out_vol((real_beta_map_vol < 0) & (real_beta_map_vol > null_beta_cutoff_map_twotail_lower)) = 0; % zero out subthreshold values.
        survivingbetavals = numel(find(two_tailed_beta_out_vol(:))); % number of voxel beta values that survive permutation at the requested p value.

        % Write voxelwise thresholded (for two tailed both tails)
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map.nii'); % in the future abs() this.
        spm_write_vol(variables.vo, two_tailed_beta_out_vol);

        % Write ABS of voxelwise thresholded (for two tailed both tails) so that we can measure clusters in one fell swoop.
        variables.vo.fname = fullfile(variables.output_folder.clusterwise,'perm_beta_2tailed_voxwise_ABS.nii');
        spm_write_vol(variables.vo, abs(two_tailed_beta_out_vol));

        % Cluster the resulting data.
        [cimg, ctab, peaks] = cluster(variables.vo.fname, 0,1); % Threshold at 0 and save outputs

        if ~isempty(ctab)
            % Calculate clusterwise P value for each cluster.
            clustertablefile = fullfile(variables.output_folder.clusterwise,'perm_beta_2tailed_voxwise_ABS_clusttab.txt');
            T=readtable(clustertablefile);
            T.clusterP = nan(size(T,1),1);
            for r = 1 : size(T,1)
                cur_nvox = T.nvox(r);   
                T.clusterP(r) = sum(cur_nvox < sorted_clusters)/numel(sorted_clusters);
            end
            U = [T(:,1) T(:,end) T(:,2:end-1)]; % Reorder columns

            % Write out reduced table with cluster p values...
            writetable(U,fullfile(variables.output_folder.clusterwise,'Table of clusters.txt')); % overwrite
            delete(clustertablefile) % clean up this intermediary file.
        else
            T=[];
        end
        delete(fullfile(variables.output_folder.clusterwise,'perm_beta_2tailed_voxwise_ABS.nii')) % and clean up this intermediary file we used to trick cluster() into working in one fell swoop.

        % Perform clusterwise thresholding (for two tailed both tails)
        clusterthresh = sorted_clusters(two_tailed_thresh_index)-1;
        out_map = remove_scatter_clusters(two_tailed_beta_out_vol, clusterthresh);

        % Write clusterwise thresholded (for two tailed both tails)
        variables.vo.fname = fullfile(variables.output_folder.clusterwise,'Thresholded by cluster size.nii');
        spm_write_vol(variables.vo, out_map);

        % Need to rename the clustered idx file to match the other tailed analyses...
        movefile(fullfile(variables.output_folder.clusterwise,'perm_beta_2tailed_voxwise_ABS_clustidx.nii'),fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded beta map_clustidx.nii'));
