function largest_cur_cluster_size = get_cur_largest_cluster_size_cfwer(parameters,variables,cur_per_voxelwise_masked_3d) % thresholded)
    % this is only for one-tailed analyses right now
    % the last input is a binary mask of the current permutation with 1's in voxels that survive the
    % thresholding whole-brain cfwer thresholding

%     permtype = parameters.tailshort;

    CC = bwconncomp(cur_per_voxelwise_masked_3d, 6);
    largest_cur_cluster_size = max(cellfun(@numel,CC.PixelIdxList(1,:))); % max val for numels in each cluster object found

    if isempty(largest_cur_cluster_size)
        largest_cur_cluster_size = 0;
%     else % threshold the volume and write it out.
%         if parameters.SavePostClusterwiseThresholdedPermutations % then save them...
%             save_post_clusterwise_thresholded_permutation(variables,parameters,thresholded,largest_cur_cluster_size,permtype);
%         end
    end