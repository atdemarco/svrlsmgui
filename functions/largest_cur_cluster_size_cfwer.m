function largest_cur_cluster_size = get_cur_cluster_size_cfwer(parameters,variables,thresholded)
    % this is only for one-tailed analyses...

    permtype = parameters.tailshort;
    thresholded.thresholded_mask=thresholded.([permtype '_threshed']); % note dynamic field reference

    thresholded.testvol_thresholded = thresholded.thresholded_mask; % now evaluate the surviving voxels for clusters...
    CC = bwconncomp(thresholded.testvol_thresholded, 6);
    largest_cur_cluster_size = max(cellfun(@numel,CC.PixelIdxList(1,:))); % max val for numels in each cluster object found

    if isempty(largest_cur_cluster_size)
        largest_cur_cluster_size = 0;
    else % threshold the volume and write it out.
        if parameters.SavePostClusterwiseThresholdedPermutations % then save them...
            save_post_clusterwise_thresholded_permutation(variables,parameters,thresholded,largest_cur_cluster_size,permtype);
        end
    end