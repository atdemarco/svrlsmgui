function largest_cur_cluster_size = get_cur_largest_cluster_wise(parameters,options,variables,thresholded)
    switch parameters.tails
        case options.hypodirection{1}
            permtype= 'pos';
            thresholded.thresholded_mask=thresholded.pos_threshed;
        case options.hypodirection{2}
            permtype= 'neg';
            thresholded.thresholded_mask=thresholded.neg_threshed;
        case options.hypodirection{3}
            permtype= 'twotail';
            thresholded.thresholded_mask=thresholded.twotail_threshed;
    end

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