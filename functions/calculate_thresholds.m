function thresholds = calculate_thresholds(parameters,variables)
    % note we're working in variables.m_idx not variables.l_idx from now on...

    % one tail pos (high scores bad)
    thresholds.voxelwise_p_value = parameters.voxelwise_p;
    thresholds.pos_thresh_index = median([1 round((1-thresholds.voxelwise_p_value) * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % row 9500 in 10000 permutations.
    thresholds.pos_beta_map_cutoff = nan(1,length(variables.m_idx));
    thresholds.one_tail_pos_alphas = nan(1,length(variables.m_idx)); %reserve space...

    % one tail neg (high scores good)
    thresholds.neg_thresh_index = median([1 round(thresholds.voxelwise_p_value * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % so row 500 in 10000 permutations
    thresholds.neg_beta_map_cutoff = nan(1,length(variables.m_idx));
    thresholds.one_tail_neg_alphas = nan(1,length(variables.m_idx)); %reserve space...

    % one-tailed catch-all index -- this works with compare_real_beta() and replaces the flipped sorting business...
    % it's the same as one tailed neg, since we switch-case sort small first no matter the tail in compare_real_beta()
    thresholds.onetail_cutoff_index =  median([1 round(thresholds.voxelwise_p_value * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % so row 500 in 10000 permutations