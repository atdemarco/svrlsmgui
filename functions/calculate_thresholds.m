function thresholds = calculate_thresholds(parameters,variables)
    thresholds.voxelwise_p_value = parameters.voxelwise_p;
    thresholds.pos_thresh_index = median([1 round((1-thresholds.voxelwise_p_value) * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % row 9500 in 10000 permutations.
    thresholds.pos_beta_map_cutoff = nan(1,length(variables.m_idx));
    thresholds.one_tail_pos_alphas = nan(1,length(variables.m_idx));

    thresholds.neg_thresh_index = median([1 round(thresholds.voxelwise_p_value * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % so row 500 in 10000 permutations
    thresholds.neg_beta_map_cutoff = nan(1,length(variables.m_idx));
    thresholds.one_tail_neg_alphas = nan(1,length(variables.m_idx));

    thresholds.two_tailed_thresh_index_neg = median([1 round(((thresholds.voxelwise_p_value/2)) * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % row 250 in 10000 permutations.
    thresholds.two_tailed_thresh_index = median([1 round((1-(thresholds.voxelwise_p_value/2)) * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % row 9750 in 10000 permutations.
    thresholds.two_tailed_beta_map_cutoff_pos = nan(1,length(variables.m_idx));
    thresholds.two_tailed_beta_map_cutoff_neg = nan(1,length(variables.m_idx));
    thresholds.twotails_alphas = nan(1,length(variables.m_idx));
