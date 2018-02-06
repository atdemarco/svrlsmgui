function thresholded = build_and_write_beta_cutoffs(options,parameters,variables,thresholds,thresholded)
    switch parameters.tails
        case options.hypodirection{1} % One-tailed positive tail...
            thresholded = write_beta_cutoff_pos_tail(variables,thresholds,thresholded);
        case options.hypodirection{2} % One-tailed negative tail...
            thresholded = write_beta_cutoff_neg_tail(variables,thresholds,thresholded);
        case options.hypodirection{3} % Both tails..
            thresholded = write_beta_cutoff_two_tailed(variables,thresholds,thresholded);
    end


function thresholded = write_beta_cutoff_pos_tail(variables,thresholds,thresholded)
    % Now write out beta cutoff map.
    thresholded.thresholded_pos = zeros(variables.vo.dim(1:3)); % make a zeros template....
    thresholded.thresholded_pos(variables.m_idx) = thresholds.pos_beta_map_cutoff; % put the 95th percentil beta values back into the lesion indices in a full volume
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (positive tail).nii');
    spm_write_vol(variables.vo, thresholded.thresholded_pos);

function thresholded = write_beta_cutoff_neg_tail(variables,thresholds,thresholded)
    % Now beta cutoff map for one-taled negative tail...
    thresholded.thresholded_neg = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    thresholded.thresholded_neg(variables.m_idx) = thresholds.neg_beta_map_cutoff; % put the 5th percentil beta values back into the lesion indices in a full volume
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (negative tail).nii');
    spm_write_vol(variables.vo, thresholded.thresholded_neg);

function thresholded = write_beta_cutoff_two_tailed(variables,thresholds,thresholded)
    % Two-tailed upper tail
    thresholded.thresholded_twotail_upper = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    thresholded.thresholded_twotail_upper(variables.m_idx) = thresholds.two_tailed_beta_map_cutoff_pos; % put the 2.5th percentil beta values back into the lesion indices in a full volume
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (two tail, upper).nii');
    spm_write_vol(variables.vo, thresholded.thresholded_twotail_upper);

    % Two-tailed lower tail
    thresholded.thresholded_twotail_lower =  zeros(variables.vo.dim(1:3)); % make a zeros template....        
    thresholded.thresholded_twotail_lower(variables.m_idx) = thresholds.two_tailed_beta_map_cutoff_neg; % put the 2.5th percentil beta values back into the lesion indices in a full volume
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (two tail, lower).nii');
    spm_write_vol(variables.vo, thresholded.thresholded_twotail_lower);