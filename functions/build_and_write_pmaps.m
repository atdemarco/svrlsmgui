function thresholded = build_and_write_pmaps(options,parameters,variables,thresholds)

    if parameters.do_CFWER
        parameters.voxelwise_p = variables.cfwer_single_pval_answer; % replace this for the purposes of these function calls.
    end

    switch parameters.tails
        case options.hypodirection{1} % One-tailed positive tail...
            thresholded = write_p_maps_pos_tail(parameters,variables,thresholds);
        case options.hypodirection{2} % One-tailed negative tail...
            thresholded = write_p_maps_neg_tail(parameters,variables,thresholds);
        case options.hypodirection{3} % Both tails..
            thresholded = write_p_maps_two_tailed(parameters,variables,thresholds);
    end

function thresholded = write_p_maps_pos_tail(parameters,variables,thresholds)
    thresholded.thresholded_pos = zeros(variables.vo.dim(1:3)); % make a zeros template....

    if parameters.invert_p_map_flag % it's already inverted
        thresholded.thresholded_pos(variables.m_idx) = thresholds.one_tail_pos_alphas;

        % Write unthresholded P-map for the positive tail
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P values (inv).nii');
        spm_write_vol(variables.vo, thresholded.thresholded_pos);
        
        % Now write out the thresholded P-map for the positive tail
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
        thresholded.thresholded_pos(thresholded.thresholded_pos < (1-parameters.voxelwise_p)) = 0; % zero out sub-threshold p value voxels (note the 1-p)
        spm_write_vol(variables.vo, thresholded.thresholded_pos);

    else % we must invert the values
        thresholded.thresholded_pos(variables.m_idx) = 1 - thresholds.one_tail_pos_alphas;
    
        % Write unthresholded P-map for the positive tail
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P values.nii');
        spm_write_vol(variables.vo, thresholded.thresholded_pos);
        
        % Now write out the thresholded P-map for the positive tail
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values.nii');
        thresholded.thresholded_pos(thresholded.thresholded_pos > parameters.voxelwise_p) = 0; % zero out voxels whose values are greater than p
        spm_write_vol(variables.vo, thresholded.thresholded_pos);
    end


function thresholded = write_p_maps_neg_tail(parameters,variables,thresholds)
    thresholded.thresholded_neg = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    if parameters.invert_p_map_flag  % it's already inverted...
        thresholded.thresholded_neg(variables.m_idx) = thresholds.one_tail_neg_alphas;
        % write out unthresholded negative p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P values (inv).nii');
        spm_write_vol(variables.vo, thresholded.thresholded_neg);
        % write out thresholded negative p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
        thresholded.thresholded_neg(thresholded.thresholded_neg < (1-parameters.voxelwise_p)) = 0; % zero out subthreshold p value voxels (note 1-p)
        spm_write_vol(variables.vo, thresholded.thresholded_neg);
    else
        thresholded.thresholded_neg(variables.m_idx) = 1 - thresholds.one_tail_neg_alphas;
        % write out unthresholded negative p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P values.nii');
        spm_write_vol(variables.vo, thresholded.thresholded_neg);
        % write out thresholded negative p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values.nii');
        thresholded.thresholded_neg(thresholded.thresholded_neg > parameters.voxelwise_p) = 0; % zero out voxels whose values are greater than p
        spm_write_vol(variables.vo, thresholded.thresholded_neg);
    end



function thresholded = write_p_maps_two_tailed(parameters,variables,thresholds)
    thresholded.thresholded_twotails = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    if parameters.invert_p_map_flag % it's already inverted...
        thresholded.thresholded_twotails(variables.m_idx) = thresholds.twotails_alphas;
        % write out unthresholded p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P values (inv).nii');
        spm_write_vol(variables.vo, thresholded.thresholded_twotails);
        % write out thresholded p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
        thresholded.thresholded_twotails(thresholded.thresholded_twotails < (1-(parameters.voxelwise_p/2))) = 0; % zero out subthreshold p value voxels (note 1-p)
        spm_write_vol(variables.vo, thresholded.thresholded_twotails);
    else
        thresholded.thresholded_twotails(variables.m_idx) = 1 - thresholds.twotails_alphas;
        % write out unthresholded p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P values.nii');
        spm_write_vol(variables.vo, thresholded.thresholded_twotails);
        % write out thresholded p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
        thresholded.thresholded_twotails(thresholded.thresholded_twotails > (parameters.voxelwise_p/2)) = 0; % zero out supra-alpha p value voxels
        spm_write_vol(variables.vo, thresholded.thresholded_twotails);
    end




