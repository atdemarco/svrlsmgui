function [thresholded,variables] = build_and_write_pmaps(options,parameters,variables,thresholds)
    switch parameters.tailshort % parameters.tails
        case 'pos' % options.hypodirection{1} % One-tailed positive tail... high scores BAD
            [thresholded,variables] = write_p_maps_pos_tail(parameters,variables,thresholds);
        case 'neg' % options.hypodirection{2} % One-tailed negative tail... high scores GOOD
            [thresholded,variables] = write_p_maps_neg_tail(parameters,variables,thresholds);
        case 'two' % options.hypodirection{3} % Both tails..
            [thresholded,variables] = write_p_maps_two_tailed(parameters,variables,thresholds);
    end

function [thresholded,variables] = write_p_maps_pos_tail(parameters,variables,thresholds)
    if parameters.do_CFWER % note that the parameters struct isn't returned by this function so nothing changes outside this function scope
        parameters.voxelwise_p = variables.cfwerinfo.cfwer_single_pval_answer;
    end

    %% Create and write the un-inverted versions...
    thresholded.thresholded_pos = zeros(variables.vo.dim(1:3)); % make a zeros template....
    thresholded.thresholded_pos(variables.m_idx) = thresholds.one_tail_pos_alphas;
    
    % Write unthresholded P-map for the positive tail
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P map.nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_pos);
    variables.files_created.unthresholded_pmap = variables.vo.fname;

    % Calculate z map
    zmap = p2z(thresholded.thresholded_pos); % note the call to p2z() here

    % write out unthresholded positive Z map
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded Z map.nii');
    svrlsmgui_write_vol(variables.vo, zmap);
    variables.files_created.unthresholded_zmap = variables.vo.fname;

    % Now write out the thresholded P-map for the positive tail
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded P map.nii');
    zero_these_vox = thresholded.thresholded_pos > parameters.voxelwise_p;
    thresholded.thresholded_pos(zero_these_vox) = 0; % zero out voxels whose values are greater than p
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_pos);
    variables.files_created.thresholded_pmap = variables.vo.fname;
    
    % write out thresholded positive Z map
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded Z map.nii');
    zmap(zero_these_vox) = 0; % apply the mask we calculated like 20 lines ago
    svrlsmgui_write_vol(variables.vo, zmap);
    variables.files_created.thresholded_zmap = variables.vo.fname;

    %% Create and write the inverted versions.
    thresholded.thresholded_pos = zeros(variables.vo.dim(1:3)); % make a zeros template....
    thresholded.thresholded_pos(variables.m_idx) = 1 - thresholds.one_tail_pos_alphas;

    % Write unthresholded P-map for the positive tail
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P map (inv).nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_pos);
    variables.files_created.unthresholded_pmap_inv = variables.vo.fname;

    % Now write out the thresholded P-map for the positive tail
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded P map (inv).nii');
    thresholded.thresholded_pos(thresholded.thresholded_pos < (1-parameters.voxelwise_p)) = 0; % zero out sub-threshold p value voxels (note the 1-p)
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_pos);
    variables.files_created.thresholded_pmap_inv = variables.vo.fname;


function [thresholded,variables] = write_p_maps_neg_tail(parameters,variables,thresholds)
    if parameters.do_CFWER % note that the parameters struct isn't returned by this function so nothing changes outside this function scope
        parameters.voxelwise_p = variables.cfwerinfo.cfwer_single_pval_answer;
    end
    
    %% Create and write the non-inverted versions.
    thresholded.thresholded_neg = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    thresholded.thresholded_neg(variables.m_idx) = thresholds.one_tail_neg_alphas;
    
    % write out unthresholded negative p map
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P map.nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_neg);
    variables.files_created.unthresholded_pmap = variables.vo.fname;

    % Calculate z map
    zmap = p2z(thresholded.thresholded_neg); % note the call to p2z() here
    
    % write out unthresholded negative Z map
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded Z map.nii');
    svrlsmgui_write_vol(variables.vo, zmap);
    variables.files_created.unthresholded_zmap = variables.vo.fname;

    % write out thresholded negative p map
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded P map.nii');
    zero_these_vox = thresholded.thresholded_neg > parameters.voxelwise_p;
    thresholded.thresholded_neg(zero_these_vox) = 0; % zero out voxels whose values are greater than p
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_neg);
    variables.files_created.thresholded_pmap = variables.vo.fname;

    % write out thresholded negative Z map
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded Z map.nii');
    zmap(zero_these_vox) = 0; % apply the mask we calculated like 20 lines ago
    svrlsmgui_write_vol(variables.vo, zmap);
    variables.files_created.thresholded_zmap = variables.vo.fname;
    
    %% Create and write the inverted versions (P maps, not Z maps)
    thresholded.thresholded_neg = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    thresholded.thresholded_neg(variables.m_idx) = 1 - thresholds.one_tail_neg_alphas;
    
    % write out unthresholded negative p map
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P map (inv).nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_neg);
    variables.files_created.unthresholded_pmap_inv = variables.vo.fname;
    
    % write out thresholded negative p map
    thresholded.thresholded_neg(thresholded.thresholded_neg < (1-parameters.voxelwise_p)) = 0; % zero out subthreshold p value voxels (note 1-p)
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Voxelwise thresholded P map (inv).nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_neg);
    variables.files_created.thresholded_pmap_inv = variables.vo.fname;

function [thresholded,variables] = write_p_maps_two_tailed(parameters,variables,thresholds)
    error('There should be an abs() around here but there isn''t, make sure this works right, especially with CFWER')
    if parameters.do_CFWER % note that the parameters struct isn't returned by this function so nothing changes outside this function scope
        parameters.voxelwise_p = variables.cfwerinfo.cfwer_single_pval_answer.onetail.pos;
    end

    thresholded.thresholded_twotails = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    if parameters.invert_p_map_flag % it's already inverted...
        thresholded.thresholded_twotails(variables.m_idx) = thresholds.twotails_alphas;
        % write out unthresholded p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P values (inv).nii');
        svrlsmgui_write_vol(variables.vo, thresholded.thresholded_twotails);
        % variables.files_created.unthresholded_pmap{1} = variables.vo.fname;
        % variables.files_created.unthresholded_pmap{2} = variables.vo.fname;
        
        % write out thresholded p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
        thresholded.thresholded_twotails(thresholded.thresholded_twotails < (1-(parameters.voxelwise_p/2))) = 0; % zero out subthreshold p value voxels (note 1-p)
        svrlsmgui_write_vol(variables.vo, thresholded.thresholded_twotails);
        % variables.files_created.unthresholded_pmap{1} = variables.vo.fname;
        % variables.files_created.unthresholded_pmap{2} = variables.vo.fname;
        
        warning('add writing z map thresholded and unthresholded here')
    else
        thresholded.thresholded_twotails(variables.m_idx) = 1 - thresholds.twotails_alphas;
        % write out unthresholded p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Unthresholded P values.nii');
        svrlsmgui_write_vol(variables.vo, thresholded.thresholded_twotails);
        % variables.files_created.unthresholded_pmap{1} = variables.vo.fname;
        % variables.files_created.unthresholded_pmap{2} = variables.vo.fname;

        warning('add writing z map thresholded and unthresholded here')
        
        % write out thresholded p map
        variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
        thresholded.thresholded_twotails(thresholded.thresholded_twotails > (parameters.voxelwise_p/2)) = 0; % zero out supra-alpha p value voxels
        svrlsmgui_write_vol(variables.vo, thresholded.thresholded_twotails);
        % variables.files_created.unthresholded_pmap{1} = variables.vo.fname;
        % variables.files_created.unthresholded_pmap{2} = variables.vo.fname;
    end