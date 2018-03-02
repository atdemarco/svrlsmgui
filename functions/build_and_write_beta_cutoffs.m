function [thresholded,variables] = build_and_write_beta_cutoffs(options,parameters,variables,thresholds,thresholded)
    switch parameters.tailshort % parameters.tails
        case 'pos' % One-tailed positive tail... high scores bad 
            [thresholded,variables] = write_beta_cutoff_pos_tail(variables,thresholds,thresholded);
        case 'neg' % One-tailed negative tail... high scores good
            [thresholded,variables] = write_beta_cutoff_neg_tail(variables,thresholds,thresholded);
        case 'two' % Both tails..
            warning('not enabled at the moment...')
            [thresholded,variables] = write_beta_cutoff_two_tailed(variables,thresholds,thresholded);
    end

function [thresholded,variables] = write_beta_cutoff_pos_tail(variables,thresholds,thresholded)
    % Now write out beta cutoff map.
    thresholded.thresholded_pos = zeros(variables.vo.dim(1:3)); % make a zeros template....
    thresholded.thresholded_pos(variables.m_idx) = thresholds.pos_beta_map_cutoff; % put the 95th percentil beta values back into the lesion indices in a full volume
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (positive tail).nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_pos);
    variables.files_created.betamask = variables.vo.fname;

function [thresholded,variables] = write_beta_cutoff_neg_tail(variables,thresholds,thresholded)
    % Now beta cutoff map for one-taled negative tail...
    thresholded.thresholded_neg = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    thresholded.thresholded_neg(variables.m_idx) = thresholds.neg_beta_map_cutoff; % put the 5th percentil beta values back into the lesion indices in a full volume
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (negative tail).nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_neg);
    variables.files_created.betamask = variables.vo.fname;

function [thresholded,variables] = write_beta_cutoff_two_tailed(variables,thresholds,thresholded)
    warning('make sure these tails are right after code refactor') % ad 2/14/18
    % Two-tailed upper tail
    thresholded.thresholded_twotail_upper = zeros(variables.vo.dim(1:3)); % make a zeros template....        
    thresholded.thresholded_twotail_upper(variables.m_idx) = thresholds.two_tailed_beta_map_cutoff_pos; % put the 2.5th percentil beta values back into the lesion indices in a full volume
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (two tail, upper).nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_twotail_upper);
    variables.files_created.betamask{1} = variables.vo.fname;
    
    % Two-tailed lower tail
    thresholded.thresholded_twotail_lower =  zeros(variables.vo.dim(1:3)); % make a zeros template....        
    thresholded.thresholded_twotail_lower(variables.m_idx) = thresholds.two_tailed_beta_map_cutoff_neg; % put the 2.5th percentil beta values back into the lesion indices in a full volume
    variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (two tail, lower).nii');
    svrlsmgui_write_vol(variables.vo, thresholded.thresholded_twotail_lower);
    variables.files_created.betamask{2} = variables.vo.fname; % note second cell index