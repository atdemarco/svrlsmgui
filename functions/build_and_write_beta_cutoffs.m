function [thresholded,variables] = build_and_write_beta_cutoffs(options,parameters,variables,thresholds,thresholded)
    switch parameters.tailshort % parameters.tails
        case 'pos' % One-tailed positive tail... high scores bad 
            [thresholded,variables] = write_beta_cutoff_pos_tail(variables,thresholds,thresholded);
        case 'neg' % One-tailed negative tail... high scores good
            [thresholded,variables] = write_beta_cutoff_neg_tail(variables,thresholds,thresholded);
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