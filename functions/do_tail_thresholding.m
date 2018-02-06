function thresholded = do_tail_thresholding(parameters,thresholded,options,templatevol)
    switch parameters.tails
        case options.hypodirection{1}
            thresholded.pos_threshed = templatevol .* (templatevol>=thresholded.thresholded_pos); % elementwise greater than operator to threshold positive tail of test betas
        case options.hypodirection{2}
            thresholded.neg_threshed = templatevol .* (templatevol<=thresholded.thresholded_neg); % elementwise less than operator to threshold negative tail of test betas.
        case options.hypodirection{3} % Now build a two-tailed thresholded version
            thresholded.threshmask = and(templatevol > 0,templatevol >= thresholded.thresholded_twotail_upper) | and(templatevol < 0,templatevol <= thresholded.thresholded_twotail_lower);
            thresholded.twotail_threshed = templatevol .* thresholded.threshmask; % mask the two-tailed beta mask for this null data...
    end
