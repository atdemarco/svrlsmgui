function thresholded = do_tail_thresholding(parameters,thresholded,templatevol)
    switch parameters.tailshort
        case 'pos' % high scores bad
            thresholded.pos_threshed = templatevol .* (templatevol>=thresholded.thresholded_pos); % elementwise greater than operator to threshold positive tail of test betas
        case 'neg' % highs scores good
            thresholded.neg_threshed = templatevol .* (templatevol<=thresholded.thresholded_neg); % elementwise less than operator to threshold negative tail of test betas.
        case 'two' % two-tailed
            warning('disabled temporarily')
            thresholded.threshmask = and(templatevol > 0,templatevol >= thresholded.thresholded_twotail_upper) | and(templatevol < 0,templatevol <= thresholded.thresholded_twotail_lower);
            thresholded.twotail_threshed = templatevol .* thresholded.threshmask; % mask the two-tailed beta mask for this null data...
    end