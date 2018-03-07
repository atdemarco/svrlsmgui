function [p_val,thresh_beta] = compare_real_beta(realbeta,nullvec,tail,cutoff_index)
    % get a p value for the real beta realtive to the null beta vector...
    betavec = [realbeta nullvec(:)']; % include the observed beta value in with the null generated ones...
    switch tail
        case 'pos' % highest (beta) value have significant p values
            [vals,p] = sort(betavec,'descend');
        case 'neg' % smallest (beta) values has significant p values
            [vals,p] = sort(betavec,'ascend'); % this is default sort() behavior.
        case 'two'
            error('not supported atm')
    end   
   
    r = 1:length(betavec);
    r(p) = r;
    all_p_vals = r/length(betavec);
    p_val = all_p_vals(1);
    thresh_beta = vals(cutoff_index); % this is the beta at the threshold! ... % thresholds.neg_beta_map_cutoff(col)