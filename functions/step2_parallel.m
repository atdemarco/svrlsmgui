function [parameters,variables,thresholds] = step2_parallel(handles,parameters,variables,thresholds,all_perm_data)
    handles = UpdateProgress(handles,'Sorting null betas for each lesioned voxel in the dataset (parallelized).',1);
    L = length(variables.m_idx);
    tails = parameters.tails; % so not a broadcast variable.
    Opt1 = handles.options.hypodirection{1};
    Opt2 = handles.options.hypodirection{2};
    Opt3 = handles.options.hypodirection{3};
    ori_beta_vals = variables.ori_beta_vals; % for parfor...
    
    
    % we need to make a billion files and combine them like in parallelized step 1...
    if parameters.do_CFWER
    end
    
    parfor col = 1 : length(variables.m_idx)
        check_for_interrupt(parameters)
        curcol = extractSlice(all_perm_data,col,L); % note this is a function at the bottom of this file..
        observed_beta = ori_beta_vals(col); % original observed beta value.
        curcol_sorted = sort(curcol); % smallest values at the top..
        
        switch tails
            case Opt1
                one_tail_pos_alphas(col) = sum(observed_beta > curcol_sorted)/numel(curcol_sorted); % percent of values observed_beta is greater than.
                pos_beta_map_cutoff(col) = curcol_sorted(thresholds.pos_thresh_index); % so the 9500th at p of 0.05 on 10,000 permutations
            case Opt2
                one_tail_neg_alphas(col) = sum(observed_beta < curcol_sorted)/numel(curcol_sorted); % percent of values observed_beta is greater than.
                neg_beta_map_cutoff(col) = curcol_sorted(thresholds.neg_thresh_index); % so the 500th at p of 0.05 on 10,000 permutations
            case Opt3
                two_tailed_beta_map_cutoff_pos(col) = curcol_sorted(thresholds.two_tailed_thresh_index); % 250...
                two_tailed_beta_map_cutoff_neg(col) = curcol_sorted(thresholds.two_tailed_thresh_index_neg); % 9750...
                twotails_alphas(col) = sum(abs(observed_beta) > abs(curcol_sorted))/numel(curcol_sorted); % percent of values observed_beta is greater than.
        end
    end
    
    % we do this now since we can't index into fields of 'thresholds' in a parfor loop.
    switch tails
        case Opt1
            thresholds.one_tail_pos_alphas = one_tail_pos_alphas;
            thresholds.pos_beta_map_cutoff = pos_beta_map_cutoff;
        case Opt2
            thresholds.one_tail_neg_alphas = one_tail_neg_alphas;
            thresholds.neg_beta_map_cutoff = neg_beta_map_cutoff;
        case Opt3
            thresholds.two_tailed_beta_map_cutoff_pos =two_tailed_beta_map_cutoff_pos;
            thresholds.two_tailed_beta_map_cutoff_neg =two_tailed_beta_map_cutoff_neg;
            thresholds.twotails_alphas = twotails_alphas;
    end
    
    
    