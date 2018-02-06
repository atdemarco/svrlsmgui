function variables = run_beta_PMU2(parameters, variables, beta_map,handles)   
    variables.ori_beta_vals = beta_map(variables.m_idx).'; % Store original observed beta values.
    
    tic;

    %% step 1: create our permutation data in a giant file we'll read back in slices
    [handles,parameters] = step1(handles,parameters,variables);

    %% Calculate the thresholds (indices, whatnot) based on user settings.
    thresholds = calculate_thresholds(parameters,variables);

    %% Read in gigantic memory mapped file... no matter whether we parallelized or not.
    all_perm_data = memmapfile(parameters.outfname_big,'Format','single');

    %% step 2: sort the betas in the huge data file data and create cutoff values
    [parameters,variables,thresholds] = step2(handles,parameters,variables,thresholds,all_perm_data);

    if parameters.do_CFWER
        [variables.cfwer_pval_null_dist,variables.cfwer_nextvaldist] = get_cfwer_dist(handles,parameters,variables);
%         assignin('base','pval_null_dist',variables.cfwer_pval_null_dist)
%         assignin('base','nextvaldist',variables.cfwer_nextvaldist)
        variables.cfwer_single_pval_answer = prctile(variables.cfwer_pval_null_dist,1-parameters.cfwer_p_value); % this is the value that is used to threshold the data...

        % save variables.cfwer_single_pval_answer
        % save the variables.cfwer_pval_null_dist
        % save the variables.cfwer_nextvaldist
        % make the summaryoutput file utilize these files.

        % we'll use the same function here for cfwer but change the p values inside and not write out the beta cutoff maps...
        thresholded = build_and_write_pmaps(handles.options,parameters,variables,thresholds);

    else
    
        %% Construct volumes of the solved p values and write them out - and write out beta cutoff maps, too
        thresholded = build_and_write_pmaps(handles.options,parameters,variables,thresholds);
        thresholded = build_and_write_beta_cutoffs(handles.options,parameters,variables,thresholds,thresholded);

        do_cluster_thresholding_of_permutations(handles,parameters,variables,all_perm_data,thresholded)
    end
    
    %% cleanup
    [handles,parameters,variables] = cleanup(handles,parameters,variables);

