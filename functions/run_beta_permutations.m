function variables = run_beta_permutations(parameters, variables, beta_map,handles)
    % used to be called 'run_beta_PMU2.m' which was the refactoring of run_beta_PMU.m - the original
    parameters = ShortenTailString(parameters); % set parameters.tailshort to 'pos','neg', or 'two'
    variables.ori_beta_vals = beta_map(variables.m_idx).'; % Store original observed beta values.

    tic;
    
    %% step 1: create our permutation data in a giant file that we'll read back in slices (can be parallelized or not)
    % this creates a giant file containing the permutation data of svr-beta values
    % note that the architecture of the giant binary file is [perm1][perm2]...[permN] 
    % importantly, the length of each permutation data is length(variables.m_idx), not length(variables.l_idx)
    % this is because m_idx contains voxels that are included in our analysis (meets minimum lesion cutoff)
    % and l_idx contains all voxels with at least one lesion -- all of those voxels (all l_idx are submitted to the svr analysis).
    % but we only save out m_idx indices, and only perform sorting and pval conversionm etc with those m_idx indices
    [handles,parameters,predAndLoss] = step1(handles,parameters,variables);

    variables.predAndLoss_perms = predAndLoss; % store so we get them back in our original caller function

    assignin('base','variables',variables)
    error('a')
    
    %% Calculate the thresholds (indices, whatnot) based on user settings.
    thresholds = calculate_thresholds(parameters,variables);

    %% Read in gigantic memory mapped file whether we are parallelizing or not
    all_perm_data = memmapfile(parameters.outfname_big,'Format','single');

    %% step 2: sort the betas in the huge data file data and create cutoff values
    % if desired, CFWER null p-maps are also created (nothing more) in this process
    [parameters,variables,thresholds] = step2(handles,parameters,variables,thresholds,all_perm_data);
    
    if parameters.do_CFWER
        variables = get_cfwer_dist(handles,parameters,variables);
        [thresholded,variables] = build_and_write_pmaps(handles.options,parameters,variables,thresholds); % we'll use the same function here for cfwer but change the p values inside and not write out the beta cutoff maps...
        variables = do_cfwer_clustering(handles,parameters,variables,all_perm_data,thresholded);
    else
        %% Construct volumes of the solved p values and write them out - and write out beta cutoff maps, too
        [thresholded,variables] = build_and_write_pmaps(handles.options,parameters,variables,thresholds);
        [thresholded,variables] = build_and_write_beta_cutoffs(handles.options,parameters,variables,thresholds,thresholded);

        variables = do_cluster_thresholding_of_permutations(handles,parameters,variables,all_perm_data,thresholded);
    end
    
    %% cleanup
    [handles,parameters,variables] = cleanup(handles,parameters,variables);