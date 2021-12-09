function variables = evaluate_clustering_results(handles,variables,parameters)
    % Read in real beta map
    real_beta_map_vol = spm_read_vols(spm_vol(variables.files_created.unthresholded_betamap));

    % Read in permutation-generated cluster sizes...
    clustervals = load(variables.files_created.largest_clusters);
    sorted_clusters = sort(clustervals.all_max_cluster_sizes);

    switch parameters.tails % we now call apply_clustering_one_tail for all one-tailed analyses.
        case handles.options.hypodirection{1} % Positive tail - high scores bad
            variables = apply_clustering_one_tail(parameters,variables,real_beta_map_vol,sorted_clusters,'pos'); 
        case handles.options.hypodirection{2} % Negative tail - high scores good
            variables = apply_clustering_one_tail(parameters,variables,real_beta_map_vol,sorted_clusters,'neg');
    end

    if ~isempty(variables.clusterresults.ctab)
        variables.clusterresults.survivingclusters = sum(variables.clusterresults.T.clusterP < parameters.clusterwise_p);
        variables.clusterresults.totalclusters = numel(variables.clusterresults.T.clusterP);
    else
        variables.clusterresults.survivingclusters = 0;
        variables.clusterresults.totalclusters = 0;
    end