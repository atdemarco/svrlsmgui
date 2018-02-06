function [survivingclusters,totalclusters,survivingbetavals,hypothdirection,clusterthresh] = evaluate_clustering_results(handles,variables,parameters)
    % Read in real beta map
    real_beta_map_vol = spm_read_vols(spm_vol(fullfile(variables.output_folder.base,'Beta map (unthresholded).nii')));

    % Read in permutation generated cluster sizes...
    clustervals = load(fullfile(variables.output_folder.clusterwise,'Largest clusters.mat'));
    sorted_clusters = sort(clustervals.all_max_cluster_sizes);
    switch parameters.tails
        case handles.options.hypodirection{1} % Positive tail
            [ctab,T,hypothdirection,survivingbetavals,clusterthresh] = apply_clustering_pos_tail(parameters,variables,real_beta_map_vol,sorted_clusters);
        case handles.options.hypodirection{2} % Negative tail
            [ctab,T,hypothdirection,survivingbetavals,clusterthresh] = apply_clustering_neg_tail(parameters,variables,real_beta_map_vol,sorted_clusters);
        case handles.options.hypodirection{3} % Two-tailed
            [ctab,T,hypothdirection,survivingbetavals,clusterthresh] = apply_clustering_two_tailed(parameters,variables,real_beta_map_vol,sorted_clusters);
    end

    if ~isempty(ctab)
        survivingclusters = sum(T.clusterP < parameters.clusterwise_p);
        totalclusters = numel(T.clusterP);
    else
        survivingclusters = 0;
        totalclusters = 0;
    end
