function [parameters,variables,thresholds] = step2_notparallel(handles,parameters,variables,thresholds,all_perm_data)
    handles = UpdateProgress(handles,'Sorting null betas for each lesioned voxel in the dataset (not parallelized).',1);
    h = waitbar(0,sprintf('Sorting null betas for each lesioned voxel in the dataset (N = %d).\n',length(variables.m_idx)),'Tag','WB');
    
    if parameters.do_CFWER
        % where we'll put our p-value converted volumetric data...
        parameters.outfname_big_p = fullfile(variables.output_folder.clusterwise,['pmu_p_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);
        fileID = fopen(parameters.outfname_big_p,'w');
    end
    
    dataRef = all_perm_data.Data; % will this eliminate some overhead
    L = length(variables.m_idx); % for each voxel...
    for col = 1 : length(variables.m_idx)
        check_for_interrupt(parameters)
        curcol = dataRef(col:L:end); % index out each column using skips the length of the data...
        observed_beta = variables.ori_beta_vals(col); % original observed beta value.
        curcol_sorted = sort(curcol); % smallest values at the top..

        if parameters.do_CFWER % Begin Mirman, 2017
            p_vec=nan(size(curcol)); % allocate space
            all_ind = 1:numel(p_vec); % we'll reuse this vector
            for i = all_ind % for each svr beta value in the vector
                ind_to_compare = setdiff(all_ind,i);
                p_vec(i) = 1 - mean(curcol(i) < curcol(ind_to_compare));
            end
            fwrite(fileID, p_vec,'single');
        end

        % Compute beta cutoff values and a pvalue map for the observed betas.
        switch parameters.tails
            case handles.options.hypodirection{1} % 'one_positive'
                thresholds.one_tail_pos_alphasone_tail_pos_alphas(col) = sum(observed_beta > curcol_sorted)/numel(curcol_sorted); % percent of values observed_beta is greater than.
                thresholds.pos_beta_map_cutoff(col) = curcol_sorted(thresholds.pos_thresh_index); % so the 9500th at p of 0.05 on 10,000 permutations
            case handles.options.hypodirection{2} %'one_negative'
                thresholds.one_tail_neg_alphas(col) = sum(observed_beta < curcol_sorted)/numel(curcol_sorted); % percent of values observed_beta is greater than.
                thresholds.neg_beta_map_cutoff(col) = curcol_sorted(thresholds.neg_thresh_index); % so the 500th at p of 0.05 on 10,000 permutations
            case handles.options.hypodirection{3} % 'two'
                thresholds.two_tailed_beta_map_cutoff_pos(col) = curcol_sorted(thresholds.two_tailed_thresh_index); % 250...
                thresholds.two_tailed_beta_map_cutoff_neg(col) = curcol_sorted(thresholds.two_tailed_thresh_index_neg); % 9750...
                thresholds.twotails_alphas(col) = sum(abs(observed_beta) > abs(curcol_sorted))/numel(curcol_sorted); % percent of values observed_beta is greater than.
        end
        waitbar(col/L,h) % show progress.
    end
    close(h)
    
    if parameters.do_CFWER
        fclose(fileID); % close the cfwer permutation data output file.
    end
