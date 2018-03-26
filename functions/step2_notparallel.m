function [parameters,variables,thresholds] = step2_notparallel(handles,parameters,variables,thresholds,all_perm_data)
    handles = UpdateProgress(handles,'Sorting null betas for each lesioned voxel in the dataset (not parallelized).',1);
    if parameters.do_CFWER
        handles = UpdateProgress(handles,'Also assembling CFWER null p-maps (not parallelized).',1);
    end

    message = sprintf('Sorting null betas for each lesioned voxel in the dataset (N = %d).',length(variables.m_idx));
    svrlsm_waitbar(parameters.waitbar,0,message);

    dataRef = all_perm_data.Data; % will this eliminate some overhead
    L = length(variables.m_idx); % for each voxel...
    for col = 1 : L
        if ~mod(col,50) % to reduce num of calls...
            check_for_interrupt(parameters)
            svrlsm_waitbar(parameters.waitbar,col/L)
        end

        curcol = dataRef(col:L:end); % index out each voxel column using skips the length of the data ---> this returns curcol, which has the length of npermutations
        observed_beta = variables.ori_beta_vals(col); % original observed beta value.
        curcol_sorted = sort(curcol); % smallest values at the top.. removed since we se compare_real_beta()
        
        if parameters.do_CFWER
            if col == 1 % open where we'll put our p-value converted volumetric data...
                parameters.outfname_big_p = fullfile(variables.output_folder.clusterwise,['pmu_p_maps_N_' num2str(length(variables.m_idx)) '.bin']);
                fileID = fopen(parameters.outfname_big_p,'w');
            end
            
            p_vec = betas2pvals(curcol,parameters.tailshort);
            
            fwrite(fileID, p_vec,'single'); % add our results from this voxel to the same giant file....

            if col == L % length(variables.m_idx) % then it's the last permutation, and our file is written.
                fclose(fileID); % close the cfwer permutation data output file
            end
        end

%         assignin('base','parameters',parameters)
%         assignin('base','observed_beta',observed_beta)
%         assignin('base','col',col)
%         assignin('base','thresholds',thresholds)
%         assignin('base','curcol_sorted',curcol_sorted)
%         assignin('base','curcol',curcol)
%         error('a')

        % Compute beta cutoff values and a pvalue map for the observed betas.
        switch parameters.tailshort
            case 'pos' % - high scores bad 
%                thresholds.one_tail_pos_alphas(col) = sum(observed_beta < curcol_sorted)/numel(curcol_sorted); % percent of null values LESS extreme than observed beta - note, not p values
                 [thresholds.one_tail_pos_alphas(col), ...
                     thresholds.pos_beta_map_cutoff(col)] = compare_real_beta(observed_beta,curcol,parameters.tailshort,thresholds.onetail_cutoff_index); % try to preclude p values of 0
%                thresholds.pos_beta_map_cutoff(col) = curcol_sorted(thresholds.pos_thresh_index); % so the 9500th at p of 0.05 on 10,000 permutations
            case 'neg' % - high scores good
%                thresholds.one_tail_neg_alphas(col) = sum(observed_beta > curcol_sorted)/numel(curcol_sorted); % percent of null values LESS extreme than observed beta - note, not p values
                 [thresholds.one_tail_neg_alphas(col), ...
                     thresholds.neg_beta_map_cutoff(col)] = compare_real_beta(observed_beta,curcol,parameters.tailshort,thresholds.onetail_cutoff_index); % try to preclude p values of 0
%                thresholds.neg_beta_map_cutoff(col) = curcol_sorted(thresholds.neg_thresh_index); % so the 500th at p of 0.05 on 10,000 permutations
%             case 'two' % 'two'
%                 warning('Check that these tails are right after code refactor') % ad 2/14/18
%                 thresholds.two_tailed_beta_map_cutoff_pos(col) = curcol_sorted(thresholds.two_tailed_thresh_index); % 250...
%                 thresholds.two_tailed_beta_map_cutoff_neg(col) = curcol_sorted(thresholds.two_tailed_thresh_index_neg); % 9750...
%                 thresholds.twotails_alphas(col) = sum(abs(observed_beta) > abs(curcol_sorted))/numel(curcol_sorted); % percent of values observed_beta is greater than.
        end
    end