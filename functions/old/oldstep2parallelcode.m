%% original parallel code V
%     %% For each voxel, calculate the beta cutoff and p-value based on our permutation data - also, convert to pvalue volumes if CFWER is requested.
%     parfor col = 1 : length(variables.m_idx) % this is for each voxel in the brain, cutting across permutations
%         check_for_interrupt(parameters)
%         curcol = extractSlice(all_perm_data,col,L);
%         observed_beta = ori_beta_vals(col); % original observed beta value.
%         curcol_sorted = sort(curcol); % smallest values at the left/top
%         
%         %% Calculate P values for the single *observed beta map* relative to the null permutation beta volumes.
%         switch tail
%             case 'pos' % high scores bad
%                 one_tail_pos_alphas(col) = sum(observed_beta < curcol_sorted)/nperms;
%                 pos_beta_map_cutoff(col) = curcol_sorted(pos_thresh_index); % so the 9500th at p of 0.05 on 10,000 permutations
%             case 'neg' % high scores good
%                 one_tail_neg_alphas(col) = sum(observed_beta > curcol_sorted)/nperms;
%                 neg_beta_map_cutoff(col) = curcol_sorted(neg_thresh_index); % so the 500th at p of 0.05 on 10,000 permutations
% %            case 'two' % two-tailed...
% %                 warning('Check that these tails are right after code refactor') % ad 2/14/18
% %                 two_tailed_beta_map_cutoff_pos(col) = curcol_sorted(two_tailed_thresh_index); % 250...
% %                 two_tailed_beta_map_cutoff_neg(col) = curcol_sorted(two_tailed_thresh_index_neg); % 9750...
% %                 twotails_alphas(col) = sum(abs(observed_beta) > abs(curcol_sorted))/numel(curcol_sorted); % percent of values observed_beta is greater than.
%         end
%        
%         % Save this permutation as p values.
%         if do_CFWER % then we need to make a billion files and combine them like in parallelized step 1...
%             p_vec = betas2pvals(curcol,tail)
%             % each of these files corresponds to the p values for all the permutations of a single voxel
%             fileID = fopen(fullfile(outpath,['pmu_p_map_' num2str(col) '_of_' num2str(total_cols) '.bin']),'w');
%             fwrite(fileID, p_vec,'single');
%             fclose(fileID);
%         end
%     end