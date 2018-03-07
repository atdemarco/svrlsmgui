function [parameters,variables,thresholds] = step2_parallel(handles,parameters,variables,thresholds,all_perm_data)
    handles = UpdateProgress(handles,'Sorting null betas for each lesioned voxel in the dataset (parallelized).',1);
    L = length(variables.m_idx);
    tail = parameters.tailshort; % so not a broadcast variable.

    ori_beta_vals = variables.ori_beta_vals; % for parfor...
    do_CFWER = parameters.do_CFWER; % for parfor...
    total_cols =  length(variables.m_idx); % note this must be m_idx since the data in our giant file is stored in frames of length m_idx not length l_idx.
    nperms = parameters.PermNumVoxelwise;
    outpath = variables.output_folder.clusterwise;

%     pos_thresh_index = thresholds.pos_thresh_index;
%     neg_thresh_index = thresholds.neg_thresh_index;
      onetail_thresh_index = thresholds.onetail_cutoff_index; % for use with compare_real_beta() 
    
    % two_tailed_thresh_index = thresholds.two_tailed_thresh_index;
    % two_tailed_thresh_index_neg = thresholds.two_tailed_thresh_index_neg;

    %% Begin parfeval code 
    batch_job_size = 500; % this is going to be optimal for different systems/#cores/jobs
    njobs = ceil(total_cols/batch_job_size); % gotta round up to capture all indices

    %% Schedule the jobs...
    p = gcp(); % get current parallel pool
    for j = 1 : njobs
        this_job_start_index = ((j-1)*batch_job_size) + 1;
        this_job_end_index = min(this_job_start_index + batch_job_size-1,total_cols); % need min so we don't go past valid indices
        job_indices = this_job_start_index:this_job_end_index;
        f(j) = parfeval(p,@parallel_step2_batch_fcn,2,job_indices,all_perm_data,total_cols,tail,ori_beta_vals,onetail_thresh_index,L,nperms,do_CFWER,outpath);
    end
    
    alphas = cell(1,njobs); %reserve space - note we want to accumulate in a row here
    betamapcutoff = cell(1,njobs); %reserve space - note we want to accumulate in a row here
    
    %% Monitor job progress...
    msg = 'Sorting null betas for each lesioned voxel in the dataset (parallelized).';
    svrlsm_waitbar(parameters.waitbar,0,msg) % update waitbar progress...
    for j = 1 : njobs
        check_for_interrupt(parameters) % allow user to interrupt
        [idx, jobalphas,jobbetamapcutoffs] = fetchNext(f);
        alphas{idx} = jobalphas; % combine these cells afterward
        betamapcutoff{idx} = jobbetamapcutoffs; % combine these cells afterward
        svrlsm_waitbar(parameters.waitbar,j/njobs) % update waitbar progress...
    end

    alphas = cell2mat(alphas); % combine afterward
    betamapcutoff = cell2mat(betamapcutoff); % combine afterward

    %% Now compute beta cutoff values and a pvalue map for the observed betas.
    % ...we do this now since we can't index into fields of the 'thresholds'variable in a parfor loop
    switch tail
        case 'pos' % pos - high scores bad
            thresholds.one_tail_pos_alphas = alphas;
            thresholds.pos_beta_map_cutoff = betamapcutoff;
        case 'neg' % neg - high scores good
            thresholds.one_tail_neg_alphas = alphas;
            thresholds.neg_beta_map_cutoff = betamapcutoff;
%         case 'two' % two-tailed
%             thresholds.two_tailed_beta_map_cutoff_pos =two_tailed_beta_map_cutoff_pos;
%             thresholds.two_tailed_beta_map_cutoff_neg =two_tailed_beta_map_cutoff_neg;
%             thresholds.twotails_alphas = twotails_alphas;
    end
    
    %% Now, since we're parallelized, like in step 1, get all those individual files into one big file that we can memmap
    if do_CFWER    
        svrlsm_waitbar(parameters.waitbar,0,'Consolidating null p-map files...');
        parameters.outfname_big_p = fullfile(variables.output_folder.clusterwise,['pmu_p_maps_N_' num2str(total_cols) '.bin']);
        fileID = fopen(parameters.outfname_big_p,'w');
        % this data is output such that each numel(variables.m_idx) contiguous sequential values refer to a voxel across all permutations:  [P1V1, P2V1, P3V1, P4V1, ..., npermutations]
        for col = 1 : total_cols % can't parallelize this since we need order to be right.
            if ~mod(500,col) % to reduce num of calls...
                check_for_interrupt(parameters)
                svrlsm_waitbar(parameters.waitbar,col/total_cols);
            end
            curpermfilepath = fullfile(outpath,['pmu_p_map_' num2str(col) '_of_' num2str(total_cols) '.bin']);
            cur_perm_data = memmapfile(curpermfilepath,'Format','single');
            fwrite(fileID, cur_perm_data.Data,'single');
            clear cur_perm_data; % remove memmap from memory.
            delete(curpermfilepath); % delete it since we don't want the data hanging around...
        end
        fclose(fileID); % close big file
        svrlsm_waitbar(parameters.waitbar,0,'');
    end
    
%% For each voxel, calculate the beta cutoff and p-value based on our permutation data - also, convert to pvalue volumes if CFWER is requested.
function [alphas,betamapcutoffs] = parallel_step2_batch_fcn(this_job_cols,all_perm_data,total_cols,tail,ori_beta_vals,onetail_thresh_index,L,nperms,do_CFWER,outpath)
    alphas = nan(1,numel(this_job_cols)); % pre-allocate spac
    betamapcutoffs = nan(1,numel(this_job_cols)); % pre-allocate spac
    for jobcolind = 1:numel(this_job_cols) % this is for each voxel in the brain, cutting across permutations
        col = this_job_cols(jobcolind);
        curcol = extractSlice(all_perm_data,col,L);
        observed_beta = ori_beta_vals(col); % original observed beta value.
        %[curcol_sorted] = sort(curcol); % Smallest values at the left/top
        
        %% Calculate P values for the single *observed beta map* relative to the null permutation beta volumes.
        switch tail
            case {'pos','neg'}
                [alphas(jobcolind),betamapcutoffs(jobcolind)] = compare_real_beta(observed_beta,curcol,tail,onetail_thresh_index); % we can do both tails without the switch as long as it's one-tailed.
%             case 'pos' % high scores bad
%                 % 'ascend' is the default assort behavior.
%                 %alphas(jobcolind) = sum(observed_beta < curcol_sorted)/nperms;
%                 [alphas(jobcolind),betamapcutoffs(jobcolind)] = compare_real_beta(observed_beta,curcol,tail,onetail_thresh_index);
%                % betamapcutoffs(jobcolind) = curcol_sorted(pos_thresh_index); % so the 9500th at p of 0.05 on 10,000 permutations
%             case 'neg' % high scores good
%                 %alphas(jobcolind) = sum(observed_beta > curcol_sorted)/nperms;
%                 [alphas(jobcolind),betamapcutoffs(jobcolind)] = compare_real_beta(observed_beta,curcol,tail,onetail_thresh_index);
%                % betamapcutoffs(jobcolind) = curcol_sorted(neg_thresh_index); % so the 500th at p of 0.05 on 10,000 permutations
            case 'two' % two-tailed...
                warning('Check that these tails are right after code refactor') % ad 2/14/18
% %                 two_tailed_beta_map_cutoff_pos(col) = curcol_sorted(two_tailed_thresh_index); % 250...
% %                 two_tailed_beta_map_cutoff_neg(col) = curcol_sorted(two_tailed_thresh_index_neg); % 9750...
% %                 twotails_alphas(col) = sum(abs(observed_beta) > abs(curcol_sorted))/numel(curcol_sorted); % percent of values observed_beta is greater than.
        end
       
        % Save this permutation as p values.
        if do_CFWER % then we need to make a billion files and combine them like in parallelized step 1...
            p_vec = betas2pvals(curcol,tail); % each of these files corresponds to the p values for all the permutations of a single voxel
            fileID = fopen(fullfile(outpath,['pmu_p_map_' num2str(col) '_of_' num2str(total_cols) '.bin']),'w');
            fwrite(fileID, p_vec,'single');
            fclose(fileID);
        end
    end