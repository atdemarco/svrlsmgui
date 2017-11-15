function variables = run_beta_PMU(parameters, variables, cmd, beta_map,handles)
    options = handles.options;
    zerostemplate = zeros(variables.vo.dim(1:3)); % make a zeros template....
    
    ori_beta_val = beta_map(variables.m_idx).'; % Original observed beta values.
    tic;
    
    %% Generate permutation data
    if parameters.parallelize % try to parfor it...
        handles = UpdateProgress(handles,'Computing beta map permutations (parallelized)...',1);
        sparseLesionData = sparse(variables.lesion_dat);
        lidx = variables.l_idx;
        midx = variables.m_idx;
        betascale = variables.beta_scale;
        
        % create permutations beforehand.
        permdata = nan(numel(variables.one_score),parameters.PermNumVoxelwise); % each COL will be a permutation.
        npermels = size(permdata,1);
        for r = 1 : size(permdata,2) % each col...
            permdata(:,r) = variables.one_score(randperm(npermels));
        end
        
        outpath = variables.output_folder.clusterwise;
        totalperms = parameters.PermNumVoxelwise;
        uselibsvm = parameters.useLibSVM;
        parfor PermIdx=1:parameters.PermNumVoxelwise
            trial_score = permdata(:,PermIdx); % extract the row of permuted data.
            if uselibsvm
                m = svmtrain(trial_score,sparseLesionData,cmd); %#ok<SVMTRAIN>
            else
                [m,~,~] = ComputeMatlabSVRLSM(parameters,variables);
            end
            
            alpha = m.sv_coef.';
            SVs = m.SVs;
            pmu_beta_map = betascale * alpha*SVs;
            tmp_map = zerostemplate; % zeros(nx, ny, nz);
            tmp_map(lidx) = pmu_beta_map;
            pmu_beta_map = tmp_map(midx).';
            
            % Save this permutation....
            fileID = fopen(fullfile(outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(totalperms) '.bin']),'w');
            fwrite(fileID, pmu_beta_map,'single');
            fclose(fileID);
        end
        
        % now get all those individual files into one big file that we can memmap to.
        outfname_big = fullfile(variables.output_folder.clusterwise,['pmu_beta_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);
        fileID = fopen(outfname_big,'w');
        
        for PermIdx=1:parameters.PermNumVoxelwise
            curpermfilepath = fullfile(outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(totalperms) '.bin']);
            cur_perm_data = memmapfile(curpermfilepath,'Format','single');
            fwrite(fileID, cur_perm_data.Data,'single');
            clear cur_perm_data; % remove memmap from memory.
            delete(curpermfilepath); % delete it since we don't want the data hanging around...
        end
        fclose(fileID); % close big file
        
    else % not parallelized
        
        handles = UpdateProgress(handles,'Computing beta map permutations (not parallelized)...',1);

        % This is where we'll save our GBs of permutation data output...
        outfname_big = fullfile(variables.output_folder.clusterwise,['pmu_beta_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);
        fileID = fopen(outfname_big,'w');

        h = waitbar(0,'Computing beta permutations...','Tag','WB');
        for PermIdx=1:parameters.PermNumVoxelwise
            % random permute subjects order
            loc = randperm(length(variables.one_score));
            trial_score = variables.one_score(loc);
            
            % Which package to use to compute SVM solution
            if parameters.useLibSVM
                m = svmtrain(trial_score,sparse(variables.lesion_dat),cmd); %#ok<SVMTRAIN>
            else
                if PermIdx == 1 
                    variables.orig_one_score = variables.one_score;
                end
                
                variables.one_score = trial_score;
                
                [m,~,~] = ComputeMatlabSVRLSM(parameters,variables);
                
                if PermIdx == parameters.PermNumVoxelwise % put it back after we've done all permutations...
                    variables.one_score = variables.orig_one_score;
                end
            end
            
            % compute the beta map
            if parameters.useLibSVM
                alpha = m.sv_coef';
                SVs = m.SVs;
                pmu_beta_map = variables.beta_scale * alpha*SVs;
            else % MATLAB's version.
                alpha = m.Alpha'; 
                SVs = m.SupportVectors;
                pmu_beta_map = variables.beta_scale * alpha*SVs;
            end
            
            
            tmp_map = zerostemplate;
            tmp_map(variables.l_idx) = pmu_beta_map;
            pmu_beta_map = tmp_map(variables.m_idx).';
            
            % Save this permutation....
            fwrite(fileID, pmu_beta_map,'single');
            
            % Display progress.
            elapsed_time = toc;
            remain_time = round(elapsed_time * (parameters.PermNumVoxelwise - PermIdx)/(PermIdx));
            remain_time_h = floor(remain_time/3600);
            remain_time_m = floor((remain_time - remain_time_h*3600)/60);
            remain_time_s = floor(remain_time - remain_time_h*3600 - remain_time_m*60);
            prompt_str = sprintf(['Permutation ', num2str(PermIdx), '/', num2str(parameters.PermNumVoxelwise), ': Est. remaining time: ', num2str(remain_time_h), ' h ', num2str(remain_time_m), ' m ' num2str(remain_time_s), 's\n']);
            waitbar(PermIdx/parameters.PermNumVoxelwise,h,prompt_str) % show progress.
            
        end
        
        close(h) % close the waitbar...
        fclose(fileID); % close the pmu data output file.
    end
    
    %% Read in gigantic memory mapped file... no matter whether we parallelized or not.
    all_perm_data = memmapfile(outfname_big,'Format','single'); % does the single precision hurt the analysis?

    %% FWE cluster correction based on permutation analysis
    voxelwise_p_value = parameters.voxelwise_p;
    pos_thresh_index = median([1 round((1-voxelwise_p_value) * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % row 9500 in 10000 permutations.
    pos_beta_map_cutoff = nan(1,length(variables.m_idx));
    one_tail_pos_alphas = nan(1,length(variables.m_idx));

    neg_thresh_index = median([1 round(voxelwise_p_value * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % so row 500 in 10000 permutations
    neg_beta_map_cutoff = nan(1,length(variables.m_idx));
    one_tail_neg_alphas = nan(1,length(variables.m_idx));

    two_tailed_thresh_index_neg = median([1 round(((voxelwise_p_value/2)) * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % row 250 in 10000 permutations.
    two_tailed_thresh_index = median([1 round((1-(voxelwise_p_value/2)) * parameters.PermNumVoxelwise) parameters.PermNumVoxelwise]); % row 9750 in 10000 permutations.
    two_tailed_beta_map_cutoff_pos = nan(1,length(variables.m_idx));
    two_tailed_beta_map_cutoff_neg = nan(1,length(variables.m_idx));
    twotails_alphas = nan(1,length(variables.m_idx));

if parameters.parallelize % try to parfor it...
    handles = UpdateProgress(handles,'Sorting null betas for each lesioned voxel in the dataset (parallelized).',1);
    L = length(variables.m_idx);
    tails = parameters.tails; % so not a broadcast variable.
    Opt1 = options.hypodirection{1};
    Opt2 = options.hypodirection{2};
    Opt3 = options.hypodirection{3};
    parfor col = 1 : length(variables.m_idx) 
        curcol = extractSlice(all_perm_data,col,L); % note this is a function at the bottom of this file..
        observed_beta = ori_beta_val(col); % original observed beta value.
        curcol_sorted = sort(curcol); % smallest values at the top..
        switch tails
            case Opt1
                one_tail_pos_alphas(col) = sum(observed_beta > curcol_sorted)/numel(curcol_sorted); % percent of values observed_beta is greater than.
                pos_beta_map_cutoff(col) = curcol_sorted(pos_thresh_index); % so the 9500th at p of 0.05 on 10,000 permutations
            case Opt2
                one_tail_neg_alphas(col) = sum(observed_beta < curcol_sorted)/numel(curcol_sorted); % percent of values observed_beta is greater than.
                neg_beta_map_cutoff(col) = curcol_sorted(neg_thresh_index); % so the 500th at p of 0.05 on 10,000 permutations
            case Opt3
                two_tailed_beta_map_cutoff_pos(col) = curcol_sorted(two_tailed_thresh_index); % 250...
                two_tailed_beta_map_cutoff_neg(col) = curcol_sorted(two_tailed_thresh_index_neg); % 9750...
                twotails_alphas(col) = sum(abs(observed_beta) > abs(curcol_sorted))/numel(curcol_sorted); % percent of values observed_beta is greater than.
        end
    end
else % otherwise, not parallelized
    
    doNewFWE = false; % concept from Mirman
    
    if doNewFWE
        outfname_big_pvals = fullfile(variables.output_folder.clusterwise,['pmu_beta_map_p_vals_N_' num2str(parameters.PermNumVoxelwise) '.bin']);
        pvalfileID = fopen(outfname_big_pvals,'w');
    end
    
    
    handles = UpdateProgress(handles,'Sorting null betas for each lesioned voxel in the dataset (not parallelized).',1);
    h = waitbar(0,sprintf('Sorting null betas for each lesioned voxel in the dataset (N = %d).\n',length(variables.m_idx)),'Tag','WB');
    dataRef = all_perm_data.Data; % will this eliminate some overhead  
    L = length(variables.m_idx);
    for col = 1 : length(variables.m_idx) 
        curcol = dataRef(col:L:end); % index out each column using skips the length of the data...
        observed_beta = ori_beta_val(col); % original observed beta value.
        curcol_sorted = sort(curcol); % smallest values at the top..
        
        if doNewFWE
            %             inputs: NNN and ALPHA  = 0.05
            %             FWER:
            %             1. generate 10,000 beta maps
            %             2. at each voxel, use leave-one-out to convert each into a p-value. now we have 10,000 p-maps.
            %             3. for each of the 10,000 p-maps, take the NNNth highest p-value.
            %             4. then take the p-value corresponding to the ALPHAth quantile of the histogram of NNNth highest pvalues
            %             5. this is a map-wide voxelwise threshold that does correction at both levels
            
            % in this loop we are going to do step 2, and create e.g. 10,000 p-vals
            
            p_vec = nan(size(curcol));
            all_ind = 1 : numel(curcol); % we'll reuse this vector
            for i = all_ind % for each svr beta value in the vector
                ind_to_compare = setdiff(all_ind,i);
                p_vec(i) = 1 - mean(curcol_sorted(i) < curcol_sorted(ind_to_compare)); % NB: IS THIS IS TAIL DEPENDENT!!!!!! DOES IT NEED TO BE DIFFERENT FOR ONE POS, ONE NEG, and TWOTAIL?
            end
            fwrite(pvalfileID, p_vec,'single');
        end
        
        % Compute beta cutoff values and a pvalue map for the observed betas.
        switch parameters.tails
            case options.hypodirection{1} % 'one_positive'
                one_tail_pos_alphas(col) = sum(observed_beta > curcol_sorted)/numel(curcol_sorted); % percent of values observed_beta is greater than.
                pos_beta_map_cutoff(col) = curcol_sorted(pos_thresh_index); % so the 9500th at p of 0.05 on 10,000 permutations
            case options.hypodirection{2} %'one_negative'
                one_tail_neg_alphas(col) = sum(observed_beta < curcol_sorted)/numel(curcol_sorted); % percent of values observed_beta is greater than.
                neg_beta_map_cutoff(col) = curcol_sorted(neg_thresh_index); % so the 500th at p of 0.05 on 10,000 permutations
            case options.hypodirection{3} % 'two'
                two_tailed_beta_map_cutoff_pos(col) = curcol_sorted(two_tailed_thresh_index); % 250...
                two_tailed_beta_map_cutoff_neg(col) = curcol_sorted(two_tailed_thresh_index_neg); % 9750...
                twotails_alphas(col) = sum(abs(observed_beta) > abs(curcol_sorted))/numel(curcol_sorted); % percent of values observed_beta is greater than.
        end        
        waitbar(col/L,h) % show progress.
    end
    close(h)
    if doNewFWE
        fclose(pvalfileID); % close the pmu pval data output file.
    end
end
    %% Construct volumes of the solved alpha values and write them out - and write out beta cutoff maps, too
    switch parameters.tails
        case options.hypodirection{1} % 'one_positive' % One-tailed positive tail...
            if parameters.invert_p_map_flag % it's already inverted
                variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
            else
                one_tail_pos_alphas = 1 - one_tail_pos_alphas;
                variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values.nii');
            end
            thresholded_pos = zerostemplate; % zeros(nx,ny,nz); % reserve space;
            thresholded_pos(variables.m_idx) = one_tail_pos_alphas;
            spm_write_vol(variables.vo, thresholded_pos);
            
            % Now write out beta cutoff map.
            thresholded_pos = zerostemplate; % zeros(nx,ny,nz); % reserve space
            thresholded_pos(variables.m_idx) = pos_beta_map_cutoff; % put the 95th percentil beta values back into the lesion indices in a full volume
            variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (positive tail).nii');
            spm_write_vol(variables.vo, thresholded_pos);
            
        case options.hypodirection{2} % 'one_negative' % One-tailed negative tail...
            
            if parameters.invert_p_map_flag  % it's already inverted...
                variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
            else
                one_tail_neg_alphas = 1 - one_tail_neg_alphas;
                variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values.nii');
            end
            thresholded_neg = zerostemplate; % zeros(nx,ny,nz); % reserve space;
            thresholded_neg(variables.m_idx) = one_tail_neg_alphas;
            spm_write_vol(variables.vo, thresholded_neg);
            
            % Now beta cutoff map.
            thresholded_neg = zerostemplate; % zeros(nx,ny,nz); % reserve space;
            thresholded_neg(variables.m_idx) = neg_beta_map_cutoff; % put the 5th percentil beta values back into the lesion indices in a full volume
            variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (negative tail).nii');
            spm_write_vol(variables.vo, thresholded_neg);
            
        case options.hypodirection{3} %'two'% Both tails..
            
            if parameters.invert_p_map_flag % it's already inverted...
                variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values (inv).nii');
            else
                twotails_alphas = 1 - twotails_alphas;
                variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Thresholded P values.nii');
            end
            
            thresholded_twotails = zerostemplate; % zeros(nx,ny,nz); % reserve space;
            thresholded_twotails(variables.m_idx) = twotails_alphas;
            spm_write_vol(variables.vo, thresholded_twotails);
            
            % Bet cutoff maps
            % Two-tailed upper tail
            thresholded_twotail_upper = zerostemplate; %zeros(nx,ny,nz); % reserve space;
            thresholded_twotail_upper(variables.m_idx) = two_tailed_beta_map_cutoff_pos; % put the 2.5th percentil beta values back into the lesion indices in a full volume
            variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (two tail, upper).nii');
            spm_write_vol(variables.vo, thresholded_twotail_upper);
            
            % Two-tailed lower tail
            thresholded_twotail_lower = zerostemplate; % zeros(nx,ny,nz); % reserve space;
            thresholded_twotail_lower(variables.m_idx) = two_tailed_beta_map_cutoff_neg; % put the 2.5th percentil beta values back into the lesion indices in a full volume
            variables.vo.fname = fullfile(variables.output_folder.voxelwise,'Beta value cutoff mask (two tail, lower).nii');
            spm_write_vol(variables.vo, thresholded_twotail_lower);
    end
    
    if doNewFWE % then let's traverse the pvalue frames we've created and get the Nth highest (Nth_highest_pval) and construct a a histogram from that.
        Nth_highest_pval = 100;
        all_perm_data_pvals = memmapfile(outfname_big_pvals,'Format','single');
        handles = UpdateProgress(handles,'Traversing null data p-value frames and noting Nth largest p-value...',1);

        for f = 1 : parameters.PermNumVoxelwise % go through each frame of generated betas in the null data...
            waitbar(f/parameters.PermNumVoxelwise,h) % show progress.
            %frame_length = length(variables.m_idx);
            %frame_start_index = 1+((f-1)*frame_length); % +1 since not zero indexing
            %frame_end_index = (frame_start_index-1)+frame_length; % -1 so we are not 1 too long.
            %relevant_data_frame = all_perm_data.Data(frame_start_index:frame_end_index); % extract the frame
            
            relevant_data_frame = parameters.PermNumVoxelwise(f:parameters.PermNumVoxelwise:end); % need to skip every e.g., 10,000 since this is organized different than the big pmu binary file
            
            %templatevol = zerostemplate;
            %templatevol(variables.m_idx) = relevant_data_frame; % put the beta values back in indices.

%              if parameters.SavePreThresholdedPermutations % then write out raw voxel NON-thresholded images for this permutation.
%                 variables.vo.fname = fullfile(variables.output_folder.clusterwise,['UNthreshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
%                 spm_write_vol(variables.vo, templatevol);
%              end

            switch parameters.tails
                case options.hypodirection{1} 
                    sorted_data_frame = sort(relevant_data_frame);
                    curmaxpvalue = sorted_data_frame(end-Nth_highest_pval);
                case options.hypodirection{2} 
                    sorted_data_frame = sort(relevant_data_frame);
                    curmaxpvalue = sorted_data_frame(Nth_highest_pval); %NB: NOT END-
                case options.hypodirection{3} % Now build a two-tailed thresholded version
                    sorted_data_frame = sort(relevant_data_frame);
                    curmaxpvalue = sorted_data_frame(end-Nth_highest_pval); % WHAT TO DO HERE???????
            end
            
            all_max_pvalue(f) = curmaxpvalue;

%             if parameters.SavePostVoxelwiseThresholdedPermutations % then write out raw voxel thresholded images for this permutation.
%                 switch parameters.tails
%                     case options.hypodirection{1} % 'one_positive'
%                         variables.vo.fname = fullfile(variables.output_folder.clusterwise,['pos_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
%                         spm_write_vol(variables.vo, pos_threshed);
%                     case options.hypodirection{2} %'one_negative'
%                         variables.vo.fname = fullfile(variables.output_folder.clusterwise,['neg_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
%                         spm_write_vol(variables.vo, neg_threshed);
%                     case options.hypodirection{3} %'two'
%                         variables.vo.fname = fullfile(variables.output_folder.clusterwise,['twotail_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
%                         spm_write_vol(variables.vo, twotail_threshed);
%                 end
%             end

%             if f <= parameters.PermNumClusterwise
%                 switch parameters.tails
%                     case options.hypodirection{1}
%                         permtype='pos';
%                         thresholded_mask=pos_threshed;
%                     case options.hypodirection{2}
%                         permtype='neg';
%                         thresholded_mask=neg_threshed;
%                     case options.hypodirection{3}
%                         permtype='twotail';
%                         thresholded_mask=twotail_threshed;
%                 end
% 
%                 testvol_thresholded = thresholded_mask; % now evaluate the surviving voxels for clusters...
%                 CC = bwconncomp(testvol_thresholded, 6);
%                 largest_cluster_size = max(cellfun(@numel,CC.PixelIdxList(1,:))); % max val for numels in each cluster object found
%                 if isempty(largest_cluster_size)
%                     largest_cluster_size = 0;
%                 else % threshold the volume and write it out.
%                     if parameters.SavePostClusterwiseThresholdedPermutations % then save them...
%                         out_map = remove_scatter_clusters(testvol_thresholded, largest_cluster_size-1);
%                         variables.vo.fname = fullfile(variables.output_folder.clusterwise,[permtype '_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumClusterwise) '_largest_cluster.nii']);
%                         spm_write_vol(variables.vo, out_map);
%                     end
%                 end
%                 all_max_cluster_sizes(f) = largest_cluster_size; % record...
%             end
        end
        sorted_all_max_pvalue = sort(all_max_pvalue);
        disp('final map-wide pvalue:')
        sorted_all_max_pvalue(pos_thresh_index) % or neg thresh ..... 
        
    else % do regular cluster correction.
        % Now for each permuted beta map, apply the beta mask and determine largest surviving cluster.
        h = waitbar(0,'Applying voxelwise beta mask to null data and noting largest null clusters...','Tag','WB');

        % Reconstruct the volumes so we can threshold and examine cluster sizes
        all_max_cluster_sizes = nan(parameters.PermNumClusterwise,1); % reserve space.

        if parameters.PermNumClusterwise > parameters.PermNumVoxelwise, error('Cannot sample more cluster permutations than have been generated by the voxelwise permutation procedure.'); end

        handles = UpdateProgress(handles,'Applying voxelwise beta mask to null data and noting largest null clusters...',1);

        for f = 1 : parameters.PermNumVoxelwise % go through each frame of generated betas in the null data...
            waitbar(f/parameters.PermNumVoxelwise,h) % show progress.

            frame_length = length(variables.m_idx);
            frame_start_index = 1+((f-1)*frame_length); % +1 since not zero indexing
            frame_end_index = (frame_start_index-1)+frame_length; % -1 so we are not 1 too long.
            relevant_data_frame = all_perm_data.Data(frame_start_index:frame_end_index); % extract the frame
            templatevol = zerostemplate;
            templatevol(variables.m_idx) = relevant_data_frame; % put the beta values back in indices.

             if parameters.SavePreThresholdedPermutations % then write out raw voxel NON-thresholded images for this permutation.
                variables.vo.fname = fullfile(variables.output_folder.clusterwise,['UNthreshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
                spm_write_vol(variables.vo, templatevol);
             end

            switch parameters.tails
                case options.hypodirection{1} 
                    pos_threshed = templatevol .* (templatevol>=thresholded_pos); % elementwise greater than operator to threshold positive tail of test betas
                case options.hypodirection{2} 
                    neg_threshed = templatevol .* (templatevol<=thresholded_neg); % elementwise less than operator to threshold negative tail of test betas.
                case options.hypodirection{3} % Now build a two-tailed thresholded version
                threshmask = and(templatevol > 0,templatevol >= thresholded_twotail_upper) | and(templatevol < 0,templatevol <= thresholded_twotail_lower);
                twotail_threshed = templatevol .* threshmask; % mask the two-tailed beta mask for this null data...
            end

            if parameters.SavePostVoxelwiseThresholdedPermutations % then write out raw voxel thresholded images for this permutation.
                switch parameters.tails
                    case options.hypodirection{1} % 'one_positive'
                        variables.vo.fname = fullfile(variables.output_folder.clusterwise,['pos_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
                        spm_write_vol(variables.vo, pos_threshed);
                    case options.hypodirection{2} %'one_negative'
                        variables.vo.fname = fullfile(variables.output_folder.clusterwise,['neg_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
                        spm_write_vol(variables.vo, neg_threshed);
                    case options.hypodirection{3} %'two'
                        variables.vo.fname = fullfile(variables.output_folder.clusterwise,['twotail_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumVoxelwise) '.nii']);
                        spm_write_vol(variables.vo, twotail_threshed);
                end
            end

            if f <= parameters.PermNumClusterwise
                switch parameters.tails
                    case options.hypodirection{1}
                        permtype='pos';
                        thresholded_mask=pos_threshed;
                    case options.hypodirection{2}
                        permtype='neg';
                        thresholded_mask=neg_threshed;
                    case options.hypodirection{3}
                        permtype='twotail';
                        thresholded_mask=twotail_threshed;
                end

                testvol_thresholded = thresholded_mask; % now evaluate the surviving voxels for clusters...
                CC = bwconncomp(testvol_thresholded, 6);
                largest_cluster_size = max(cellfun(@numel,CC.PixelIdxList(1,:))); % max val for numels in each cluster object found
                if isempty(largest_cluster_size)
                    largest_cluster_size = 0;
                else % threshold the volume and write it out.
                    if parameters.SavePostClusterwiseThresholdedPermutations % then save them...
                        out_map = remove_scatter_clusters(testvol_thresholded, largest_cluster_size-1);
                        variables.vo.fname = fullfile(variables.output_folder.clusterwise,[permtype '_threshed_perm_' num2str(f) '_of_' num2str(parameters.PermNumClusterwise) '_largest_cluster.nii']);
                        spm_write_vol(variables.vo, out_map);
                    end
                end
                all_max_cluster_sizes(f) = largest_cluster_size; % record...
            end
        end
        close(h)

        % Save the resulting cluster lists
        fname = 'Largest clusters.mat';
        save(fullfile(variables.output_folder.clusterwise,fname),'all_max_cluster_sizes')

        handles = UpdateProgress(handles,'Cleaning up null data...',1);

        % Clean up as necessary
        if ~parameters.SavePermutationData
            delete(outfname_big) % delete the monster bin file with raw permutation data in it.
        end
    end % end doNewFWE
    
% for parallelization to eliminate large overhead transfering to and from workers
function sliceData = extractSlice(all_perm_data,col,L)
    sliceData = all_perm_data.Data(col:L:end);