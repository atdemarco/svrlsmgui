function handles = createPvalsForDissociation(handles)
     handles = UpdateProgress(handles,'Dissocation analysis - computing pvals of delta-beta maps...',1);
     if handles.parameters.runfromgui, handles.parameters.waitbar = [handles.progressaxes_rectangle handles.progressaxes_text]; end % so we can actually update the progress bar
     
     types = {'disjunction','conjunction'};
     for T = 1 : numel(types)
         handles = UpdateProgress(handles,['Dissocation analysis - computing pvals of for ' types{T} '...'],1);
         handles = computeThisDissociationType(handles,types{T});
     end
     
function handles = computeThisDissociationType(handles,dissoctype)
    % We pass handles so we can update the wait bar etc...
     dissociation = handles.dissociation; % for convenience
     parameters = handles.parameters; % for convenience - we can use the most recent waitbar, etc.
    
    %% Ok now read back in our delta svrb file to convert to p values - dynamically pull the conjunction or dissociation null permutation data
    all_perm_data = memmapfile(dissociation.([dissoctype '_svrb_file']),'Format','single');

    dataRef = all_perm_data.Data; % will this eliminate some overhead
    L = length(dissociation.m_idx); % for each voxel feature  - this m_idx only refers to voxel features overlapping in both main effect analyses
    for col = 1 : L % For each voxel feature (each svrb-delta)
        if ~mod(col,50) % to reduce num of calls...
            check_for_interrupt(parameters)
            svrlsm_waitbar(parameters.waitbar,col/L)
        end
        
        curcol = dataRef(col:L:end); % index out each voxel column using skips the length of the data ---> this returns curcol, which has the length of npermutations
        observed_beta = dissociation.(['orig_betas_' dissoctype])(col); % original observed beta value - note dynamic field name, are we pulling from dissociation or conjunction?
        
        if parameters.do_CFWER
            if col == 1 % open where we'll put our p-value converted volumetric data...
                pvalfile = ['outfname_big_p_' dissoctype]; % dynamic 
                dissociation.(pvalfile) = fullfile(dissociation.output_folder.base,['pmu_p_maps_N_' num2str(length(dissociation.m_idx)) '_' dissoctype '.bin']); % note dynamic field name and file name
                fileID = fopen(dissociation.(pvalfile),'w'); % note dynamic field name 
            end
            
            p_vec = betas2pvals(curcol,'pos'); %parameters.tailshort);
            fwrite(fileID, p_vec,'single'); % add our results from this voxel to the same giant file....
            
            if col == L, fclose(fileID); end % then it's the last permutation, and our file is written - close the cfwer permutation data output file
        end
        
        % Compute beta cutoff values and a pvalue map for the real observed betas
        onetail_cutoff_index = median([1 round(dissociation.voxp * dissociation.nperms) dissociation.nperms]);
        [one_tail_pos_alphas(col), pos_beta_map_cutoff(col)] = compare_real_beta(observed_beta,curcol,'pos',onetail_cutoff_index);
        [one_tail_neg_alphas(col), neg_beta_map_cutoff(col)] = compare_real_beta(observed_beta,curcol,'neg',onetail_cutoff_index);
    end
    
    handles.dissociation = dissociation; % so we get the pval file names we created....
    handles.dissociation.(['one_tail_pos_alphas_' dissoctype]) = one_tail_pos_alphas; % note dynamic field name
    handles.dissociation.(['pos_beta_map_cutoff_' dissoctype]) = pos_beta_map_cutoff; % note dynamic field name    
    handles.dissociation.(['one_tail_neg_alphas_' dissoctype]) = one_tail_neg_alphas; % note dynamic field name
    handles.dissociation.(['neg_beta_map_cutoff_' dissoctype]) = neg_beta_map_cutoff; % note dynamic field name