function [success,handles] = runDissociation(hObject,eventdata,handles)
    success = 1; % by default
    %% First, run the two main analyses we'll draw on
    handles.dissociation = [];
    handles.dissociation.starttime = now;

%     try
        handles = runMainDissociationAnalyses(hObject,eventdata,handles);

        %% Construct some general info we'll use to do the analysis...
        first_analysis = handles.dissociation.maineffects{1};
        handles.dissociation.output_folder.base = fullfile(first_analysis.baseoutputdir,'dissociation');
        handles.dissociation.voxp = first_analysis.voxelwise_p;
        handles.dissociation.nperms = first_analysis.PermNumVoxelwise;
        handles.dissociation.vo = first_analysis.vo; % we reuse this for a template
        handles.dissociation.cfwer_p_value = first_analysis.cfwer_p_value;
        handles.dissociation.cfwer_v_value = first_analysis.cfwer_v_value;

        %% Now compute the combined beta maps to compare against to determine conjunction/disjunction
        handles = createBetaDeltaMapForDissociation(handles); % create new beta map (deltas)

        %% Now convert those beta maps into p maps for all permutatoins
        handles = createPvalsForDissociation(handles); % convert those betas into p values

        %% We run this once for disjunction and once for conjunction
        dissoctype = {'conjunction','disjunction'};
        for d = 1 : numel(dissoctype)
            handles.dissociation.current_dissoctype = dissoctype{d};
            % disp(['Saving results for ' handles.dissociation.current_dissoctype])
            handles = computeClusterResultsDissociation(handles);
        end

        %% OK, make the summary output file
        tosave = [];
        tosave.variables = handles.variables;
        tosave.dissociation = handles.dissociation;
        tosave.parmsfile = fullfile(handles.dissociation.output_folder.base,'Dissociation Parameters');% fullfile('Dissociation Parameters');

        handles.parameters.parmsfile = tosave.parmsfile;
        handles.parameters.output_folders = handles.variables.output_folder;
        save(tosave.parmsfile,'tosave')

        WriteDissociationSummary(tosave.parmsfile)

        %     %% cleanup
        %     [handles,parameters,variables] = cleanup(handles,parameters,variables);
%     catch
%         success = 0;   
%     end
    
function handles = computeClusterResultsDissociation(handles)
    %% Run this for each tail - pos and neg...
    tails = {'pos','neg'};
    for T = 1 : numel(tails)
        if strcmp(tails{T},'pos') && strcmp(handles.dissociation.current_dissoctype,'conjunction')
            continue; % skip the positive tail of the conjunction, since that would identify areas involved in *NEITHER* behavior...
        end 
        handles.parameters.tailshort = tails{T};
        handles = writeEachTailOut(handles);
    end

function handles = writeEachTailOut(handles)
    %% Create and fill out the thresholds struct (we use this for applying thresholding and clustering
    thresholds = calculate_thresholds(handles.parameters,handles.dissociation); % these should pass all the info we need to create the 'thresholds' struct
    dissoctype = handles.dissociation.current_dissoctype;
    
    % Transfer this information from the dissociation struct in handles to the thresholds struct 
    fnames = {'one_tail_pos_alphas','one_tail_neg_alphas','pos_beta_map_cutoff','neg_beta_map_cutoff'};
    for f = 1 : numel(fnames) % set these all depending on the type of dissocation...
        thresholds.(fnames{f}) = handles.dissociation.([fnames{f} '_' dissoctype]);
    end

    %% Conditionally set the output folders... by tail, and dissociation type.
    % note here we suffix with the tail direction as in --- [folderstring '_' handles.parameters.tailshort]
    folderstring = sprintf('Voxwise p%s (%d perms)',strrep(num2str(handles.parameters.voxelwise_p),'0.',''),handles.parameters.PermNumVoxelwise);
    handles.dissociation.output_folder.voxelwise = fullfile(handles.dissociation.output_folder.base,dissoctype,[folderstring '_' handles.parameters.tailshort]);
    
    folderstring = sprintf('Clustwise p%s (%d perms)',strrep(num2str(handles.parameters.clusterwise_p),'0.',''),handles.parameters.PermNumClusterwise);
    handles.dissociation.output_folder.clusterwise = fullfile(handles.dissociation.output_folder.voxelwise,[folderstring '_' handles.parameters.tailshort]);
    
    folderstring = sprintf('%s CFWER p%s at v%d (%d perms)',handles.dissociation.current_dissoctype,strrep(num2str(handles.parameters.cfwer_p_value),'0.',''),handles.parameters.cfwer_v_value,handles.parameters.PermNumClusterwise);
    handles.dissociation.output_folder.cfwer = fullfile(handles.dissociation.output_folder.base,dissoctype,[folderstring '_' handles.parameters.tailshort]);
    
    
    %% Now apply the multiple comparisons correction (cfwer or, maybe, regular cluster extent)
    if handles.parameters.do_CFWER
        % Overwrite voxelwise and clusterwise directory names used for regular non-cfwer output
        handles.dissociation.output_folder.voxelwise = handles.dissociation.output_folder.cfwer;
        handles.dissociation.output_folder.clusterwise = handles.dissociation.output_folder.cfwer;
        
        % This function is equivalent to the single get_cfwer_dist() for dissociations
        handles = makeCfwerDistributionForDissociation(handles); % Construct the info we'll need to do the cfwer thresholding using the regular machinery

        %% set relevant files on the basis of the current dissociation type
        handles.dissociation.cfwerinfo = handles.dissociation.(['cfwerinfo_' handles.dissociation.current_dissoctype]);
        handles.dissociation.files_created.largest_clusters = handles.dissociation.files_created.(['largest_clusters_' handles.dissociation.current_dissoctype]);

        options = []; % not used..
        [thresholded,variables] = build_and_write_pmaps(options,handles.parameters,handles.dissociation,thresholds); %#ok<ASGLU>
        variables = do_cfwer_clustering([],handles.parameters,variables,[],[]);

        %% Let's store files_created.unthresholded_betamap like when we do cluster-extent thresholding, so that we can refer to it in the output.
        if strcmp(dissoctype,'disjunction'), unthresholded_betamap = fullfile(handles.dissociation.output_folder.base,'A_disjunct_B_svrb.nii');
        else, unthresholded_betamap = fullfile(handles.dissociation.output_folder.base,'A_conjunct_B_svrb.nii');
        end
        variables.files_created.unthresholded_betamap = unthresholded_betamap;
    else
        %% Do regular cluster extent thresholding...
        mkdir(handles.dissociation.output_folder.voxelwise)
        mkdir(handles.dissociation.output_folder.clusterwise)
        %% Construct volumes of the solved p values and write them out - and write out beta cutoff maps, too
        options = []; % not used..
        [thresholded,variables] = build_and_write_pmaps(options,handles.parameters,handles.dissociation,thresholds);
        [thresholded,variables] = build_and_write_beta_cutoffs(options,handles.parameters,variables,thresholds,thresholded);
        
        all_perm_data = memmapfile(handles.dissociation.([dissoctype '_svrb_file']),'Format','single'); % Get a handle to the file with all null svrb data
        variables = do_cluster_thresholding_of_permutations(handles,handles.parameters,variables,all_perm_data,thresholded);
       
        if strcmp(dissoctype,'disjunction'), unthresholded_betamap = fullfile(handles.dissociation.output_folder.base,'A_disjunct_B_svrb.nii');
        else, unthresholded_betamap = fullfile(handles.dissociation.output_folder.base,'A_conjunct_B_svrb.nii');
        end
        variables.files_created.unthresholded_betamap = unthresholded_betamap;
        
        variables = evaluate_clustering_results(handles,variables,handles.parameters);
    end
    handles.variables = variables; % so we can save in the output (we only need one version...)