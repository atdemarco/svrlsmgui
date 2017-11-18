function handles = PopulateGUIFromParameters(handles)
    [~,path_to_show] = fileparts(handles.parameters.lesion_img_folder); % only show final directory.
    set(handles.lesionfoldereditbox,'string',path_to_show,'tooltipstring',handles.parameters.lesion_img_folder)

    [~,file,ext] = fileparts(handles.parameters.score_file); % only show final directory.
    set(handles.scorefileeditbox,'string',[file ext],'tooltipstring',handles.parameters.score_file)

    [~,path_to_show] = fileparts(handles.parameters.analysis_out_path); % only show final directory.
    set(handles.outputfoldereditbox,'string',path_to_show,'tooltipstring',handles.parameters.analysis_out_path)

    warning('off', 'MATLAB:table:ModifiedVarnames');
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames') % MATLAB 2017a
    
    %% Read in the score file.
    handles.scorefiledata = readtable(handles.parameters.score_file); 

    %% Find the column for lesion files: 'RegistryCode'
    handles.scorefile.subjectrows = find(~cellfun(@isempty,handles.scorefiledata.RegistryCode));
    handles.scorefile.nsubs_in_scorefile = numel(handles.scorefile.subjectrows);

    check_variables = [handles.parameters.score_name handles.parameters.control_variable_names];

    % add - don't show registry code as one_score by default.
    
    for c = 1 : numel(check_variables)
        if ~any(strcmp(check_variables{c},handles.scorefiledata.Properties.VariableNames)) % then variable doesn't exist in file
            if c == 1
                handles.parameters.score_name = handles.scorefiledata.Properties.VariableNames{1}; % replace with the First col in this file...
            else
                handles.parameters.control_variable_names(c-1) = []; % remove from the list...
            end
        end
    end

    check_variables = [handles.parameters.score_name handles.parameters.control_variable_names];
    track_variables = zeros(handles.scorefile.nsubs_in_scorefile,numel(check_variables));

    %% Try to discover who has the currently specified behaviors necessary for the currently specified analysis    
    for c = 1 : numel(check_variables)
        coldata = handles.scorefiledata.(check_variables{c})(handles.scorefile.subjectrows);
        if isa(coldata,'cell')
            track_variables(:,c) = ~cellfun(@isempty,coldata);
        else % it's numbers...
            track_variables(:,c) = ~isnan(coldata);
        end
    end
    num_subs_that_have_all_variables = sum(all(track_variables,2));
    
    %% Try to discover the actual lesion files specified in the design file
    handles.scorefile.haslesion_file = zeros(numel(handles.scorefile.subjectrows),1);
    for p = handles.scorefile.subjectrows'
        lesion_found = 0; % default to *not found*
        if exist(fullfile(handles.parameters.lesion_img_folder,[handles.scorefiledata.RegistryCode{handles.scorefile.subjectrows(p)} '.nii']),'file')
            lesion_found = 1;
        end
        handles.scorefile.haslesion_file(p) = lesion_found;
    end
    handles.scorefile.total_lesions_found = sum(handles.scorefile.haslesion_file);

    %% Display who has all the required lesion files and scores for the current configuration...

    set(handles.lesionfilepresenttextbox,'String',sprintf('%d/%d subjects have lesion files.',handles.scorefile.total_lesions_found,handles.scorefile.nsubs_in_scorefile))
    set(handles.behavioralscorepresenttextbox,'String',sprintf('%d/%d subjects have required scores.',num_subs_that_have_all_variables,handles.scorefile.nsubs_in_scorefile))

    tmp = handles.scorefiledata.Properties.VariableNames;

    set([handles.scorenamepopupmenu handles.realcovariateslistbox handles.potentialcovariateslist handles.addcovariate handles.removecovariate],'enable','on')
    set([handles.scorenamepopupmenu handles.potentialcovariateslist],'String',tmp);

    desired_score_name_index = find(strcmp(handles.parameters.score_name,tmp));

    if isempty(desired_score_name_index)
        handles.parameters.score_name = tmp{1}; % default to 1....
        desired_score_name_index = find(strcmp(handles.parameters.score_name,tmp));
    end
    
    %% Paint the options to the GUI controls.
    
    set(handles.scorenamepopupmenu,'Value',desired_score_name_index)

    if ~isempty(handles.parameters.control_variable_names)
        set(handles.realcovariateslistbox,'String',handles.parameters.control_variable_names,'Value',min([get(handles.realcovariateslistbox,'Value') numel(handles.parameters.control_variable_names)]))
    else
        set(handles.realcovariateslistbox,'String',[],'Value',[])
    end

    if isempty(handles.parameters.control_variable_names)
        set(handles.removecovariate,'enable','off')
    else
        set(handles.removecovariate,'enable','on')
    end

    % Configure the run/view analysis button enabled or not
    if handles.parameters.analysis_is_completed == 0 % hasn't been run yet.
        set(handles.runanalysisbutton,'enable','on')
        set(handles.viewresultsbutton,'enable','off')
    elseif handles.parameters.analysis_is_completed == 1 % successfully completed.
        set(handles.runanalysisbutton,'enable','off')
        set(handles.viewresultsbutton,'enable','on')
    elseif handles.parameters.analysis_is_completed == 2 % there was an error...
        set(handles.runanalysisbutton,'enable','on')
        set(handles.viewresultsbutton,'enable','off')
    end

    set(handles.analysisnameeditbox,'String',handles.parameters.analysis_name)
    set(handles.lesionthresholdeditbox,'String',num2str(handles.parameters.lesion_thresh))
    set(handles.invertpmapcheckbox,'value',handles.parameters.invert_p_map_flag)

    set([handles.applycovariatestobehaviorcheckbox handles.applycovariatestolesioncheckbox],'enable','on'); % so we can modify them
    set(handles.applycovariatestobehaviorcheckbox,'value',handles.parameters.apply_covariates_to_behavior)
    set(handles.applycovariatestolesioncheckbox,'value',handles.parameters.apply_covariates_to_lesion)

    if numel(handles.parameters.control_variable_names) == 0 % if necessary, disable the checkboxes (if there's no covariates added)
        set([handles.applycovariatestobehaviorcheckbox handles.applycovariatestolesioncheckbox],'enable','off');
    end

    set(handles.npermutationseditbox,'string',num2str(handles.parameters.PermNumVoxelwise))

    set(handles.cluster_voxelwise_p_editbox,'string',strrep(num2str(handles.parameters.voxelwise_p),'0.','.')) % remove leading zero at PT request 5/2/17
    set(handles.clusterwisepeditbox,'string',strrep(num2str(handles.parameters.clusterwise_p),'0.','.')) % remove leading zero at PT request 5/2/17

    set(handles.permutationtestingcheckbox,'value',handles.parameters.DoPerformPermutationTesting)

    if ~handles.parameters.DoPerformPermutationTesting
        set(get(handles.permutationtestingpanel,'children'),'enable','off')
        set(handles.permutationtestingcheckbox,'enable','on')
    else
        set(get(handles.permutationtestingpanel,'children'),'enable','on')
    end
    
    set(handles.lesionvolumecorrectiondropdown,'Value',find(strcmp(handles.parameters.lesionvolcorrection,get(handles.lesionvolumecorrectiondropdown,'String'))))
    set(handles.hypodirectiondropdown,'Value',         find(strcmp(handles.parameters.tails              ,get(handles.hypodirectiondropdown,'String'))))

    UpdateTitleBar(handles);
