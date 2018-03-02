function handles = PopulateGUIFromParameters(handles)
    [~,path_to_show] = fileparts(handles.parameters.lesion_img_folder); % only show final directory.
    set(handles.lesionfoldereditbox,'string',path_to_show,'tooltipstring',handles.parameters.lesion_img_folder)

    [~,file,ext] = fileparts(handles.parameters.score_file); % only show final directory.
    set(handles.scorefileeditbox,'string',[file ext],'tooltipstring',handles.parameters.score_file)

    [~,path_to_show] = fileparts(handles.parameters.analysis_out_path); % only show final directory.
    set(handles.outputfoldereditbox,'string',path_to_show,'tooltipstring',handles.parameters.analysis_out_path)

    warning('off', 'MATLAB:table:ModifiedVarnames');
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames') % MATLAB 2017a
    
    opts = detectImportOptions(handles.parameters.score_file);
    handles.scorefiledata = readtable(handles.parameters.score_file,opts); % added to support e.g., MAC CSV files 1/31/18

    % tell user if registrycode is not a field
    
    fieldnames = handles.scorefiledata.Properties.VariableNames;
    if ~any(strcmp(fieldnames,'RegistryCode'))
        msgbox('Note: The selected CSV file does not appear to have a "RegistryCode" column listing subject lesion files. Please correct this and reload this file, otherwise the analysis will not run.')
        return
        % in the future, try to guess which is the subject list volumn...
    end
        
    handles.scorefile.subjectrows = find(~cellfun(@isempty,handles.scorefiledata.RegistryCode));
    handles.scorefile.nsubs_in_scorefile = numel(handles.scorefile.subjectrows);

    check_variables = [handles.parameters.score_name handles.parameters.control_variable_names];

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

    for c = 1 : numel(check_variables)
        coldata = handles.scorefiledata.(check_variables{c})(handles.scorefile.subjectrows);
        if isa(coldata,'cell')
            track_variables(:,c) = ~cellfun(@isempty,coldata);
        else % it's numbers...
            track_variables(:,c) = ~isnan(coldata);
        end
    end
    
    num_subs_that_have_all_variables = sum(all(track_variables,2));

    handles.scorefile.haslesion_file = zeros(numel(handles.scorefile.subjectrows),1);
    for p = handles.scorefile.subjectrows'
        lesion_found = 0; % default to *not found*
        if exist(fullfile(handles.parameters.lesion_img_folder,[handles.scorefiledata.RegistryCode{handles.scorefile.subjectrows(p)} '.nii']),'file')
            lesion_found = 1;
        end
        handles.scorefile.haslesion_file(p) = lesion_found;
    end
    handles.scorefile.total_lesions_found = sum(handles.scorefile.haslesion_file);

    %%

    set(handles.lesionfilepresenttextbox,'String',sprintf('%d/%d subjects have lesion files.',handles.scorefile.total_lesions_found,handles.scorefile.nsubs_in_scorefile))
    set(handles.behavioralscorepresenttextbox,'String',sprintf('%d/%d subjects have required scores.',num_subs_that_have_all_variables,handles.scorefile.nsubs_in_scorefile))

    tmp = handles.scorefiledata.Properties.VariableNames;
    
    % new 1/23/18
    registrycode_column_name = 'RegistryCode';
    registrycode_column = strcmp(registrycode_column_name,tmp);
    tmp(registrycode_column) = []; % remove.
    % %
    
    set([handles.scorenamepopupmenu handles.realcovariateslistbox handles.potentialcovariateslist handles.addcovariate handles.removecovariate],'enable','on')
    set([handles.scorenamepopupmenu handles.potentialcovariateslist],'String',tmp);

    desired_score_name_index = find(strcmp(handles.parameters.score_name,tmp));

    if isempty(desired_score_name_index)
        handles.parameters.score_name = tmp{1}; % default to 1....
        desired_score_name_index = find(strcmp(handles.parameters.score_name,tmp));
    end

    set(handles.scorenamepopupmenu,'Value',desired_score_name_index)

    if ~isempty(handles.parameters.control_variable_names)
        set(handles.realcovariateslistbox,'String',handles.parameters.control_variable_names,'Value',min([get(handles.realcovariateslistbox,'Value') numel(handles.parameters.control_variable_names)]))
    else
        set(handles.realcovariateslistbox,'String',[],'Value',[])
    end
    
    if numel(handles.parameters.control_variable_names) == 0
        set([handles.removecovariate handles.realcovariateslistbox],'enable','off')
    else
        set([handles.removecovariate handles.realcovariateslistbox],'enable','on')
    end

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
    set([handles.applycovariatestobehaviorcheckbox handles.applycovariatestolesioncheckbox],'enable','on'); % so we can modify them
    set(handles.applycovariatestobehaviorcheckbox,'value',handles.parameters.apply_covariates_to_behavior)
    set(handles.applycovariatestolesioncheckbox,'value',handles.parameters.apply_covariates_to_lesion)

    if numel(handles.parameters.control_variable_names) == 0 % if necessary, disable the checkboxes (if there's no covariates added)
        set([handles.applycovariatestobehaviorcheckbox handles.applycovariatestolesioncheckbox],'enable','off');
    end
    
    %% display permutation testing pane choices

    % set all the children of the pane to visible -- at the end of this code block we'll change the visibility.
    set(get(handles.permutationtestingpanel,'children'),'visible','on')

    % checkboxes...
    set(handles.permutationtestingcheckbox,'value',handles.parameters.DoPerformPermutationTesting)
    set(handles.do_cfwer_checkbox,'value',handles.parameters.do_CFWER)
        
    % fill in all the values in the pane
    set(handles.npermutationseditbox,'string',num2str(handles.parameters.PermNumVoxelwise))
    set(handles.cluster_voxelwise_p_editbox,'string',strrep(num2str(handles.parameters.voxelwise_p),'0.','.')) % remove leading zero - at PT request 5/2/17
    set(handles.clusterwisepeditbox,'string',strrep(num2str(handles.parameters.clusterwise_p),'0.','.')) % remove leading zero - at PT request 5/2/17
    set(handles.cfwer_v_value_editbox,'string',num2str(handles.parameters.cfwer_v_value));
    set(handles.cfwer_p_value_editbox,'string',strrep(num2str(handles.parameters.cfwer_p_value),'0.','.')) % remove leading zero

    if ~handles.parameters.DoPerformPermutationTesting
        set(get(handles.permutationtestingpanel,'children'),'enable','off')
        set(handles.permutationtestingcheckbox,'enable','on')
    else
        set(get(handles.permutationtestingpanel,'children'),'enable','on')
    end
    
    if handles.parameters.do_CFWER
        set([handles.cluster_voxelwise_p_editbox handles.clusterwisepeditbox handles.text18 handles.text13],'visible','off'); % hide the normal options
    else % hide the cwfer options...
        set([handles.cfwer_v_value_editbox handles.cfwer_p_value_editbox handles.text28 handles.text29],'visible','off'); % hide the cwfer options
    end
    
    %%

    set(handles.interrupt_button,'Enable','off') % only enable during an analysis....
    
    set(handles.lesionvolumecorrectiondropdown,'Value',find(strcmp(handles.parameters.lesionvolcorrection,get(handles.lesionvolumecorrectiondropdown,'String'))))
    

    % Accommodate reading in old config files with "one tail pos", "one tail neg" and "two tails" old labeling for hypothesis directionality
    if ~any(strcmp(handles.parameters.tails,handles.options.hypodirection)) % the options are labeled in the old way...
        if strfind(handles.parameters.tails,'pos') % one tail pos.
            handles.parameters.tails = handles.options.hypodirection{1}; % first option is on tailed pos no matter the label.
        elseif strfind(handles.parameters.tails,'neg') % one tail neg
            handles.parameters.tails = handles.options.hypodirection{2}; % second option is one tailed neg no matter the label.
        else % two-tailed
             handles.parameters.tails = handles.options.hypodirection{3}; % third option is two tailed no matter the label.
        end
    end
    
    set(handles.hypodirectiondropdown,'Value',         find(strcmp(handles.parameters.tails              ,get(handles.hypodirectiondropdown,'String'))))

    UpdateTitleBar(handles);
