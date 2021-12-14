function handles = runMainDissociationAnalyses(hObject,eventdata,handles)
    % handles.parameters.orig_name = handles.parameters.analysis_name
    %% Run the two main analyses for dissociations - we'll use the output to compute the difference, etc
    ORIG_CONTROL_VARS = handles.parameters.control_variable_names; % store;
    COV_OTHER_TAIL  = any(strcmp('Dissociation: Opposite behavior',ORIG_CONTROL_VARS));

    for B = 1 : 2 % Run the full analysis for each main effect
        handles.parameters.SavePermutationData = true; % we'll need this data, so keep it around
        curbehav = handles.parameters.double_dissociation_behaviors{B};
        handles = UpdateProgress(handles,['Dissocation analysis ' num2str(B) ': ' curbehav],1);
        handles.parameters.score_name = curbehav; % update score_name so we analyze the right behavior
        
        %% control for other tail?
        if COV_OTHER_TAIL 
            if B == 1, A = 2; else, A = 1; end
            otherbehav = handles.parameters.double_dissociation_behaviors{A}; 
            % replace the covariate called 'Dissociation: Opposite behavior' with the other behavior for this iteration            
            handles.parameters.control_variable_names{strcmp('Dissociation: Opposite behavior',ORIG_CONTROL_VARS)} = otherbehav; % slot it in there.
        end

        [success,handles] = RunSingleAnalysis(hObject,eventdata,handles);
        handles.dissociation.maineffects{B} = handles.parameters; % so we can figure out where we should read the data from for each main effect svrb map...
    end
    handles.parameters.control_variable_names = ORIG_CONTROL_VARS; % put back so the GUI displays the right information