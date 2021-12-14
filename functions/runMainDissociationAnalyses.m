function handles = runMainDissociationAnalyses(hObject,eventdata,handles)
    % handles.parameters.orig_name = handles.parameters.analysis_name
    %% Run the two main analyses for dissociations - we'll use the output to compute the difference, etc
    for B = 1 : 2 % Run the full analysis for each main effect
        handles.parameters.SavePermutationData = true; % we'll need this data, so keep it around
        curbehav = handles.parameters.double_dissociation_behaviors{B};
        handles = UpdateProgress(handles,['Dissocation analysis ' num2str(B) ': ' curbehav],1);
        handles.parameters.score_name = curbehav; % update score_name so we analyze the right behavior
        
        %% control for other tail?
        COV_OTHER_TAIL  = false;
        if COV_OTHER_TAIL 
            if B == 1, A = 2; else, A = 1; end
            otherbehav = handles.parameters.double_dissociation_behaviors{A}; 
            handles.parameters.control_variable_names = {otherbehav};
            handles.parameters.apply_covariates_to_behavior = 0;
            handles.parameters.apply_covariates_to_lesion = 0;
        end
        
        [success,handles] = RunSingleAnalysis(hObject,eventdata,handles);
        handles.dissociation.maineffects{B} = handles.parameters; % so we can figure out where we should read the data from for each main effect svrb map...
    end