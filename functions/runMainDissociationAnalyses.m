function handles = runMainDissociationAnalyses(hObject,eventdata,handles)
    % handles.parameters.orig_name = handles.parameters.analysis_name
    %% Run the two main analyses for dissociations - we'll use the output to compute the difference, etc
    for B = 1 : 2 % Run the full analysis for each main effect
        handles.parameters.SavePermutationData = true; % we'll need this data, so keep it around
        curbehav = handles.parameters.double_dissociation_behaviors{B};
        handles = UpdateProgress(handles,['Dissocation analysis ' num2str(B) ': ' curbehav],1);
        handles.parameters.score_name = curbehav; % update score_name so we analyze the right behavior
        % dissoclabel = ['dissociation_mainpart_' num2str(B) 'of2'];
        % handles.parameters.analysis_name = dissoclabel; % dynamic directory field name - so we can reference from interaction module
        [success,handles] = RunSingleAnalysis(hObject,eventdata,handles);
        handles.dissociation.maineffects{B} = handles.parameters; % so we can figure out where we should read the data from for each main effect svrb map...
    end