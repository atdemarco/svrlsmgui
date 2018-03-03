function handles = LoadParametersFromSVRLSMFile(handles,hObject,filepath)
    tmp = load(filepath); 
    handles.parameters = tmp.tosave; % tosave is the name of the variable we save; it's arbitrary but invariant.
    handles = PopulateGUIFromParameters(handles); % show what we've just loaded.    
    guidata(hObject, handles); % Update handles structure
    handles = UpdateProgress(handles,['Loaded ' handles.parameters.parameter_file_name],0);
