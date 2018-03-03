function handles = SaveSVRLSMGUIFile(handles,hObject)
    handles.parameters.is_saved = 1;
    handles.parameters.datetime_save = datetime;
    tosave = handles.parameters; % arbitrary, but I think it was giving me issues saving the field parameters from struct handles...
    save(handles.parameters.parameter_file_name,'tosave') % write the file
    handles = UpdateProgress(handles,['Saved ' handles.parameters.parameter_file_name],1);
    handles = LoadParametersFromSVRLSMFile(handles,hObject,handles.parameters.parameter_file_name);
