function yesno = OutputDirectoryAlreadyExists(handles)
    yesno = 0;
    if exist(handles.parameters.baseoutputdir,'dir')
        yesno = 1;
    end
