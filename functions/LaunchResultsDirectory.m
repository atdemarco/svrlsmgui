function LaunchResultsDirectory(hObject,eventdata,handles)
     fulloutdir = fullfile(handles.parameters.analysis_out_path,handles.parameters.analysis_name);
     OpenDirectoryInNativeWindow(fulloutdir)
    
    % Open summary file in MATLAB's web browser.
    fpath = fileparts(handles.parameters.parmsfile);
    overviewhtmlfile = fullfile(fpath,'overview.html');
    if exist(overviewhtmlfile,'file')
        web(overviewhtmlfile) % Launch in MATLAB's "web browser"
    end