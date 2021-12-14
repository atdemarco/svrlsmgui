function LaunchResultsDirectory(hObject,eventdata,handles)
     fulloutdir = fullfile(handles.parameters.analysis_out_path,handles.parameters.analysis_name);     
     fpath = fileparts(handles.parameters.parmsfile);
     overviewhtmlfile = fullfile(fpath,'overview.html');
    
    % Open summary file in MATLAB's web browser.
    OpenDirectoryInNativeWindow(fulloutdir)
    if exist(overviewhtmlfile,'file'), web(overviewhtmlfile); end % Launch in MATLAB's "web browser"