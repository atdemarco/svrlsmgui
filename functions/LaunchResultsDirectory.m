function LaunchResultsDirectory(hObject,eventdata,handles)
    fulloutdir = fullfile(handles.parameters.analysis_out_path,handles.parameters.analysis_name);
    if ismac || isunix % use open ...
        [~] = system(['open "' fulloutdir '"']);
    elseif ispc % use winopen
        winopen(fulloutdir)
    else
        warndlg('Cannot open output directory because I cannot determine the OS you are using.')
    end
    
    % Open overview file... %dev1
    fpath = fileparts(handles.parameters.parmsfile);
    overviewhtmlfile = fullfile(fpath,'overview.html');
    if exist(overviewhtmlfile,'file')
        web(overviewhtmlfile) % Launch in MATLAB's "web browser"
    end