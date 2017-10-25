function UpdateTitleBar(handles)
    if ~handles.parameters.is_saved, saved_tag = '*';
    else saved_tag = '';
    end
    
    if ~exist(handles.parameters.parameter_file_name,'file'), fname = ['Unsaved analysis named ''' handles.parameters.analysis_name ''''];
    else, fname = handles.parameters.parameter_file_name;
    end
    
    set(gcf, 'Name', [fname saved_tag]);