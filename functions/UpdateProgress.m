function handles = UpdateProgress(handles,message,add_to_listbox)
    if isfield(handles,'progresslistbox') % don't do this if the window's not open...
        olddata = get(handles.progresslistbox,'string');
        message = [datestr(now,13) ' ' message]; % HH:MM:SS
        if isempty(olddata)
            newdata = {message};
        else
            if ~add_to_listbox, newdata = {message};
            else newdata = [olddata ; message];
            end
        end
        numcells = numel(newdata);
        set(handles.progresslistbox,'string',newdata,'value',numcells);
        drawnow;
    end