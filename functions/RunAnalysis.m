function [success,handles] = RunAnalysis(hObject,eventdata,handles)
    if isempty(gcbo) || isa(gcbo,'matlab.ui.container.Menu')
        handles.parameters.runfromgui = 0;
    else
       handles.parameters.runfromgui = 1;
    end
    % Are we doing a single analysis, or a dissocation analysis?
    [success,handles] = singleOrDouble(hObject,eventdata,handles);

function [success,handles] = singleOrDouble(hObject,eventdata,handles)
    if handles.parameters.PERMIT_DOUBLE_DISSOCIATIONS && handles.parameters.run_double_dissociation % Then run the new procedure
        [success,handles] = runDissociation(hObject,eventdata,handles);
    else
        [success,handles] = RunSingleAnalysis(hObject,eventdata,handles);
    end