function [success,handles] = RunAnalysisNoGUI(varargin) % hObject,eventdata,handles)
    input = varargin{1};
    if nargin > 1, handles = varargin{2}; end % if batch was invoked from the GUI, we get the handle to the listbox to display progress etc...
    
    switch class(input)
        case 'char' % then it's a path
            handles.parameters = load(input); 
            handles.parameters = handles.parameters.tosave;
        case 'struct' % then it's a raw struct
            handles.parameters = input;
    end

    [success,handles] = RunAnalysis([],[],handles); %success = RunAnalysis(hObject,eventdata,handles)
    