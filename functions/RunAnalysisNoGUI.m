function success = RunAnalysisNoGUI(input) % hObject,eventdata,handles)
switch class(input)
    case 'char' % then it's a path
        error('double check i''m working right.')
        handles.parameters = load(input); 
    case 'struct' % then it's a raw struct
        handles.parameters = input;
end

success = RunAnalysis([],[],handles); %success = RunAnalysis(hObject,eventdata,handles)