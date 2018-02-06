function [handles,parameters] = step1(handles,parameters,variables)
    if parameters.parallelize
        [handles,parameters] = step1_parallel(handles,parameters,variables);
    else
        [handles,parameters] = step1_notparallel(handles,parameters,variables);
    end
