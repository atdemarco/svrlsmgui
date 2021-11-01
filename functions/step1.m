function [handles,parameters,predAndLoss] = step1(handles,parameters,variables)
    if parameters.parallelize
        [handles,parameters,predAndLoss] = step1_parallel(handles,parameters,variables);
    else
        [handles,parameters,predAndLoss] = step1_notparallel(handles,parameters,variables);
    end