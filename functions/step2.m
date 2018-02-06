function [parameters,variables,thresholds] = step2(handles,parameters,variables,thresholds,all_perm_data)
if parameters.parallelize
    [parameters,variables,thresholds] = step2_parallel(handles,parameters,variables,thresholds,all_perm_data);
else
    [parameters,variables,thresholds] = step2_notparallel(handles,parameters,variables,thresholds,all_perm_data);
end
