function [handles,parameters] = step1(handles,parameters,variables)
    if parameters.method.mass_univariate
        if parameters.parallelize
            error('implement the massunivariate step 1 -- parallelize')
        else
            error('implement the massunivariate step 1 -- non parallelize')
        end
    else
        if parameters.parallelize
            [handles,parameters] = step1_parallel(handles,parameters,variables);
        else
            [handles,parameters] = step1_notparallel(handles,parameters,variables);
        end
    end