function optim_string = CurrentOptimString(parameters)
    % assemble a string summarizing the current optimization choices
    % obj function (search strategy, # iters/divs)
    switch parameters.optimization.search_strategy 
        case 'Bayes Optimization', strategy_string = 'bayopt';
        case 'Grid Search', strategy_string = 'grid';
        case 'Random Search', strategy_string = 'rand';
        otherwise, error('Unknown search strategy string.')
    end

    switch parameters.optimization.search_strategy 
        case 'Grid Search', auxstr = ['(' strategy_string ', ' num2str(parameters.optimization.grid_divisions) ' divs)'];
        otherwise, auxstr = ['(' strategy_string ', ' num2str(parameters.optimization.iterations) ' iters)'];
    end

    switch parameters.optimization.objective_function
        case {'Predict Behavior','Maximum Correlation','Resubstitution Loss'}
            objfctn = parameters.optimization.objective_function;
            objfctn = strrep(objfctn,'Predict Behavior','Pred Behav'); % shorten for menu display niceness
            objfctn = strrep(objfctn,'Maximum Correlation','Max Corr'); % shorten for menu display niceness
            objfctn = strrep(objfctn,'Resubstitution Loss','Resub Loss'); % shorten for menu display niceness
        otherwise
            error(['Unknown objective function string: ' parameters.optimization.objective_function])
    end
    
    optim_string = sprintf('%s %s',objfctn,auxstr);
    
