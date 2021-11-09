function parameters = svrlsm_optimizehyperparameters(parameters,variables)
    % https://stackoverflow.com/questions/39636898/box-constraint-in-libsvm-package-compare-matlab-fitcsvm-and-libsvm-options
     svrlsm_waitbar(parameters.waitbar,0,['Hyperparameter optimization (' CurrentOptimString(parameters) ')'])
 
     %% configure and run the optimization
     results = generic_hyperopts(parameters,variables);
     switch parameters.optimization.search_strategy
         case 'Bayes Optimization'
             %% store the optimization results
             for f = results.bestPoint.Properties.VariableNames % save the best field values...
                 curname = f{1};
                 curval = results.bestPoint.(curname);
                 if strcmp(curname,'KernelScale'), parameters.optimization.best.sigma = curval; end
                 if strcmp(curname,'BoxConstraint'), parameters.optimization.best.cost = curval; end
                 if strcmp(curname,'Epsilon'), parameters.optimization.best.epsilon = curval; end
                 if strcmp(curname,'Standardize'), parameters.optimization.best.standardize = curval == categorical(true); end % need it as a Boolean...
             end
         case {'Grid Search','Random Search'} % then we have a rank!
             bestRow = find(results.Rank == 1);
             tmpresults=removevars(results,{'Rank','Objective'});
             for f = tmpresults.Properties.VariableNames
                 curname = f{1};
                 curval = results.(curname)(bestRow); % extract the best row...
                 if strcmp(curname,'KernelScale'), parameters.optimization.best.sigma = curval; end
                 if strcmp(curname,'BoxConstraint'), parameters.optimization.best.cost = curval; end
                 if strcmp(curname,'Epsilon'), parameters.optimization.best.epsilon = curval; end
                 if strcmp(curname,'Standardize'), parameters.optimization.best.standardize = curval == categorical(true); end % need it as a Boolean...
             end
         otherwise
             error(['Unknown search strategy ...' parameters.optimization.search_strategy '?'])
     end
     
     parameters.optimization.results = results;
     svrlsm_waitbar(parameters.waitbar,0,'')