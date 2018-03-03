function parameters = svrlsm_optimizehyperparameters(parameters,variables)
% https://stackoverflow.com/questions/39636898/box-constraint-in-libsvm-package-compare-matlab-fitcsvm-and-libsvm-options

%% Which search strategy to use
 doGeneric = true;
 if ~doGeneric 
     switch parameters.optimization.search_strategy
         case 'Bayes Optimization'
             results = svrlsm_bayesopt(parameters,variables);
    %          for f = results.bestPoint.Properties.VariableNames % save the best field values...
    %              parameters.optimization.best.(f{1}) = results.bestPoint.(f{1}); % dynamic field names.
    %          end
             parameters.optimization.results = results;
         case 'Grid Search' 
    %         parameters.optimization.results = svrlsm_gridsearch(parameters,variables);
         case 'Random Search'
    %         parameters.optimization.results = svrlsm_randomsearch(parameters,variables);
     end

     for f = results.bestPoint.Properties.VariableNames % save the best field values...
         curname = f{1};
         curval = results.bestPoint.(curname);
         if strcmp(curname,'sigma'), parameters.optimization.best.sigma = curval; end
         if strcmp(curname,'box'), parameters.optimization.best.cost = curval; end
         if strcmp(curname,'epsilon'), parameters.optimization.best.epsilon = curval; end
         if strcmp(curname,'standardize'), parameters.optimization.best.standardize = curval == categorical(true); end
     end

 elseif doGeneric
    svrlsm_waitbar(parameters.waitbar,0,['Hyperparameter optimization (' CurrentOptimString(parameters) ')'])

     results = generic_hyperopts(parameters,variables);
     
     for f = results.bestPoint.Properties.VariableNames % save the best field values...
         curname = f{1};
         curval = results.bestPoint.(curname);
         
         if strcmp(curname,'KernelScale')
             parameters.optimization.best.sigma = curval;
         end
         
         if strcmp(curname,'BoxConstraint')
             parameters.optimization.best.cost = curval;
         end
         
         if strcmp(curname,'Epsilon')
             parameters.optimization.best.epsilon = curval;
         end
         
         if strcmp(curname,'Standardize')
             parameters.optimization.best.standardize = curval == categorical(true); % need it as a boolean...
         end
     end
     
     parameters.optimization.results = results;
    svrlsm_waitbar(parameters.waitbar,0,'')

 end
 
  