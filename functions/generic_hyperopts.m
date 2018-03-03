function results = generic_hyperopts(parameters,variables)
    % This attempts to optimize using matlab's default hyperparameter optimization
    % options ... 
    
    % because it doesn't use bayesopt, it cannot update the progress bar (not output function)
    
    % For clarity...
   lesiondata = variables.lesion_dat;
   behavdata = variables.one_score;

   assignin('base','lesiondata',lesiondata)
   assignin('base','behavdata',behavdata)
   
   error('a')
   
   params = hyperparameters('fitrsvm',lesiondata,behavdata); 
   
   % turn optimize off for all hyperparams by default
   for N = 1:numel(params)
       params(N).Optimize = false;
   end
    
   if parameters.optimization.params_to_optimize.sigma
       params(strcmp({params.Name},'KernelScale')).Optimize = true;
       params(strcmp({params.Name},'KernelScale')).Range = [parameters.optimization.params_to_optimize.sigma_range];
   end
   
   if parameters.optimization.params_to_optimize.cost
       params(strcmp({params.Name},'BoxConstraint')).Optimize = true;
       params(strcmp({params.Name},'BoxConstraint')).Range = [parameters.optimization.params_to_optimize.cost_range];
   end
   
   if parameters.optimization.params_to_optimize.epsilon
       params(strcmp({params.Name},'Epsilon')).Optimize = true;
       params(strcmp({params.Name},'Epsilon')).Range = [parameters.optimization.params_to_optimize.epsilon_range];
   end
   
   if parameters.optimization.params_to_optimize.standardize
       params(strcmp({params.Name},'Standardize')).Optimize = true;
   end
   
   
   % parameters.optimization.objective_function 
   
   svrlsm_waitbar(parameters.waitbar,0,['Hyperparameter optimization (generic bayesopt: <add details> ...'])
   
   hyperoptoptions = struct('AcquisitionFunctionName','expected-improvement-plus', ...
       'MaxObjectiveEvaluations',parameters.optimization.iterations, ...
       'UseParallel',parameters.parallelize, ...
       'Repartition',parameters.optimization.crossval.repartition, ...
       'KFold',parameters.optimization.crossval.nfolds, ... %'PlotFcn',[], ...
       'Verbose',myif(parameters.optimization.verbose_during_optimization,2,0));
       % 'OutputFcn',@optim_outputfun); % verbose is either 0 or 2...
   
   Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf', 'OptimizeHyperparameters',params, ...
       'HyperparameterOptimizationOptions', hyperoptoptions); % default is KFold 5!

   results = Mdl.HyperparameterOptimizationResults;
   
   assignin('base','results',results)
   assignin('base','Mdl',Mdl)