function results = generic_hyperopts(parameters,variables)
    % This attempts to optimize using matlab's default hyperparameter optimization options ... 
    % because it doesn't use bayesopt, it cannot update the progress bar (not output function)
    
   % For clarity...
   lesiondata = variables.lesion_dat;
   behavdata = variables.one_score;

   %% Configure hyperparam options
   params = hyperparameters('fitrsvm',lesiondata,behavdata); 
   
   % Turn optimize off for all hyperparams by default
   for N = 1:numel(params)
       params(N).Optimize = false;
   end
   
   % Then turn them on as necessary, and set their ranges from the user analysis config
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
   
   optimizeropts = resolveoptimizeropts(parameters);
   
   hyperoptoptions = struct('AcquisitionFunctionName','expected-improvement-plus', optimizeropts{:});
%        'Optimizer', resolveoptimizer(parameters)
%        'MaxObjectiveEvaluations',parameters.optimization.iterations, ...
%        'UseParallel',parameters.parallelize, ...
%        'Repartition',parameters.optimization.crossval.repartition, ...
%        'KFold',parameters.optimization.crossval.nfolds, ... %'PlotFcn',[], ...
%        'Verbose',myif(parameters.optimization.verbose_during_optimization,2,0));
   
       % 'OutputFcn',@optim_outputfun); % verbose is either 0 or 2...
   
   Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf', 'OptimizeHyperparameters',params, ...
       'HyperparameterOptimizationOptions', hyperoptoptions);

   results = Mdl.HyperparameterOptimizationResults;
%    
%    assignin('base','results',results)
%    assignin('base','Mdl',Mdl)
%    
%    
function optimizeropts = resolveoptimizeropts(parameters)

    switch parameters.optimization.search_strategy 
        case 'Bayes Optimization'
            optimchoice = 'bayesopt';
        case 'Grid Search' 
            optimchoice = 'gridsearch';
        case 'Random Search'
            optimchoice = 'randomsearch';
    end
    
    repartitionopt = myif(parameters.optimization.crossval.do_crossval, ...
        {'Repartition',parameters.optimization.crossval.repartition},{});
    itersornumdivsstr = myif(strcmp(optimchoice,'gridsearch'),'NumGridDivisions','MaxObjectiveEvaluations'); % do we need grid divs or do we need max objective evaluations...?
    itersornumdivs = myif(strcmp(optimchoice,'gridsearch'),parameters.optimization.grid_divisions,parameters.optimization.iterations);
    parameters.optimization.crossval.do_crossval = true;
    warning('Crossvalidation set to on for hyperparam opt.')
    nfolds = myif(parameters.optimization.crossval.do_crossval,parameters.optimization.crossval.nfolds,1);
    
    optimizeropts = {'Optimizer', optimchoice, ...
        itersornumdivsstr,itersornumdivs, ...
       'UseParallel',parameters.parallelize, ...
        repartitionopt{:}, ...
       'KFold', nfolds, ...
       'Verbose',myif(parameters.optimization.verbose_during_optimization,2,0), ...
       'ShowPlots',false};

   %optimizeropts{:}