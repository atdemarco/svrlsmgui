function results = svrlsm_bayesopt(parameters,variables)
    % For clarity...
    lesiondata = variables.lesion_dat;
    behavdata = variables.one_score;

   %% Which parameters to optimize and ranges to optim over
   params = hyperparameters('fitrsvm',lesiondata,behavdata);
   standrange = params(strcmp('Standardize',{params.Name})).Range;
   
   sigma = optimizableVariable('sigma',parameters.optimization.params_to_optimize.sigma_range,'Transform','log','optimize',parameters.optimization.params_to_optimize.sigma); % gamma ~ sigma...
   box = optimizableVariable('box',parameters.optimization.params_to_optimize.cost_range,'Transform','log','optimize',parameters.optimization.params_to_optimize.cost); % cost  == boxconstraint
   epsilon = optimizableVariable('epsilon',parameters.optimization.params_to_optimize.epsilon_range,'Transform','none','optimize',parameters.optimization.params_to_optimize.epsilon);
   %standardize = optimizableVariable('standardize',standrange,'Transform','none','Type','categorical','optimize',parameters.optimization.params_to_optimize.standardize);
   %standrange = [0 1]; % [true false]; % {categorical(true), categorical(false)}; % ?!
   standardize = optimizableVariable('standardize',standrange,'Transform','none','Type','categorical','optimize',parameters.optimization.params_to_optimize.standardize);
   svrlsm_waitbar(parameters.waitbar,0,['Hyperparameter optimization (bayesopt: ' parameters.optimization.objective_function ')...'])
   
     % Has the user requested cross-validation or not?
%      if parameters.optimization.crossval.do_crossval
%          % Switch definition of anonymous function for objective function for optimization
%          switch parameters.optimization.objective_function
%              case 'Predict Behavior'
%                  minfn = @(x) predict_behavior_w_crossval(x,lesiondata,behavdata,parameters);
%              case 'Maximum Correlation'
%                  minfn = @(x) max_correlation_w_crossval(x,lesiondata,behavdata,parameters);
%              case 'Resubstitution Loss'
%                  minfn = @(x) resubstitution_loss_w_crossval(x,lesiondata,behavdata,parameters);
%              otherwise
%                  error(['Unknown objective function string: ' parameters.optimization.objective_function])
%          end
% 
%          % Start with this partition -- it will be replaced each call to the objective function if Repartition is 'true'
%          parameters.initialPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds);
% 
%          %results = bayesopt(minfn,[sigma,box], ...
%          results = bayesopt(minfn,[gamma,box,epsilon,standardize], ...
%              'IsObjectiveDeterministic',false,'AcquisitionFunctionName','expected-improvement-plus', ...
%              'MaxObjectiveEvaluations',parameters.optimization.iterations, 'UseParallel',parameters.parallelize, ...
%              'PlotFcn',[],'Verbose',0);
% 
%      else % No cross-validation was requested, so run these without cross-validation...
%          
%          % Switch definition of anonymous function for objective function for optimization
%          switch parameters.optimization.objective_function
%              case 'Predict Behavior'
%                  minfn = @(x) predict_behavior(x,lesiondata,behavdata,parameters);
%                 %minfn = @(x)sqrt(sum((behavdata - resubPredict(fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',x.box,'KernelScale',x.sigma,'Epsilon',parameters.epsilon,'Standardize',false,'true')))).^2);
%              case 'Maximum Correlation'
%                  minfn = @(x) max_correlation(x,lesiondata,behavdata,parameters);
%              case 'Resubstitution Loss'
%                  minfn = @(x) resubstitution_loss(x,lesiondata,behavdata,parameters);
%              otherwise 
%                  error(['Unknown objective function string: ' parameters.optimization.objective_function])
%          end
%          
%         results = bayesopt(minfn,[sigma,box], ... % ,epsilon,standardize], ... % the parameters we're going to try to optimize
%             'IsObjectiveDeterministic',true,'AcquisitionFunctionName','expected-improvement-plus', ...
%             'MaxObjectiveEvaluations',parameters.optimization.iterations, 'UseParallel',parameters.parallelize, ...
%             'PlotFcn',[],'Verbose',0,'OutputFcn',@optim_outputfun);
% 
%      end

%%

%     % Start with this partition -- it will be replaced each call to the objective function if Repartition is 'true'
    parameters.initialPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds);
% assignin('base','lesiondata',lesiondata)
% error('s')
    minfn = @(x)combined_bayesopt_objective_functions3(x,sparse(lesiondata),behavdata,parameters);

    results = bayesopt(minfn,[sigma box epsilon standardize], ... 
        'IsObjectiveDeterministic',false, ... % parameters.optimization.crossval.do_crossval, ...
        'AcquisitionFunctionName','expected-improvement-plus', ... % 'ExplorationRatio', 0.7, ... % 0.5 default
        'MaxObjectiveEvaluations',parameters.optimization.iterations, ...
        'UseParallel',parameters.parallelize, ...
        'NumCoupledConstraints', 1, 'AreCoupledConstraintsDeterministic', false, ... % ~parameters.optimization.crossval.do_crossval, ... % if crossval, then it's not deterministic.
        'PlotFcn',[],'OutputFcn',@optim_outputfun, ...
        'Verbose',myif(parameters.optimization.verbose_during_optimization,2,0));  % verbose is either 0 or 2...

    svrlsm_waitbar(parameters.waitbar,0,''); % clear this.

% %% Objective function: Maximum correlation - no cross-validation
% function objective = max_correlation(x,lesiondata,behavdata,parameters)
%     svrtype = myif(parameters.useLibSVM,'libSVM','MATLAB');
%     switch svrtype
%         case 'libSVM'
%             error('transform gamma to sigma here?')
%             m = svmtrain(behavdata,sparse(lesiondata),get_libsvm_spec(x.cost,x.gamma,parameters.epsilon));% note get_libsvm_spec() function
%             predicted = svmpredict(behavdata, sparse(lesiondata), m, '-q');
%         case 'MATLAB' % since we take care of standardization in preprocessing
%             Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',x.box,'KernelScale',x.sigma,'Epsilon',parameters.epsilon,'Standardize',false);
%             predicted = predict(Mdl,lesiondata);
%     end
%     corrvals = corrcoef(predicted,behavdata); % compute the correlation...
%     objective = -1 * corrvals(2,1); % grab the r. --- and multiply by negative one since we want to maximize the value...
% 
% %% Objective function: Predict behavior - no cross-validation
function objective = predict_behavior(x,lesiondata,behavdata,parameters)
%     svrtype = myif(parameters.useLibSVM,'libSVM','MATLAB');
%      switch svrtype
%          case 'libSVM'
%              error('transform gamma to sigma here?')
%               m = svmtrain(behavdata,sparse(lesiondata),get_libsvm_spec(x.cost,x.gamma,parameters.epsilon));% note get_libsvm_spec() function
%               predicted = svmpredict(behavdata, sparse(lesiondata), m, '-q');
%          case 'MATLAB' % since we take care of standardization in preprocessing
%             touse = get_cur_optim_iter_parms(x,parameters); % returns constants for parms we're not optimizing on
%            Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',touse.cost,'KernelScale',touse.sigma,'Standardize',touse.standardize,'Epsilon',touse.epsilon);
             Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',x.cost,'KernelScale',x.sigma,'Standardize',false,'Epsilon',parameters.epsilon);
            predicted = predict(Mdl,lesiondata);
%      end
    difference =  behavdata(:) - predicted(:);
    mean_abs_diff = mean(abs(difference));
    objective = mean_abs_diff;  % we want to minimize this value.
%     
% %% Objective function: Resubstitution loss - no cross-validation
% function objective = resubstitution_loss(x,lesiondata,behavdata,parameters)
% 
%      touse = get_cur_optim_iter_parms(x,parameters); % returns constants for parms we're not optimizing on
%      Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',touse.cost,'KernelScale',touse.gamma,'Standardize',touse.standardize,'Epsilon',touse.epsilon);
% %	 Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',x.box,'KernelScale',x.sigma,'Epsilon',parameters.epsilon,'Standardize',false);
%      objective = resubLoss(Mdl);
% 
% %% Objective function: Resubstitution loss - with K-fold cross-validation
% function objective = resubstitution_loss_w_crossval(x,lesiondata,behavdata,parameters)
%      if parameters.optimization.crossval.repartition
%          curPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds);
%      else
%          curPartition = parameters.initialPartition;
%      end
%      
%      touse = get_cur_optim_iter_parms(x,parameters); % returns constants for parms we're not optimizing on
%      Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',touse.cost,'KernelScale',touse.gamma,'Standardize',touse.standardize,'Epsilon',touse.epsilon,'CVPartition',curPartition);
%               
% %     Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',x.box,'KernelScale',x.sigma,'Standardize',false,'Epsilon',parameters.epsilon,'CVPartition', curPartition);
% %        myif(parameters.optimization.crossval.repartition,partition,partition.repartition)); 
%     
%      objective = kfoldLoss(Mdl);
% 
% %% Objective function: Predict behavior - with K-fold cross-validation
% function objective = predict_behavior_w_crossval(x,lesiondata,behavdata,parameters)
% %     svrtype = myif(parameters.useLibSVM,'libSVM','MATLAB');
% %      switch svrtype
% %          case 'libSVM'
% %              error('transform gamma to sigma here?')
% %               m = svmtrain(behavdata,sparse(lesiondata),get_libsvm_spec(x.cost,x.gamma,parameters.epsilon));% note get_libsvm_spec() function
% %               predicted = svmpredict(behavdata, sparse(lesiondata), m, '-q');
% %          case 'MATLAB' % since we take care of standardization in preprocessing
%              
%               if parameters.optimization.crossval.repartition
%                   curPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds);
%               else
%                   curPartition = parameters.initialPartition;
%               end
%               
%               touse = get_cur_optim_iter_parms(x,parameters); % returns constants for parms we're not optimizing on
%               Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',touse.cost,'KernelScale',touse.gamma,'Standardize',touse.standardize,'Epsilon',touse.epsilon,'CVPartition',curPartition);
%               predicted = kfoldPredict(Mdl);
% %      end
%     difference =  behavdata(:) - predicted(:);
%     mean_abs_diff = mean(abs(difference));
%     objective = mean_abs_diff;  % we want to minimize this value.
% 
% 
% %% Objective function: Maximum correlation - with K-fold cross-validation
% function objective = max_correlation_w_crossval(x,lesiondata,behavdata,parameters)
% %     svrtype = myif(parameters.useLibSVM,'libSVM','MATLAB');
% %     switch svrtype
% %         case 'libSVM'
% %             error('transform gamma to sigma here?')
% %             m = svmtrain(behavdata,sparse(lesiondata),get_libsvm_spec(x.cost,x.gamma,parameters.epsilon));% note get_libsvm_spec() function
% %             predicted = svmpredict(behavdata, sparse(lesiondata), m, '-q');
% %         case 'MATLAB' % since we take care of standardization in preprocessing
%              if parameters.optimization.crossval.repartition
%                  curPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds);
%              else
%                  curPartition = parameters.initialPartition.repartition;
%              end
%              touse = get_cur_optim_iter_parms(x,parameters); % returns constants for parms we're not optimizing on
%              Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',touse.cost,'KernelScale',touse.gamma,'Standardize',touse.standardize,'Epsilon',touse.epsilon,'CVPartition',curPartition);
%              
%             predicted = kfoldPredict(Mdl);
%     %end
%     corrvals = corrcoef(predicted,behavdata); % compute the correlation...
%     objective = -1 * corrvals(2,1); % grab the r. --- and multiply by negative one since we want to maximize the value...
% 

% function objective = combined_bayesopt_objective_functions(x,lesiondata,behavdata,parameters)
% 
%     % hyperparameter values to use for this iteration
%     touse = get_cur_optim_iter_parms(x,parameters); % returns constants for parms we're not optimizing on
%     
%     % add in cross-validation parameters dynamically as necessary...
%     crossvalparms = {}; % default
%     if parameters.optimization.crossval.do_crossval
%         if parameters.optimization.crossval.repartition
%             curPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds);
%         else
%             curPartition = parameters.initialPartition.repartition;
%         end
%         crossvalparms ={'CVPartition',curPartition};
%     end
%     
%     Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',touse.cost,'KernelScale',touse.sigma,'Standardize',touse.standardize,'Epsilon',touse.epsilon,crossvalparms{:}); % 'CVPartition',curPartition);
%     
%     switch parameters.optimization.objective_function
%         case 'Predict Behavior' % embedded conditional of whether we predict or kfoldpredict
%             if parameters.optimization.crossval.do_crossval
%                 predicted = kfoldPredict(Mdl);
%             else
%                 predicted = predict(Mdl,lesiondata);
%             end
%             
%             difference =  behavdata(:) - predicted(:);
%             mean_abs_diff = mean(abs(difference));
%             objective = mean_abs_diff;  % we want to minimize this value.
% 
%         case 'Maximum Correlation' % embedded conditional of whether we predict or kfoldpredict
%             if parameters.optimization.crossval.do_crossval
%                 predicted = kfoldPredict(Mdl);
%             else
%                 predicted = predict(Mdl,lesiondata);
%             end
%             
%             corrvals = corrcoef(predicted,behavdata); % compute the correlation...
%             objective = -1 * corrvals(2,1); % grab the r. --- and multiply by negative one since we want to maximize the value...
%             
%         case 'Resubstitution Loss'
%             if parameters.optimization.crossval.do_crossval
%                 objective = kfoldLoss(Mdl);
%             else
%                 objective = resubLoss(Mdl);
%             end
%         otherwise
%             error(['Unknown objective function string: ' parameters.optimization.objective_function])
%     end
%                  
% function touse = get_cur_optim_iter_parms(x,parameters)
%     parms = {'cost','sigma','epsilon','standardize'};
%     for f = parms
%         if isfield(x,f{1}) % then we're asked to optimize it so return:
%             touse.(f{1}) = x.(f{1}); % the dynamic value for each iteration
%         else
%             touse.(f{1}) = parameters.(f{1}); % or the constant supplied by the user
%         end
%     end