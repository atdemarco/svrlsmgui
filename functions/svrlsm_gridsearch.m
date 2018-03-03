function results = svrlsm_gridsearch(parameters,variables)
%      if parameters.optimization.crossval.do_crossval
%          % Switch definition of anonymous function for objective function for optimization
%          switch parameters.optimization.objective_function
%              case 'Predict Behavior'
%                  results = fitrsvm(
%                  minfn = @(x) predict_behavior_w_crossval(x,lesiondata,behavdata,parameters);
%              case 'Maximum Correlation'
%                  minfn = @(x) max_correlation_w_crossval(x,lesiondata,behavdata,parameters);
%              case 'Resubstitution Loss'
%                  minfn = @(x) resubstitution_loss_w_crossval(x,lesiondata,behavdata,parameters);
%              otherwise
%                  error(['Unknown objective function string: ' parameters.optimization.objective_function])
%          end
%          % Start with this partition -- it will be replaced each call to the objective function if Repartition is 'true'
%          parameters.initialPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds);
% % 
% %          results = bayesopt(minfn,[sigma,box], ...
% %              'IsObjectiveDeterministic',false,'AcquisitionFunctionName','expected-improvement-plus', ...
% %              'MaxObjectiveEvaluations',parameters.optimization.iterations, 'UseParallel',parameters.parallelize, ...
% %              'PlotFcn',[],'Verbose',0);
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
%     svrlsm_waitbar(parameters.waitbar,0,''); % clear this.
