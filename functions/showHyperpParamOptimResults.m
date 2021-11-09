function showHyperpParamOptimResults(parms)
error('i am defunct. delete me.')
    return 
    %% This function plots the visual results of hyperparameter optimization (if requested) into the output summary.
%     hyperoptoptions.NumGridDivisions=8
%     parms(1).Range(2) = 2000000; % kernel scsale to bigger range! 
%% if using grid search OR randomsearch, then we have a Rank column in the .results output table
%% if we use bayesopt we have...Mdl.HyperparameterOptimizationResults.bestPoint

%  hyperoptoptions.Optimizer='bayesopt'
%  % hyperoptoptions.Optimizer='randomsearch'
%  % hyperoptoptions.Optimizer='gridsearch'
% 
%      Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf', 'OptimizeHyperparameters',parms, 'HyperparameterOptimizationOptions', hyperoptoptions)
% 

    %%
    
    disp('Each hyperparameter vs the objective')

    % 2d plots vs the objective 
    ff=figure;
    subplot(2,2,1)
    scatter(results.BoxConstraint,results.Objective)
    xlabel('Box Constraint')
    ylabel('Objective')
    axis square

    subplot(2,2,2)
    scatter(results.Epsilon,results.Objective)
    xlabel('Epsilon')
    ylabel('Objective')
    axis square

    subplot(2,2,3)
    scatter(results.KernelScale,results.Objective)
    xlabel('Kernel Scale')
    ylabel('Objective')
    axis square

    subplot(2,2,4)
    scatter(results.Standardize,results.Objective)
    xlabel('Standardize')
    ylabel('Objective')
    axis square

    % 

    ff=figure;
    subplot(2,3,1)
    histogram(results.Objective)
    xlabel('Objective')
    ylabel('Count of param combos')

    subplot(2,3,2)

