function [objective,constraint] = combined_bayesopt_objective_functions(x,lesiondata,behavdata,parameters)
    % retrieve hyperparameter values to use for this iteration (see function for explanation of why we actually hafta do it this way)

    %touse = get_cur_optim_iter_parms(x,parameters); % returns constants for parms we're not optimizing on
    
    %% if we try to optimize epsilon, we end up getting a billion infeasible points.... what if we try to optimize standardize?
    % Do not crossval initially, so we can count the IsSupportVector field (which crossvalidated model objects do not have)
    %Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',touse.cost,'KernelScale',touse.sigma,'Standardize',touse.standardize,'Epsilon',touse.epsilon);
    %Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',x.box,'KernelScale',x.sigma,'Standardize',strcmp(x.standardize,'true'),'Epsilon',x.epsilon);
    %Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','BoxConstraint',x.box,'KernelScale',x.sigma,'Standardize',strcmp(x.standardize,'true'),'Epsilon',x.epsilon);
    
    Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','KernelScale',x.sigma,'BoxConstraint',x.box,'Epsilon',x.epsilon);

    % Do cross-validation as necessary...
    if parameters.optimization.crossval.do_crossval
        if parameters.optimization.crossval.repartition % create a partition each iteration
            curPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds);
        else % use the same initial one...
            curPartition = parameters.initialPartition.repartition;
        end
        CVMdl = crossval(Mdl,'CVPartition',curPartition); % Now crossval if desired.
    end
     
    switch parameters.optimization.objective_function
        case 'Predict Behavior' % embedded conditional of whether we predict or kfoldpredict
            if parameters.optimization.crossval.do_crossval
                predicted = kfoldPredict(CVMdl);
            else
                predicted = predict(Mdl,lesiondata);
            end
            
            difference =  behavdata(:) - predicted(:);
            mean_abs_diff = mean(abs(difference));
            objective = mean_abs_diff;  % we want to minimize this value.

        case 'Maximum Correlation' % embedded conditional of whether we predict or kfoldpredict
            if parameters.optimization.crossval.do_crossval
                predicted = kfoldPredict(CVMdl);
            else
                predicted = predict(Mdl,lesiondata);
            end
            
            corrvals = corrcoef(predicted,behavdata); % compute the correlation...
            objective = -1 * corrvals(2,1); % grab the r. --- and multiply by negative one since we want to maximize the value...
            
        case 'Resubstitution Loss'
            if parameters.optimization.crossval.do_crossval
                objective = kfoldLoss(CVMdl);
            else
                objective = resubLoss(Mdl);
            end
        otherwise
            error(['Unknown objective function string: ' parameters.optimization.objective_function])
    end
    
    min_pct_is_svs = .8; % at least 80% of the subjects should be used as a support vector.
    nsubs = numel(behavdata);
    min_n_support_vectors = round(min_pct_is_svs * nsubs);
    constraint = min_n_support_vectors - sum(Mdl.IsSupportVector); % negative means constraint is SATISFIED
    