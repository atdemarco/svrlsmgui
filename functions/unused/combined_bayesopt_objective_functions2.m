function [objective,constraint] = combined_bayesopt_objective_functions2(x,lesiondata,behavdata,parameters)
    % this replicates Zhang et al 2014 with a bayes optimization search by maximizing the 
    n_crossval = 30; % 40 from Zhang et al...
    vals = nan(1,n_crossval);
    nsvs = nan(1,n_crossval);

    required_percent_svs = .80; %
    nsubs = numel(behavdata);
    required_svs = round(required_percent_svs * nsubs);

    for N = 1 : numel(vals)
%          if strcmp(x.standardize,'true') % this probably does NOT work and is screwing up hyperparameter optimization altogether!
%              dostand = true;
%          else
%              dostand = false;
%          end
%         
        Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','KernelScale',x.sigma,'BoxConstraint',x.box,'Epsilon',x.epsilon,'Standardize',x.standardize); % 'Standardize',dostand);
        nsvs(N) = (required_svs + .5) - sum(Mdl.IsSupportVector); % negative means constraint is SATISFIED
        
        curPartition = cvpartition(behavdata,'k',parameters.optimization.crossval.nfolds); % repartition the model
        CVMdl = crossval(Mdl,'CVPartition',curPartition); % crossvalidate with k folds...
        
        switch parameters.optimization.objective_function
            case 'Predict Behavior'
                predicted = kfoldPredict(CVMdl);
                difference =  behavdata(:) - predicted(:);
                vals(N) = mean(abs(difference)); % minimize the mean absolute difference
            case 'Maximum Correlation' % embedded conditional of whether we predict or kfoldpredict
                predicted = kfoldPredict(CVMdl);
                corrvals = corrcoef(predicted,behavdata); % compute the correlation...
                vals(N) = -1 * corrvals(2,1); % grab the r. --- and multiply by negative one since we want to maximize the value...
                %vals(N) = corrvals(2,1); % grab the r. --- and multiply by negative one since we want to maximize the value...
            case 'Resubstitution Loss'
                vals(N) = kfoldLoss(CVMdl);
            otherwise
                error(['Unknown objective function string: ' parameters.optimization.objective_function])
        end
    end
    objective = mean(vals);
    constraint = mean(nsvs);