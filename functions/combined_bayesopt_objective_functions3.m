function [objective,constraint] = combined_bayesopt_objective_functions3(x,lesiondata,behavdata,parameters)
    lesiondata=full(lesiondata); % reconstruct the full volume
    
    % this replicates Zhang et al 2014 with a bayes optimization search by maximizing the 
    n_crossval = 5; % % Zhang et al. used 40.
    vals = nan(1,n_crossval);
    nsvs = nan(1,n_crossval);

    required_percent_svs = .80; %
    nsubs = numel(behavdata);
    required_svs = round(required_percent_svs * nsubs);

%     z.sigma = parameters.sigma;
%     z.box = parameters.cost;
%     z.epsilon = parameters.epsilon;
%     z.standardize = parameters.standardize;
%     
%     fieldnames_to_optimize = fieldnames(x);
%     for f = 1 : numel(fieldnames_to_optimize)
%         switch fieldnames_to_optimize{f}
%             case 'box' % cost
%                 z.box = x.box;
%             case 'sigma'
%                 z.sigma = x.sigma;
%             case 'epsilon'
%                 z.epsilon = x.epsilon;
%             case 'standardize'
%                 z.standardize = strcmp(x.standardize,'true');
%         end
%     end

    for N = 1 : numel(vals)
        Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','KernelScale',x.sigma,'BoxConstraint',x.box,'Epsilon',x.epsilon,'Standardize',strcmp(x.standardize,'true'));
        %Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','KernelScale',x.sigma,'BoxConstraint',x.box,'Epsilon',x.epsilon); % ,'Standardize',strcmp(x.standardize,'true'));
        %Mdl = fitrsvm(lesiondata,behavdata,'KernelFunction','rbf','KernelScale',z.sigma,'BoxConstraint',z.box,'Epsilon',z.epsilon,'Standardize',z.standardize);
        
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
                
                %vals(N) = -1 * corrvals(2,1); % grab the r. --- and multiply by negative one since we want to maximize the value...
                vals(N) = corrvals(2,1); % grab the r. --- and multiply by negative one since we want to maximize the value...
            case 'Resubstitution Loss'
                vals(N) = kfoldLoss(CVMdl);
            otherwise
                error(['Unknown objective function string: ' parameters.optimization.objective_function])
        end
    end
    objective = mean(vals);
    constraint = mean(nsvs);