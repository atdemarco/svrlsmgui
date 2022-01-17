function [beta_map, variables] = get_beta_map(parameters, variables)
    variables.one_score = variables.one_score(:);
    variables.predict=[]; % placeholder
    
    % Decide whether we will do mass univariate or multivariate...
    if parameters.method.mass_univariate
        [beta_map,variables] = get_beta_map_mu(parameters,variables); % retrieve mass univariate beta map.
    else % do multivariate
        [beta_map,variables,Mdl] = get_beta_map_svr(parameters,variables); % retrieve svrlsm beta map.

        %% Compute and save the Model predictions for these real data
        predict.realData = variables.one_score;
        switch class(Mdl)
            case 'RegressionSVM'
                predict.resubPredict = Mdl.resubPredict;
                predict.resubLossMSE = Mdl.resubLoss('LossF','mse');
                predict.resubLossEps = Mdl.resubLoss('LossF','eps');
            case 'classreg.learning.partition.RegressionPartitionedSVM'
                predict.resubPredict = Mdl.kfoldPredict;
                predict.resubLossMSE = Mdl.kfoldLoss('LossF','mse');
                predict.resubLossEps = Mdl.kfoldLoss('LossF','eps');
        end
        variables.predict = predict; % store for return
        % compute the r2... - R-Squared is the ratio of Sum of Squares Regression (SSR) and Sum of Squares Total (SST).
        SSR = sum((predict.resubPredict-mean(predict.realData)).^2); % Sum of Squares Regression (SSR)
        SST = sum((predict.realData-mean(predict.realData)).^2); % Sum of Squares Total (SST).
        variables.predict.r2 = 1- (SSR/SST); % R2 = 1 - (SSR /SST) % adjusted 1/17/2021
    end
end
