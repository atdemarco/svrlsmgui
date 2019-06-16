function [Mdl,w,variables] = ComputeMatlabSVRLSM(parameters,variables)

    % Decide whether we'll use optimized parameters or not...
    box = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
    sigma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, parameters.optimization.best.sigma, parameters.sigma); % no longer derive from sigma
    standardize = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.standardize, parameters.optimization.best.standardize, parameters.standardize); 
    epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);
    variables.one_score = variables.one_score(:);
    
        Mdl = fitrsvm(variables.lesion_dat,variables.one_score,'ObservationsIn','rows', 'KernelFunction','rbf', 'KernelScale',sigma,'BoxConstraint',box,'Standardize',standardize,'Epsilon',epsilon);
        w = Mdl.Alpha.'*Mdl.SupportVectors;

        % as of v0.8 9/29/17 we have customized scaling available in parameters.svscaling
        variables.beta_scale = 10 / prctile(abs(w),parameters.svscaling); % parameters.svscaling is e.g, 100 or 99 or 95 % 10/max(abs(w));
