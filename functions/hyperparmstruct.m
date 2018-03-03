function hyperparms = hyperparmstruct(parameters)
    % retrieve the hyperparameters that will be used as a struct.
    hyperparms.cost = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
    hyperparms.sigma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, parameters.optimization.best.sigma, parameters.sigma);
    hyperparms.gamma = sigma2gamma(hyperparms.sigma); % for libsvm
    hyperparms.standardize = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.standardize, parameters.optimization.best.standardize, parameters.standardize);
    hyperparms.epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);
    
    % sometimes we end up with .standardize as a categorical which is annoying, so let's replace with a boolean
    if iscategorical(hyperparms.standardize)
        hyperparms.standardize = myif(hyperparms.standardize=='true',true,false);
    end
    
    