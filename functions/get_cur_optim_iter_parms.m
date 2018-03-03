function touse = get_cur_optim_iter_parms(x,parameters)
% helper function to retrieve the parameters we'll use for a given
% iteration of an objective function of hyperparameter optimization...
% we do this because you can't optimize a subset of parameters because if
% their field name occurs in the objective function it throws and error --
% so we do it all dynamically..
% ad - 2-23-18

    for f = {'cost','sigma','epsilon','standardize'} % the potential params to optimize
        if isfield(x,f{1}) % then we're asked to optimize it so return:
            touse.(f{1}) = x.(f{1}); % the dynamic value for each iteration
        else
            touse.(f{1}) = parameters.(f{1}); % or the constant supplied by the user
        end
    end