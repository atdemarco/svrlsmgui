function [beta_map, variables] = get_beta_map(parameters, variables)
    variables.one_score = variables.one_score(:);
    % Decide whether we will do mass univariate or multivariate...
    if parameters.method.mass_univariate
        [beta_map,variables] = get_beta_map_mu(parameters,variables); % retrieve mass univariate beta map.
    else % do multivariate
        [beta_map,variables] = get_beta_map_svr(parameters,variables); % retrieve svrlsm beta map.
    end
end
