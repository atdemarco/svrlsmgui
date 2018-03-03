function variables = optimalParameterReport(parameters,variables)
    mkdir(variables.output_folder.hyperparameterinfo)
    %variables = optimalParameterReport(parameters,variables);
    %variables.files_created.cfwerinfo = fullcfwerout;

    % based on Zhang et al (2014)
    % returns results for two measures of hyper parameter optimality:

    % 1. reproducibility index and
    % 2. prediction accuracy

    if parameters.useLibSVM
        error('Currently, optimalParameterReport() only supports MATLAB''s SVR machinery.')
    end
    
    hyperparameter_quality.pred_accuracy = computePredictionAccuracy(parameters,variables);
    hyperparameter_quality.repro_index = computeReproducibilityIndex(parameters,variables);
    
    
    %% save the results so we can use them later
    fname = 'hyperparameter_quality.mat';
    fpath = fullfile(variables.output_folder.hyperparameterinfo,fname);   
    save(fpath,'hyperparameter_quality')
    
    variables.files_created.hyperparameter_quality = fpath;
    variables.hyperparameter_quality = hyperparameter_quality;
    
% 1. reproducibility index
%    - solve the SVRLSM N times for a random subset of 80% of subjects (Zhang et al used 85 of 106 pts for each rerun)
%    - compute pairwise correlations between all the resulting SVR-B maps
%    - mean of the correlation coefficients is the  REPRODUCIBILITY INDEX
function repro_index = computeReproducibilityIndex(parameters,variables)
    % for clarity extract these.
    behavdata = variables.one_score;
    lesiondata = variables.lesion_dat;
    
    N_subsets_to_perform = 40; % Zhang et al., 2015
    subset_include_percentage = .8; % Zhang et al., 2015
    
    nsubs = numel(behavdata);
    nholdin = round(subset_include_percentage*nsubs); % how many subjects should be in each subset?
    %subset_indices = cell2mat(cellfun(@(x) randperm(nsubs,x)',repmat({nholdin},1,N_subsets_to_perform),'uni',false)); % each column contains a unique permutation

    % Decide whether we'll use optimized parameters or not...
    cost = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
    sigma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, parameters.optimization.best.sigma, parameters.sigma); % no longer derive from sigma
    standardize = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.standardize, parameters.optimization.best.standardize, parameters.standardize);
    epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);

    N_subset_results = cell(1,N_subsets_to_perform); % save the correlation coefficients in this vector, then average them

    for N = 1 : N_subsets_to_perform
        includesubs = randperm(nsubs,nholdin); % each column contains a unique permutation
        Mdl = fitrsvm(lesiondata(includesubs,:),behavdata(includesubs,:),'ObservationsIn','rows','KernelFunction','rbf', 'KernelScale',sigma,'BoxConstraint',cost,'Standardize',standardize,'Epsilon',epsilon);
        w = Mdl.Alpha.'*Mdl.SupportVectors;
        if N == 1  % compute initial beta_scale to reuse...
            beta_scale = 10 / prctile(abs(w),parameters.svscaling); % parameters.svscaling is e.g, 100 or 99 or 95
        end 
        w = w'*beta_scale; % the beta here is kind of irrelevant because it just scales the data, and correlation is insensitive to scaling...
        N_subset_results{N} = w; % save the w-map...
    end
    
    % now for each pair of w maps in N_subset_results{:}, compute a correlation coefficient...
    C = nchoosek(1:numel(N_subset_results),2); % rows of pairwise indices to correlate...
    repro_index_correlation_results = nan(size(C,1),1); % reseve space...
    for row = 1 : size(C,1)
        corrvals = corrcoef(N_subset_results{C(row,1)},N_subset_results{C(row,2)}); % compute the correlation...
        r = corrvals(2,1); % grab the r
        repro_index_correlation_results(row) = r; % store.
    end
    
    repro_index.data = repro_index_correlation_results;
    repro_index.mean = mean(repro_index_correlation_results);
    repro_index.std = std(repro_index_correlation_results);
    
% 2. prediction accuracy
%    - mean correlation coefficient between predixcted scores and testing
%    scores with 40 5-fold cross-validations.
function pred_accuracy = computePredictionAccuracy(parameters,variables)
    behavdata = variables.one_score; % for clarity extract these.
    lesiondata = variables.lesion_dat;
    
    N_crossvals_to_perform = 40; % Zhang et al., 2015
    nfolds = 5; % Zhang et al., 2014
%    partition = cvpartition(behavdata,'k',nfolds); % parameters.optimization.crossval.nfolds); % if 'repartition' is enabled, then repartition at each iteration.

    % Decide whether we'll use optimized parameters or not...
    cost = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
    
    cost = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost)
    
    parameters.optimization.params_to_optimize.cost
    parameters.optimization.best.cost
    parameters.cost
    
    sigma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, parameters.optimization.best.sigma, parameters.sigma); % no longer derive from sigma
    epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);
    standardize = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.standardize, parameters.optimization.best.standardize, parameters.standardize); 

    N_crossval_results = nan(1,N_crossvals_to_perform); % save the correlation coefficients in this vector, then average them

    for N = 1 : N_crossvals_to_perform
        Mdl = fitrsvm(lesiondata,behavdata,'ObservationsIn','rows','KernelFunction','rbf', 'KernelScale',sigma,'BoxConstraint',cost,'Standardize',standardize,'Epsilon',epsilon,'Kfold',nfolds); % this is a 5-fold
        predicted = kfoldPredict(Mdl);
        corrvals = corrcoef(predicted,behavdata); % compute the correlation...
        r = corrvals(2,1); % grab the r
        N_crossval_results(N) = r; % save.
    end
    
    pred_accuracy.data = N_crossval_results; % can plot a distribution if we like...
    pred_accuracy.mean = mean(N_crossval_results);
    pred_accuracy.std = std(N_crossval_results);
    
    