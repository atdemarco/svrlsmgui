function variables = optimalParameterReport(parameters,variables)
    mkdir(variables.output_folder.hyperparameterinfo)
    %variables = optimalParameterReport(parameters,variables);
    %variables.files_created.cfwerinfo = fullcfwerout;

    % based on Zhang et al (2014)
    % returns results for two measures of hyper parameter optimality:

    % 1. reproducibility index and
    % 2. prediction accuracy
    % 3. mean absolute difference...

    [hyperparameter_quality.pred_accuracy, ...
     hyperparameter_quality.mean_abs_diff] = computePredictionAccuracy(parameters,variables);
 
    hyperparameter_quality.repro_index = computeReproducibilityIndex(parameters,variables);
    hyperparameter_quality.behavioral_predictions = computerBehavioralPrediction(parameters,variables); % this will be used for our WritePredictBehaviorReport in summary output.
    
    %% Save the results so we can use them later
    fname = 'hyperparameter_quality.mat';
    fpath = fullfile(variables.output_folder.hyperparameterinfo,fname);   
    save(fpath,'hyperparameter_quality','-v7.3') % if it's bigger than 2 GB we need v7.3 anyway...
    variables.files_created.hyperparameter_quality = fpath;
    variables.hyperparameter_quality = hyperparameter_quality;
    
function behavioral_predictions = computerBehavioralPrediction(parameters,variables)
    behavdata = variables.one_score; % for clarity extract these.
    lesiondata = variables.lesion_dat;

    %nfolds = 5; % Zhang et al., 2014
    nfolds = parameters.hyperparameter_quality.report.nfolds;

    hyperparms = hyperparmstruct(parameters);
    
    svrlsm_waitbar(parameters.waitbar,0,'Hyperparameter optimization: Storing behavioral predictions.');
    if parameters.useLibSVM
        libsvmstring = get_libsvm_spec(hyperparms.cost,hyperparms.gamma,hyperparms.epsilon); % we may need to re-apply standardization...?
        crossvalinds = crossvalind('KFold', variables.SubNum, nfolds);
        lesiondata=sparse(lesiondata);
        m = svmtrain(behavdata, lesiondata, libsvmstring); %#ok<SVMTRAIN> -- use all data
        behavioral_predictions.Mdl = m;
        for curfold=1:nfolds
            infold = crossvalinds == curfold;
            outoffold = crossvalinds ~= curfold;
            m = svmtrain(behavdata(outoffold), lesiondata(outoffold,:), libsvmstring); %#ok<SVMTRAIN>
            behavioral_predictions.XVMdl_predicted(infold) = svmpredict(behavdata(infold), lesiondata(infold,:), m, '-q'); % predict OOF observations based on cur fold data.
        end
        behavioral_predictions.XVMdl = [];
        
    else % Run with MATLAB machinery.
        Mdl = fitrsvm(lesiondata,behavdata,'ObservationsIn','rows','KernelFunction','rbf', 'KernelScale',hyperparms.sigma,'BoxConstraint',hyperparms.cost,'Standardize',hyperparms.standardize,'Epsilon',hyperparms.epsilon);
        XVMdl = crossval(Mdl,'KFold',nfolds); % this is a 5-fold

        behavioral_predictions.Mdl = Mdl;
        behavioral_predictions.XVMdl = XVMdl;
        behavioral_predictions.XVMdl_predicted = kfoldPredict(XVMdl);
    end
    svrlsm_waitbar(parameters.waitbar,0,'');
    
% 1. reproducibility index
%    - solve the SVRLSM N times for a random subset of 80% of subjects (Zhang et al used 85 of 106 pts for each rerun)
%    - compute pairwise correlations between all the resulting SVR-B maps
%    - mean of the correlation coefficients is the  REPRODUCIBILITY INDEX

function repro_index = computeReproducibilityIndex(parameters,variables)
    % For clarity extract these.
    behavdata = variables.one_score;
    lesiondata = variables.lesion_dat;
    
    % restore to 40 crossval steps once we introduce parallelization -- maybe make it user-configurable
    %N_subsets_to_perform = 2; % 10; % Zhang et al., 2015
    %subset_include_percentage = .8; % Zhang et al., 2015
    subset_include_percentage = parameters.hyperparameter_quality.report.repro_ind_subset_pct;
    N_subsets_to_perform = parameters.hyperparameter_quality.report.n_replications;

    nsubs = numel(behavdata);
    nholdin = round(subset_include_percentage*nsubs); % how many subjects should be in each subset?

    hyperparms = hyperparmstruct(parameters);
    
    N_subset_results = cell(1,N_subsets_to_perform); % save the correlation coefficients in this vector, then average them
    percent_obs_are_SVs = nan(1,N_subsets_to_perform); % count the number of SVs...
    
    svrlsm_waitbar(parameters.waitbar,0,'Hyperparameter optimization: Measuring reproducibility index.');
    for N = 1 : N_subsets_to_perform
        svrlsm_waitbar(parameters.waitbar,N/N_subsets_to_perform);
        includesubs = randperm(nsubs,nholdin); % each column contains a unique permutation
        
        if parameters.useLibSVM % do libsvm...
            libsvmstring = get_libsvm_spec(hyperparms.cost,hyperparms.gamma,hyperparms.epsilon); % we may need to re-apply standardization...?
            lesiondata=sparse(lesiondata);
            m = svmtrain(behavdata(includesubs), lesiondata(includesubs,:), libsvmstring); %#ok<SVMTRAIN>
            percent_obs_are_SVs(N) = m.totalSV / variables.SubNum; % should we instead use: sum(m.sv_coef == parameters.cost)) % (this is the number of bounded SVs)
            w = m.sv_coef'*m.SVs;
        else % do with MATLAB
            Mdl = fitrsvm(lesiondata(includesubs,:),behavdata(includesubs,:),'ObservationsIn','rows','KernelFunction','rbf', 'KernelScale',hyperparms.sigma,'BoxConstraint',hyperparms.cost,'Standardize',hyperparms.standardize,'Epsilon',hyperparms.epsilon);
            percent_obs_are_SVs(N) = sum(Mdl.IsSupportVector) / numel(Mdl.IsSupportVector);
            w = Mdl.Alpha.'*Mdl.SupportVectors;
        end
        
        if N == 1  % compute initial beta_scale to reuse...
            beta_scale = 10 / prctile(abs(w),parameters.svscaling); % parameters.svscaling is e.g, 100 or 99 or 95
        end
        
        w = w'*beta_scale; % the beta here is kind of irrelevant because it just scales the data, and correlation is insensitive to scaling...
        N_subset_results{N} = w; % save the w-map...
    end
    
    % Now for each pair of w maps in N_subset_results{:}, compute a correlation coefficient...
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
    repro_index.percent_obs_are_SVs = percent_obs_are_SVs; % number of support vectors...
    
% 2. prediction accuracy
%    - mean correlation coefficient between predixcted scores and testing scores with 40 5-fold cross-validations.
function [pred_accuracy,mean_abs_diff] = computePredictionAccuracy(parameters,variables)
    behavdata = variables.one_score; % for clarity extract these.
    lesiondata = variables.lesion_dat;
    
    % restore to 40 crossval steps once we introduce parallelization -- maybe make it user-configurable
    %N_crossvals_to_perform = 2; % 10; % Zhang et al., 2015
    % nfolds = 5; % Zhang et al., 2014

    nfolds  = parameters.hyperparameter_quality.report.nfolds;
    N_crossvals_to_perform = parameters.hyperparameter_quality.report.n_replications;
    
    hyperparms = hyperparmstruct(parameters);

    N_crossval_correl_results = nan(1,N_crossvals_to_perform); % save the correlation coefficients in this vector, then average them
    N_crossval_MAD_results = nan(1,N_crossvals_to_perform); % average mean absolute difference between predicted and real.
    percent_obs_are_SVs = nan(1,N_crossvals_to_perform); % count the number of SVs...
    
    if parameters.useLibSVM % try to save time by only doing this once...
        lesiondata=sparse(lesiondata);
    end
    
    svrlsm_waitbar(parameters.waitbar,0,'Hyperparameter optimization: Measuring prediction accuracy.');
    
    predicted = nan(variables.SubNum,1); % reserve space.
    for N = 1 : N_crossvals_to_perform
        svrlsm_waitbar(parameters.waitbar,N/N_crossvals_to_perform);
        if parameters.useLibSVM
            libsvmstring = get_libsvm_spec(hyperparms.cost,hyperparms.gamma,hyperparms.epsilon); % we may need to re-apply standardization...?
            crossvalinds = crossvalind('KFold', variables.SubNum, nfolds);
            for curfold=1:nfolds
                infold = crossvalinds == curfold;
                outoffold = crossvalinds ~= curfold;
                m = svmtrain(behavdata(outoffold), lesiondata(outoffold,:), libsvmstring); %#ok<SVMTRAIN>
                predicted(infold) = svmpredict(behavdata(infold), lesiondata(infold,:), m, '-q'); % predict OOF observations based on cur fold data.
            end
        else
            Mdl = fitrsvm(lesiondata,behavdata,'ObservationsIn','rows','KernelFunction','rbf', 'KernelScale',hyperparms.sigma,'BoxConstraint',hyperparms.cost,'Standardize',hyperparms.standardize,'Epsilon',hyperparms.epsilon);
            percent_obs_are_SVs(N) = sum(Mdl.IsSupportVector) / numel(Mdl.IsSupportVector);
            XVMdl = crossval(Mdl,'Kfold',nfolds); % this is a 5-fold
            predicted = kfoldPredict(XVMdl);
        end
        corrvals = corrcoef(predicted,behavdata); % compute the correlation...
        r = corrvals(2,1); % grab the r
        N_crossval_correl_results(N) = r; % save.
        N_crossval_MAD_results(N) = mean(abs(predicted-behavdata));
    end
    svrlsm_waitbar(parameters.waitbar,0,'');
    pred_accuracy.data = N_crossval_correl_results; % can plot a distribution if we like...
    pred_accuracy.mean = mean(N_crossval_correl_results);
    pred_accuracy.std = std(N_crossval_correl_results);
    pred_accuracy.percent_obs_are_SVs = percent_obs_are_SVs; % number of support vectors...
    
    mean_abs_diff.data = N_crossval_MAD_results; % can plot a distribution if we like... redundant w pred_accuracy
    mean_abs_diff.mean = mean(N_crossval_MAD_results);
    mean_abs_diff.std = std(N_crossval_MAD_results);
    mean_abs_diff.percent_obs_are_SVs = percent_obs_are_SVs; % number of support vectors... redundant w pred_accuracy