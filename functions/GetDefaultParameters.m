function parameters = GetDefaultParameters(varargin) % (handles)
    % To add: make this a "default" .mat file, configurable in preferences...
    parameters.gui_version = 0.0;  % default.
    parameters.populated_to_gui = false; % default to false.
    if nargin > 0
%        if ishandle(varargin{1})
            parameters.populated_to_gui = true; %handles = varargin{1}; % called figure gui.
%        else
%            error('unknown input argument... should be figure handle or nothing')
%        end
    end
    
    parameters.method.mass_univariate = false; % if true, use mass univariate LSM, otherwise use SVR
    
    parameters.useLibSVM = 0; % now try to default to MATLAB ...
%     if handles.details.libsvm
%         parameters.useLibSVM = 1; % 1 = libSVM and 0 = MATLAB
%     else 
%         parameters.useLibSVM = 0; % 1 = libSVM and 0 = MATLAB
%     end

    parameters.svscaling = 100; % defaults to max (100th percentile) - added 9/29/17

    options = lsmtb_options; % should be made obsolete some time soon!
    %parameters.tails = handles.options.hypodirection{2}; % high scores good.
    parameters.tails = options.hypodirection{2}; % high scores good.
    parameters.datetime_run = []; % when the analysis was run.
    parameters.datetime_save = []; % when the analysis was run.
    parameters.analysis_is_completed = 0;
    parameters.parameter_file_name = ''; % identity of file if parameters were opened from a file.
    parameters.parallelize = true; % do not parallelize...
	
    mypath = fileparts(which('svrlsmgui'));
    parameters.analysis_name = 'Unnamed';
    parameters.is_saved = 0;
    
    %parameters.analysis_root = fullfile(mypath,'output');
    parameters.analysis_out_path = fullfile(mypath,'output'); % parameters.analysis_root; % is this a good default?
    parameters.score_file = fullfile(mypath,'default','PNT.csv');
    parameters.score_name = 'Sim_ROI_123';
    parameters.lesion_img_folder = fullfile(mypath,'default','lesion_imgs');
    
    parameters.DoPerformPermutationTesting = true;
    parameters.do_CFWER = false;
    parameters.cfwer_v_value = 1; % do we want to stick with this?
    parameters.cfwer_p_value = .05; % do we want to stick with tihs as default?

    % defaults -- can be set manually, but will be overridden if optimization is on.
    parameters.gamma = 5; % from Zhang et al., 2014 -- libsvm uses this
    parameters.sigma = gamma2sigma(parameters.gamma); % matlab uses this
    parameters.cost = 30; % from Zhang et al., 2014
    parameters.epsilon = .1; % default from libSVM for epsilon.
    parameters.standardize = true; % default behavior from Zhang et al, 2014

    % when current setting equals this, a (default) string will appear in the menu.
    parameters.svr_defaults.gamma = parameters.gamma;
    parameters.svr_defaults.sigma = parameters.sigma;
    parameters.svr_defaults.cost = parameters.cost;
    parameters.svr_defaults.epsilon = parameters.epsilon;
    parameters.svr_defaults.standardize = parameters.standardize;

    %parameters.lesionvolcorrection = handles.options.lesionvolumecorrection{3}; % both
    parameters.lesionvolcorrection = options.lesionvolumecorrection{3}; % both

    %parameters.beta_map = 1; 
    %parameters.sensitivity_map = 1; 

    %parameters.invert_p_map_flag = 1; % Inverse p-map, i.e., use 1-p instead of p for display on MRIcron.
    parameters.control_variable_names = {}; % None at first ...  

    parameters.lesion_thresh = 10; % The least lesion subject number for a voxel to be considered in the following analysis. 
    parameters.voxelwise_p = 0.005;
    parameters.clusterwise_p = 0.05;
    parameters.PermNumVoxelwise = 10000;
    parameters.PermNumClusterwise = 10000;

    % todo: make these apply to specific covariates (i.e. a boolean vector instead of the scalar for each)
    parameters.apply_covariates_to_behavior = 0;
    parameters.apply_covariates_to_lesion = 0;

    % output summary options
    parameters.do_make_summary = 1; % handles.summary_create_summary
    parameters.summary.narrative = true; % handles.summary_narrative_summary
    parameters.summary.beta_map = true;  % handles.summary_svrbetamap
    parameters.summary.voxelwise_thresholded = true; % handles.summary_thresholded_pmap
    parameters.summary.clusterwise_thresholded = true; % handles.summary_clusterthresholded
    parameters.summary.cfwer_diagnostics = true; % handles.summary_cfwerdiagnostics
    parameters.summary.variable_diagnostics = true; % handles.model_variablediagnostics
    parameters.summary.cluster_stability = true; % handles.summary_clusterstability
    parameters.summary.hyperparameter_optimization_record = true; % handles.summary_paramoptimization
    parameters.summary.parameter_assessment = true; % handles.summary_parameterassessment
    parameters.summary.lesion_overlap = true; %handles.summary_lesionoverlap
    parameters.summary.predictions = false; % write out predictions with r squared and whatnot

    % debug save output options
    parameters.SavePreThresholdedPermutations = 0;
    parameters.SavePostVoxelwiseThresholdedPermutations = 0;
    parameters.SavePostClusterwiseThresholdedPermutations = 0;
    parameters.SavePermutationData = 0; % leave the giant bin file?
    parameters.do_use_cache_when_available = false; % don't regenerate the big binary files if we have them around.
    % now the permutation data for the cfwer...
    parameters.SavePermutationPData = 0; % giant bin file for the cfwer
    parameters.SaveNullPMapsPreThresholding = 0;
    parameters.SaveNullPMapsPostThresholding = 0;

    % Hyperparameter Optimization - new as of Feb 2018
    parameters.optimization.do_optimize = false; % default.
    parameters.optimization.verbose_during_optimization = false; % spit out the progress?
    parameters.optimization.search_strategy = 'Bayes Optimization'; % or 'Grid Search' or 'Random Search'
    parameters.optimization.objective_function = 'Resubstitution Loss'; % 'Predict Behavior'; % or 'Maximum Correlation'
    parameters.optimization.iterations = 200; % by default allow 200 iterations for optimization.
    parameters.optimization.grid_divisions = 10; % by default, use 10 grid divisions...
    
    % if optimization is turned on, then do optimize these:
    parameters.optimization.params_to_optimize.cost = true; % optimize by default
    parameters.optimization.params_to_optimize.cost_range = [.001 300];
    parameters.optimization.params_to_optimize.cost_range_default = parameters.optimization.params_to_optimize.cost_range;
    parameters.optimization.params_to_optimize.sigma = true; % optimize by default  --- in the context of optimization this parameter is referred to as SIGMA not GAMMA and is in units of SIGMA.
    parameters.optimization.params_to_optimize.sigma_range = [.001 300]; 
    parameters.optimization.params_to_optimize.sigma_range_default = parameters.optimization.params_to_optimize.sigma_range;
    parameters.optimization.params_to_optimize.epsilon = true ; % don't optimize by default
    parameters.optimization.params_to_optimize.epsilon_range = [.005 3];
    parameters.optimization.params_to_optimize.epsilon_range_default = parameters.optimization.params_to_optimize.epsilon_range;
    parameters.optimization.params_to_optimize.standardize = true; % don't optimize by default
    parameters.optimization.params_to_optimize.standardize_range = [true false];
    
    % Cross-validation for hyperparameter optimization
    parameters.optimization.crossval.do_crossval = true;
    parameters.optimization.crossval.method = 'kfold'; % may be others in the future?
    parameters.optimization.crossval.nfolds = 5;
    parameters.optimization.crossval.nfolds_default = parameters.optimization.crossval.nfolds;
    parameters.optimization.crossval.repartition = true; % repartition at each iteration
    
    % Hyperparameter quality report
    parameters.hyperparameter_quality.report.nfolds = 5;
    parameters.hyperparameter_quality.report.repro_ind_subset_pct = .8; % 80%
    parameters.hyperparameter_quality.report.n_replications = 10; % how many times to repeat ...
    
    % What to do with the lesion data when we read it in.
    parameters.imagedata.do_binarize = true; % binarize image data in read_lesion_imgs (should be off for vbmish things)
    parameters.imagedata.do_resample = false;
    parameters.imagedata.resample_to = 2; % 2mm iso if do_resample is true;
    
    % ICA Beta stuff...
    parameters.beta.do_ica_on_lesiondata = false;
    
    if parameters.populated_to_gui % note that varargin{1} == "handles"
        UpdateProgress(varargin{1},'Retrieved default parameters...',1); %#ok<*NASGU> %handles = UpdateProgress(varargin{1},'Retrieved default parameters...',1); %#ok<*NASGU>
    end