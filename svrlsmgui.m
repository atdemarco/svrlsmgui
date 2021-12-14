function varargout = svrlsmgui(varargin)
% SVRLSMGUI MATLAB code for svrlsmgui.fig
%      SVRLSMGUI, by itself, creates a new SVRLSMGUI or raises the existing
%      singleton*.
%
%      H = SVRLSMGUI returns the handle to a new SVRLSMGUI or the handle to
%      the existing singleton*.
%
%      SVRLSMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SVRLSMGUI.M with the given input arguments.
%
%      SVRLSMGUI('Property','Value',...) creates a new SVRLSMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before svrlsmgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to svrlsmgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help svrlsmgui

% Last Modified by GUIDE v2.5 02-Nov-2021 10:41:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename,'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @svrlsmgui_OpeningFcn, 'gui_OutputFcn',  @svrlsmgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , 'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1}); 
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:}); 
end

% End initialization code - DO NOT EDIT

% --- Executes just before svrlsmgui is made visible.
function svrlsmgui_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject; % Choose default command line output for svrlsmgui

    addPathsAsNecessary; % Add paths as necessary...
    
    handles = ConfigureSVRLSMGUIOptions(handles);
    handles.details = CheckIfNecessaryFilesAreInstalled(handles);

    handles = populateGUI(handles);

    % 0.02 - trying to clean it up to run on a variety of systems - 4/24/17
    % 0.03 - first version to be used by other individuals officially - 5/1/17
    % 0.04 - 5/2/17
    % 0.05 - 7/26/17 - moved to the linux machine, working on improving stability and fixing bugs
    % 0.06 - August 2017 - added parallelization, continued development with paper
    % 0.07 - September 2017 - summary output, fixed bug in two tail thresholding, public release...
    % 0.08 - added diagnostic plot for behavioral nuisance model in summary
    %        output file (corrplot); added custom support vector scaling
    %        other than max of the map when backprojecting the analysis
    %        hyperplane
    % 0.10 - January 2018 - fixing reported bugs
    % 0.15 - massive code refactoring and implementation of CFWER
    % 2.0 - revision at end of 2021
    %       -- no libsvm support
    %       -- no sv scaling
    
    handles.parameters.gui_version = 2; % 0.15; % version of the the gui
    guidata(hObject, handles); % Update handles structure

function addPathsAsNecessary
    % Do we need to add the functions subdirectory to the path?
    pathCell = regexp(path, pathsep, 'split');
    myPath = fileparts(mfilename('fullpath'));
    if ~any(strcmp(myPath,pathCell)), addpath(myPath); end
    
    functionsPath = fullfile(myPath,'functions');
    if ~any(strcmp(functionsPath,pathCell)), addpath(functionsPath); end
    

function handles = DisableAll(handles)
    set(get(handles.analysispreferencespanel,'children'),'enable','off')
    set(get(handles.covariatespanel,'children'),'enable','off')
    set(get(handles.permutationtestingpanel,'children'),'enable','off')
    set([handles.viewresultsbutton handles.interrupt_button handles.runanalysisbutton],'Enable','off') % handles.cancelanalysisbutton 
    set(handles.optionsmenu,'enable','off') % since viewing this menu references parameters that may not be loaded.  
    msgbox('One or more necessary component is missing from MATLAB''s path. Address the message in the SVRLSMgui window and restart this gui.')

function handles = populateGUI(handles)
    if handles.details.stats_toolbox && handles.details.spm
        handles = UpdateProgress(handles,'All necessary functions are available...',1);
        handles.parameters = GetDefaultParameters(handles);
        handles = PopulateGUIFromParameters(handles);
    elseif ~handles.details.spm 
        handles = UpdateProgress(handles,'SPM12 functions not available. Download and/or add SPM12 to the MATLAB path and relaunch the SVRLSMGUI.',1);
        handles = DisableAll(handles);
    elseif ~handles.details.stats_toolbox
         handles = UpdateProgress(handles,'No SVR algorithm available. Install Statistics Toolbox in MATLAB and then relaunch the GUI.',1);
         handles = DisableAll(handles);
    end

    handles.parameters.parallelize = handles.details.can_parallelize; % override default

function handles = UpdateCurrentAnalysis(handles,hObject)
    changemade = true; % default
    
    switch get(gcbo,'tag') % use gcbo to see what the cbo is and determine what field it goes to -- and to validate
        case 'no_map_crossvalidation' % turns off crossvalidation...
            handles.parameters.crossval.do_crossval = false;
        case 'kfold_map_crossvalidation'
            answer = inputdlg(sprintf('Enter the numbers of folds for crossvalidation.'), 'Number of folds', 1,{num2str(handles.parameters.crossval.nfolds)});
             if isempty(answer), return; end % cancel pressed
             str = str2num(answer{1});
             if isempty(str) || str <= 0 || ~isint(str)
                 changemade=false;
                 warndlg('Input must be a positive integer.');
             else % update the parameter value.
                handles.parameters.crossval.do_crossval = true;
                handles.parameters.crossval.method = 'kfold';
                handles.parameters.crossval.nfolds = str;
             end
        case 'hyperparm_quality_report_options'
            set(handles.hyperparm_qual_n_folds,'Label',['Folds: ' num2str(handles.parameters.hyperparameter_quality.report.nfolds)])
            set(handles.repro_index_subset_percentage,'Label',['Subset %: ' num2str(handles.parameters.hyperparameter_quality.report.repro_ind_subset_pct)])
            set(handles.hyperparm_qual_n_replications,'Label',['Replications: ' num2str(handles.parameters.hyperparameter_quality.report.n_replications)])
        case 'hyperparm_qual_n_folds'
            answer = inputdlg(sprintf('Enter the number of folds for hyperparameter quality testing:'), ...
                'Number of folds', 1,{num2str(handles.parameters.hyperparameter_quality.report.nfolds)});
            if isempty(answer), return; end % cancel pressed
            str = str2num(answer{1});
            if isempty(str) || str <= 1 || ~isint(str)
                changemade=false;
                warndlg('Input must be a positive integer greater than 1.');
            else % update the parameter value.
                handles.parameters.hyperparameter_quality.report.nfolds = str;
            end
        case 'repro_index_subset_percentage'
            answer = inputdlg('Enter the % of sample to use for computing w-map reproducibility index [0-1]:', ...
                'Sample percent', 1,{num2str(handles.parameters.hyperparameter_quality.report.repro_ind_subset_pct)});
            if isempty(answer), return; end % cancel pressed
            str = str2num(answer{1});
            if isempty(str) || str <= 0 || str >= 1 % ~isint(str)
                changemade=false;
                warndlg('Input must be a value between 0 and 1 (i.e. 0 and 100%).');
            else % update the parameter value.
                handles.parameters.hyperparameter_quality.report.repro_ind_subset_pct = str;
            end
        case 'hyperparm_qual_n_replications'
            answer = inputdlg(sprintf('Enter the number of replications for hyperparameter quality testing:'), ...
                'Number of replications', 1,{num2str(handles.parameters.hyperparameter_quality.report.n_replications)});
            if isempty(answer), return; end % cancel pressed
            str = str2num(answer{1});
            if isempty(str) || str <= 1 || ~isint(str)
                changemade=false;
                warndlg('Input must be a positive integer greater than 1.');
            else % update the parameter value.
                handles.parameters.hyperparameter_quality.report.n_replications = str;
            end
        case 'image_data_options_parent_menu'
            set(handles.do_binarize_data_menu_option,'Checked',myif(handles.parameters.imagedata.do_binarize,'on','off'))
        case 'do_binarize_data_menu_option'
            handles.parameters.imagedata.do_binarize = ~handles.parameters.imagedata.do_binarize;
        case 'set_resolution_parent_menu_option'
            set(handles.manual_analysis_resolution_menu,'Label',['Manual: ' num2str(handles.parameters.imagedata.resample_to) ' mm'], ...
                'checked',myif(handles.parameters.imagedata.do_resample,'on','off'));
            set(handles.do_not_resample_images_menu,'Checked',myif(handles.parameters.imagedata.do_resample,'off','on'));
            %set(handles.manual_analysis_resolution_menu,'Enable','on') % Until we finish implementation.
        case 'do_not_resample_images_menu'
            handles.parameters.imagedata.do_resample = false; % turn resampling off.
        case 'manual_analysis_resolution_menu'
            answer = inputdlg(sprintf('Enter the size in mm for voxels to be resampled to.'), ...
                'Resample size (mm)', 1,{num2str(handles.parameters.imagedata.resample_to)});
            if isempty(answer), return; end % cancel pressed
            str = str2num(answer{1});
            if isempty(str) || str <= 0 || ~isint(str)
                changemade=false;
                warndlg('Input must be a positive integer in millimeters.');
            else % update the parameter value.
                handles.parameters.imagedata.do_resample = true;
                handles.parameters.imagedata.resample_to = str;
            end
        case 'open_lesion_folder_button'
            OpenDirectoryInNativeWindow(handles.parameters.lesion_img_folder)
        case 'open_score_file_button'
            openFileInSystemViewer(handles.parameters.score_file)
        case 'open_output_folder_button'
            OpenDirectoryInNativeWindow(handles.parameters.analysis_out_path)
        case 'ica_lesion_decompose_option'
            handles.parameters.beta.do_ica_on_lesiondata = ~handles.parameters.beta.do_ica_on_lesiondata;
        case 'requirements_menu'
        case 'search_strategy_options'
            set(handles.optimization_iterations_menu_option,'Label',['Iterations: ' num2str(handles.parameters.optimization.iterations)])
            set(handles.griddivs_optimization_menu_option,'Label',['Grid Divs: ' num2str(handles.parameters.optimization.grid_divisions)])
            return
        case 'summary_prediction_menu'
            handles.parameters.summary.predictions = ~handles.parameters.summary.predictions;
        case 'lsm_method_parent_menu'
            set(get(handles.lsm_method_parent_menu,'children'),'checked','off')
            if handles.parameters.method.mass_univariate
                set(handles.mass_univariate_menu_option,'checked','on')
                set(handles.svr_parent_menu,'enable','off')
            else
                set(handles.multivariate_lsm_option,'checked','on')
                set(handles.svr_parent_menu,'enable','on')
            end

            if handles.details.stats_toolbox
                set(handles.multivariate_lsm_option,'enable','on')
            else
                handles.parameters.method.mass_univariate = true; % at least...
                set(handles.multivariate_lsm_option,'enable','off')            
            end
        case 'mass_univariate_menu_option'
            handles.parameters.method.mass_univariate = true;
        case 'multivariate_lsm_option'
            handles.parameters.method.mass_univariate = false;
        case 'svr_parent_menu'
        case 'optimization_is_verbose_menu'
            handles.parameters.optimization.verbose_during_optimization = ~handles.parameters.optimization.verbose_during_optimization;
        case 'summary_lesionoverlap'
            handles.parameters.summary.lesion_overlap = ~handles.parameters.summary.lesion_overlap;
        case 'summary_paramoptimization'
            handles.parameters.summary.hyperparameter_optimization_record = ~handles.parameters.summary.hyperparameter_optimization_record;
        case 'summary_create_summary'
            handles.parameters.do_make_summary = ~handles.parameters.do_make_summary;
        case 'summary_narrative_summary'
            handles.parameters.summary.narrative = ~handles.parameters.summary.narrative;
        case 'summary_svrbetamap'
            handles.parameters.summary.beta_map = ~handles.parameters.summary.beta_map;
        case 'summary_voxelwise_thresholded'
            handles.parameters.summary.voxelwise_thresholded = ~handles.parameters.summary.voxelwise_thresholded;
        case 'summary_clusterwise_thresholded'
            handles.parameters.summary.clusterwise_thresholded = ~handles.parameters.summary.clusterwise_thresholded;
        case 'summary_cfwerdiagnostics'
            handles.parameters.summary.cfwer_diagnostics = ~handles.parameters.summary.cfwer_diagnostics;
        case 'model_variablediagnostics'
            handles.parameters.summary.variable_diagnostics = ~handles.parameters.summary.variable_diagnostics;
        case 'summary_clusterstability'
            handles.parameters.summary.cluster_stability = ~handles.parameters.summary.cluster_stability;
        case 'summary_parameterassessment'
            handles.parameters.summary.parameter_assessment = ~handles.parameters.summary.parameter_assessment;
        case 'do_use_cache_menu'
            handles.parameters.do_use_cache_when_available = ~handles.parameters.do_use_cache_when_available;
        case 'crossval_menu_option_none'
            handles.parameters.optimization.crossval.do_crossval = false; % disable crossvalidation.
        case 'crossval_menu_option_kfold'
            msg = sprintf('Enter the number of folds for cross-validation (default = %d):',handles.parameters.optimization.crossval.nfolds_default);
            answer = inputdlg(msg, ...
                'Number of folds', 1,{num2str(handles.parameters.optimization.crossval.nfolds)});
            warning('bug here - if you click cancel it causes an error - in future, embed the str<=0 and isint within the isempty... and do that with other str2num(answer{1})''s as well')
            if isempty(answer), return; end % cancel pressed
            str = str2num(answer{1});
            if isempty(str) || str <= 0 || ~isint(str)
                changemade=false;
                warndlg('Input must be a positive integer.');
            else % update the parameter value.
                handles.parameters.optimization.crossval.do_crossval = true; % enable crossvalidation.
                handles.parameters.optimization.crossval.nfolds = str;
            end

        case 'do_repartition_menu_option' % flip this choice.
            handles.parameters.optimization.crossval.repartition = ~handles.parameters.optimization.crossval.repartition;
        case 'standardize_menu' % set the value to use if not optimization
            vals = {'true','false'};
            msg = sprintf('Standardize value (default = %s):',myif(handles.parameters.svr_defaults.standardize,'true','false'));
            s = listdlg('PromptString',msg,'SelectionMode','single', ...
                'ListString',vals,'InitialValue',myif(handles.parameters.standardize,1,2), ...
                'Name','Standardize Parameter','ListSize',[250 80]);

            if isempty(s)
                changemade=false; % cancelled...
            else
                handles.parameters.standardize = str2num(vals{s}); % myif(v==1,true,false); % new value.
            end

        case 'epsilon_menu' % set the value to use if not optimization
            msg = sprintf('Enter new value for epsilon (default = %0.2f):',handles.parameters.svr_defaults.epsilon);
            % add min and max range - dev1
            answer = inputdlg(msg,'Epsilon Parameter',1,{num2str(handles.parameters.epsilon)});
            if isempty(answer), return; end % cancel pressed
            numval = str2num(answer{1});
            if isnumeric(numval) && ~isempty(numval)
                handles.parameters.epsilon = numval;
            end

        % Whether to allow optimization of given hyperparameter
        case 'do_optimize_cost_menu'
            if ~handles.parameters.optimization.params_to_optimize.cost % then enabling it will prompt for new values...
                minmsg = sprintf('Minimum (default = %0.3f):',handles.parameters.optimization.params_to_optimize.cost_range_default(1));
                maxmsg = sprintf('Maximum (default = %0.2f):',handles.parameters.optimization.params_to_optimize.cost_range_default(2));
                msg = {minmsg,maxmsg};
                defaultans = {num2str(handles.parameters.optimization.params_to_optimize.cost_range(1)),num2str(handles.parameters.optimization.params_to_optimize.cost_range(2))};
                answer = inputdlg(msg,'Cost Range',1,defaultans);
                if isempty(answer), return; end % cancel pressed
                if any(cellfun(@isempty,answer)), warndlg('Invalid Cost range.'); return; end
                minval = str2num(answer{1}); maxval = str2num(answer{2});
                if ~all([isnumeric(minval) isnumeric(maxval)]), warndlg('Cost range values must be numbers.'); return; end
                if ~all([minval maxval] > 0), warndlg('Cost range values must be positive.'); return; end
                handles.parameters.optimization.params_to_optimize.cost_range = sort([minval maxval]);
            end

            handles.parameters.optimization.params_to_optimize.cost = ~handles.parameters.optimization.params_to_optimize.cost;
        case 'do_optimize_gamma_menu'
            if ~handles.parameters.optimization.params_to_optimize.sigma % then enabling it will prompt for new values...
                defaultmin = handles.parameters.optimization.params_to_optimize.sigma_range_default(1);
                defaultmax = handles.parameters.optimization.params_to_optimize.sigma_range_default(2);
                minmsg = sprintf('Minimum (default = %0.3f):',defaultmin); % convert to gamma if necessary
                maxmsg = sprintf('Maximum (default = %0.2f):',defaultmax); % convert to gamma if necessary
                msg = {minmsg,maxmsg};
                oldmin = handles.parameters.optimization.params_to_optimize.sigma_range(1);
                oldmax = handles.parameters.optimization.params_to_optimize.sigma_range(2);
                defaultans = {num2str(oldmin),num2str(oldmax)}; % convert to gamma if necessary
                answer = inputdlg(msg,'Sigma Range',1,defaultans);  % display as gamma if necessary
                if isempty(answer), return; end % cancel pressed
                if any(cellfun(@isempty,answer)), warndlg('Invalid sigma range.'); return; end
                minval = str2num(answer{1}); maxval = str2num(answer{2});
                if ~all([isnumeric(minval) isnumeric(maxval)]), warndlg('Sigma range values must be numbers.'); return; end % display as gamma if necessary
                if ~all([minval maxval] > 0), warndlg('Sigma range values must be positive.'); return; end % display as gamma if necessary
                handles.parameters.optimization.params_to_optimize.sigma_range = sort([minval maxval]); % convert gamma back to sigma if necessary.
            end

            handles.parameters.optimization.params_to_optimize.sigma = ~handles.parameters.optimization.params_to_optimize.sigma;
        case 'do_optimize_epsilon_menu'
            if ~handles.parameters.optimization.params_to_optimize.epsilon % then enabling it will prompt for new values...
                minmsg = sprintf('Minimum (default = %0.3f):',handles.parameters.optimization.params_to_optimize.epsilon_range_default(1));
                maxmsg = sprintf('Maximum (default = %0.2f):',handles.parameters.optimization.params_to_optimize.epsilon_range_default(2));
                msg = {minmsg,maxmsg};
                defaultans = {num2str(handles.parameters.optimization.params_to_optimize.epsilon_range(1)),num2str(handles.parameters.optimization.params_to_optimize.epsilon_range(2))};
                answer = inputdlg(msg,'Epsilon Range',1,defaultans);
                if isempty(answer), return; end % cancel pressed
                if any(cellfun(@isempty,answer)), warndlg('Invalid Epsilon range.'); return; end
                minval = str2num(answer{1}); maxval = str2num(answer{2});
                if ~all([isnumeric(minval) isnumeric(maxval)]), warndlg('Epsilon range values must be numbers.'); return; end
                if ~all([minval maxval] > 0), warndlg('Epsilon range values must be positive.'); return; end
                handles.parameters.optimization.params_to_optimize.epsilon_range = sort([minval maxval]);
            end
            handles.parameters.optimization.params_to_optimize.epsilon = ~handles.parameters.optimization.params_to_optimize.epsilon;
        case 'do_optimize_standardize_menu'
            handles.parameters.optimization.params_to_optimize.standardize = ~handles.parameters.optimization.params_to_optimize.standardize;

        % Optimization misc choices
        case 'optimization_iterations_menu_option'
            answer = inputdlg('Enter a new number of iterations:','Number of optimization iterations',1,{num2str(handles.parameters.optimization.iterations)});
            if isempty(answer), return; end % cancel pressed
            str = str2num(answer{1});
            if isempty(str) || str <= 0 || ~isint(str)
                changemade=false;
                warndlg('Input must be a positive integer.');
            else % update the parameter value.
                handles.parameters.optimization.iterations = str;
            end
        case 'griddivs_optimization_menu_option'
            answer = inputdlg('Enter a new number of grid divisions:','Number of optimization grid divisions',1,{num2str(handles.parameters.optimization.grid_divisions)});
            if isempty(answer), return; end % cancel pressed
            str = str2num(answer{1});
            if isempty(str) || str <= 0 || ~isint(str)
                changemade=false;
                warndlg('Input must be a positive integer.');
            else % update the parameter value.
                handles.parameters.optimization.grid_divisions = str;
            end    

        % Optimization search strategy choice
        case 'random_search_menu_option'
            handles.parameters.optimization.search_strategy = 'Random Search';
        case 'bayes_optimization_menu_choice'
            handles.parameters.optimization.search_strategy = 'Bayes Optimization';
        case 'gridsearch_option'
            handles.parameters.optimization.search_strategy = 'Grid Search';        
        % Optimization objective function choice
        case 'predictbehavior_optimize_menu_choice' % bayes opt predict behavior
            handles.parameters.optimization.objective_function = 'Predict Behavior';
        case 'correlation_optimize_menu_choice' % maximum correlation w behavior
            handles.parameters.optimization.objective_function = 'Maximum Correlation';
        case 'resubloss_optimize_menu_choice'
            handles.parameters.optimization.objective_function = 'Resubstitution Loss';
        % Whether or not to optimize hyperparameters
        case 'no_optimize_menu_choice'
            handles.parameters.optimization.do_optimize = false;
        case 'current_optimization_menu_option'
            handles.parameters.optimization.do_optimize = true; % will use whatever setting is configured.
        case 'save_pre_thresh'
            handles.parameters.SavePreThresholdedPermutations = ~handles.parameters.SavePreThresholdedPermutations;
        case 'retain_big_binary_file'
            handles.parameters.SavePermutationData = ~handles.parameters.SavePermutationData;
        case 'save_post_vox_thresh'
            handles.parameters.SavePostVoxelwiseThresholdedPermutations = ~handles.parameters.SavePostVoxelwiseThresholdedPermutations;
        case 'save_post_clusterwise_thresholded'
            handles.parameters.SavePostClusterwiseThresholdedPermutations= ~handles.parameters.SavePostClusterwiseThresholdedPermutations;
        case 'save_unthresholded_pmaps_cfwer'
            handles.parameters.SaveNullPMapsPreThresholding = ~handles.parameters.SaveNullPMapsPreThresholding;
        case 'save_thresholded_pmaps_cfwer'
            handles.parameters.SaveNullPMapsPostThresholding = ~handles.parameters.SaveNullPMapsPostThresholding;
        case 'retain_big_binary_pval_file'
            handles.parameters.SavePermutationPData = ~handles.parameters.SavePermutationPData;
        case 'parallelizemenu'
            handles.parameters.parallelize = ~handles.parameters.parallelize;
        case 'applycovariatestobehaviorcheckbox'
            handles.parameters.apply_covariates_to_behavior = get(gcbo,'value');
        case 'applycovariatestolesioncheckbox'
            handles.parameters.apply_covariates_to_lesion = get(gcbo,'value');
        case 'lesionvolumecorrectiondropdown'
            contents = get(handles.lesionvolumecorrectiondropdown,'string');
            newval = contents{get(handles.lesionvolumecorrectiondropdown,'value')};
            handles.parameters.lesionvolcorrection = newval;
        case 'hypodirectiondropdown'
            contents = get(handles.hypodirectiondropdown,'string');
            newval = contents{get(handles.hypodirectiondropdown,'value')};
            if find(strcmp(newval,handles.options.hypodirection)) == 3 % Disable two-tails from the GUI
                warndlg('Two-tailed hypothesis tests are not available in this version of SVRLSMGUI.')
                changemade = false;
            else
                handles.parameters.tails = newval;
            end
        case 'addcovariate'
            % first, what are we trying to add?
            contents = get(handles.potentialcovariateslist,'String'); 
            newcovariate = contents{get(handles.potentialcovariateslist,'Value')};
            if strcmp(newcovariate,handles.parameters.score_name) % check if it's our main score name...
                warndlg('This variable is already chosen as the main outcome in the analysis. If you''d like to add it as a covariate, remove it from "Score Name."')
                changemade = false;
            elseif any(strcmp(newcovariate,handles.parameters.control_variable_names)) % only if it's not on our list of covariates already.
                warndlg('This variable is already on the list of covariates. You may not add it twice.')
                changemade = false;
            else % it's new
                handles.parameters.control_variable_names{end+1} = newcovariate;
            end
        case 'removecovariate'
           contents = get(handles.potentialcovariateslist,'String'); % what are we trying to remove? 
           newcovariate = contents{get(handles.potentialcovariateslist,'Value')};
           if any(strcmp(newcovariate,handles.parameters.control_variable_names)) % only if it IS on our list!
               index_to_remove = strcmp(newcovariate,handles.parameters.control_variable_names);
               handles.parameters.control_variable_names(index_to_remove) = [];
           else
                changemade = false;
           end
        case 'chooselesionfolderbutton'
            folder_name = uigetdir(handles.parameters.lesion_img_folder,'Choose a folder containing lesion files for this analysis.');
            if folder_name % if folder_name == 0 then cancel was clicked.
                [~,attribs] = fileattrib(folder_name); % we need read access from here.
                if attribs.UserRead
                    handles.parameters.lesion_img_folder = folder_name;
                else
                    warndlg('You do not have read access to the directory you selected for the lesion files. Adjust the permissions and try again.')
                    changemade = false;                
                end
            else
                changemade = false;
            end
        case 'choosescorefilebutton'
            [FileName,PathName] = uigetfile(fullfile(fileparts(handles.parameters.score_file),'*.csv'),'Select a file with behavioral scores.');
            if FileName
                scorefile_name =  fullfile(PathName,FileName);
                handles.parameters.score_file = scorefile_name;
                handles.parameters.control_variable_names = {}; % also clear covariates...  
            else % cancel was clicked.
                changemade = false;
            end
        case 'chooseoutputfolderbutton'
            folder_name = uigetdir(handles.parameters.analysis_out_path,'Choose a folder in which to save this analysis.');
            if folder_name
                [~,attribs] = fileattrib(folder_name); % we need read/write access from here.
                if attribs.UserRead && attribs.UserWrite
                    handles.parameters.analysis_out_path = folder_name;
                else
                    warndlg('You do not have read and write access to the directory you selected to save your output. Adjust the permissions and try again.')
                    changemade = false;                                
                end
            else
                changemade = false;
            end
        case 'scorenamepopupmenu' % User has changed the one_score in question...
            contents = get(gcbo,'string');
            newval = contents{get(gcbo,'value')};
            handles.parameters.run_double_dissociation = false; % by default...
            changemade = true; % ?
            %% Special handling for dissociations by detecting "Dissociation" in the name. This initiates a 
            if contains(newval,'Dissociation')
                % ok now get the user defined dissociation...
                othercontents = contents(1:end-1); % the last one is  Dissocation...

                % remove previously selected covariates
                selected_covariates = get(handles.realcovariateslistbox,'string');
                selected_covariates(get(handles.realcovariateslistbox,'value'));
                othercontents = setdiff(othercontents,selected_covariates); 

                [indx1,tf1] = listdlg('PromptString',{'Select Behavior 1:'},'SelectionMode','single','ListString',othercontents);
                if isempty(indx1), return; end % make sure this was filled out
                behav1 = othercontents{indx1};
                residualcontents = setdiff(othercontents,behav1);
                [indx2,tf2] = listdlg('PromptString',{'Select Behavior 2:'},'SelectionMode','single','ListString',residualcontents);
                if isempty(indx2), return; end % make sure both were filled out.
                behav2 = residualcontents{indx2};
                %['Dissocate: ' behav1 ' & ' behav2]
                handles.parameters.double_dissociation_behaviors = {behav1,behav2};
                handles.parameters.run_double_dissociation = true;
            end

            if any(strcmp(newval,handles.parameters.control_variable_names))
                warndlg('This variable is already chosen as a covariate. If you''d like to use it as the outcome of interest, remove it as a covariate.')
                changemade = false;
            else
                handles.parameters.score_name = newval;
            end
        case 'analysisnameeditbox'
            handles.parameters.analysis_name = get(gcbo,'string');
        case 'lesionthresholdeditbox'
            str = str2num(get(gcbo,'string')); %TO ADD: also make sure this doesn''t exceed the number of lesions available in the data
            if isempty(str) || ~isint(str) 
                warndlg('Input must be a positive integer.');
                changemade = false;
            else % update the parameter value.
                handles.parameters.lesion_thresh = str;
            end
        case 'computebetamapcheckbox'
            handles.parameters.beta_map  = get(gcbo,'value');
        case 'computesensitivitymapcheckbox'
            handles.parameters.sensitivity_map = get(gcbo,'value');
        case 'cluster_voxelwise_p_editbox' 
            str = str2num(get(gcbo,'string'));
            if isempty(str) || str <= 0 || str >= 1 
                changemade = false;
                warndlg('Input must be a number between 0 and 1.');
            else % update the parameter value.
                handles.parameters.voxelwise_p = str;
            end
        case 'clusterwisepeditbox'
            str = str2num(get(gcbo,'string'));
            if isempty(str) || str <= 0 || str >= 1
                warndlg('Input must be a number between 0 and 1.');
                changemade = false;
            else % update the parameter value.
                handles.parameters.clusterwise_p = str;
            end
        case 'npermutationseditbox' % This is voxelwise permutations
            str = str2num(get(gcbo,'string'));
            if isempty(str) || str<=0 || ~isint(str)
                changemade = false;
                warndlg('Input must be a positive integer.');
            else % update the parameter value.
                handles.parameters.PermNumVoxelwise = str;
            end
        case 'permutationtestingcheckbox' % enable/disable permutation testing.
            handles.parameters.DoPerformPermutationTesting = get(gcbo,'value');
        case 'do_cfwer_checkbox'
            handles.parameters.do_CFWER = get(gcbo,'value');
        case 'cfwer_v_value_editbox'
            str = str2num(get(gcbo,'string'));
            if isempty(str) || str <= 0 || ~isint(str)
                changemade=false;
                warndlg('Input must be a positive integer.');
            else % update the parameter value.
                handles.parameters.cfwer_v_value = str; % note that this is cubic millimeters and will later be converted into # voxels (in read_lesioned_images)
            end
        case 'cfwer_p_value_editbox'
            str = str2num(get(gcbo,'string'));
            if isempty(str) || str<=0 || str >= 1 % not a valid p value...
                changemade=false;
                warndlg('Input must be a positive number less than 1.');
            else % update the parameter value.
                handles.parameters.cfwer_p_value = str;
            end
         case 'atlasparcellation_menu_option_parent' % 2021 
             if ~isfield(handles.parameters,'use_atlas_parcellation') % for earlier version compatibliity that didn't have this field
                 handles.parameters.use_atlas_parcellation = false;
                 handles.parameters.analysis_parcellation_file = '';
             end

            if strcmp(handles.parameters.analysis_parcellation_file,'')
                maskstring = 'Select Parcellation...';
            else
                maskstring = handles.parameters.analysis_parcellation_file;
            end
            
            handles.dataparcellation_to_use_menu_option.Text = maskstring; % what mask are we using?
            if handles.parameters.use_atlas_parcellation % then populate the menu item and set the checkbox...
                handles.do_not_apply_parcellation_menu_option.Checked = false; % "Voxelwise"
                handles.dataparcellation_to_use_menu_option.Checked = true;
            else % don't use parcellation
                handles.do_not_apply_parcellation_menu_option.Checked = true; % "Voxelwise"
                handles.dataparcellation_to_use_menu_option.Checked = false;
            end
         case 'analysismask_menu_option_parent' % 2021
             if strcmp(handles.parameters.analysis_mask_file,'')
                 maskstring = 'Select mask...';
             else
                 maskstring = handles.parameters.analysis_mask_file;
             end
             handles.datamask_to_use_menu_option.Text = maskstring;

             if handles.parameters.use_analysis_mask % then populate the menu item and set the checkbox...
                 handles.do_not_apply_datamask_menu_option.Checked = false;
                 handles.datamask_to_use_menu_option.Checked = true;
             else % don't use mask.
                 handles.do_not_apply_datamask_menu_option.Checked = true;
                 handles.datamask_to_use_menu_option.Checked = false;
             end
        case 'do_not_apply_parcellation_menu_option'
             handles.parameters.use_atlas_parcellation = false;
             changemade = true;
        case 'dataparcellation_to_use_menu_option'
            handles.parameters.use_atlas_parcellation = true;
            [FileName,PathName] = uigetfile('*.nii','Select an atlas parcellation file within which to run the analysis.');
            filepath = fullfile(PathName,FileName);
            if ~exist(filepath,'file'), return; end
            handles.parameters.analysis_parcellation_file = filepath; 
            changemade = true;
        case 'do_not_apply_datamask_menu_option' 
             handles.parameters.use_analysis_mask = false;
             changemade = true;
        case 'datamask_to_use_menu_option'
            handles.parameters.use_analysis_mask = true;
            [FileName,PathName] = uigetfile('*.nii','Select a mask file within which to run the analysis.');
            filepath = fullfile(PathName,FileName);
            if ~exist(filepath,'file'), return; end
            handles.parameters.analysis_mask_file = filepath; 
            changemade = true;
        otherwise
            warndlg(['Unknown callback object ' get(gcbo,'tag') ' - has someone modified the code?'])
    end

    if changemade % then set to not saved...
        handles.parameters.is_saved = 0;
    end

    if ~handles.parameters.is_saved % then the analysis configuration was modified, so it hasn't been completed in its current state.
        handles.parameters.analysis_is_completed = 0; % set "is completed" to 0 (in its current state) so the user can click run button, and cannot click show output button.
        handles = PopulateGUIFromParameters(handles);
    end

    guidata(hObject, handles); % Update handles structure
    UpdateTitleBar(handles); % update title bar to show if we have any changes made.

function doignore = IgnoreUnsavedChanges(handles)
    if ~isfield(handles,'parameters') % something bad probably happened with gui initiation.
    	doignore = 1; 
    elseif isfield(handles.parameters,'is_saved') && ~handles.parameters.is_saved % then prompt if the user wants to continue or cancel.
        choice = questdlg('If you continue you will lose unsaved changes to this analysis configuration.', 'Unsaved Changes', 'Continue Anyway','Cancel','Cancel');
        switch choice
            case 'Continue Anyway', doignore = 1;
            case 'Cancel', doignore = 0;
        end
    else % don't hang 
        doignore = 1;
    end

% --- Outputs from this function are returned to the command line.
function varargout = svrlsmgui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function newmenu_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
    if IgnoreUnsavedChanges(handles)
        handles.parameters = GetDefaultParameters(handles);
        handles = PopulateGUIFromParameters(handles);
    end

function openmenu_Callback(hObject, eventdata, handles)
    if IgnoreUnsavedChanges(handles)
        [FileName,PathName] = uigetfile('*.mat','Select an SVRLSMGUI parameters file.');
        filepath = fullfile(PathName,FileName);
        if ~exist(filepath,'file'), return; end
        handles = LoadParametersFromSVRLSMFile(handles,hObject,filepath);
    end

function closemenu_Callback(hObject, eventdata, handles)
    if IgnoreUnsavedChanges(handles)
        delete(gcf)
    end

function savemenu_Callback(hObject, eventdata, handles)
    if exist(handles.parameters.parameter_file_name,'file')
        handles = SaveSVRLSMGUIFile(handles,hObject); % do the actual save.
    else
        saveasmenu_Callback(hObject, eventdata, handles)
    end
   
function saveasmenu_Callback(hObject, eventdata, handles)
    if exist(handles.parameters.parameter_file_name,'file')
        defaultsavename = handles.parameters.parameter_file_name; % fileparts(handles.parameters.parameter_file_name);
    else
        defaultsavename = fullfile(pwd,'Unnamed.mat');
    end
    [file,path] = uiputfile('*.mat','Save SVRLSM GUI parameter file as...',defaultsavename);
    if file == 0, return; end % then cancel was pressed

    handles.parameters.parameter_file_name = fullfile(path,file);
    handles = SaveSVRLSMGUIFile(handles,hObject); % do the actual save.

function quitmenu_Callback(hObject, eventdata, handles)
    close(gcf) % to trigger close request fcn which handles unsaved changes...
    
function onlinehelpmenu_Callback(hObject, eventdata, handles)
    web('https://github.com/atdemarco/svrlsmgui/wiki')

function figure1_CreateFcn(hObject, eventdata, handles)

function aboutmenu_Callback(hObject, eventdata, handles)
    helpstr = ['SVRLSM GUI ' num2str(handles.parameters.gui_version) ', Andrew DeMarco 2017-2021, based on Zhang et al. (2014)'];
    helpdlg(helpstr,'About');

function runanalysisbutton_Callback(hObject, eventdata, handles)
    [success,handles] = RunAnalysis(hObject,eventdata,handles); % now returns handles 10/26/17

    set(handles.runanalysisbutton,'Enable','on'); %set(handles.cancelanalysisbutton,'visible','off')

    % Re-enable interface...
    set(get(handles.permutationtestingpanel,'children'),'enable','on');
    set(get(handles.analysispreferencespanel,'children'),'enable','on');
    set(get(handles.covariatespanel,'children'),'enable','on');
    
    switch success
        case 1 % success
            handles.parameters.analysis_is_completed = 1; % Completed...
            handles = UpdateProgress(handles,'Analysis has completed successfully.',1);
        case 0 % failure
            handles.parameters.analysis_is_completed = 2; % Error...
            handles = UpdateProgress(handles,'Analysis encountered an error and did not complete...',1);
            set(handles.interrupt_button,'enable','off') % disable.
            rethrow(handles.error)
        case 2 % interrupted
            handles.parameters.analysis_is_completed = 2; % Error...
            handles = UpdateProgress(handles,'Analysis was interrupted by user...',1);            
    end
    
    guidata(hObject, handles); % Update handles structure
    handles = PopulateGUIFromParameters(handles); % refresh gui so we can enable/disable control variable as necessary.

% Select in the dropdown list the selected item.
function realcovariateslistbox_Callback(hObject, eventdata, handles)
    if isempty(get(hObject,'String')), return; end % hack for error
    contents = cellstr(get(hObject,'String'));
    val = contents{get(hObject,'Value')};
    dropdownoptions = get(handles.potentialcovariateslist,'string');
    set(handles.potentialcovariateslist,'value',find(strcmp(val,dropdownoptions)))

function optionsmenu_Callback(hObject, eventdata, handles)
    yn = {'off','on'};
    set(handles.parallelizemenu,'Checked',yn{1+handles.parameters.parallelize}) % is parallelization selected by user?
    set(handles.parallelizemenu,'Enable',yn{1+handles.details.can_parallelize}) % can we parallelize on this platform?

function save_perm_data_Callback(hObject, eventdata, handles) % update the subitems with checkboxes
    yn = {'off','on'};
    set(handles.save_post_clusterwise_thresholded,'Checked',yn{1+handles.parameters.SavePostClusterwiseThresholdedPermutations})
    set(handles.save_post_vox_thresh,'Checked',yn{1+handles.parameters.SavePostVoxelwiseThresholdedPermutations})
    set(handles.save_pre_thresh,'Checked',yn{1+handles.parameters.SavePreThresholdedPermutations})
    % the cfwer files...
    set(handles.save_unthresholded_pmaps_cfwer,'Checked',yn{1+handles.parameters.SaveNullPMapsPreThresholding})
    set(handles.save_thresholded_pmaps_cfwer,'Checked',yn{1+handles.parameters.SaveNullPMapsPostThresholding})

function svrmenu_Callback(hObject, eventdata, handles)
        set(handles.use_lib_svm,'checked','off')
        set(handles.use_matlab_svr,'checked','on')

    if handles.details.stats_toolbox
        set(handles.use_matlab_svr,'enable','on')
    else
        set(handles.use_matlab_svr,'enable','off')
    end
    set(handles.use_lib_svm,'enable','off')

function cost_menu_Callback(hObject, eventdata, handles)
    msg = sprintf('Enter new parameter value for Cost/BoxConstraint (default = %0.1f)',handles.parameters.svr_defaults.cost);
    answer = inputdlg(msg,'Cost Parameter',1,{num2str(handles.parameters.cost)});
    if isempty(answer), return; end % cancel pressed
    numval = str2num(answer{1}); %#ok<*ST2NM>
    if isnumeric(numval) && ~isempty(numval)
        handles.parameters.cost = numval;
        handles.parameters.is_saved = 0;
        guidata(hObject, handles);
        handles = PopulateGUIFromParameters(handles);
    end
   
function gamma_menu_Callback(hObject, eventdata, handles)
    % This function used to switch out the name Gamma and Sigma depending on SVR implementation (matlab vs libsvm) - but in v2 it does not support libsvm.
    defaultval = handles.parameters.svr_defaults.sigma; 
    msg = sprintf('Enter new parameter value for %s (default = %0.1f)','Sigma',defaultval);
    oldval = handles.parameters.sigma;
    answer = inputdlg(msg,'Sigma Parameter',1,{num2str(oldval)}); 
    if isempty(answer), return; end % cancel pressed
    numval = str2num(answer{1});
    if isnumeric(numval) && ~isempty(numval)
        handles.parameters.sigma = numval;
        handles.parameters.is_saved = 0;
        guidata(hObject, handles);
        handles = PopulateGUIFromParameters(handles);
    end
    
function open_batch_job_Callback(hObject, eventdata, handles)
    folder_name = uigetdir(pwd,'Choose a folder containing .mat config files of your analyses.');
    if ~folder_name, return; end % if folder_name == 0 then cancel was clicked.
    files = dir(fullfile(folder_name,'*.mat'));
    fname = {files.name};
    [s,v] = listdlg('PromptString','Choose the analyses to run:','SelectionMode','multi','ListString',fname);
    if ~v, return; end % cancelled..
    for f = 1:numel(s)
        curs=s(f);
        curfname = fname{curs};
        curfile = fullfile(folder_name,curfname);
%         try % so one or more can fail without stopping them all.
            handles = UpdateProgress(handles,['Batch: Starting file ' curfname '...'],1);
            [success,handles] = RunAnalysisNoGUI(curfile,handles); % second parm, handles, allows access to gui elements, etc.
            handles = UpdateProgress(handles,['Batch: Finished file ' curfname '.'],1);
            
            switch success
                case 1 % success
                    handles.parameters.analysis_is_completed = 1; % Completed...
                    handles = UpdateProgress(handles,['Batch: Finished file successfully ' curfname '.'],1);
                case 0 % failure
                    handles.parameters.analysis_is_completed = 2; % Error...
                    %handles = UpdateProgress(handles,'Analysis encountered an error and did not complete...',1);
                    handles = UpdateProgress(handles,['Batch: Analysis encountered an error and did not complete: ' curfname '.'],1);
                    if isfield(handles,'interrupt_button')
                        set(handles.interrupt_button,'enable','off') % disable.
                    end
                    rethrow(handles.error)
            end

%         catch
%             disp('return from open_batch_job in svrlsmgui() with error')
%             msg = ['A batch job specified by file ' fname{curs} ' encountered an error and was aborted.'];
%             warning(msg)
%         end
    end
    
    handles = UpdateProgress(handles,'Batch: All batch jobs done.',1);            

function figure1_CloseRequestFcn(hObject, eventdata, handles)
    if IgnoreUnsavedChanges(handles), delete(hObject); end

function checkforupdates_Callback(hObject, eventdata, handles)
    disp('Option no longer supported.')
    
%% callbacks -- replace these in the future.
function controlvariablepopupmenu_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function computebetamapcheckbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);
function computesensitivitymapcheckbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);
function analysisnameeditbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function permutation_unthresholded_checkbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function permutation_voxelwise_checkbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function permutation_largest_cluster_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function progresslistbox_Callback(hObject, eventdata, handles)

function potentialcovariateslist_Callback(hObject, eventdata, handles)

function viewresultsbutton_Callback(hObject, eventdata, handles)
    LaunchResultsDirectory(hObject,eventdata,handles);
    
function npermutationsclustereditbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function parallelizemenu_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);
    
function debug_menu_Callback(hObject, eventdata, handles)

function use_matlab_svr_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function interrupt_button_Callback(hObject, eventdata, handles) % attempt to interrupt an ongoing analysis
    set(hObject,'enable','off') % so user doesn't click a bunch...
%     poolobj = gcp('nocreate'); % < this doesn't stop, it just relaunches the parpool...
%     delete(poolobj); % try to stop what's happening if there's parallelization happening
    set(gcf,'userdata','cancel')
    guidata(hObject, handles); % Update handles structure so it saves...

function no_optimize_menu_choice_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function correlation_optimize_menu_choice_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);
    
function resubloss_optimize_menu_choice_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function predictbehavior_optimize_menu_choice_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function optimize_menu_Callback(hObject, eventdata, handles)
    set(get(hObject,'Children'),'Checked','off') % uncheck all child menus
    cur_optim_string = CurrentOptimString(handles.parameters);
    set(handles.current_optimization_menu_option,'Label',cur_optim_string)
    opts = [handles.parameters_to_optimize_menu handles.search_strategy_menu_option handles.objective_function_menu_option handles.crossvalidation_parent_menu handles.optimization_is_verbose_menu];
    if ~handles.parameters.optimization.do_optimize % then no optimization...
        set(handles.no_optimize_menu_choice,'Checked','on')
        set(opts,'enable','off') % Visual cue that these don't apply when optimization is off
    else 
        set(handles.current_optimization_menu_option,'Checked','on')
        set(opts,'enable','on') % Visual cue that these don't apply when optimization is off
    end
     if handles.parameters.optimization.verbose_during_optimization 
         set(handles.optimization_is_verbose_menu,'Checked','on')
     end

function current_optimization_menu_option_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function parameters_menu_Callback(hObject, eventdata, handles)
    do_opt = handles.parameters.optimization.do_optimize; % for convenience

    if do_opt && handles.parameters.optimization.params_to_optimize.cost
        label = ['Cost: optimize (' num2str(handles.parameters.optimization.params_to_optimize.cost_range(1)) ' - ' num2str(handles.parameters.optimization.params_to_optimize.cost_range(2)) ')'];
        set(handles.cost_menu,'label',label,'enable','off');
    else
        set(handles.cost_menu,'label',['Cost: ' num2str(handles.parameters.cost) myif(handles.parameters.svr_defaults.cost == handles.parameters.cost,' (default)','')],'enable','on');
    end

    if do_opt && handles.parameters.optimization.params_to_optimize.sigma
        oldmin = handles.parameters.optimization.params_to_optimize.sigma_range(1);
        oldmax = handles.parameters.optimization.params_to_optimize.sigma_range(2);
        label = ['Sigma: optimize (' num2str(oldmin) ' - ' num2str(oldmax) ')'];
        set(handles.gamma_menu,'label',label,'enable','off');
    else 
        cursigma = handles.parameters.sigma; 
        set(handles.gamma_menu,'label',['Sigma : ' num2str(cursigma) myif(handles.parameters.svr_defaults.sigma == handles.parameters.sigma,' (default)','')],'enable','on');
    end
    
    if do_opt && handles.parameters.optimization.params_to_optimize.epsilon
        label = ['Epsilon: optimize (' num2str(handles.parameters.optimization.params_to_optimize.epsilon_range(1)) ' - ' num2str(handles.parameters.optimization.params_to_optimize.epsilon_range(2)) ')'];
        set(handles.epsilon_menu,'label',label,'enable','off');
    else
        set(handles.epsilon_menu,'label',['Epsilon: ' num2str(handles.parameters.epsilon) myif(handles.parameters.svr_defaults.epsilon == handles.parameters.epsilon,' (default)','')],'enable','on');       
    end

    if do_opt && handles.parameters.optimization.params_to_optimize.standardize
        set(handles.standardize_menu,'label','Standardize: optimize (yes/no)','enable','off');
    else
        set(handles.standardize_menu,'label',['Standardize: ' myif(handles.parameters.standardize,'true','false') myif(handles.parameters.svr_defaults.standardize == handles.parameters.standardize,' (default)','')],'enable','on');
    end
    
function epsilon_menu_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);
    
function standardize_menu_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);
    
%% Objective function...
function objective_function_menu_option_Callback(hObject, eventdata, handles)
    set(get(hObject,'Children'),'Checked','off') % uncheck all child menus
    set([handles.predictbehavior_optimize_menu_choice handles.correlation_optimize_menu_choice],'enable','off')
    switch handles.parameters.optimization.objective_function
        case 'Predict Behavior'
            set(handles.predictbehavior_optimize_menu_choice,'Checked','on')
        case 'Maximum Correlation'
            set(handles.correlation_optimize_menu_choice,'Checked','on')
        case 'Resubstitution Loss'
            set(handles.resubloss_optimize_menu_choice,'Checked','on')
        otherwise
            error('Unknown optimization objective function.')
    end
    
%% Search strategy ...
function search_strategy_menu_option_Callback(hObject, eventdata, handles)
    set(get(hObject,'Children'),'Checked','off') % uncheck all child menus
    switch handles.parameters.optimization.search_strategy
        case 'Bayes Optimization'
            set(handles.bayes_optimization_menu_choice,'Checked','on')
        case 'Grid Search'
            set(handles.gridsearch_option,'Checked','on')
        case 'Random Search'
            set(handles.random_search_menu_option,'Checked','on')
        otherwise
            error('Unknown optimization objective function.')
    end

function parameters_to_optimize_menu_Callback(hObject, eventdata, handles)
    set(get(hObject,'children'),'checked','off')
    
    if handles.parameters.optimization.params_to_optimize.cost
        label = ['Cost (' num2str(handles.parameters.optimization.params_to_optimize.cost_range(1)) ' - ' num2str(handles.parameters.optimization.params_to_optimize.cost_range(2)) ')'];
        set(handles.do_optimize_cost_menu,'label',label,'checked','on')
    end
    
    curmin = handles.parameters.optimization.params_to_optimize.sigma_range(1); % convert to sigma if necessary...
    curmax = handles.parameters.optimization.params_to_optimize.sigma_range(2); % convert to sigma if necessary...
    label = ['Sigma (' num2str(curmin) ' - ' num2str(curmax) ')'];  % convert to sigma if necessary...
    set(handles.do_optimize_gamma_menu,'label',label,'checked',myif(handles.parameters.optimization.params_to_optimize.sigma,'on','off'))
    
    if handles.parameters.optimization.params_to_optimize.epsilon
        label = ['Epsilon (' num2str(handles.parameters.optimization.params_to_optimize.epsilon_range(1)) ' - ' num2str(handles.parameters.optimization.params_to_optimize.epsilon_range(2)) ')'];
        set(handles.do_optimize_epsilon_menu,'label',label,'checked','on')
    end
    
    if handles.parameters.optimization.params_to_optimize.standardize
        label = 'Standardize (yes/no)';
        set(handles.do_optimize_standardize_menu,'label',label,'checked','on')
    end
    
function crossvalidation_parent_menu_Callback(hObject, eventdata, handles)
    set(get(hObject,'children'),'checked','off')
    set(handles.crossval_menu_option_kfold,'label',['K-Fold: ' num2str(handles.parameters.optimization.crossval.nfolds) ' folds' myif(handles.parameters.optimization.crossval.nfolds == handles.parameters.optimization.crossval.nfolds_default,' (default)','')])
    if ~handles.parameters.optimization.crossval.do_crossval
        set(handles.crossval_menu_option_none,'checked','on')
        set(handles.do_repartition_menu_option,'enable','off')
    else
        if handles.parameters.optimization.crossval.repartition 
            set(handles.do_repartition_menu_option,'checked','on','enable','on')
        end
        if strcmp(handles.parameters.optimization.crossval.method,'kfold')
            set(handles.crossval_menu_option_kfold,'checked','on')
        else
            error('Unknown crossvalidation option string.')
        end
    end
    
function parent_cache_menu_Callback(hObject, eventdata, handles)
    yn = {'off','on'};
    set(handles.do_use_cache_menu,'Checked',yn{1+handles.parameters.do_use_cache_when_available})
    set(handles.retain_big_binary_file,'Checked',yn{1+handles.parameters.SavePermutationData})
    set(handles.retain_big_binary_pval_file,'Checked',yn{1+handles.parameters.SavePermutationPData})

function output_summary_menu_Callback(hObject, eventdata, handles)
    yn = {'off','on'};
    set(handles.summary_create_summary,'checked',yn{1+handles.parameters.do_make_summary});
    set(handles.summary_narrative_summary,'checked',yn{1+handles.parameters.summary.narrative});
    set(handles.summary_svrbetamap,'checked',yn{1+handles.parameters.summary.beta_map});
    set(handles.summary_voxelwise_thresholded,'checked',yn{1+handles.parameters.summary.voxelwise_thresholded});
    set(handles.summary_clusterwise_thresholded,'checked',yn{1+handles.parameters.summary.clusterwise_thresholded});
    set(handles.summary_cfwerdiagnostics,'checked',yn{1+handles.parameters.summary.cfwer_diagnostics});
    set(handles.model_variablediagnostics,'checked',yn{1+handles.parameters.summary.variable_diagnostics});
    set(handles.summary_clusterstability,'checked',yn{1+handles.parameters.summary.cluster_stability});
    set(handles.summary_parameterassessment,'checked',yn{1+handles.parameters.summary.parameter_assessment});
    set(handles.summary_paramoptimization,'checked',yn{1+handles.parameters.summary.hyperparameter_optimization_record});
    set(handles.summary_lesionoverlap,'checked',yn{1+handles.parameters.summary.lesion_overlap});
    set(handles.summary_prediction_menu,'checked',yn{1+handles.parameters.summary.predictions});
    if handles.parameters.do_make_summary
        set(get(hObject,'children'),'enable','on')
    else
        set(get(hObject,'children'),'enable','off')
        set(handles.summary_create_summary,'enable','on')
    end

function requirements_menu_Callback(hObject, eventdata, handles)
    set(handles.spm12_installed_menu,'checked',myif(handles.details.spm,'on','off')) % this will not specifically detect spm12 though!
    set(handles.parcomp_toolbox_installed_menu,'checked',myif(handles.details.can_parallelize,'on','off'))
    set(handles.stats_toolbox_installed_menu,'checked',myif(handles.details.stats_toolbox,'on','off'))
    set(handles.matlab_version_installed_menu,'checked','on') % what's the requirement?
    set(get(handles.requirements_menu,'children'),'enable','off')

function beta_options_menu_Callback(hObject, eventdata, handles)
    set(handles.ica_lesion_decompose_option,'checked',myif(handles.parameters.beta.do_ica_on_lesiondata,'on','off'))

function crossvalidate_map_parent_menu_Callback(hObject, eventdata, handles)
    set(get(hObject,'children'),'checked','off')
%     parameters.crossval.do_crossval = false; % by default do not do crossvalidated output...
%     parameters.crossval.method = 'kfold';
%     parameters.crossval.nfolds = 5;
%     parameters.crossval.nfolds_default = parameters.crossval.nfolds;
    set(handles.kfold_map_crossvalidation,'label',['K-Fold: ' num2str(handles.parameters.crossval.nfolds) ' folds' myif(handles.parameters.optimization.crossval.nfolds == handles.parameters.crossval.nfolds_default,' (default)','')])
    if ~handles.parameters.crossval.do_crossval
        set(handles.no_map_crossvalidation,'checked','on')
        %set(handles.do_repartition_menu_option,'enable','off')
    else
%         if handles.parameters.optimization.crossval.repartition 
%             set(handles.do_repartition_menu_option,'checked','on','enable','on')
%         end
        if strcmp(handles.parameters.crossval.method,'kfold')
            set(handles.kfold_map_crossvalidation,'checked','on')
        else
            error('Unknown crossvalidation option string.')
        end
    end

function no_map_crossvalidation_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function kfold_map_crossvalidation_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);
