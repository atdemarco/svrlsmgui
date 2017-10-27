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

% Last Modified by GUIDE v2.5 29-Sep-2017 09:52:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @svrlsmgui_OpeningFcn, ...
                   'gui_OutputFcn',  @svrlsmgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function handles = ConfigureSVRLSMGUIOptions(handles)
    handles.options.lesionvolumecorrection = {'Regress on Behavior','Regress on Lesion','Regress on Both','DTLVC','None'};
    handles.options.hypodirection = {'One-tailed (positive)','One-tailed (negative)','Two-tailed'};
    set(handles.lesionvolumecorrectiondropdown,'String',handles.options.lesionvolumecorrection)
    set(handles.hypodirectiondropdown,'String',handles.options.hypodirection)

% --- Executes just before svrlsmgui is made visible.
function svrlsmgui_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject; % Choose default command line output for svrlsmgui
    
    % Hide some debugging controls...
    tohide = [handles.permutation_unthresholded_checkbox handles.permutation_voxelwise_checkbox handles.permutation_largest_cluster handles.save_raw_permutation_data handles.text17 handles.text21 handles.npermutationsclustereditbox];
    set(tohide,'visible','off')
    
    handles = ConfigureSVRLSMGUIOptions(handles);
    
    handles.details = CheckIfNecessaryFilesAreInstalled(handles);
    
    if handles.details.stats_toolbox && handles.details.spm && handles.details.libsvm
        handles = UpdateProgress(handles,'All necessary functions are available...',1);
        handles.parameters = GetDefaultParameters(handles);
        handles = PopulateGUIFromParameters(handles);
    elseif ~handles.details.spm 
        handles = UpdateProgress(handles,'SPM functions not available. Download and/or add SPM to the MATLAB path and relaunch the SVRLSMGUI.',1);
        handles = DisableAll(handles);
    elseif ~handles.details.stats_toolbox && ~handles.details.libsvm
        handles = UpdateProgress(handles,'No SVR algorithm available. Install Statistics Toolbox in MATLAB or compile and install libSVM and relaunch the GUI.',1);
        handles = DisableAll(handles);
    elseif ~handles.details.stats_toolbox && handles.details.libsvm % yes stats toolbox, no libsvm
        handles = UpdateProgress(handles,'Only libSVM is available to compute images. MATLAB''s SVM will not be available.',1);
        handles.parameters = GetDefaultParameters(handles);
        handles = PopulateGUIFromParameters(handles);
    elseif handles.details.stats_toolbox && ~handles.details.libsvm % no stats toolbox, yes libsvm
        handles = UpdateProgress(handles,'Only MATLAB''s Stats Toolbox is available to compute images. libSVM will not be available.',1);
        handles.parameters = GetDefaultParameters(handles);
        handles = PopulateGUIFromParameters(handles);
    end
    
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
    
    handles.parameters.gui_version = 0.08; % version of the the gui

    guidata(hObject, handles); % Update handles structure


function handles = DisableAll(handles)
    set(get(handles.analysispreferencespanel,'children'),'enable','off')
    set(get(handles.covariatespanel,'children'),'enable','off')
    set(get(handles.permutationtestingpanel,'children'),'enable','off')
    msgbox('One or more necessary component is missing from MATLAB''s path. Address the message in the SVRLSMgui window and restart this gui.')

function handles = LoadParametersFromSVRLSMFile(handles,hObject,filepath)
    tmp = load(filepath); 
    handles.parameters = tmp.tosave; % tosave is the name of the variable we save; it's arbitrary but invariant.
    handles = PopulateGUIFromParameters(handles); % show what we've just loaded.    
    guidata(hObject, handles); % Update handles structure
    handles = UpdateProgress(handles,['Loaded ' handles.parameters.parameter_file_name],0);
    
function parameters = GetDefaultParameters(handles)
    % To add: make this a "default" .mat file, configurable in preferences...
    if handles.details.libsvm
        parameters.useLibSVM = 1; % 1 = libSVM and 0 = MATLAB
    else %if handles.details.stats_toolbox 
        parameters.useLibSVM = 0; % 1 = libSVM and 0 = MATLAB
    end

    parameters.svscaling = 100; % defaults to max (100th percentile) - added 9/29/17
    
    parameters.control_variable_name = []; % vestigial.
    parameters.tails = handles.options.hypodirection{1};
    parameters.datetime_run = []; % when the analysis was run.
    parameters.datetime_save = []; % when the analysis was run.
    parameters.analysis_is_completed = 0;
    parameters.parameter_file_name = ''; % identity of file if parameters were opened from a file.
    parameters.parallelize = false; % do not parallelize...
	
    mypath = fileparts(which('svrlsmgui'));
    
    parameters.analysis_root = fullfile(mypath,'output'); % fullfil(mypath,''/Users/LauraSkipper/Documents/MATLAB/andrew/svrlsm/v01';
    parameters.score_file = fullfile(mypath,'default','PNT.csv');
    parameters.score_name = 'Sim_ROI_123';
    parameters.lesion_img_folder = fullfile(mypath,'default','lesion_imgs'); %parameters.analysis_root
    
    parameters.DoPerformPermutationTesting = 1;

    parameters.gamma = 5;
    parameters.cost = 30;

    parameters.lesionvolcorrection = handles.options.lesionvolumecorrection{2}; % DTLVC

    parameters.beta_map = 1; 
    parameters.sensitivity_map = 1; 

    parameters.invert_p_map_flag = 1; % Inverse p-map, i.e., use 1-p instead of p for display on MRIcron.
    parameters.control_variable_names = {}; % None at first ...  

    parameters.lesion_thresh = 10; % The least lesion subject number for a voxel to be considered in the following analysis. 
    parameters.voxelwise_p = 0.005;
    parameters.clusterwise_p = 0.05;
    parameters.PermNumVoxelwise = 10000;
    parameters.PermNumClusterwise = 10000;

    parameters.apply_covariates_to_behavior = 0;
    parameters.apply_covariates_to_lesion = 0;

    parameters.analysis_name = 'Unnamed';
    parameters.analysis_out_path = parameters.analysis_root; % is this a good default?
    parameters.is_saved = 0;
    
    parameters.SavePreThresholdedPermutations = 0;
    parameters.SavePostVoxelwiseThresholdedPermutations = 0;
    parameters.SavePostClusterwiseThresholdedPermutations = 0;
    parameters.SavePermutationData = 0; % leave the giant bin file?
    
    handles = UpdateProgress(handles,'Retrieved default parameters...',1); %#ok<*NASGU>

function handles = UpdateCurrentAnalysis(handles,hObject)
changemade = true; % default
switch get(gcbo,'tag') % use gcbo to see what the cbo is and determine what field it goes to -- and to validate
    case 'parallelizemenu'
        handles.parameters.parallelize = ~handles.parameters.parallelize;
    case 'applycovariatestobehaviorcheckbox'
        handles.parameters.apply_covariates_to_behavior = get(gcbo,'value');
    case 'applycovariatestolesioncheckbox'
        handles.parameters.apply_covariates_to_lesion = get(gcbo,'value');
        handles.parameters.is_saved = 0;
    case 'lesionvolumecorrectiondropdown'
        contents = get(handles.lesionvolumecorrectiondropdown,'string');
        newval = contents{get(handles.lesionvolumecorrectiondropdown,'value')};
        handles.parameters.lesionvolcorrection = newval;
    case 'hypodirectiondropdown'
        contents = get(handles.hypodirectiondropdown,'string');
        newval = contents{get(handles.hypodirectiondropdown,'value')};
        handles.parameters.tails = newval;
    case 'addcovariate'
        % first, what are we trying to add?
        contents = get(handles.potentialcovariateslist,'String'); 
        newcovariate = contents{get(handles.potentialcovariateslist,'Value')};
        if strcmp(newcovariate,handles.parameters.score_name) %dev1 % check if it's our main score name...
            warndlg('This variable is already chosen as the main outcome in the analysis. If you''d like to add it as a covariate, remove it from "Score Name."')
            changemade = false;
        elseif any(strcmp(newcovariate,handles.parameters.control_variable_names)) % only if it's not on our list of covariates already.
            warndlg('This variable is already on the list of covariates. You may not add it twice.')
            changemade = false;
        else % it's new
            handles.parameters.control_variable_names{end+1} = newcovariate;
        end
    case 'removecovariate'
       % what are we trying to remove? 
       contents = get(handles.potentialcovariateslist,'String');
       newcovariate = contents{get(handles.potentialcovariateslist,'Value')};
       if any(strcmp(newcovariate,handles.parameters.control_variable_names)) % only if it IS on our list!
           index_to_remove = strcmp(newcovariate,handles.parameters.control_variable_names);
           handles.parameters.control_variable_names(index_to_remove) = [];
       else
            changemade = false;
       end
    case 'chooselesionfolderbutton'
        folder_name = uigetdir(pwd,'Choose a folder containing lesion files for this analysis.');
        if folder_name % if folder_name == 0 then cancel was clicked.
            handles.parameters.lesion_img_folder = folder_name;
        else
            changemade = false;
        end
    case 'choosescorefilebutton'
        [FileName,PathName] = uigetfile('*.csv','Select a file with behavioral scores.');
        if FileName % if FileName == 0 then cancel was clicked.
            scorefile_name =  fullfile(PathName,FileName);
            handles.parameters.score_file = scorefile_name;
        else
            changemade = false;
        end
    case 'chooseoutputfolderbutton'
        folder_name = uigetdir(pwd,'Choose a folder in which to save this analysis.');
        if folder_name % if folder_name == 0 then cancel was clicked.
            handles.parameters.analysis_out_path = folder_name;
        else
            changemade = false;
        end
    case 'scorenamepopupmenu' % User has changed the one_score in question...
        contents = get(gcbo,'string');
        newval = contents{get(gcbo,'value')};
        if any(strcmp(newval,handles.parameters.control_variable_names))
            warndlg('This variable is already chosen as a covariate. If you''d like to use it as the outcome of interest, remove it as a covariate.')
            changemade = false;
        else
            handles.parameters.score_name = newval;
        end
    case 'analysisnameeditbox'
        handles.parameters.analysis_name = get(gcbo,'string');
    case 'lesionthresholdeditbox'
        str = get(gcbo,'string'); %TO ADD: also make sure this doesn''t exceed the number of lesions available in the data
        if isempty(str2num(str)) % then it's not a valid number...
            set(gcbo,'string','10'); % change to default....
            warndlg('Input must be numerical');
            changemade = false;
        else % update the parameter value.
            handles.parameters.lesion_thresh = str2num(str);
        end
    case 'computebetamapcheckbox'
        handles.parameters.beta_map  = get(gcbo,'value');
    case 'computesensitivitymapcheckbox'
        handles.parameters.sensitivity_map = get(gcbo,'value');
    case 'invertpmapcheckbox'
        handles.parameters.invert_p_map_flag = get(gcbo,'value');
    case 'cluster_voxelwise_p_editbox' 
        str = get(gcbo,'string');
        newval = str2num(str);
        if isempty(newval) || newval <= 0 || newval >= 1 % then it's not a valid number...
            %set(gcbo,'string','.005');
            warndlg('Input must be a number between 0 and 1.');
            %handles = PopulateGUIFromParameters(handles);           
            changemade = false;
        else % update the parameter value.
            handles.parameters.voxelwise_p = str2num(str);
        end
    case 'clusterwisepeditbox'
        str = get(gcbo,'string');
        newval = str2num(str);
        if isempty(newval) || newval <= 0 || newval >= 1 % then it's not a valid value...
            %set(gcbo,'string','.05');
            warndlg('Input must be a number between 0 and 1.');
            %handles = PopulateGUIFromParameters(handles);           
            changemade = false;
        else % update the parameter value.
            %handles.parameters.is_saved = 0;
            handles.parameters.clusterwisepeditbox = str2num(str);
        end
    case 'sv_scaling_95th_percentile'
        handles.parameters.svscaling = 95;
    case 'sv_scaling_99th_percentile'
        handles.parameters.svscaling = 99;
    case 'maxsvscaling'
        handles.parameters.svscaling = 100; % default
    case 'npermutationseditbox' % This is voxelwise permutations
        str = get(gcbo,'string');
        if isempty(str2num(str)) % then it's not a valid number...
            set(gcbo,'string','10000');
            warndlg('Input must be numerical');
%         elseif str2num(str) < handles.parameters.PermNumClusterwise
%             warndlg('Voxelwise permutations cannot be fewer than clusterwise permutations.');
%             % handles = PopulateGUIFromParameters(handles);           
%             changemade = false;
        else % update the parameter value.
            %handles.parameters.is_saved = 0;
            handles.parameters.PermNumVoxelwise = str2num(str);
            % ClusterWise is overridden and set at analysis runtime.
        end
    case 'npermutationsclustereditbox' % This is clusterwise permutations
        warning('This control is disabled and should not be used.')
%         str = get(gcbo,'string');
%         if isempty(str2num(str)) % then it's not a valid number...
%             warndlg('Input must be numerical');
%             %handles = PopulateGUIFromParameters(handles);
%             set(gcbo,'string',num2str(handles.parameters.PermNumClusterwise));
%             changemade = false;
%         elseif str2num(str) > handles.parameters.PermNumVoxelwise %#ok<*ST2NM>
%             warndlg('Cannot have more clusterwise permutations than voxelwise permutations.');
%             %handles = PopulateGUIFromParameters(handles);
%             set(gcbo,'string',num2str(handles.parameters.PermNumClusterwise));
%             changemade = false;
%         else % update the parameter value.
%             handles.parameters.PermNumClusterwise = str2num(str);
%         end
%         %handles.parameters.is_saved = 0;
    case 'permutation_unthresholded_checkbox'
        handles.parameters.SavePreThresholdedPermutations = get(gcbo,'value');
        %handles.parameters.is_saved = 0;
    case 'permutation_voxelwise_checkbox'
        handles.parameters.SavePostVoxelwiseThresholdedPermutations = get(gcbo,'value');
        %handles.parameters.is_saved = 0;
    case 'permutation_largest_cluster'
        handles.parameters.SavePostClusterwiseThresholdedPermutations = get(gcbo,'value');
        %handles.parameters.is_saved = 0;
    case 'save_raw_permutation_data'
        handles.parameters.SavePermutationData = get(gcbo,'value'); % leave the giant bin file?
        %handles.parameters.is_saved = 0;
    case 'permutationtestingcheckbox' % enable/disable permutation testing.
        handles.parameters.DoPerformPermutationTesting = get(gcbo,'value'); % leave the giant bin file?
        %handles.parameters.is_saved = 0;
    otherwise
        warndlg('Unknown callback object - has someone modified the code?')
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

function details = CheckIfNecessaryFilesAreInstalled(handles)
    % Try to add private functions to path....
    mypath = fileparts(which('svrlsmgui'));
    addpath(fullfile(mypath,'functions')) % when called 'private' this was unnecessary, but let's add it.
    % make sure nifti toolbox is on path.
    addpath(fullfile(mypath,'functions','nifti')) % nifti toolbox
    addpath(fullfile(mypath,'functions','libsvm-3.18')) % libsvm...
    addpath(fullfile(mypath,'functions','libsvm-3.18','matlab')) % libsvm's matlab directory...
    
    handles = UpdateProgress(handles,'Checking if necessary files are installed...',1);

    spm_found = which('spm','-all');
    if isempty(spm_found)
        handles = UpdateProgress(handles,'Warning: SPM is not installed and/or visible on MATLAB''s path.',1);
        details.spm = 0;
    else
        handles = UpdateProgress(handles,'SPM is installed and visible on MATLAB''s path.',1);
        details.spm = 1;
    end
    
    svmtrain_found = which('svmtrain','-all'); % nb: svmtrain is also the name of a statistics toolbox function.
    correct_svmtrains = cellfun(@(x) strfind(x,'libsvm'),svmtrain_found,'Uni',false);
    
    if ~isempty(correct_svmtrains{1})
        handles = UpdateProgress(handles,'libsvm is installed and visible on MATLAB''s path.',1);
        details.libsvm = 1;
    elseif ~all(isempty(correct_svmtrains)) % one of the function is right, but not at top of path
        handles = UpdateProgress(handles,'libsvm may be installed but it appears to be overloaded MATLAB''s path (svmtrain.m?).',1);
        details.libsvm = 0;
    else
        handles = UpdateProgress(handles,'libsvm is either not installed or not visible on MATLAB''s path.',1);
        details.libsvm = 0;
    end
    
    % Now for the MATLAB statistics svr functions
    
    matlab_stats_found = license('checkout','statistics_toolbox');
    if matlab_stats_found
        handles = UpdateProgress(handles,'MATLAB Statistics Toolbox license found.',1);
        details.stats_toolbox = 1;
    else
        handles = UpdateProgress(handles,'MATLAB Statistics Toolbox license not found.',1);
        details.stats_toolbox = 0;
    end
    
function doignore = IgnoreUnsavedChanges(handles)
    if isfield(handles.parameters,'is_saved') && ~handles.parameters.is_saved % then prompt if the user wants to continue or cancel.
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

function filemenu_Callback(hObject, eventdata, handles)
% nothing goes here.

function helpmenu_Callback(hObject, eventdata, handles)
% nothing goes here.

function newmenu_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
    if IgnoreUnsavedChanges(handles)
        handles.parameters = GetDefaultParameters(handles);
        handles = PopulateGUIFromParameters(handles);
    end

function openmenu_Callback(hObject, eventdata, handles)
    if IgnoreUnsavedChanges(handles)
        [FileName,PathName] = uigetfile('*.mat','Select an SVRLSMGUI parameters file.');
        filepath = fullfile(PathName,FileName);
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
        defaultsavename = fileparts(handles.parameters.parameter_file_name);
    else
        defaultsavename = [handles.parameters.analysis_name '.mat'];
    end
    [file,path] = uiputfile(defaultsavename,'Save SVRLSM GUI parameter file as...');
    if file == 0 % then cancel was pressed
        return;
    end

    handles.parameters.parameter_file_name = fullfile(path,file);
    handles.parameters.parameter_file_name
    handles = SaveSVRLSMGUIFile(handles,hObject); % do the actual save.
    
function handles = SaveSVRLSMGUIFile(handles,hObject)
    handles.parameters.is_saved = 1;
    handles.parameters.datetime_save = datetime;
    tosave = handles.parameters; % arbitrary, but I think it was giving me issues saving the field parameters from struct handles...
    save(handles.parameters.parameter_file_name,'tosave') % write the file
    handles = UpdateProgress(handles,['Saved ' handles.parameters.parameter_file_name],1);
    handles = LoadParametersFromSVRLSMFile(handles,hObject,handles.parameters.parameter_file_name);

function preferencesmenu_Callback(hObject, eventdata, handles)
    prefs_to_update = preferences(handles.parameters,handles.details);
    if ~isempty(prefs_to_update) % then prefs window "cancel" button was not clicked - update vals.
        handles.parameters.gamma = prefs_to_update.gamma;
        handles.parameters.cost = prefs_to_update.cost;
        handles.parameters.useLibSVM = prefs_to_update.useLibSVM;
        handles.parameters.is_saved = 0;
        handles = PopulateGUIFromParameters(handles);
    end
    guidata(hObject, handles); % Update handles structure

function quitmenu_Callback(hObject, eventdata, handles)
    if IgnoreUnsavedChanges(handles), close(gcf); end
    
function scorefileeditbox_Callback(hObject, eventdata, handles)
% This should never trigger.

function scorefileeditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%ok
function choosescorefilebutton_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

% ok
function scorenamepopupmenu_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function scorenamepopupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ok
function controlvariablepopupmenu_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function controlvariablepopupmenu_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lesionfoldereditbox_Callback(hObject, eventdata, handles)
% this should never trigger

function lesionfoldereditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ok
function chooselesionfolderbutton_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

% ok
function computebetamapcheckbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

% ok
function computesensitivitymapcheckbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

% ok
function invertpmapcheckbox_Callback(hObject, eventdata, handles) %#ok<*INUSL>
handles = UpdateCurrentAnalysis(handles,hObject);

% ok
function analysisnameeditbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function analysisnameeditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ok
function onlinehelpmenu_Callback(hObject, eventdata, handles)
%web('https://docs.google.com/a/email.arizona.edu/document/d/1McqkxVTmhzkE2tesAk_d4zp_qS_xt8WWGQlZ1lUONmY/edit?usp=sharing','-browser')
web('https://github.com/atdemarco/svrlsmgui/wiki')

function figure1_CreateFcn(hObject, eventdata, handles)

% ok
function clusterwisepeditbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function clusterwisepeditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ok
function lesionthresholdeditbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function lesionthresholdeditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function npermutationseditbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function save_raw_permutation_data_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function permutation_unthresholded_checkbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function permutation_voxelwise_checkbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function permutation_largest_cluster_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function npermutationseditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function progresslistbox_Callback(hObject, eventdata, handles)

function generictext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function progresslistbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function aboutmenu_Callback(hObject, eventdata, handles)
    helpstr = ['SVRLSM GUI ' num2str(handles.parameters.gui_version) ', Andrew DeMarco 2017, based on Zhang et al. (2014)'];
    helpdlg(helpstr,'About');

function runanalysisbutton_Callback(hObject, eventdata, handles)
    handles.parameters.datetime_run = date; % when the analysis was run.
    handles.parameters.PermNumClusterwise = handles.parameters.PermNumVoxelwise; % override the user so that these two values are the same.
    
    [success,handles] = RunAnalysis(hObject,eventdata,handles); % now returns handles 10/26/17
    
    if success
       handles.parameters.analysis_is_completed = 1; % Completed...
       handles = UpdateProgress(handles,'Analysis has completed successfully.',1);
    else
        handles.parameters.analysis_is_completed = 2; % Error...
        handles = UpdateProgress(handles,'Analysis encountered an error and did not complete...',1);
    end
    
    guidata(hObject, handles); % Update handles structure
    handles = PopulateGUIFromParameters(handles); % refresh gui so we can enable/disable control variable as necessary.

function permutationtestingcheckbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function cluster_voxelwise_p_editbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function cluster_voxelwise_p_editbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function outputfoldereditbox_Callback(hObject, eventdata, handles)

function outputfoldereditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function chooseoutputfolderbutton_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function viewresultsbutton_Callback(hObject, eventdata, handles)
    LaunchResultsDirectory(hObject,eventdata,handles);

function LaunchResultsDirectory(hObject,eventdata,handles)
    fulloutdir = fullfile(handles.parameters.analysis_out_path,handles.parameters.analysis_name);
    if ismac || isunix % use open ...
        [~] = system(['open "' fulloutdir '"']);
    elseif ispc % use winopen
        winopen(fulloutdir)
    else
        warndlg('Cannot open output directory because I cannot determine the OS you are using.')
    end
    
    % Open overview file... %dev1
    fpath = fileparts(handles.parameters.parmsfile);
    overviewhtmlfile = fullfile(fpath,'overview.html');
    if exist(overviewhtmlfile,'file')
        web(overviewhtmlfile) % Launch in MATLAB's "web browser"
    end
    
function npermutationsclustereditbox_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

function npermutationsclustereditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Select in the dropdown list the selected item.
function realcovariateslistbox_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
val = contents{get(hObject,'Value')};
dropdownoptions = get(handles.potentialcovariateslist,'string');
set(handles.potentialcovariateslist,'value',find(strcmp(val,dropdownoptions)))

% --- Executes during object creation, after setting all properties.
function realcovariateslistbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function addcovariate_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function potentialcovariateslist_Callback(hObject, eventdata, handles)

function potentialcovariateslist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function removecovariate_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function lesionvolumecorrectiondropdown_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function lesionvolumecorrectiondropdown_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function hypodirectiondropdown_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function hypodirectiondropdown_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function applycovariatestobehaviorcheckbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function applycovariatestolesioncheckbox_Callback(hObject, eventdata, handles)
    handles = UpdateCurrentAnalysis(handles,hObject);

function optionsmenu_Callback(hObject, eventdata, handles)
if handles.parameters.parallelize
    set(handles.parallelizemenu,'Checked','on')
else
    set(handles.parallelizemenu,'Checked','off')
end

function parallelizemenu_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

% --------------------------------------------------------------------
function SVscalingmenu_Callback(hObject, eventdata, handles)
    children = get(handles.SVscalingmenu,'children');
    set(children,'Checked','off')

    switch handles.parameters.svscaling
        case 100 % these are ordered "backward" seeming, so 3rd is top in list
            set(children(3),'Checked','on')
        case 99
            set(children(2),'Checked','on')
        case 95
            set(children(1),'Checked','on')
    end
% --------------------------------------------------------------------
function maxsvscaling_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

% --------------------------------------------------------------------
function sv_scaling_99th_percentile_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);

% --------------------------------------------------------------------
function sv_scaling_95th_percentile_Callback(hObject, eventdata, handles)
handles = UpdateCurrentAnalysis(handles,hObject);
