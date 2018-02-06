function [success,handles] = RunAnalysis(hObject,eventdata,handles)
if isempty(gcbo) || isa(gcbo,'matlab.ui.container.Menu')
    handles.parameters.runfromgui = 0;
else
    handles.parameters.runfromgui = 1;
end
try

handles.parameters = ValidateSVRLSMParameters(handles.parameters);  % fill in any missing parms

%% Validate parameters - in future move to other function
handles.parameters.time = [];
handles.parameters.time.starttime = datestr(now);

handles.parameters.datetime_run = date; % when the analysis was run. 
handles.parameters.PermNumClusterwise = handles.parameters.PermNumVoxelwise; % override the user so that these two values are the same.
    
tic; % this is what we'll use to time the execution of the analysis...

handles.parameters.baseoutputdir = fullfile(handles.parameters.analysis_out_path,handles.parameters.analysis_name,handles.parameters.datetime_run);

if OutputDirectoryAlreadyExists(handles)
    prompt={'An output directory with this name already exists. To avoid overwriting that output, change the name of your analysis. Or continue to overwrite.'};
    name = 'Rename analysis';
    defaultans = {handles.parameters.analysis_name};
    answer = inputdlg(prompt,name,[1 40],defaultans);
    if isempty(answer)
        success= 1;
        return
    end
    handles.parameters.analysis_name = answer{1};
    if handles.parameters.runfromgui 
        handles = PopulateGUIFromParameters(handles); 
    end
end

handles = UpdateProgress(handles,'Beginning analysis...',1);

% Expose the Cancel Analysis button...
if handles.parameters.runfromgui
    set(handles.runanalysisbutton,'visible','off')
    set(handles.cancelanalysisbutton,'visible','on')
    set(handles.cancelanalysisbutton,'string','Cancel Analysis')
    set(gcf,'userdata',[]); % in case it was set to 'cancel' previously we don't want to re-trigger
    
    % disable analysis controls while running.
    set(get(handles.permutationtestingpanel,'children'),'enable','off')
    set(get(handles.analysispreferencespanel,'children'),'enable','off')
    set(get(handles.covariatespanel,'children'),'enable','off')
end

parameters = handles.parameters; % Make a local copy of parameters struct ... for convenience.

% Tell which SVM will be used...
if parameters.useLibSVM, svmmethod = 'libSVM';
else svmmethod = 'MATLAB''s SVM';
end

handles = UpdateProgress(handles,['Analyses will be computed with ' svmmethod '.'],1);


%% Read behavior score and/or make sure all the requested behavioral scores are present for the patients in the file.
handles = UpdateProgress(handles,'Reading behavioral scores...',1);

variables = read_behavior_score(parameters);

%% Setup our output directory.
variables.output_folder.analysis = fullfile(handles.parameters.analysis_out_path,handles.parameters.analysis_name);

folderstring = [parameters.score_name ', ' lower(parameters.tails)]; % Sim_ROI_123, One Tailed (Negative)
variables.output_folder.base = fullfile(variables.output_folder.analysis,handles.parameters.datetime_run,folderstring);
variables.output_folder.covariates = fullfile(variables.output_folder.base,'Covariates');
folderstring = sprintf('Voxelwise p%s based on %d permutations',strrep(num2str(parameters.voxelwise_p),'0.',''),parameters.PermNumVoxelwise);
variables.output_folder.voxelwise = fullfile(variables.output_folder.base,folderstring);
folderstring = sprintf('Clusterwise p%s based on %d permutations',strrep(num2str(parameters.clusterwise_p),'0.',''),parameters.PermNumClusterwise);
variables.output_folder.clusterwise = fullfile(variables.output_folder.voxelwise,folderstring);

success = CreateDirectory(variables.output_folder.base); %#ok<NASGU> % this is the date subdirectory.

handles.parameters.parmsfile = fullfile(variables.output_folder.base,'Analysis Parameters');

% Record subjects in the analysis, and those who were excluded due to missi
handles.parameters.excluded_subjects = variables.removed_subjects; % subjects that were removed from analysis because of missing 1 or more variables...
handles.parameters.nsubjects = variables.SubNum;

tosave = handles.parameters;
save(tosave.parmsfile,'tosave') % write the file

%% Read lesion images...
handles = UpdateProgress(handles,'Reading lesion images...',1);

% remove people with non-existent lesion data
has_no_lesion = cellfun(@(x) exist(fullfile(parameters.lesion_img_folder,[x '.nii']),'file'),variables.SubjectID) == 0;
n_without_lesions=sum(has_no_lesion);
if n_without_lesions > 0
    handles = UpdateProgress(handles,[num2str(n_without_lesions) ' lesion file(s) not found, so subject(s) removed from analysis.'],1);
    variables.SubNum = variables.SubNum - n_without_lesions;
    variables.scorefiledata(has_no_lesion,:) = []; % remove rows
    variables.one_score(has_no_lesion) = []; % remove rows
    variables.SubjectID(has_no_lesion) = []; % remove rows
    if isempty(variables.SubjectID), error('No subjects are included in this analysis.'); end
else
    handles = UpdateProgress(handles,'All lesion files found.',1);
end

variables = read_lesion_imgs(parameters, variables);


handles = UpdateProgress(handles,'Successfully read behavioral scores and lesion images...',1);
handles = UpdateProgress(handles,sprintf('Running analysis for ''%s''...', parameters.score_name),1);

%% Decide what we have to covary out of the behavioral data (one score) and the lesion data
switch handles.parameters.lesionvolcorrection
    case 'Regress on Behavior'
        include_lesionvol_in_behavioral_nuisance_model = 1;
        include_lesionvol_in_brain_nuisance_model = 0;
    case 'Regress on Both'
        include_lesionvol_in_behavioral_nuisance_model = 1;
        include_lesionvol_in_brain_nuisance_model = 1;
    case 'Regress on Lesion'
        include_lesionvol_in_behavioral_nuisance_model = 0;
        include_lesionvol_in_brain_nuisance_model = 1;
    case 'None'
        include_lesionvol_in_behavioral_nuisance_model = 0;
        include_lesionvol_in_brain_nuisance_model = 0;
    case 'DTLVC'
        include_lesionvol_in_behavioral_nuisance_model = 0;
        include_lesionvol_in_brain_nuisance_model = 0;
        handles = UpdateProgress(handles,sprintf('Applying DTLVC transformation to lesion data voxel values.'),1);
        variables.lesion_dat = variables.lesion_dat ./ repmat(sqrt(variables.lesion_vol),1,size(variables.lesion_dat,2));
    otherwise
        error('Unknown lesion vol correction string.')
end

%% Store values so later we can determine if one_score is correlated with lesion volume prior to any correction
handles.parameters.one_score = variables.one_score;
handles.parameters.lesion_vol = variables.lesion_vol;

%% Construct and run behavioral nuisance model
handles.parameters.behavioralmodeldata = []; % empty unless we have behavioral model to run - for summary diagnostics

if isempty(parameters.control_variable_names) % then make sure we won't try to apply any covariates...
    handles.parameters.apply_covariates_to_behavior = 0;
    handles.parameters.apply_covariates_to_lesion = 0;
end

behavioral_nuisance_model_options = [handles.parameters.apply_covariates_to_behavior include_lesionvol_in_behavioral_nuisance_model];
if any(behavioral_nuisance_model_options)
    modelspec = 'one_score ~';
    tmp = [];
    tmp.one_score = variables.one_score;

    switch num2str(behavioral_nuisance_model_options)
        case '1  0'
            handles = UpdateProgress(handles,sprintf('Behavior nuisance model will include behavioral covariate(s) but not lesion size.'),1);
            for c = 1 : numel(handles.parameters.control_variable_names)
                curcovariate = handles.parameters.control_variable_names{c};
                modelspec = [modelspec '+' curcovariate]; %#ok<AGROW>
                tmp.(curcovariate) = variables.scorefiledata.(curcovariate);
            end
        case '0  1'
            handles = UpdateProgress(handles,sprintf('Behavior nuisance model will include lesion size but not behavioral covariate(s).'),1);
            modelspec = [modelspec ' LesionVolInternal'];
            tmp.LesionVolInternal = variables.lesion_vol;
        case '1  1'
            handles = UpdateProgress(handles,sprintf('Behavior nuisance model will include both behavioral covariate(s) and lesion size.'),1);
            for c = 1 : numel(handles.parameters.control_variable_names)
                curcovariate = handles.parameters.control_variable_names{c};
                modelspec = [modelspec '+' curcovariate]; %#ok<AGROW>
                tmp.(curcovariate) = variables.scorefiledata.(curcovariate);
            end
            modelspec = [modelspec '+ LesionVolInternal'];
            tmp.LesionVolInternal = variables.lesion_vol;
    end

    % Clean up and actually run the model. Save the results back to one_score field.
    modelspec = strrep(modelspec,'~+','~'); % remove leading plus sign of there is one
    t = struct2table(tmp);
    handles = UpdateProgress(handles,sprintf(['Behavioral nuisance model spec: ' modelspec]),1);
    if strcmp(modelspec,'one_score ~')
        handles = UpdateProgress(handles,sprintf('Skipping running an empty behavioral nuisance model.'),1);
    else
        mdl = fitlm(t,modelspec,'intercept',true);
        
        variables.one_score = mdl.Residuals.Raw(:); % save raw residuals (residualized behavioral data) for analysis.
        variables.one_score = variables.one_score + repmat(mdl.Coefficients.Estimate(1),size(variables.one_score)); % 8/7/17 -- add estimated intercept back in.
        handles = UpdateProgress(handles,sprintf('Behavior nuisance model complete.'),1);

        % Now collect behavioral nuisance model data for diagnostics in summary - 9/25/17
        data = mdl.Variables;
        data.Properties.VariableNames{1} = parameters.score_name;
        tmp = table(variables.one_score,'VariableNames',{[parameters.score_name '_corrected']});
        handles.parameters.behavioralmodeldata = [data tmp]; % save so we can do model diagnostics in summary
    end
else
    handles = UpdateProgress(handles,sprintf('No behavior nuisance model will be constructed.'),1);
end

check_for_interrupt(parameters)
    
%% Construct and run voxelwise lesion data nuisance model
brain_nuisance_model_options = [handles.parameters.apply_covariates_to_lesion include_lesionvol_in_brain_nuisance_model];
if any(brain_nuisance_model_options)
    tmp = [];
    switch num2str(brain_nuisance_model_options)
        case '1  0'
            handles = UpdateProgress(handles,sprintf('Lesion data nuisance model will include behavioral covariates but not lesion size.'),1);
            for c = 1 : numel(handles.parameters.control_variable_names)
                curcovariate = handles.parameters.control_variable_names{c};
                tmp.(curcovariate) = variables.scorefiledata.(curcovariate);
            end
        case '0  1'
            handles = UpdateProgress(handles,sprintf('Lesion data nuisance model will include lesion size but not behavioral covariates.'),1);
            tmp.LesionVolInternal = variables.lesion_vol(:);
        case '1  1'
            handles = UpdateProgress(handles,sprintf('Lesion data nuisance model will include both behavioral covariate(s) and lesion size.'),1);
            for c = 1 : numel(handles.parameters.control_variable_names)
                curcovariate = handles.parameters.control_variable_names{c};
                tmp.(curcovariate) = variables.scorefiledata.(curcovariate);
            end
            tmp.LesionVolInternal = variables.lesion_vol(:);
    end
    
    % Run the model and save the results back to the lesion_dat matrix
    variables.lesion_nuisance_model = struct2array(tmp); % make into regular ol' data columns...
    handles = UpdateProgress(handles,sprintf('Beginning to run lesion nuisance model with these covariates: %s\n',strjoin(fieldnames(tmp)')),1);
    variables = continuize_lesions(variables,handles.parameters); % will automatically use the field .lesion_nuisance_model
    variables.lesion_dat = variables.lesion_dat2;
    handles = UpdateProgress(handles,sprintf('Lesion nuisance model complete.'),1);
else
    handles = UpdateProgress(handles,sprintf('No lesion data nuisance model will be employed.'),1);
end

if ~isfield(handles,'options')
    handles.options.lesionvolumecorrection = {'Regress on Behavior','Regress on Lesion','Regress on Both','DTLVC','None'};
    handles.options.hypodirection = {'One-tailed (positive)','One-tailed (negative)','Two-tailed'};
end

check_for_interrupt(parameters)

%% Compute the real beta map
handles = UpdateProgress(handles,'Computing beta map...',1);
variables.one_score = variables.one_score*100/max(abs(variables.one_score));
[beta_map, variables] = get_beta_map(parameters, variables);

check_for_interrupt(parameters)

%% Permutation test
if parameters.DoPerformPermutationTesting
    handles = UpdateProgress(handles,'Creating output directories...',1);
    success = CreateDirectory(variables.output_folder.voxelwise); %#ok<NASGU>
    success = CreateDirectory(variables.output_folder.clusterwise); %#ok<NASGU>

    [variables] = run_beta_PMU2(parameters, variables, beta_map,handles);
    
    if parameters.do_CFWER
        warning('add me')
        % Evaluate CFWER results...
    else    
        % Evaluate clustering results
        [survivingclusters,totalclusters,survivingbetavals,hypothdirection,clusterthresh] = evaluate_clustering_results(handles,variables,parameters);
        
        handles = UpdateProgress(handles,sprintf('Results of analysis (%s):',hypothdirection),1);
        handles = UpdateProgress(handles,sprintf('%d voxels survive voxelwise threshold (P < %g, %d perms).',survivingbetavals,parameters.voxelwise_p,parameters.PermNumVoxelwise),1);
        handles = UpdateProgress(handles,sprintf('%d of %d clusters survive clusterwise threshold (P < %g, k > %d voxels, %d perms).',survivingclusters,totalclusters,parameters.clusterwise_p,clusterthresh,parameters.PermNumClusterwise),1);
    end
end

handles.parameters.time.endtime = datestr(now);
handles.parameters.time.runduration = toc;
handles.parameters.analysis_is_completed = 1;
tosave=handles.parameters;

%% Finish the analysis...
try delete(parmsfile); end %#ok<TRYNC>
save(tosave.parmsfile,'tosave') % write the file again so we know the analysis completed
success = 1;

check_for_interrupt(parameters)

handles = UpdateProgress(handles,sprintf('Starting summary file...'),1);
htmlout = SummarizeAnalysis(tosave.parmsfile); % if desired...

if isempty(htmlout) % then no summary requsted
    handles = UpdateProgress(handles,sprintf('Done, no summary file requested.'),1);
else
    handles = UpdateProgress(handles,sprintf('Done writing summary file...'),1);    
end

catch ME % If the analysis encounters an error of some sort...
      success = 0; % failure.
      all_waitbars = findobj(allchild(0), 'flat', 'Tag', 'WB');
      close([all_waitbars]); % clean up.

      if strcmp(get(gcf,'userdata'),'cancel') % if the running was canceled, then let's make sure to clean up what happened...
          set(handles.cancelanalysisbutton,'string','Cancel Analysis')
          set(gcf,'userdata',[]); % in case it was set to 'cancel' previously we don't want to re-trigger
          success = 2; % cancelled...
          % don't rethrow here -- stop gracefully.
      else
        handles.error = ME;
%         rethrow(ME)
      end
end