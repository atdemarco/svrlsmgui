function variables = read_behavior_score(parameters)
    opts = detectImportOptions(parameters.score_file); % addded to support e.g., MAC CSV files 1/31/18
    variables.scorefiledata = readtable(parameters.score_file,opts); 
    variables.covariates=[]; % placeholder.
    data = variables.scorefiledata; % shorthand...
    varsToCheck = [parameters.score_name parameters.control_variable_names];
    for v = 1 : numel(varsToCheck)
        curVarName = varsToCheck{v};
        if ~any(strcmp(curVarName,data.Properties.VariableNames))
            error(['Unable to find column named ' curVarName ' in your input CSV file, halting the analysis.'])
        end
        if v == 1 % then it's the outcome measure, the main value... save it.
            variables.one_score = data.(curVarName); 
            variables.SubjectID = data.RegistryCode; % we just need to do this once obviously
        else % collect the covariate values...
            variables.covariates = [variables.covariates data.(curVarName)]; % append a column
        end
    end

    parameters.control_variable_name=[]; % remove
    
    % remove subjects without behavior score of the main variable or any covariate
    remove_idx = find(any(isnan([variables.one_score variables.covariates]),2)); % any nan the rows and convert to indices
    
    if numel(remove_idx) == 1
        warning('Subject %s was not included in the analysis because lack of one or more behavior score.\n', variables.SubjectID{remove_idx});
    elseif numel(remove_idx) > 1
        warning('The following %d subjects have been removed from the analysis because they are missing one or more behavioral score.\n', numel(remove_idx));
        fprintf('%s\n', variables.SubjectID{remove_idx});
    end
    
    variables.excluded.no_behavior = variables.SubjectID(remove_idx);
    
    variables.one_score(remove_idx) = [];
    if ~isempty(variables.covariates)
        variables.covariates(remove_idx,:) = []; % remove whole rows.
    end
    variables.SubjectID(remove_idx) = [];
    variables.scorefiledata(remove_idx,:) = [];

    variables.SubNum = length(variables.SubjectID); % number of subject used in the analysis