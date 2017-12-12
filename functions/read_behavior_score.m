function variables = read_behavior_score(parameters)
    variables.scorefiledata = readtable(parameters.score_file);
    variables.covariates=[]; % placeholder.
    variables.covariates_class = {};
    data = variables.scorefiledata; % shorthand...
    varsToCheck = [parameters.score_name parameters.control_variable_names];
    for v = 1 : numel(varsToCheck)
        curVarName = varsToCheck{v};
        if ~any(strcmp(curVarName,data.Properties.VariableNames))
            error(['Unable to find column named ' curVarName ' in your input CSV file, halting the analysis.'])
        end
        
        if v == 1 % then it's the outcome measure, the main value... save it.
            variables.one_score = data.(curVarName); 
            variables.one_score_class = class(variables.one_score); % for supporting categoricals...
            variables.SubjectID = data.RegistryCode; % we just need to do this once obviously
        else % collect the covariate values...
            this_covariate_data = data.(curVarName);
            this_covariate_class = class(this_covariate_data);
            variables.covariates_class{end+1} = this_covariate_class;
            switch this_covariate_class
                case 'cell'
                    empties = cellfun(@isempty, this_covariate_data); % to replace with nans.
                    unique_categories = unique(this_covariate_data);
                    unique_categories(cellfun(@isempty,unique_categories)) = [];
%                     n_unique_categories = numel(unique_categories);
%                     if n_unique_categories > 2
%                         error(['Categorical variables are only supported up to two categories, but found ' num2str(n_unique_categories) '.'])
%                     end
                     ref_name = unique_categories{1};
%                     disp(['Reference value for covariate named ' curVarName ' defaults to ' ref_name '.'])
                    this_covariate_data = double(strcmp(this_covariate_data,ref_name));
                    this_covariate_data(empties) = NaN; % so we remove the subject 
            end
                
            variables.covariates = [variables.covariates this_covariate_data]; % append a column
            
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
    
    variables.removed_subjects = variables.SubjectID(remove_idx);
    
    variables.one_score(remove_idx) = [];
    if ~isempty(variables.covariates)
        variables.covariates(remove_idx,:) = []; % remove whole rows.
    end
    variables.SubjectID(remove_idx) = [];
    variables.scorefiledata(remove_idx,:) = [];

    variables.SubNum = length(variables.SubjectID); % number of subject used in the analysis