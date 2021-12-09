function WriteSummaryNarrative(parms)
    if parms.useLibSVM, svmtype = 'libSVM'; else svmtype = 'MATLAB'; end  % What software was used for SVR?
    if parms.parallelize,  parall = ''; else parall = 'not '; end  % Was the analysis parallelized?
    if parms.optimization.do_optimize, optim= ''; else optim = 'not '; end  % Add whether hyperparameter optimization was used
    if parms.runfromgui, runfromgui = ''; else runfromgui = ' not '; end  % Was the analysis run from the gui?

    % parms.excluded_subjects - has three fields depending on why subject was excluded
    % added 4/18/18 - forced concatenation of excluded_names on cell array rows with (:)'
    excluded_names = [parms.excluded_subjects.no_behavior(:)' parms.excluded_subjects.no_lesion(:)' parms.excluded_subjects.novoxels(:)'];
    excluded_names(cellfun(@isempty,excluded_names)) = []; % remove empties...
    nexcluded = numel(excluded_names);
    excluded_names = strjoin(excluded_names,', ');
    nsubs = parms.nsubjects; 
    
    % write the final n subjects.

    % is lesion volume and one_score correlated prior to correction?
    [rho,pval] = corrcoef(parms.one_score(:),parms.lesion_vol(:));
    rho = rho(2,1);
    pval = pval(2,1);

    %% Assemble narrative summary...
    % Hypothesis direction and behavior variable name
    descr = ['This analysis named ''' parms.analysis_name ''' tested the hypothesis (' lower(parms.tails) ') that there is a relationship between lesion status and the behavior score ''' parms.score_name '''.'];
    descr = [ descr ' ' num2str(nsubs) ' subjects were listed for inclusion, and ' num2str(nexcluded) ' were excluded due to missing behavioral data, lesion data, or for having no voxels inside the minimum lesion cutoff mask (' excluded_names ').'];

    % add note about data being binarized, resampled, etc, as necessary.
    if parms.imagedata.do_binarize, descr = [descr ' Each lesion volume was binarized prior to the analysis.']; end
    if parms.imagedata.do_resample, descr = [descr ' Each lesion volume was resampled to ' num2str(parms.imagedata.resample_to) ' mm³ voxels prior to analysis.']; end
    if parms.use_analysis_mask, descr = [descr ' This analysis was carried out within a specified inclusion mask.']; end % Dec 2021
    
    % define some convenient anonymous functions
    nice_p = @(pval) myif(pval < .001,'p < .001',strrep(sprintf('p = %.2f', pval),' 0.',' .')); % enforce no leading unnecessary 0 
    nice_r = @(rho) strrep(sprintf('r = %.2f', rho),' 0.',' .'); % enforce no leading unnecessary 0 
    
    % Was it corrected for lesion volume?
    if strcmp(parms.lesionvolcorrection,'None'), descr = [descr ' Data was not corrected for lesion volume.'];
    else, descr = [descr ' Data was corrected for lesion volume via ' parms.lesionvolcorrection '.'];
    end

    % What were the covariates - were they used?
    n_behav_covariates = numel(parms.control_variable_names);
    if n_behav_covariates > 0
        if n_behav_covariates > 1
            waswere = 'covariates were'; itthey = 'they';
        else
            waswere = 'covariate was'; itthey = 'it';
        end

        % What text do we write for the covariate inclusion info...
        if ~parms.apply_covariates_to_behavior && ~parms.apply_covariates_to_lesion
            descr = [descr ' Although ' num2str(n_behav_covariates) ' behavioral ' waswere ' specified (' strjoin(parms.control_variable_names,', ') '), ' itthey ' were not included in a nuisance model for either the behavioral score or lesion data prior to SVR.'];
        elseif parms.apply_covariates_to_behavior && ~parms.apply_covariates_to_lesion
            descr = [descr ' ' num2str(n_behav_covariates) ' behavioral ' waswere ' specified (' strjoin(parms.control_variable_names,', ') ') and covaried out of the behavioral outcome variable ''' parms.score_name ''' in a nuisance model. No nuisance model was applied to the lesion data.'];
        elseif parms.apply_covariates_to_behavior && parms.apply_covariates_to_lesion
            descr = [descr ' ' num2str(n_behav_covariates) ' behavioral ' waswere ' specified (' strjoin(parms.control_variable_names,', ') ') and covaried out of both the behavioral outcome variable ''' parms.score_name ''' and the lesion data in nuisance models.'];
        elseif ~parms.apply_covariates_to_behavior && parms.apply_covariates_to_lesion
            descr = [descr ' ' num2str(n_behav_covariates) ' behavioral ' waswere ' specified (' strjoin(parms.control_variable_names,', ') ') and covaried out of the lesion data in a nuisance model. No nuisance model was applied to the behavioral outcome variable ''' parms.score_name '''.'];
        end

        % When no lesion size correction is included but there is at least one behavioral covariate, report the post-nuisance-model r between lesion size and primary behavior
        if strcmp(parms.lesionvolcorrection,'None') && isfield(parms,'behavioralmodeldata') && ~isempty(parms.behavioralmodeldata) % earlier versions don't have this (added 0.08, 9/25/17)
            post_correction_one_score = parms.behavioralmodeldata.([parms.score_name '_corrected']);
            % is lesion volume and one_score correlated after correction?
            [rho,pval] = corrcoef(post_correction_one_score(:),parms.lesion_vol(:));
            rho = rho(2,1);
            pval = pval(2,1);

            descr = [ descr ' Following correction with the behavioral nuisance model, the behavioral score under investigation is correlated with lesion volume across the patient group at ' nice_r(rho) '.']; % ', ' nice_p(pval) '.'];
        end

    else % no covariates.
        descr = [descr ' No other covariates were specified.']; % , so no nuisance model was applied to the behavioral score or the lesion data prior to SVR
    end

    % Was permutation testing performed?
    if parms.DoPerformPermutationTesting, descr = [descr ' The resulting SVR-&beta; values were thresholded at p < ' strrep(num2str(parms.voxelwise_p),'0.','.') ' and corrected for cluster size at p < ' strrep(num2str(parms.clusterwise_p),'0.','.') ', both based on ' num2str(parms.PermNumVoxelwise) ' permutations.'];
    else descr = [ descr ' Permutation testing was not performed at either the voxel level or cluster level, meaning the results should be interpreted as tentative at very best.']; 
    end

    % What was the minimum lesion overlap set to?
    descr = [descr ' The analysis was restricted to voxels with at least ' num2str(parms.lesion_thresh) ' overlapping lesions.'];

    % Date run, what svr software was used, what gamma and cost was- -- and parallelization
    descr = [descr [' The analysis was run on ' parms.datetime_run ' using ' svmtype '''s SVR procedures (' parall 'parallelized, ' runfromgui 'run from the GUI), with parameters gamma = ' num2str(parms.gamma) ' and cost = ' num2str(parms.cost) '.' ]];

    %% Was crossvalidation conducted on the svr map output?
    if parms.crossval.do_crossval
        cxvalstring = [' SVR-&beta; maps were generated using ' num2str(parms.crossval.nfolds) '-fold cross-validation.'];
    else
        cxvalstring = ' SVR-&beta; maps were not generated using any cross-validation.';
    end
    descr = [descr cxvalstring];
        
    % Hyperparameter report.
    hyperparm_report = ['Hyperparameter optimization was ' optim ' utilized.'];
    if parms.optimization.do_optimize
        hyperparm_report = [hyperparm_report ' The method that was employed was ' CurrentOptimString(parms) '.'];
    end
    descr  = [descr hyperparm_report];
    
    % How long did it take to run?
    one_hour = 60*60; % seconds.
    if parms.time.runduration > one_hour % then report in hours and minutes...
        hours = parms.time.runduration/one_hour;
        minutes = round(rem(hours,1) * 60); % what percent of 60 minutes?
        if floor(hours) > 1, hplural = 's'; else hplural=''; end  %#ok<*SEPEX>
        if minutes > 1, mplural = 's'; else mplural=''; end 
        durationstring = [num2str(floor(hours)) ' hour' hplural ' and ' num2str(minutes) ' minute' mplural '.'];
    else % report in minutes
        minutes = parms.time.runduration/60; % divide by 60 secs.
        seconds = round(rem(minutes,1) * 60); % what percent of 60 secs?
        if floor(minutes) > 1, mplural = 's'; else mplural=''; end
        if seconds > 1, splural = 's'; else splural=''; end 
        durationstring = [num2str(floor(minutes)) ' minute' mplural ' and ' num2str(seconds) ' second' splural '.'];
    end
    descr = [descr ' The analysis completed in ' durationstring];

    % Print the narrative summary to the file.
    fprintf(parms.fileID,'<h2>Analysis summary</h2>');
    fprintf(parms.fileID,['<p>' descr '</p>']);