function [handles,parameters,predAndLoss] = step1_notparallel(handles,parameters,variables)
    % This is where we'll save our GBs of permutation data output...
    parameters.outfname_big = fullfile(variables.output_folder.clusterwise,['pmu_beta_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);
    predAndLoss= {}; % reserve so in case we're not using svr, we have something to return...
    
    %% Try to use cache to skip this step by relying on cached permutation data
    if can_skip_generating_beta_perms(parameters,variables)
        warndlg('Relying on cached data not fully supported/reliable yet. Aborting.')
        error('This is not supported')
        handles = UpdateProgress(handles,'Using cached beta map permutation data...',1);
        return
    else
        handles = UpdateProgress(handles,'Computing beta map permutations (not parallelized)...',1);
        svrlsm_waitbar(parameters.waitbar,0,['Computing beta permutations (' myif(parameters.method.mass_univariate,'mass univariate','multivariate') ')...']);
    end
    
    %% If we got here then we need to generate the permutation data
    fileID = fopen(parameters.outfname_big,'w');
    
    for PermIdx=1:parameters.PermNumVoxelwise
        % Random permute subjects order for this permutation.
        trial_score = variables.one_score(randperm(length(variables.one_score)));
        if parameters.method.mass_univariate % Estimate via mass-univariate
            pmu_beta_map = nan(size(variables.lesion_dat,2),1); % reserve space
            for vox = 1 : size(variables.lesion_dat,2)
                 [Q, R] = qr(trial_score, 0); % use the householder transformations to compute the qr factorization of an n by p matrix x.
                 y = double(variables.lesion_dat(:,vox));% / 10000; % why divide by 10,000?
                 pmu_beta_map(vox) = R \ (Q' * y);  % equivalent to fitlm's output: lm.Coefficients.Estimate
            end
        else % Estimate via multivariate svr
            if PermIdx == 1, variables.orig_one_score = variables.one_score; end % store this
            variables.one_score = trial_score; % this is so we can use the same ComputeMatlabSVRLSM function :)
            [m,w,~,predAndLoss] = ComputeMatlabSVRLSM(parameters,variables); % ComputeMatlabSVRLSM will utilize our optimized parameters if available...
            
            predAndLoss{PermIdx} = predAndLoss; % store ... we accumulate..
            
            if PermIdx == parameters.PermNumVoxelwise, variables.one_score = variables.orig_one_score; end % restore this once we're done all our permutations

            % Compute the beta map here (but only if we didn't already compute it already by necessity via crossvalidation)
            if ~parameters.crossval.do_crossval % conditional added to support crossvalidated betamap option in June 2019
                pmu_beta_map = variables.beta_scale * m.Alpha * m.SupportVectors;
            else % use pre-computed (and averaged) beta map(s) -- this should only be available with svr in matlab specifically (not mass univariate, and not libsvm right now)
                pmu_beta_map = w; % here contains an average of the crossvalidated fold models' beta values, so we don't have to scale or do anything here, it's already all done in the ComputeMatlabSVRLSM() function
            end
        end

        tmp_map = zeros(variables.vo.dim(1:3)); % make a zeros template....        
        tmp_map(variables.l_idx) = pmu_beta_map;
        pmu_beta_map = tmp_map(variables.m_idx).';

        % Save this permutation to our single giant file.
        fwrite(fileID, pmu_beta_map,'single'); 

        check_for_interrupt(parameters)
        if ~mod(PermIdx,20) % update every 20 indices - Display progress.
            prompt_str = get_step1_prog_string(PermIdx,parameters);
            svrlsm_waitbar(parameters.waitbar,PermIdx/parameters.PermNumVoxelwise,prompt_str);
        end
    end
    svrlsm_waitbar(parameters.waitbar,0,''); % clear
    fclose(fileID); % close the pmu data output file.

% The string we print saying the progress.
function prompt_str = get_step1_prog_string(PermIdx,parameters)
    elapsed_time = toc;
    remain_time = round(elapsed_time * (parameters.PermNumVoxelwise - PermIdx)/(PermIdx));
    remain_time_h = floor(remain_time/3600);
    remain_time_m = floor((remain_time - remain_time_h*3600)/60);
    remain_time_s = floor(remain_time - remain_time_h*3600 - remain_time_m*60);
    prompt_str = sprintf('Permutation %d/%d: Est. remaining time: %dh %dm %ds', PermIdx, parameters.PermNumVoxelwise, remain_time_h, remain_time_m,remain_time_s);