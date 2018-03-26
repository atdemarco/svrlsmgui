function [handles,parameters] = step1_notparallel(handles,parameters,variables)
    % This is where we'll save our GBs of permutation data output...
    parameters.outfname_big = fullfile(variables.output_folder.clusterwise,['pmu_beta_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);
    
    %% Try to use cache to skip this step by relying on cached permutation data
    if can_skip_generating_beta_perms(parameters,variables)
        warndlg('Relying on cached data not fully supported/reliable yet.')
        error('This is not supported')
        handles = UpdateProgress(handles,'Using cached beta map permutation data...',1);
        return
    else
        handles = UpdateProgress(handles,'Computing beta map permutations (not parallelized)...',1);
        svrlsm_waitbar(parameters.waitbar,0,'Computing beta permutations...');
    end
    
    %% If we got here then we need to generate the permutation data
    
    fileID = fopen(parameters.outfname_big,'w');
    
    for PermIdx=1:parameters.PermNumVoxelwise
        check_for_interrupt(parameters)

        % random permute subjects order for this permutation.
        trial_score = variables.one_score(randperm(length(variables.one_score)));

        % Which package to use to compute SVM solution
        if parameters.useLibSVM % then use libSVM ...
%             box = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
%             gamma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, sigma2gamma(parameters.optimization.best.sigma), sigma2gamma(parameters.sigma)); % now derive from sigma...
%             epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);
%             libsvmstring = get_libsvm_spec(box,gamma,epsilon);
            hyperparms = hyperparmstruct(parameters);
            libsvmstring = get_libsvm_spec(hyperparms.cost,hyperparms.gamma,hyperparms.epsilon); % Standardization is already applied.
            m = svmtrain(trial_score,sparse(variables.lesion_dat),libsvmstring); %#ok<SVMTRAIN>
        else % otherwise use MATLAB...
            if PermIdx == 1
                variables.orig_one_score = variables.one_score;  % store this.
            end
            
            variables.one_score = trial_score; % this is so we can use the same ComputeMatlabSVRLSM function :)
            [m,~,~] = ComputeMatlabSVRLSM(parameters,variables); % ComputeMatlabSVRLSM will utilize our optimized parameters if available...
            
            if PermIdx == parameters.PermNumVoxelwise % put it back after we've done all permutations...
                variables.one_score = variables.orig_one_score; % restore this.
            end
        end

        % Compute the beta map
%         if parameters.useLibSVM
%             alpha = m.sv_coef';
%             SVs = m.SVs;
%         else % MATLAB's version.
%             alpha = m.Alpha';
%             SVs = m.SupportVectors;
%         end
%         
        % Compute the beta map
        alpha = m.(myif(parameters.useLibSVM,'sv_coef','Alpha'))'; % note dynamic field reference
        SVs = m.(myif(parameters.useLibSVM,'SVs','SupportVectors')); % note dynamic field reference
        
        pmu_beta_map = variables.beta_scale * alpha * SVs;

        tmp_map = zeros(variables.vo.dim(1:3)); % make a zeros template....        
        tmp_map(variables.l_idx) = pmu_beta_map;
        pmu_beta_map = tmp_map(variables.m_idx).';

        % Save this permutation to our single giant file.
        fwrite(fileID, pmu_beta_map,'single'); 

        if ~mod(PermIdx,20) % update every 20 indices...                
            % Display progress.
%             elapsed_time = toc;
%             remain_time = round(elapsed_time * (parameters.PermNumVoxelwise - PermIdx)/(PermIdx));
%             remain_time_h = floor(remain_time/3600);
%             remain_time_m = floor((remain_time - remain_time_h*3600)/60);
%             remain_time_s = floor(remain_time - remain_time_h*3600 - remain_time_m*60);
%             prompt_str = sprintf('Permutation %d/%d: Est. remaining time: %dh %dh %ds', PermIdx, parameters.PermNumVoxelwise, remain_time_h, remain_time_m,remain_time_s);
            prompt_str = get_step1_prog_string(PermIdx,parameters);
            svrlsm_waitbar(parameters.waitbar,PermIdx/parameters.PermNumVoxelwise,prompt_str);
        end
    end
    svrlsm_waitbar(parameters.waitbar,0,''); % clear
    fclose(fileID); % close the pmu data output file.

function prompt_str = get_step1_prog_string(PermIdx,parameters)
    elapsed_time = toc;
    remain_time = round(elapsed_time * (parameters.PermNumVoxelwise - PermIdx)/(PermIdx));
    remain_time_h = floor(remain_time/3600);
    remain_time_m = floor((remain_time - remain_time_h*3600)/60);
    remain_time_s = floor(remain_time - remain_time_h*3600 - remain_time_m*60);
    prompt_str = sprintf('Permutation %d/%d: Est. remaining time: %dh %dh %ds', PermIdx, parameters.PermNumVoxelwise, remain_time_h, remain_time_m,remain_time_s);