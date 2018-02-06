function [handles,parameters] = step1_notparallel(handles,parameters,variables)
    handles = UpdateProgress(handles,'Computing beta map permutations (not parallelized)...',1);
    % This is where we'll save our GBs of permutation data output...
    parameters.outfname_big = fullfile(variables.output_folder.clusterwise,['pmu_beta_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);
    fileID = fopen(parameters.outfname_big,'w');

    h = waitbar(0,'Computing beta permutations...','Tag','WB');
    for PermIdx=1:parameters.PermNumVoxelwise
        check_for_interrupt(parameters)

        % random permute subjects order
        loc = randperm(length(variables.one_score));
        trial_score = variables.one_score(loc);

        % Which package to use to compute SVM solution
        if parameters.useLibSVM
            m = svmtrain(trial_score,sparse(variables.lesion_dat),['-s 3 -t 2 -c ', num2str(parameters.cost), ' -g ', num2str(parameters.gamma), ' -q']); %#ok<SVMTRAIN>
        else
            if PermIdx == 1
                variables.orig_one_score = variables.one_score;
            end

            variables.one_score = trial_score;

            [m,~,~] = ComputeMatlabSVRLSM(parameters,variables);

            if PermIdx == parameters.PermNumVoxelwise % put it back after we've done all permutations...
                variables.one_score = variables.orig_one_score;
            end
        end

        % compute the beta map
        if parameters.useLibSVM
            alpha = m.sv_coef';
            SVs = m.SVs;
        else % MATLAB's version.
            alpha = m.Alpha';
            SVs = m.SupportVectors;
        end
        
        pmu_beta_map = variables.beta_scale * alpha*SVs;

        tmp_map = zeros(variables.vo.dim(1:3)); % make a zeros template....        
        tmp_map(variables.l_idx) = pmu_beta_map;
        pmu_beta_map = tmp_map(variables.m_idx).';

        % Save this permutation....
        fwrite(fileID, pmu_beta_map,'single');

        % Display progress.
        elapsed_time = toc;
        remain_time = round(elapsed_time * (parameters.PermNumVoxelwise - PermIdx)/(PermIdx));
        remain_time_h = floor(remain_time/3600);
        remain_time_m = floor((remain_time - remain_time_h*3600)/60);
        remain_time_s = floor(remain_time - remain_time_h*3600 - remain_time_m*60);
        prompt_str = sprintf(['Permutation ', num2str(PermIdx), '/', num2str(parameters.PermNumVoxelwise), ': Est. remaining time: ', num2str(remain_time_h), ' h ', num2str(remain_time_m), ' m ' num2str(remain_time_s), 's\n']);
        waitbar(PermIdx/parameters.PermNumVoxelwise,h,prompt_str) % show progress.

    end

    close(h) % close the waitbar...
    fclose(fileID); % close the pmu data output file.
