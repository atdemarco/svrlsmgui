function [handles,parameters] = step1_parallel(handles,parameters,variables)
    % This is where we'll save our GBs of permutation data output...
    parameters.outfname_big = fullfile(variables.output_folder.clusterwise,['pmu_beta_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);

    %% Try to use cache to skip this step by relying on cached permutation data
    if can_skip_generating_beta_perms(parameters,variables)
        handles = UpdateProgress(handles,'Using cached beta map permutation data...',1);
        return
    else
        handles = UpdateProgress(handles,'Computing beta map permutations (parallelized)...',1);
        svrlsm_waitbar(parameters.waitbar,0,'Computing beta permutations (parallelized)...');
    end
    
    %% If we got here then we need to generate the permutation data
    
    %% Cute way to create the permuted behavioral data... each column is one permutation.
    permdata = variables.one_score(cell2mat(cellfun(@(x) randperm(x)',repmat({numel(variables.one_score)},1,parameters.PermNumVoxelwise),'uni',false)));
    outpath = variables.output_folder.clusterwise;
    totalperms = parameters.PermNumVoxelwise;
    
    %% Figure out what parameters we'll be using:
        if parameters.useLibSVM
            box = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
            gamma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, sigma2gamma(parameters.optimization.best.sigma), sigma2gamma(parameters.sigma)); % now derive from sigma...
            epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);
            parameters.step1.libsvmstring = get_libsvm_spec(box,gamma,epsilon); % we standardize lesion cols already in RunAnalysis for libSVM
        else % use matlab's -- note the cell array we're creating - it's fancy since later we'll parameters.step1.matlab_svr_parms{:}
            parameters.step1.matlab_svr_parms = {'BoxConstraint', myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost), ...
                'KernelScale', myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, parameters.optimization.best.sigma, parameters.sigma), ...
                'Standardize', myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.standardize, parameters.optimization.best.standardize, parameters.standardize), ...
                'Epsilon', myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon)};
        end

    %% parfeval code
    batch_job_size = 100; % this is going to be optimal for different systems/#cores/jobs
    nperms = parameters.PermNumVoxelwise; % size(lesiondata,2);
    njobs = ceil(nperms/batch_job_size); % gotta round up to capture all indices

    %% to reduce overhead transfering data to workers...
    tmp.dims = variables.vo.dim(1:3);
    tmp.lesiondata = sparse(variables.lesion_dat); % full() it on the other end.
    tmp.outpath = variables.output_folder.clusterwise;
    tmp.beta_scale = variables.beta_scale;
    tmp.m_idx = variables.m_idx;
    tmp.l_idx = variables.l_idx;
    tmp.totalperms = totalperms;
    tmp.permdata = permdata;
    tmp.useLibSVM = parameters.useLibSVM;
    
    if tmp.useLibSVM
        tmp.libsvmstring = parameters.step1.libsvmstring;
    else
        tmp.matlab_svr_parms = parameters.step1.matlab_svr_parms;
    end
    
    %% Schedule the jobs...
    p = gcp(); % get current parallel pool
    for j = 1 : njobs
        this_job_start_index = ((j-1)*batch_job_size) + 1;
        this_job_end_index = min(this_job_start_index + batch_job_size-1,nperms); % need min so we don't go past valid indices
        this_job_perm_indices = this_job_start_index:this_job_end_index;
        tmp.this_job_perm_indices = this_job_perm_indices; % update for each set of jobs...
        f(j) = parfeval(p,@parallel_step1_batch_fcn_lessoverhead,0,tmp); %f(j) = parfeval(p,@parallel_step1_batch_fcn,0,variables,permdata,this_job_perm_indices,totalperms,parameters);
    end
    
    %% Monitor job progress...
    for j = 1 : njobs
        check_for_interrupt(parameters) % allow user to interrupt
        idx = fetchNext(f);
        svrlsm_waitbar(parameters.waitbar,j/njobs) % update waitbar progress...
    end
    
    %% Now assemble all those individual files from each parfored permutation into one big file that we can memmap
    handles = UpdateProgress(handles,'Consolidating beta map permutation data...',1);
    svrlsm_waitbar(parameters.waitbar,0,'Consolidating beta map permutation data...');
    fileID = fopen(parameters.outfname_big,'w');
    for PermIdx=1:parameters.PermNumVoxelwise
        if ~mod(PermIdx,100)
            svrlsm_waitbar(parameters.waitbar,PermIdx/parameters.PermNumVoxelwise); % update user on progress
            check_for_interrupt(parameters)
        end 
        curpermfilepath = fullfile(outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(totalperms) '.bin']);
        cur_perm_data = memmapfile(curpermfilepath,'Format','single');
        fwrite(fileID, cur_perm_data.Data,'single');
        clear cur_perm_data; % remove memmap from memory.
        delete(curpermfilepath); % delete it since we don't want the data hanging around...
    end
    svrlsm_waitbar(parameters.waitbar,0,''); % reset.
    fclose(fileID); % close big file


function parallel_step1_batch_fcn_lessoverhead(tmp)
    if ~tmp.useLibSVM % we transfer it as sparse... so we only need to full() it for matlab.
        tmp.lesiondata = full(tmp.lesiondata); 
    end
    
    for PermIdx = tmp.this_job_perm_indices % not parfor here... we're already parforring.
        trial_score = tmp.permdata(:,PermIdx); % extract the row of permuted data.
        if tmp.useLibSVM
            m = svmtrain(trial_score,tmp.lesiondata,tmp.libsvmstring); %#ok<SVMTRAIN>
            alpha = m.sv_coef';
            SVs = m.SVs;
        else % use matlab's...
            m = fitrsvm(tmp.lesiondata,trial_score,'ObservationsIn','rows', 'KernelFunction','rbf', tmp.matlab_svr_parms{:});
            alpha = m.Alpha';
            SVs = m.SupportVectors;
        end

        pmu_beta_map = tmp.beta_scale * alpha * SVs; % these contain results for this permutation for all l_idx indices 
        tmp_map = zeros(tmp.dims); % make a zeros template....        
        tmp_map(tmp.l_idx) = pmu_beta_map; % return the lesion data to their lidx indices...
        pmu_beta_map = tmp_map(tmp.m_idx).'; % extract only the midx indices, since these are the only voxels that will be output in the results -- midx contains only voxels that exceed the lesion threshold

        % Save this permutation....
        fileID = fopen(fullfile(tmp.outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(tmp.totalperms) '.bin']),'w');
        fwrite(fileID, pmu_beta_map,'single');
        fclose(fileID);
    end

% function parallel_step1_batch_fcn(variables,parameters)
%     if tmp.useLibSVM
%         sparseLesionData = sparse(variables.lesion_dat);
%     %else % we need sigma for matlab...
%         %sigma = sqrt((1/parameters.gamma)/2); % sigma derived from gamma - for matlab's svr
%     end
%     
%     dims = variables.vo.dim(1:3);
%     outpath = variables.output_folder.clusterwise;
%     
%     for PermIdx = tmp.this_job_perm_indices % not parfor here... we're already parforring.
%         trial_score = tmp.permdata(:,PermIdx); % extract the row of permuted data.
%         if parameters.useLibSVM
%             m = svmtrain(trial_score,sparse(variables.lesion_dat),parameters.step1.libsvmstring); %#ok<SVMTRAIN>
%             alpha = m.sv_coef';
%             SVs = m.SVs;
%         else % use matlab's...
%             m = fitrsvm(variables.lesion_dat,trial_score,'ObservationsIn','rows', 'KernelFunction','rbf', parameters.step1.matlab_svr_parms{:});
%             alpha = m.Alpha';
%             SVs = m.SupportVectors;
%         end
% 
%         pmu_beta_map = variables.beta_scale * alpha * SVs; % these contain results for this permutation for all l_idx indices 
%         tmp_map = zeros(dims); % make a zeros template....        
%         tmp_map(variables.l_idx) = pmu_beta_map; % return the lesion data to their lidx indices...
%         pmu_beta_map = tmp_map(variables.m_idx).'; % extract only the midx indices, since these are the only voxels that will be output in the results -- midx contains only voxels that exceed the lesion threshold
% 
%         % Save this permutation....
%         fileID = fopen(fullfile(outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(tmp.totalperms) '.bin']),'w');
%         fwrite(fileID, pmu_beta_map,'single');
%         fclose(fileID);
%     end
