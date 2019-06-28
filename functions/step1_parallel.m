function [handles,parameters] = step1_parallel(handles,parameters,variables)
    % This is where we'll save our GBs of permutation data output...
    parameters.outfname_big = fullfile(variables.output_folder.clusterwise,['pmu_beta_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);

    %% Try to use cache to skip this step by relying on cached permutation data
    if can_skip_generating_beta_perms(parameters,variables)
        error('caching disabled, cause it''s not completely supported')
        handles = UpdateProgress(handles,'Using cached beta map permutation data...',1); return
    else
        handles = UpdateProgress(handles,'Computing beta map permutations (parallelized)...',1);
        svrlsm_waitbar(parameters.waitbar,0,'Computing beta permutations (parallelized)...');
    end
    
    %% If we got here then we need to generate the permutation data
    
    %% Cute way to create the permuted behavioral data... each column is one permutation.
    permdata = variables.one_score(cell2mat(cellfun(@(x) randperm(x)',repmat({numel(variables.one_score)},1,parameters.PermNumVoxelwise),'uni',false)));
    outpath = variables.output_folder.clusterwise;
    totalperms = parameters.PermNumVoxelwise;
    hyperparms = hyperparmstruct(parameters);
    
    %% Figure out what parameters we'll be using:
    if parameters.useLibSVM
        parameters.step1.libsvmstring = get_libsvm_spec(hyperparms.cost,hyperparms.gamma,hyperparms.epsilon); % Standardization is already applied.
        tmp.libsvmstring = parameters.step1.libsvmstring;
    else % use matlab's -- note the cell array we're creating - it's fancy since later we'll parameters.step1.matlab_svr_parms{:}
        parameters.step1.matlab_svr_parms = [{'BoxConstraint', hyperparms.cost, ...
            'KernelScale', hyperparms.sigma, ...
            'Standardize', hyperparms.standardize, ...
            'Epsilon', hyperparms.epsilon} ...
            myif(parameters.crossval.do_crossval,{'KFold',parameters.crossval.nfolds},[])]; % this is for the crossvalidation option...
        tmp.matlab_svr_parms = parameters.step1.matlab_svr_parms;
        tmp.do_crossval = parameters.crossval.do_crossval; % also save this here too.
    end

    %% parfeval code
    batch_job_size = 100; % this is going to be optimal for different systems/#cores/jobs - set this through gui?
    nperms = parameters.PermNumVoxelwise;
    njobs = ceil(nperms/batch_job_size); % gotta round up to capture all indices

    %% to reduce overhead transfering data to workers...
    tmp.dims = variables.vo.dim(1:3);
    tmp.lesiondata = sparse(variables.lesion_dat); % full() it on the other end - does that save time with transfer to worker overhead?!
    tmp.outpath = variables.output_folder.clusterwise;
    
    if ~parameters.method.mass_univariate % then we need to save beta_scale as well...
        tmp.beta_scale = variables.beta_scale;
    end
    
    tmp.m_idx = variables.m_idx;
    tmp.l_idx = variables.l_idx;
    tmp.totalperms = totalperms;
    tmp.permdata = permdata;
    tmp.useLibSVM = parameters.useLibSVM;
    tmp.use_mass_univariate = parameters.method.mass_univariate;

    %% Schedule the jobs...
    p = gcp(); % get current parallel pool
    for j = 1 : njobs
        this_job_start_index = ((j-1)*batch_job_size) + 1;
        this_job_end_index = min(this_job_start_index + batch_job_size-1,nperms); % need min so we don't go past valid indices
        this_job_perm_indices = this_job_start_index:this_job_end_index;
        tmp.this_job_perm_indices = this_job_perm_indices; % update for each set of jobs...
        f(j) = parfeval(p,@parallel_step1_batch_fcn_lessoverhead,0,tmp);
    end
    
    %% Monitor job progress and allow user to bail, hopefully...
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
    if ~tmp.useLibSVM || tmp.use_mass_univariate, tmp.lesiondata = full(tmp.lesiondata); end % we transfer it as sparse... so we need to full() it for all non-libSVM methods
    
    for PermIdx = tmp.this_job_perm_indices % each loop iteration will compute one whole-brain permutation result (regardless of LSM method)
        trial_score = tmp.permdata(:,PermIdx); % extract the row of permuted data.
        if tmp.use_mass_univariate % solve whole-brain permutation PermIdx on a voxel-by-voxel basis.
            pmu_beta_map = nan(size(tmp.lesiondata,2),1); % reserve space -- are these dims right?
            for vox = 1 : size(tmp.lesiondata,2)
                 [Q, R] = qr(trial_score, 0); % use the householder transformations to compute the qr factorization of an n by p matrix x.
                 y = double(tmp.lesiondata(:,vox));% / 10000; % why divide by 10,000? %betas(vox) = R \ (Q' * y);  % equivalent to fitlm's output: lm.Coefficients.Estimate
                 pmu_beta_map(vox) = R \ (Q' * y);  % equivalent to fitlm's output: lm.Coefficients.Estimate
            end
        else % use an SVR method
            if tmp.useLibSVM
                m = svmtrain(trial_score,tmp.lesiondata,tmp.libsvmstring); %#ok<SVMTRAIN> % alpha = m.sv_coef'; % SVs = m.SVs;
            else % use matlab's...
                m = fitrsvm(tmp.lesiondata,trial_score,'KernelFunction','rbf', tmp.matlab_svr_parms{:});
            end
            
            if ~tmp.do_crossval % then compute the beta map as usual...
                pmu_beta_map = tmp.beta_scale * m.(myif(tmp.useLibSVM,'sv_coef','Alpha'))' * m.(myif(tmp.useLibSVM,'SVs','SupportVectors'));
            else % we don't need to do any computations since w should already contain a scaled/averaged model
                ws = []; % we'll accumulate in here
                for mm = 1 : numel(m.Trained)
                    curMdl = m.Trained{mm};
                    w = curMdl.Alpha.'*curMdl.SupportVectors;
                    beta_scale = 10/max(abs(w));
                    w = w.'*beta_scale;
                    ws(1:numel(w),mm) = w; % accumulate here...
                end
                w = mean(ws,2);
                pmu_beta_map = w; % here contains an average of the crossvalidated fold models' beta values
            end
        end
        
        tmp_map = zeros(tmp.dims); % make a zeros template....        
        tmp_map(tmp.l_idx) = pmu_beta_map; % return the lesion data to their lidx indices...
        pmu_beta_map = tmp_map(tmp.m_idx).'; % extract only the midx indices, since these are the only voxels that will be output in the results -- midx contains only voxels that exceed the lesion threshold

        % Save this permutation (PermIdx)....
        fileID = fopen(fullfile(tmp.outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(tmp.totalperms) '.bin']),'w');
        fwrite(fileID, pmu_beta_map,'single');
        fclose(fileID);
    end
    
function [Bouts,lesiondataout] = parallel_step1_batch_muvlsm(lesiondata,modelspec)
    for vox = 1 : size(lesiondata,2) % run of the mill loop now...
        [B,~,resids] = regress(lesiondata(:,vox),modelspec);
        lesiondataout(:,vox) = resids + repmat(B(1),size(resids));
        Bouts(:,vox) = B;
    end
