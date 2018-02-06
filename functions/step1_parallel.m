function [handles,parameters] = step1_parallel(handles,parameters,variables)
    handles = UpdateProgress(handles,'Computing beta map permutations (parallelized)...',1);
    sparseLesionData = sparse(variables.lesion_dat);
    lidx = variables.l_idx;
    midx = variables.m_idx;
    betascale = variables.beta_scale;

    % create permutations beforehand.
    permdata = nan(numel(variables.one_score),parameters.PermNumVoxelwise); % each COL will be a permutation.
    npermels = size(permdata,1);
    for r = 1 : size(permdata,2) % each col...
        permdata(:,r) = variables.one_score(randperm(npermels));
    end

    outpath = variables.output_folder.clusterwise;
    totalperms = parameters.PermNumVoxelwise;
    uselibsvm = parameters.useLibSVM;
    parfor PermIdx=1:parameters.PermNumVoxelwise
        check_for_interrupt(parameters)
        trial_score = permdata(:,PermIdx); % extract the row of permuted data.
        if uselibsvm
            m = svmtrain(trial_score,sparseLesionData,['-s 3 -t 2 -c ', num2str(parameters.cost), ' -g ', num2str(parameters.gamma), ' -q']); %#ok<SVMTRAIN>
            alpha = m.sv_coef';
            SVs = m.SVs;
        else
            [m,~,~] = ComputeMatlabSVRLSM(parameters,variables);
            alpha = m.Alpha';
            SVs = m.SupportVectors;
        end

        pmu_beta_map = betascale * alpha * SVs;
        tmp_map = zeros(variables.vo.dim(1:3)); % make a zeros template....        
        tmp_map(lidx) = pmu_beta_map;
        pmu_beta_map = tmp_map(midx).';

        % Save this permutation....
        fileID = fopen(fullfile(outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(totalperms) '.bin']),'w');
        fwrite(fileID, pmu_beta_map,'single');
        fclose(fileID);
    end

    % now get all those individual files into one big file that we can memmap to.
    parameters.outfname_big = fullfile(variables.output_folder.clusterwise,['pmu_beta_maps_N_' num2str(parameters.PermNumVoxelwise) '.bin']);
    fileID = fopen(parameters.outfname_big,'w');

    for PermIdx=1:parameters.PermNumVoxelwise
        check_for_interrupt(parameters)
        curpermfilepath = fullfile(outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(totalperms) '.bin']);
        cur_perm_data = memmapfile(curpermfilepath,'Format','single');
        fwrite(fileID, cur_perm_data.Data,'single');
        clear cur_perm_data; % remove memmap from memory.
        delete(curpermfilepath); % delete it since we don't want the data hanging around...
    end
    fclose(fileID); % close big file
