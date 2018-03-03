        
%% Old paralellized code
%     lidx = variables.l_idx;
%     midx = variables.m_idx;
%     betascale = variables.beta_scale;
%     uselibsvm = parameters.useLibSVM;
%     dims = variables.vo.dim(1:3);
%     sparseLesionData = sparse(variables.lesion_dat); % matlab svr doesn't support sparse matrices
%     sigma = sqrt((1/parameters.gamma)/2); % sigma derived from gamma - for matlab's svr
%     lesion_dat = variables.lesion_dat; % for parfor.
%     

    
%     parfor PermIdx=1:parameters.PermNumVoxelwise
%         check_for_interrupt(parameters)
%         
%         trial_score = permdata(:,PermIdx); % extract the row of permuted data.
%         if uselibsvm
%             m = svmtrain(trial_score,sparseLesionData,['-s 3 -t 2 -c ', num2str(parameters.cost), ' -g ', num2str(parameters.gamma), ' -q']); %#ok<SVMTRAIN>
%             alpha = m.sv_coef';
%             SVs = m.SVs;
%         else
%             m = fitrsvm(lesion_dat,trial_score,'ObservationsIn','rows', 'KernelFunction','rbf', 'KernelScale',sigma,'BoxConstraint',parameters.cost,'Standardize',false);
%             alpha = m.Alpha';
%             SVs = m.SupportVectors;
%         end
% 
%         pmu_beta_map = betascale * alpha * SVs;
%         tmp_map = zeros(dims); % make a zeros template....        
%         tmp_map(lidx) = pmu_beta_map;
%         pmu_beta_map = tmp_map(midx).';
% 
%         % Save this permutation....
%         fileID = fopen(fullfile(outpath,['pmu_beta_map_' num2str(PermIdx) '_of_' num2str(totalperms) '.bin']),'w');
%         fwrite(fileID, pmu_beta_map,'single');
%         fclose(fileID);
%     end