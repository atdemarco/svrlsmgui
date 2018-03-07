function [beta_map, variables] = get_beta_map(parameters, variables)
    variables.one_score = variables.one_score(:);

    % Decide whether we will do mass univariate or multivariate...
    if parameters.method.mass_univariate
        betas = nan(size(variables.l_idx)); % reserve space
        onescore = variables.one_score;
        lesiondata = variables.lesion_dat;
        if parameters.parallelize
            parfor vox = 1 : size(variables.lesion_dat,2)
                 [Q, R] = qr(onescore, 0); % use the householder transformations to compute the qr factorization of an n by p matrix x.
                 y = double(lesiondata(:,vox));% / 10000; % why divide by 10,000?
                 beta = R \ (Q' * y);  % equivalent to fitlm's output: lm.Coefficients.Estimate
                 betas(vox) = beta;
            end
        else
            lesiondata = variables.lesion_dat;
            for vox = 1 : size(variables.lesion_dat,2)
                 [Q, R] = qr(variables.one_score, 0); % use the householder transformations to compute the qr factorization of an n by p matrix x.
                 y = double(lesiondata(:,vox));% / 10000; % why divide by 10,000?
                 betas(vox) = R \ (Q' * y);  % equivalent to fitlm's output: lm.Coefficients.Estimate
            end
        end
        
        tmp = zeros(variables.vo.dim(1:3));
        beta_map = tmp; 
        tmp(variables.l_idx) = betas; % put betas back in the matrix, not w...
        beta_map(variables.m_idx) = tmp(variables.m_idx); % m_idx -> m_idx

        variables.vo.fname = fullfile(variables.output_folder.base,'Beta map (unthresholded, mass univariate).nii'); 
        svrlsmgui_write_vol(variables.vo, beta_map);
        variables.files_created.unthresholded_betamap = variables.vo.fname; % save for later reference.
        
    else % do multivariate
        if parameters.useLibSVM % using libSVM
            box = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
            gamma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, sigma2gamma(parameters.optimization.best.sigma), sigma2gamma(parameters.sigma)); % now derive from sigma...
            epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);
            libsvmstring = get_libsvm_spec(box,gamma,epsilon);
            % Standardization is already applied.
            m = svmtrain(variables.one_score,sparse(variables.lesion_dat),libsvmstring);
            w = m.sv_coef'*m.SVs;
            variables.beta_scale = 10 / prctile(abs(w),parameters.svscaling);
        else % using MATLAB
            [m,w,variables] = ComputeMatlabSVRLSM(parameters,variables);
        end

        tmp = zeros(variables.vo.dim(1:3));
        beta_map = tmp;
        
        tmp(variables.l_idx) = w'*variables.beta_scale; % return all lesion data to its l_idx indices.
        beta_map(variables.m_idx) = tmp(variables.m_idx); % m_idx -> m_idx
        
        variables.vo.fname = fullfile(variables.output_folder.base,'Beta map (unthresholded).nii'); 
        
        if ~parameters.beta.do_ica_on_lesiondata
            svrlsmgui_write_vol(variables.vo, beta_map);
        else % gotta dump the lesion parcels back into whole-brain space.
            icdata = w'*variables.beta_scale; % we'll pull all the stuff out...
            svrlsmgui_write_vol_from_icdata(parameters,variables,icdata)
        end
        
        variables.files_created.unthresholded_betamap = variables.vo.fname; % save for later reference.

        %variables.pos_idx = find(beta_map>0);
        %variables.neg_idx = find(beta_map<0);
    end
    
end
