function [beta_map, variables] = get_beta_map(parameters, variables)

    variables.one_score = variables.one_score(:);
        
    if parameters.useLibSVM % using libSVM
        m = svmtrain(variables.one_score,sparse(variables.lesion_dat),['-s 3 -t 2 -c ', num2str(parameters.cost), ' -g ', num2str(parameters.gamma), ' -q']);
        w = m.sv_coef'*m.SVs;
        % as of v0.8 9/29/17 we have customized scaling available in parameters.svscaling
        variables.beta_scale = 10 / prctile(abs(w),parameters.svscaling); % parameters.svscaling is e.g, 100 or 99 or 95
    else % using MATLAB
        [m,w,variables] = ComputeMatlabSVRLSM(parameters,variables);
    end
    
    tmp = zeros(variables.vo.dim(1:3));
    beta_map = tmp; 
    tmp(variables.l_idx) = w'*variables.beta_scale;
    beta_map(variables.m_idx) = tmp(variables.m_idx);
    variables.vo.fname = fullfile(variables.output_folder.base,'Beta map (unthresholded).nii'); 
    spm_write_vol(variables.vo, beta_map);

    variables.pos_idx = find(beta_map>0);
    variables.neg_idx = find(beta_map<0);
end
