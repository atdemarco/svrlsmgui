function [beta_map,variables] = get_beta_map_svr(parameters,variables)
% retrieve svrlsm beta map.

if parameters.useLibSVM % using libSVM
    %             box = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
    %             gamma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, sigma2gamma(parameters.optimization.best.sigma), sigma2gamma(parameters.sigma)); % now derive from sigma...
    %             epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);
    %            libsvmstring = get_libsvm_spec(box,gamma,epsilon);

    hyperparms = hyperparmstruct(parameters);
    libsvmstring = get_libsvm_spec(hyperparms.cost,hyperparms.gamma,hyperparms.epsilon); % Standardization is already applied.
    m = svmtrain(variables.one_score,sparse(variables.lesion_dat),libsvmstring); %#ok<SVMTRAIN>
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

if ~parameters.beta.do_ica_on_lesiondata % then use standard procedure...
    svrlsmgui_write_vol(variables.vo, beta_map);
else % gotta dump the lesion parcels back into whole-brain space.
    icdata = w'*variables.beta_scale; % we'll pull all the stuff out...
    svrlsmgui_write_vol_from_icdata(parameters,variables,icdata)
end

variables.files_created.unthresholded_betamap = variables.vo.fname; % save for later reference.

%variables.pos_idx = find(beta_map>0);
%variables.neg_idx = find(beta_map<0);