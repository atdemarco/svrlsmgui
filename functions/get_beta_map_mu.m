function [beta_map,variables] = get_beta_map_mu(parameters,variables)
% Calculate beta map using mass univariate implementation of VLSM
% code based on implementation from SM Wilson's vlsm2 (Bates et al., 2003)

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
    for vox = 1 : size(variables.lesion_dat,2)
         [Q, R] = qr(onescore, 0); % use the householder transformations to compute the qr factorization of an n by p matrix x.
         y = double(lesiondata(:,vox));% / 10000; % why divide by 10,000?
         betas(vox) = R \ (Q' * y);  % equivalent to fitlm's output: lm.Coefficients.Estimate
         % if mod.... check for interrupt
    end
end

tmp = zeros(variables.vo.dim(1:3));
beta_map = tmp; 
tmp(variables.l_idx) = betas; % put betas back in the matrix, not w...
beta_map(variables.m_idx) = tmp(variables.m_idx); % m_idx -> m_idx

variables.vo.fname = fullfile(variables.output_folder.base,'Beta map (unthresholded, mass univariate).nii'); 
svrlsmgui_write_vol(variables.vo, beta_map);
variables.files_created.unthresholded_betamap = variables.vo.fname; % save for later reference.

