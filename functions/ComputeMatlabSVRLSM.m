function [Mdl,w,variables] = ComputeMatlabSVRLSM(parameters,variables)
    variables.one_score = variables.one_score(:);
    sigma = sqrt((1/parameters.gamma)/2); % sigma derived from gamma
    Mdl = fitrsvm(variables.lesion_dat,variables.one_score,'ObservationsIn','rows', 'KernelFunction','rbf', 'KernelScale',sigma,'BoxConstraint',parameters.cost,'Standardize',false);
    w = Mdl.Alpha.'*Mdl.SupportVectors;
%    assignin('base','Mdl',Mdl)

    % as of v0.8 9/29/17 we have customized scaling available in parameters.svscaling
    variables.beta_scale = 10 / prctile(abs(w),parameters.svscaling); % parameters.svscaling is e.g, 100 or 99 or 95
    %variables.beta_scale = 10/max(abs(w));