function [Mdl,w,variables,predAndLoss] = ComputeMatlabSVRLSM(parameters,variables)
    % Decide whether we'll use optimized parameters or not...
    box = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.cost, parameters.optimization.best.cost, parameters.cost);
    sigma = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.sigma, parameters.optimization.best.sigma, parameters.sigma); % no longer derive from sigma
    standardize = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.standardize, parameters.optimization.best.standardize, parameters.standardize); 
    epsilon = myif(parameters.optimization.do_optimize & parameters.optimization.params_to_optimize.epsilon, parameters.optimization.best.epsilon, parameters.epsilon);
    variables.one_score = variables.one_score(:);
    
    if ~parameters.crossval.do_crossval % then business as usual.
        Mdl = fitrsvm(variables.lesion_dat,variables.one_score,'ObservationsIn','rows', 'KernelFunction','rbf', 'KernelScale',sigma,'BoxConstraint',box,'Standardize',standardize,'Epsilon',epsilon);
        w = Mdl.Alpha.'*Mdl.SupportVectors;

        % ok let's just make sure the beta scaling is calculated, and applied before returning w...
        variables.beta_scale = 10 / max(abs(w)); 
        w = w.'*variables.beta_scale;
        
        if parameters.summary.predictions % if requested...
            predAndLoss.resubPredict = Mdl.resubPredict;
            predAndLoss.resubLossMSE = Mdl.resubLoss('LossF','mse');
            predAndLoss.resubLossEps = Mdl.resubLoss('LossF','eps');
        else % not requested - don't do this, to save time
            predAndLoss.resubPredict = [];
            predAndLoss.resubLossMSE = [];
            predAndLoss.resubLossEps = [];
        end
    else % estimate a crossvalidated model and average the resulting folds... this is new as of June 2019
        Mdl = fitrsvm(variables.lesion_dat,variables.one_score,'ObservationsIn','rows', 'KernelFunction','rbf', 'KernelScale',sigma,'BoxConstraint',box,'Standardize',standardize,'Epsilon',epsilon,'KFold',parameters.crossval.nfolds);
        ws = []; % we'll accumiulate in here.
        for mm = 1 : numel(Mdl.Trained)
            curMdl = Mdl.Trained{mm};
            w = curMdl.Alpha.'*curMdl.SupportVectors;
            
            %disp('Using dynamic beta-scale in permutations for k-fold crossvalidation.')
            if ~isfield(variables,'k_fold_orig_beta_scales')
                variables.k_fold_orig_beta_scales = nan(1,numel(Mdl.Trained)); % nans will trigger into filling in
            end
            if isnan(variables.k_fold_orig_beta_scales(mm)) % then fill it in - this is an original, real data model, and crossvalid is enabled, we need to save this fold's beta_scale!
                variables.k_fold_orig_beta_scales(mm) = 10/max(abs(w)); % replace the nan with the real answer for this real, original fold (won't be recomputed for the permutations)
            end

            cur_orig_beta_scale = variables.k_fold_orig_beta_scales(mm); % for kfold xval comparison (i.e comparing this against a null hypothesis of the null dist of the same randomized *FOLD*
            beta_scale = cur_orig_beta_scale; % = 10/max(abs(w));
            
            w = w.'*beta_scale;
            ws(1:numel(w),mm) = w; % accumulate here...
        end

        w = mean(ws,2);

        variables.beta_scale = 1; % since we don't need to do additional scaling - we've already scaled... and this won't be used anyway, we will use the fold-specific beta_scales.

        if parameters.summary.predictions % if requested...
            predAndLoss.resubPredict = Mdl.kfoldPredict;
            predAndLoss.resubLossMSE = Mdl.kfoldLoss('LossF','mse');
            predAndLoss.resubLossEps = Mdl.kfoldLoss('LossF','eps');
        else % not requested - don't do this, to save time
            predAndLoss.resubPredict = [];
            predAndLoss.resubLossMSE = [];
            predAndLoss.resubLossEps = [];
        end
        
    end
    
    predAndLoss.permData = variables.one_score; % really we want to be able to know what permutatoin was used here.    