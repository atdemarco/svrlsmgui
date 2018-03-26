function WriteOptimalParameterReport(parms)
    %if parms.DoPerformPermutationTesting % && parms.do_CFWER
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Hyperparameter Quality Report</h2>');
    
    %% Is there a record of saving the hyperparameter quality info?
    if ~isfield(parms.files_created,'hyperparameter_quality')
        fprintf(parms.fileID,'It appears that hyperparameter quality was not measured.');
        return
    end
    
    %% Load data
    optim_data = load(parms.files_created.hyperparameter_quality);

%     optim_data.hyperparameter_quality.pred_accuracy
%     optim_data.hyperparameter_quality.repro_index.data
    
    if parms.useLibSVM
        didconverge = '<unclear with libSVM>';
        allparms = optim_data.hyperparameter_quality.behavioral_predictions.Mdl.Parameters;
 %         mdlparms.KernelScale
%         mdlparms.BoxConstraint
%         mdlparms.Epsilon = 
    else
        converged = optim_data.hyperparameter_quality.behavioral_predictions.Mdl.ConvergenceInfo.Converged;
        didconverge = myif(converged,'did','did not');
        mdlparms = optim_data.hyperparameter_quality.behavioral_predictions.Mdl.ModelParameters;
    end
    

    if ~parms.useLibSVM

        fprintf(parms.fileID,'The final SVR model %s converge. The model had parameters Cost/BoxConstraint = %.2f, %s = %.2f, Epsilon = %.2f, and Standardize = %s.', ...
            didconverge, ...
            mdlparms.BoxConstraint, ...
            myif(parms.useLibSVM,'Gamma','Sigma/KernelScale'), ...
            mdlparms.KernelScale, ... % mdlparms.KernelScale is Gamma or Sigma etc...
            mdlparms.Epsilon, ...
            myif(mdlparms.StandardizeData,'true','false'));
    
        nsvs = sum(optim_data.hyperparameter_quality.behavioral_predictions.Mdl.IsSupportVector);
        nsubs = optim_data.hyperparameter_quality.behavioral_predictions.Mdl.NumObservations;
        modelbias = optim_data.hyperparameter_quality.behavioral_predictions.Mdl.Bias;
        fprintf(parms.fileID,' The final model utilized data points from %d of the %d observations (individual subjects) as support vectors. The final model had a bias of %.2f. ',nsvs,nsubs,modelbias);
    else
        fprintf(parms.fileID,'Skipping certain output because use of libSVM. ')
    end
    
    fprintf(parms.fileID,'Average prediction accuracy is %.2f (SD = %.2f) and average reproducibility index is r =  %.2f (SD = %.2f). Average mean absolute difference between real and predicted behaviors is %.2f (SD = %.2f).<br><br>', ...
        optim_data.hyperparameter_quality.pred_accuracy.mean, ...
        optim_data.hyperparameter_quality.pred_accuracy.std, ...
        optim_data.hyperparameter_quality.repro_index.mean, ...
        optim_data.hyperparameter_quality.repro_index.std, ...
        optim_data.hyperparameter_quality.mean_abs_diff.mean, ...
        optim_data.hyperparameter_quality.mean_abs_diff.std);
        
    

    %% Plot prediction accuracy data
    f=figure('visible','off');
    a = subplot(1,6,1:2,'parent',f); % left pane.

    [fdata,x]=ksdensity(optim_data.hyperparameter_quality.pred_accuracy.data,'Support',[-1.01 1.01]);
    plot(x,fdata,'k','parent',a);
    hold on;
    set(a,'xlim',[-1 1]) % since it's a corr coef...
    
    %line(cfwerdata.cfwerinfo.cfwer_single_pval_answer*ones(1,2),get(a,'ylim'),'color','r','parent',a)
    ylabel('Density')
    xlabel('Correlation Coefficient (r)')
    %legend({'Null FWE distribution',sprintf('Whole-brain threshold (FWER = %.2f)',cfwerdata.cfwerinfo.cfwer_p_value)})
    title(sprintf('Prediction accuracy'))

    % label the result.
    %    txt = sprintf('FWER(%.2f,v=%d) = %.5f',cfwerdata.cfwerinfo.cfwer_p_value,cfwerdata.cfwerinfo.v_value_in_mm3_actual_floor,cfwerdata.cfwerinfo.cfwer_single_pval_answer);
    %    xlim = get(a,'xlim'); ylim = get(a,'ylim');
    %    text(min(xlim)+range(xlim)*.10,min(ylim)+range(ylim)*.10,txt,'Parent',a)

    %% Plot reproduciblity index
    a = subplot(1,6,3:4,'parent',f); % middle pane.
    [fdata,x]=ksdensity(optim_data.hyperparameter_quality.repro_index.data,'Support',[-1.01 1.01]);
    plot(x,fdata,'k','parent',a);
    hold on;
    set(a,'xlim',[-1 1]) % since it's a corr coef...
    ylabel('Density')
    xlabel('Correlation Coefficient (r)')
    %legend({'Null FWE distribution',sprintf('Whole-brain threshold (FWER = %.2f)',cfwerdata.cfwerinfo.cfwer_p_value)})
    title(sprintf('Reproducibility index'))

    %% Plot mean absolute difference...
    a = subplot(1,6,5:6,'parent',f); % right pane.
    [fdata,x]=ksdensity(optim_data.hyperparameter_quality.mean_abs_diff.data);
    plot(x,fdata,'k','parent',a);
    hold on;
    %set(a,'xlim',[-1 1]) % since it's a corr coef...
    ylabel('Density')
    xlabel('Mean abs. diff.')
    %legend({'Null FWE distribution',sprintf('Whole-brain threshold (FWER = %.2f)',cfwerdata.cfwerinfo.cfwer_p_value)})
    title(sprintf('Mean abs. diff.(obs vs pred)'))
    
    % also plot here a line showing the mean absolute difference from the average of the ORIGINAL data?
    %origdata = optim_data.hyperparameter_quality.behavioral_predictions.Mdl.Y;
    %average_diff_from_mean_real_data = mean(abs(origdata - mean(origdata)));
    
    %% Save the resultant figure and write to summary file
    figsnap = getframe(f); % capture whole figure.
    close(f); % close the fig
    fname = 'param_quality.png';
    imwrite(figsnap.cdata,fullfile(parms.picturedir,fname));
    imstr = 'Parameter quality plots';
    imtxt = ['<img src="images/' fname '" alt="' imstr '">'];
    fprintf(parms.fileID,'%s<br><br>',imtxt);
