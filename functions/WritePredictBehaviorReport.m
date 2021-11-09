function [] = WritePredictBehaviorReport(parms)
    %% Write out the top of the section
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Behavioral predictions</h2>');

    %% This will plot all 4 plots.
    parms = plotModelPredictionAndLosses(parms); 
    
    %% Finish section
    fprintf(parms.fileID,'<br><br>');

function parms = plotModelPredictionAndLosses(parms)
    fontopts={'FontSize',8};
    predictions_narrative_text = 'Two types of model losses are calculated for the real data model, including mean square error and epsilon-insensitive loss.';
    
    %% Now plot 2 types of losses for real model vs the null models - print out the rank if permutations have been done
    fields = {'MSE','Eps'};
    fig=figure('visible','off','color','white'); % we'll plot a 2x2 set of results panes
    for f = 1 : numel(fields)
        sps(f) = subplot(2,2,f);
        lossname = fields{f};
        hold on;
        xlabel(['Loss (' lossname ')'])
        ylabel('Density over perms')
        axis square
        grid on

        curfield = ['resubLoss' lossname];
        realdata = parms.predAndLoss.(curfield);
        
        predictions_narrative_text = [predictions_narrative_text ' For ' lossname ' loss, the real model achieves a value of ' num2str(round(realdata,2)) ', ']; % append information  1/2
        
        if parms.DoPerformPermutationTesting
            permdata = cell2mat(cellfun(@(x) x.(curfield),parms.predAndLoss_perms,'uni',false));
            [n,edges] = histcounts(permdata,'binwidth',.1,'norm','pdf');
            plot(edges(1:end-1),n,'k','linewidth',2);
            p = plot([realdata realdata],get(gca,'ylim'),'r-','linewidth',2);
            leg = legend({'Perm Mdls','Real Mdl'},fontopts{:});
            % Low numbers are good, so we want to figure out the rank of the real data relative to the permutations
            realdatarank = 1 + numel(permdata) - sum(realdata<permdata); % is this right?
            predictions_narrative_text = [predictions_narrative_text ' which ranks ' num2str(realdatarank) '/' num2str(num2str(parms.PermNumVoxelwise)) ' relative to the randomly permuted models.'];
        else % no permutations done
            text(mean(get(gca,'xlim')),mean(get(gca,'ylim')),'No perms.','horiz','center','vert','mid','FontSize',8) % put label that no permutations have been done, so not putting anyhting in plot
            predictions_narrative_text = [predictions_narrative_text 'but no permutation testing was performed so cannot calculate rank of real model loss, relative to random variability.'];
        end
        set(sps(f),'fontsize',7)
    end
    
    %% Now plot correlation of real predictions vs null model predictions - print out the rank if permutations have been done.
    sp = subplot(2,2,3);
    
    [realcorr,realcorrp] = corr(parms.one_score,parms.predAndLoss.resubPredict);
    predictions_narrative_text = [predictions_narrative_text ' Predicted scores correlate with the real scores at r = ' num2str(round(realcorr,3)) ', p = ' num2str(round(realcorrp,3))];
 
    if parms.DoPerformPermutationTesting % then we can plot and calculate everything.
        permcorr = cell2mat(cellfun(@(x) corr(parms.one_score,x.resubPredict),parms.predAndLoss_perms,'uni',false));
        
        % Here, high numbers (high correlations) are good - figure out the rank of the real data relative to the permutations
        realdatarank = 1 + numel(permcorr) - sum(realcorr>permcorr); % is this right?
        predictions_narrative_text =[predictions_narrative_text ', which ranks ' num2str(realdatarank) '/' num2str(num2str(parms.PermNumVoxelwise)) ' relative to the randomly permuted models.'];
       
        [n,edges] = histcounts(permcorr,'binwidth',.05,'norm','pdf');
        plot(edges(1:end-1),n,'k','linewidth',2);
        hold on;
        plot([realcorr realcorr],get(gca,'ylim'),'r-','linewidth',2)
        leg = legend({'Perm Mdls','Real Mdl'});
        set(leg,fontopts{:})
        set(gca,'xlim',[-1.01 1.01])
    else % no permutations, so don't plot much...
        text(mean(get(gca,'xlim')),mean(get(gca,'ylim')),'No perms.','horiz','center','vert','mid','FontSize',8) % put label that no permutations have been done, so not putting anyhting in plot
        predictions_narrative_text =[predictions_narrative_text ', but no permutation testing performed so cannot rank the correlation relative to random variability.'];
    end
    
    hold on;
    xlabel('Score correl r(predicted,real)')
    ylabel('Density over perms')
    set(sp,'fontsize',7)
    axis square
    grid on
    
    %%  Plot real data vs. predicted plot - plot this no matter what, since we don't need permutations.
    sp = subplot(2,2,4);
    realdata = parms.one_score;
    predicted = parms.predAndLoss.resubPredict;
    scatter(realdata,predicted,10,'k') % ,'Parent',ax)
    xlabel('Real')
    if ~parms.crossval.do_crossval
        crossvalstring = '(no crossval)';
    else % e.g. "kfold  (folds=5)"
        crossvalstring =  [parms.crossval.method ,' (folds=' num2str(parms.crossval.nfolds) ')'];
    end
    
    ylabel(['Predicted ' crossvalstring])
    axis square
    grid on;
    set(sp,'fontsize',7)
    
    set(fig,'Position',[0 0 800 400])
    
    %% Ok grab the resulting image, save it, and write the results to the output html overview file.
    imdata = frame2im(getframe(fig));
    close(fig)
    thisfname = 'prediction_model_losses.png';
    imwrite(imdata,fullfile(parms.picturedir,thisfname)); % save the image file
    fprintf(parms.fileID,'<br><br>');

    imstr = [predictions_narrative_text ' The figures below show the plots related to the model predictions,' ...
        ' including two model loss metrics, and the correlation between the predicted data and real data. If ' ...
        'permutation testing has been conducted, then model loss metrics and correlatedness are ranked for the ' ... 
        'real model results relative to randomly permuted models.'];
    fprintf(parms.fileID,'%s<br>',imstr);
    cur_alttext = 'Plots related to behavioral predictions.';
    imtxt = ['<img src="images/' thisfname '" alt="' cur_alttext '" width="70%" height="70%">'];
    fprintf(parms.fileID,'%s',imtxt);