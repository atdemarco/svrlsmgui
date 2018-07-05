function [] = WritePredictBehaviorReport(parms)
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Behavioral predictions</h2>');

    %% Is there a record of saving this info?
    if ~isfield(parms.files_created,'hyperparameter_quality')
        fprintf(parms.fileID,'It appears that behavioral predictions were not computed.');
        return
    end

    %% Load the fitted model
    tmp = load(parms.files_created.hyperparameter_quality);
    models = tmp.hyperparameter_quality.behavioral_predictions;
    %nSVs = models.Mdl.IsSupportVector;
    
    %XVMdl = models.XVMdl;
    %try
        behavior_name = parms.behavioralmodeldata.Properties.VariableNames{1};

        predicted_field = [behavior_name '_predicted'];

        corrected_field = [behavior_name '_corrected']; % lesion vol corrected.
        was_lesion_vol_corrected = any(strcmp(corrected_field,parms.behavioralmodeldata.Properties.VariableNames));

        %% Build a new table we'll use to run out linear models on.
        t.(behavior_name) = parms.behavioralmodeldata.(behavior_name); % original behavior...
        t.lesion_vol = parms.behavioralmodeldata.LesionVolInternal; % this should always be calculated...?

        %% calculate predicted behavior and transform it back into the original units...
        %predicted_behavior= kfoldPredict(XVMdl);
        predicted_behavior = models.XVMdl_predicted(:); % turn into column vector
        t.(predicted_field) = (predicted_behavior / parms.original_behavior_transformation.maxmultiplier) - parms.original_behavior_transformation.minoffset;

        %if was_lesion_vol_corrected % we can look at it with and without lesion volume correction...
        %    t.(corrected_field) = parms.behavioralmodeldata.(corrected_field);
        %end

        t=struct2table(t);

        lm_lesonly = fitlm(t,[behavior_name '~1+lesion_vol']); % lesion vol only
        lm_predonly = fitlm(t,[behavior_name '~1+' predicted_field]); % pred behavior only
        lm_both = fitlm(t,[behavior_name '~1+lesion_vol+' predicted_field]); % both

        %assignin('base','lm_lesonly',lm_lesonly)

        %T = evalc('disp(lm_lesonly)'); %evalc (c=capture)
        %fprintf(parms.fileID,'<br>Model with lesion volume only<br>'); % break before next section

        % style="text-align:right"

        fprintf(parms.fileID,'<br>%s<br>',lm2table(lm_lesonly,'Behavioral Performance Predicted by Lesion Volume Only'));
        fprintf(parms.fileID,'<br>%s<br>',lm2table(lm_predonly,'Behavioral Performance Predicted by SVR-Predicted Behavior Only'));
        fprintf(parms.fileID,'<br>%s<br>',lm2table(lm_both,'Behavioral Performance Predicted Predicted by Lesion Volume and SVR-Predicted Behavior'));
        fprintf(parms.fileID,'<br>');

        [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit([t.lesion_vol t.(predicted_field)],t.(behavior_name),'display','off');

        if history.in(end) 
            fprintf(parms.fileID,'A stepwise model (penter = .05, premove = .10) would include the predicted behavior.');
        else
            fprintf(parms.fileID,'A stepwise model (penter = .05, premove = .10) would exclude the predicted behavior.');
        end
        fprintf(parms.fileID,'<br>');


        %% Predicted behavior plot
        f = figure('visible','off');
        a=axes(f);
        hold on;

        % plot(t.(behavior_name),'ko-','parent',a)
        % hold on;
        % plot(t.(predicted_field),'ro-','parent',a)
        % xlabel(a,'Subject')
        % ylabel(a,'Score')
        % legend(a,{'Actual','Predicted'})
        % title('Behavior predicted by SVR Model')
        
        % This was modified on 5/4/18 to be scatter-plot with 
        % predicted score on the x and actual score on the y. and a best fit line.

%         x = t.(predicted_field);
%         y = t.(behavior_name);
%         scatter(x,y,'ko');
%         hold on;
%         coef_fit = polyfit(x,y,1);
%         y_fit = polyval(coef_fit,get(a,'xlim'));
%         plot(xlim,y_fit,'r');
%         
%         xlabel('Predicted score')
%         ylabel('Actual score')
%         title('Predicted vs observed behavior and best fit line')
%         pred_im = getframe(f); % capture whole figure.
%         close(f); % close the fig
% 
%         predfname = 'predicted_behav_im.png';
%         imwrite(pred_im.cdata,fullfile(parms.picturedir,predfname));
%         imstr = 'Predicted behavior and observed behavior.';
%         imtxt = ['<img src="images/' predfname '" alt="' imstr '">'];
        
        % This replaces the scatterplot code...
        add_scatterplot(parms,t.(predicted_field),'SVR Predicted score',t.(behavior_name), ...
            'Raw Score (uncorr)','SVR Predicted vs Raw Score (uncorr)','predicted_behav_im.png')
        
        %% lesion vol vs. raw score (uncorrected) - ok
        add_scatterplot(parms,t.lesion_vol,'Lesion Volume',t.(behavior_name), ...
            'Raw Score (uncorr)','Lesion Volume vs Raw Score (uncorr)','lesion_vol_vs_raw_score.png')
        
        %% SVR prediction vs. lesion-vol regressed score - ?
        if was_lesion_vol_corrected
            add_scatterplot(parms,t.(predicted_field),'SVR Predicted score',parms.behavioralmodeldata.(corrected_field), ...
            'Lesion Vol Corrected Score','SVR Predicted vs Lesion Vol Corrected Score','predicted_behav_lesvol_im.png')
        end
        
        %% combined regression prediction (lesion vol + SVR prediction) vs. raw score
        % What does this mean?
%         add_scatterplot(t.(predicted_field),'SVR Predicted score',t.(corrected_field), ...
%             'Lesion Vol Corrected Score','SVR Predicted vs Lesion Vol Corrected Score','predicted_behav_lesvol_im.png')

% 
%     catch
%         imtxt = 'ERROR HERE - no lesion correction info stored (Cause it wasnt conducted)';
%         fprintf(parms.fileID,'%s',imtxt);
%         fprintf(parms.fileID,'<br><br>'); % break before next section
%     end

    function add_scatterplot(parms,x,xlab,y,ylab,plottitle,fname)
        f = figure('visible','off');
        a=axes(f);
        hold on;
        scatter(x,y,'ko');
        hold on;
        coef_fit = polyfit(x,y,1);
        y_fit = polyval(coef_fit,get(a,'xlim'));
        plot(xlim,y_fit,'r');
        curtitle = [plottitle ' and best fit line'];
        title(curtitle);
        xlabel(xlab);
        ylabel(ylab);
        im = getframe(f); % capture whole figure.
        close(f); % close the fig
        imwrite(im.cdata,fullfile(parms.picturedir,fname));
        imtxt = ['<img src="images/' fname '" alt="' curtitle '">'];
        fprintf(parms.fileID,'%s<br><br>',imtxt);
        
%     subplot(2,1,2)
%     for c = 1 : size(data,1)
%         h1 = plot(1:size(data,1),data.(['holdOut_' num2str(c)]),['b-']); % last will overwrite handle
%         hold on;
%     end
%    h2 = plot(1:size(t,1),t.(behavior_name),'ko-');
%    h3 = plot(1:size(t,1),t.(predicted_field),'ro-');
%     xlabel(a,'Subject')
%     ylabel(a,'Behavior')
%     title('Behaviors predicted by SVR Model')
%     legend([h2 h3],{'Actual','Predicted'})

    %% Write results
%     analysistime = strrep(strrep(char(datetime),' ','_'),':','-');
%     data.Properties.Description = ['CRL''s SVR predict (version ' script_version ') run using design file ' design_file ', lesion directory ' lesion_dir ' on ' char(datetime) ' with behavior called ' behavior_name ' using minimum overlap of ' num2str(lesion_min_mask) '. Optimization was set to ' num2str(doOptimize) '.'];
%     outfname = ['svr_predict_out_' behavior_name '_' analysistime '.xls'];
%     outfname = fullfile(output_dir,outfname);
%     if exist(output_dir,'dir')
%         writetable(data,outfname) % save output...
%     else
%         warning('The specified output directory does not exist, so no output was saved.')
%     end

% function prepped = prep4output(mdl)
%     prepped = evalc('disp(mdl)'); %evalc (c=capture)
%     prepped = strrep(prepped,newline,'<br>'); % == char(10)
%     
