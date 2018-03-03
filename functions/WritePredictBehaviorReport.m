function [] = WritePredictBehaviorReport(parms)
    warning('add me')
    return

%     if mean(abs(data.(fieldname_randdelta))) < .01
%         error('Model sensitivity to actual data is very low for these hyperparameters, so predictions are invalid.')
%     end

    %% print models:
    
    %if parameters.lesionvolcorrection ~= none

    % lesion vol only
    lm_lesonly = fitlm(data,[behavior_name '~1+lesion_vol']);
    lm_lesonly.disp % show results

    % pred behavior only
    lm_predonly = fitlm(data,[behavior_name '~1+' fieldname_pred]);
    lm_predonly.disp % show results

    % both
    lm_both = fitlm(data,[behavior_name '~1+lesion_vol+' fieldname_pred]);
    lm_both.disp % show results

    % random predictor
    lm_rand = fitlm(data,[behavior_name '~1+' fieldname_predrand]);
    lm_rand.disp

    % rand both
    lm_randboth = fitlm(data,[behavior_name '~1+lesion_vol+' fieldname_predrand]);
    lm_randboth.disp % show results

    % stepwise regression and look at the .in field of the history output.
    [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit([data.lesion_vol data.(fieldname_pred)],data.(behavior_name),'display','off');

    if history.in(end) 
        disp('A stepwise model (penter = .05, premove = .10) would include the predicted behavior.')
    else
        disp('A stepwise model (penter = .05, premove = .10) would exclude the predicted behavior.')
    end
% 
%     % for RANDOM predicted behavior -- stepwise regression and look at the .in field of the history output.
%     [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit([data.lesion_vol data.(fieldname_predrand)],data.(behavior_name),'display','off');
% 
%     if history.in(end) 
%         disp('A stepwise model (penter = .05, premove = .10) would include the random predicted behavior.')
%     else
%         disp('A stepwise model (penter = .05, premove = .10) would exclude the random predicted behavior.')
%     end

    %% make figure...
    close all;

    f = figure;
    subplot(2,1,1)
    hold on;
    plot(data.(behavior_name),'ko-')
    hold on;
    plot(data.(fieldname_pred),'ro-')
    xlabel('Subject')
    ylabel('Score')
    legend({'Actual','Predicted'})

    subplot(2,1,2)
    for c = 1 : size(data,1)
        h1 = plot(1:size(data,1),data.(['holdOut_' num2str(c)]),['b-']); % last will overwrite handle
        hold on;
    end
    h2 = plot(1:size(data,1),data.(behavior_name),'ko-');
    h3 = plot(1:size(data,1),data.(fieldname_pred),'ro-');
    xlabel('Subject')
    ylabel('Behavior')
    title('Behaviors predicted in hold-out models')
    legend([h1 h2 h3],{'Holdouts','Actual','Predicted'})

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