function [cost_best, gamma_best, acc] = optimize_parameters(handles,parameters,variables, svr_cost, svr_gamma) 
% AD: parameters so we know what svr machinery to use
%[cost_best, gamma_best, acc] = optimize_parameters(variables, svr_cost, svr_gamma)
%Run 5-fold cross-validation to find best cost and gamma parameters.
%
%Inputs are the variables structure and vectors specifying the cost and
%gamma values to be considered. Zhang et al. recommended:
%   svr_cost = 1:50
%   svr_gamma = 0.5:10
%
%Returns cost and gamma parameters that produce the highest prediction
%accuracy and the full accuracy matrix. Accuracy is measured as correlation
%between predicted and actual scores.

if handles.parameters.runfromgui
    scatter(handles.paramplot,kron(svr_cost,ones(1,numel(svr_gamma))),repmat(svr_gamma,1,numel(svr_cost)),'k.');
    xlabel('Cost');
    ylabel('Gamma');
    drawnow;
end

k = 5; %Zhang et al. used 5-fold cross-validation
cross_ind = crossvalind('Kfold', variables.SubNum , k);
acc = nan(length(svr_cost), length(svr_gamma)); 
tic
for c = 1:length(svr_cost)
    for g=1:length(svr_gamma)
        %set parmaters (cost and gamma) for SVR-LSM fitting command
        cmd = ['-s 3 -t 2 -c ', num2str(svr_cost(c)), ' -g ', num2str(svr_gamma(g)), ' -q'];
        %set up vector for predicted scores
        predicted_score = nan(size(cross_ind));
        for i=1:k
            if parameters.useLibSVM % then use libSVM...
                %fit SVR-LSM model withholding k-th subset
                m = svmtrain(variables.one_score(cross_ind ~= i), sparse(variables.lesion_dat(cross_ind ~= i, :)), cmd);
                %use that model to predict scores for the k-th subset
                predicted_score(cross_ind == i) = svmpredict(variables.one_score(cross_ind == i), sparse(variables.lesion_dat(cross_ind == i, :)), m, '-q');
            else % use MATLAB's svr functionality
                disp(['Fold ' num2str(i) '....'])
                %fit SVR-LSM model withholding k-th subset
                sigma = sqrt((1/svr_gamma(g))/2); % sigma derived from gamma
                lesions_not_in_sample = variables.lesion_dat(cross_ind ~= i,:);
                behavs_not_in_sample = variables.one_score(cross_ind ~= i)';
                Mdl = fitrsvm(lesions_not_in_sample,behavs_not_in_sample,'ObservationsIn','rows', 'KernelFunction','rbf', 'KernelScale',sigma,'BoxConstraint',svr_cost(c),'Standardize',false);
                %use that model to predict scores for the k-th subset
                %behavs_in_sample = variables.one_score(cross_ind == i);
                lesions_in_sample = variables.lesion_dat(cross_ind == i,:);
                P = predict(Mdl,lesions_in_sample); % behavs_in_sample);
                predicted_score(cross_ind == i) = P; % (cross_ind == i);
            end
        end
        disp('done 5 folds:')
        %predicted_score'
        %variables.one_score
        %plot(variables.one_score, predicted_score, '.k')
        x = corrcoef(variables.one_score, predicted_score);
        acc(c, g) = x(1,2)
        
        if handles.parameters.runfromgui
%             cdata = get(handles.paramplot,'cdata');
%             set(handles.paramplot,'cdata',cdata);
            % change color of dot....
            drawnow;
        end
    end
    %time keeping
    prop_done = (c*g)/numel(acc);
    elapsed_time = toc;
    remain_time = round((elapsed_time / prop_done) - elapsed_time);
    disp(['Proportion complete: ' num2str(prop_done, 3) ', Estimated remaining time: ' num2str(remain_time/60, 3) ' minutes.'])
end
%surf(svr_gamma, svr_cost, acc)
%find max accuracy, use that to choose optimal cost and gamma values
[~, cost_best_ind] = max(max(acc,[],2)); %find row with max correlation
[~, gamma_best_ind] = max(max(acc,[],1)); %find column with max correlation
cost_best = svr_cost(cost_best_ind);
gamma_best = svr_gamma(gamma_best_ind);