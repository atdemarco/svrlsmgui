function WriteHyperParamOptimReport(parms)
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Hyperparameter optimization information</h2>');

    if ~parms.optimization.do_optimize
        fprintf(parms.fileID,'Hyperparameter optimization was not utilized, so there is no information to display.');
    else
        
    %% What parameters were optimized and range of each
    n_parms_optimized = [parms.optimization.params_to_optimize.cost + parms.optimization.params_to_optimize.sigma + parms.optimization.params_to_optimize.epsilon + parms.optimization.params_to_optimize.standardize];
    ndone = 0;
    narr = 'Hyperparameter optimization was utilized to choose parameters values for ';
    if parms.optimization.params_to_optimize.cost
        ndone = ndone + 1;
        narr = [narr sprintf('Cost/BoxConstraint (range %.3f - %.2f)',parms.optimization.params_to_optimize.cost_range(1),parms.optimization.params_to_optimize.cost_range(2)) myif(ndone == n_parms_optimized,'.','')];
    end
    
    if parms.optimization.params_to_optimize.sigma
        ndone = ndone + 1;
        narr = [narr myif(ndone > 1,myif(ndone == n_parms_optimized,', and ',', '),'') sprintf('Sigma (range %.3f - %.2f)',parms.optimization.params_to_optimize.sigma_range(1),parms.optimization.params_to_optimize.sigma_range(2)) myif(ndone == n_parms_optimized,'.','')];
    end
    
    if parms.optimization.params_to_optimize.epsilon
        ndone = ndone + 1;
        narr = [narr myif(ndone > 1,myif(ndone == n_parms_optimized,', and ',', '),'') sprintf('Epsilon (range %.3f - %.2f)',parms.optimization.params_to_optimize.epsilon_range(1),parms.optimization.params_to_optimize.epsilon_range(2)) myif(ndone == n_parms_optimized,'.','')];
        %,sprintf('epsilon (range %.3f - %.2f), ',parms.optimization.params_to_optimize.epsilon_range(1),parms.optimization.params_to_optimize.epsilon_range(2)),''), ...
    end
    
    if parms.optimization.params_to_optimize.standardize
        ndone = ndone + 1;
        narr = [narr myif(ndone > 1,myif(ndone == n_parms_optimized,', and ',', '),'') 'standardize (range true/false)' myif(ndone == n_parms_optimized,'.','')];
    end
    
    fprintf(parms.fileID,narr);
         
%         fprintf(parms.fileID,'Hyperparameter optimization was utilized to choose %s %s %s %s', ...
%             myif(parms.optimization.params_to_optimize.cost,sprintf('cost/boxconstraint (range %.3f - %.2f), ',parms.optimization.params_to_optimize.cost_range(1),parms.optimization.params_to_optimize.cost_range(2)),''), ...
%             myif(parms.optimization.params_to_optimize.sigma,sprintf('sigma (range %.3f - %.2f), ',parms.optimization.params_to_optimize.sigma_range(1),parms.optimization.params_to_optimize.sigma_range(2)),''), ...
%             myif(parms.optimization.params_to_optimize.epsilon,sprintf('epsilon (range %.3f - %.2f), ',parms.optimization.params_to_optimize.epsilon_range(1),parms.optimization.params_to_optimize.epsilon_range(2)),''), ...
%             myif(parms.optimization.params_to_optimize.standardize,'standardize (range true/false)',''));
        
    % Optimization found the resulting optimal values for the optimized parameters:
    
    % Ultimately, the parameters that were used for the analysis were....

         % The search strategy and the number of iterations/divisions used and objective function
         fprintf(parms.fileID,' Optimization was conducting using a %s search strategy with %d %s and %s as the objective function.', ...
            parms.optimization.search_strategy, ...
            myif(strcmp(parms.optimization.search_strategy,'Grid Search'),parms.optimization.grid_divisions,parms.optimization.iterations), ...
            myif(strcmp(parms.optimization.search_strategy,'Grid Search'),'grid divisions','iterations'), ...
            parms.optimization.objective_function);

         % the objective function that was used 
%          fprintf(parms.fileID,', using %s as the objective function.',);
        
%% Cross-validation for hyperparameter optimization
if ~parms.optimization.crossval.do_crossval
    fprintf(parms.fileID,' Cross-validation was not used during hyperparameter optimization.');
else
    fprintf(parms.fileID,' A %d-fold cross-validation scheme was used during hyperparameter optimization in which the folds were%srepartitioned each iteration.', ...
        parms.optimization.crossval.nfolds, ...
        myif(parms.optimization.crossval.repartition,' ',' not '));
end

%% The resulting parameters

% parms.optimization.crossval.method = 'kfold'; % may be others in the future?




    narr = 'parameter values of ';
    ndone = 0;
    if parms.optimization.params_to_optimize.cost
        ndone = ndone + 1;
        narr = [narr sprintf('Cost/BoxConstraint = %.2f',parms.optimization.best.cost) myif(ndone == n_parms_optimized,'.','')];
    end
    
    if parms.optimization.params_to_optimize.sigma
        ndone = ndone + 1;
        narr = [narr myif(ndone > 1,myif(ndone == n_parms_optimized,', and ',', '),'') sprintf('Sigma = %.2f',parms.optimization.best.sigma) myif(ndone == n_parms_optimized,'.','')];
    end
    
    if parms.optimization.params_to_optimize.epsilon
        ndone = ndone + 1;
        narr = [narr myif(ndone > 1,myif(ndone == n_parms_optimized,', and ',', '),'') sprintf('Epsilon = %.2f',parms.optimization.best.epsilon) myif(ndone == n_parms_optimized,'.','')];
    end
    
    if parms.optimization.params_to_optimize.standardize
        ndone = ndone + 1;
        narr = [narr myif(ndone > 1,myif(ndone == n_parms_optimized,', and ',', '),'') sprintf('Standardize = %s',myif(parms.optimization.best.standardize,'true','false')) myif(ndone == n_parms_optimized,'.','')];
    end
    
    fprintf(parms.fileID,'The best point observed during hyperparameter optimization evaluated the objective function to %s which corresponded to %s','<fill in>',narr);

          

         
        %% plots of the searches???
        
%         if parms.nclusters == 0
%             fprintf(parms.fileID,'%s','Although permutation testing was conducted, there were zero voxels that passed voxelwise thresolding, and thus no clusters of any size to plot.');
%         else
%             assess_interval=100;
%             clusters_to_show = 5; %clusters_to_show = min(last_significant_cluster,5); % if there's less than 5, then 
%             %cluster_stability_im = plotClusterPermStability(clusterwisedir,assess_interval,clusters_to_show);
%             cluster_stability_im = plotClusterPermStability(parms,assess_interval,clusters_to_show);
% 
%             imwrite(cluster_stability_im,fullfile(parms.picturedir,'cluster_stabil.png'));
% 
%             imstr = ['Cluster size stability over ' num2str(parms.PermNumVoxelwise) ' permutations, assessed every ' num2str(assess_interval) ' permutations. Also plotted are up to ' num2str(clusters_to_show) ' clusters regardless of significance for comparison to critical threshold.'];
%             fprintf(parms.fileID,'%s<br>',imstr);
%             cur_alttext = imstr;
%             imtxt = ['<img src="images/cluster_stabil.png" alt="' cur_alttext '">'];
%             fprintf(parms.fileID,'%s',imtxt);
%         end
    end

    fprintf(parms.fileID,'<br><br>\n');
