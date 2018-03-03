function WriteHyperParamOptimReport(parms)
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Hyperparameter optimization information</h2>');

    if ~parms.optimization.do_optimize
        fprintf(parms.fileID,'Hyperparamter optimization was not utilized, so there is no information to display.');
    else
         % what parameters were optimized and range of each
        fprintf(parms.fileID,'Hyperparameters optimization was utilized to choose %s %s %s %s', ...
            myif(parms.optimization.params_to_optimize.cost,sprintf('cost/boxconstraint (%.2f - %.2f), ',parms.optimization.params_to_optimize.cost_range(1),parms.optimization.params_to_optimize.cost_range(2)),''), ...
            myif(parms.optimization.params_to_optimize.sigma,sprintf('sigma (%.2f - %.2f), ',parms.optimization.params_to_optimize.sigma_range(1),parms.optimization.params_to_optimize.sigma_range(2)),''), ...
            myif(parms.optimization.params_to_optimize.sigma,sprintf('epsilon (%.2f - %.2f), ',parms.optimization.params_to_optimize.epsilon_range(1),parms.optimization.params_to_optimize.epsilon_range(2)),''), ...
            myif(parms.optimization.params_to_optimize.standardize,'standardize (true/false)',''));

%          parms.optimization.params_to_optimize.sigma = true; % optimize by default  --- in the context of optimization this parameter is referred to as SIGMA not GAMMA and is in units of SIGMA.
%          parms.optimization.params_to_optimize.sigma_range = [.1 100]; % < update this...
%          parms.optimization.params_to_optimize.epsilon = false; % don't optimize by default
%          parms.optimization.params_to_optimize.epsilon_range = [.1 100]; % < update this...
%          parms.optimization.params_to_optimize.standardize = false; % don't optimize by default
%          parms.optimization.params_to_optimize.standardize_range = [true false];

            % the search strategy and the number of iterations/divisions used
         fprintf(parms.fileID,'Optimization was conducting using a %s search strategy using %d %s', ...
            parms.optimization.search_strategy, ...
            myif(strcmp(parms.optimization.search_strategy,'Grid Search'),parms.optimization.iterations,parms.optimization.grid_divisions), ...
            myif(strcmp(parms.optimization.search_strategy,'Grid Search'),'iterations','grid divisions'));

         % the objective function that was used 
         fprintf(parms.fileID,'The objective function used was %s',parms.optimization.objective_function);
        
         %% crossvalidations...?
%          % Cross-validation for hyperparameter optimization
%          parms.optimization.crossval.do_crossval = true;
%          parms.optimization.crossval.method = 'kfold'; % may be others in the future?
%          parms.optimization.crossval.nfolds = 5;
%          parms.optimization.crossval.nfolds_default = parms.optimization.crossval.nfolds;
%          parms.optimization.crossval.repartition = true; % repartition at each iteration
          

         %% the resulting parameters
         
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
