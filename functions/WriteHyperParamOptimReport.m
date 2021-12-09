function WriteHyperParamOptimReport(parms)
    %% Title of section
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Hyperparameter optimization information</h2>');

    if ~parms.optimization.do_optimize
        fprintf(parms.fileID,'Hyperparameter optimization was not utilized, so there is no information to display.');
    else
        printHyperParamOptimResults(parms)
    end

    fprintf(parms.fileID,'<br><br>\n');

function printHyperParamOptimResults(parms)
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
    end
    
    if parms.optimization.params_to_optimize.standardize
        ndone = ndone + 1;
        narr = [narr myif(ndone > 1,myif(ndone == n_parms_optimized,', and ',', '),'') 'standardize (range true/false)' myif(ndone == n_parms_optimized,'.','')];
    end
    
    fprintf(parms.fileID,narr);
    
    % The search strategy and the number of iterations/divisions used and objective function
    fprintf(parms.fileID,' Optimization was conducted using a %s search strategy with %d %s and %s as the objective function.', ...
        parms.optimization.search_strategy, ...
        myif(strcmp(parms.optimization.search_strategy,'Grid Search'),parms.optimization.grid_divisions,parms.optimization.iterations), ...
        myif(strcmp(parms.optimization.search_strategy,'Grid Search'),'grid divisions','iterations'), ...
        parms.optimization.objective_function);
    
    %% Cross-validation for hyperparameter optimization
    if ~parms.optimization.crossval.do_crossval
        fprintf(parms.fileID,' Cross-validation was not utilized during hyperparameter optimization.');
    else
        fprintf(parms.fileID,' A %d-fold cross-validation scheme was used during hyperparameter optimization in which the folds were%srepartitioned at each iteration.', ...
            parms.optimization.crossval.nfolds, myif(parms.optimization.crossval.repartition,' ',' not '));
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
    
    bestval = parms.optimization.results.Objective(parms.optimization.results.Rank==1); % find the top ranked one...
    fprintf(parms.fileID,' The best point observed during hyperparameter optimization evaluated the objective function to %.2f which corresponded to %s',bestval,narr);

    showGraphics(parms)
    
function showGraphics(parms)
    results = parms.optimization.results; % here we assume grid search was done.
    
    bestrow = find(results.Rank==1);
    ff=figure('Visible','off','Position',[0 0 1200 800],'Color','White');
    markersize = 25;
    jit = @(x) x .* (1 + randn(size(x)) ./ 5);

    standardizes = {'true','false'};
    IND = 0;
    for s = 1: numel(standardizes)
        for v = 1 : 3
            IND = IND + 1;
            sp = subplot(2,3,IND);
            relstandardize = results.Standardize==standardizes{s};

            markers = markersize.*ones(size(results.Rank(relstandardize)));
            scatter3(jit(results.Epsilon(relstandardize)), jit(results.BoxConstraint(relstandardize)), jit(results.KernelScale(relstandardize)), markers, results.Rank(relstandardize));
            hold on;
            
            % plot the location of the winner with a dark thick + sign...
            % plot3(results.Epsilon(bestrow), results.BoxConstraint(bestrow), results.KernelScale(bestrow),'k+','linewidth',2)

            colormap jet
            title(['Standardize = ' standardizes{s}])
            xlabel('Epsilon'); 
            ylabel('Box Constraint'); 
            zlabel('Kernel Scale'); 
            set(gca,'xscale','log','yscale','log','zscale','log'); 
            c = colorbar; 
            c.Label.String= 'Rank';
            axis square;
            grid off;
            if v  == 1, view(0,0)
            elseif v == 2, view(90,0)
            else, view(90,90)
            end
            set(sp,'fontsize',7)
        end
    end
    
    imdata = frame2im(getframe(ff));
    close(ff)
    thisfname = 'hyperparms_rank.png';
    imwrite(imdata,fullfile(parms.picturedir,thisfname)); % save the image file

    fprintf(parms.fileID,'<br><br>');

    cur_alttext = 'Hyperparameters ranked';
    imtxt = ['<img src="images/' thisfname '" alt="' cur_alttext '" width="50%" height="50%">'];
    fprintf(parms.fileID,'%s',imtxt);
    