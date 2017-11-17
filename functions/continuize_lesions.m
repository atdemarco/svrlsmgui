function variables = continuize_lesions(variables,parameters)
% added check_for_interrupt(parameters) - 11/16/17
% support parallelization - 8/7/17 
% fixed bug - include intercept term and adding its beta back into residuals 8/7/17

lesiondataout = nan(size(variables.lesion_dat));
modelcols = variables.lesion_nuisance_model;
lesiondata = variables.lesion_dat;
const = ones(size(modelcols,1),1);
modelspec = [const modelcols];
Bouts = nan(size(modelspec,2),size(variables.lesion_dat,2));
    
if parameters.parallelize
    parfor vox = 1 : size(variables.lesion_dat,2)
        check_for_interrupt(parameters)
        [B,~,resids] = regress(lesiondata(:,vox),modelspec);
        lesiondataout(:,vox) = resids + repmat(B(1),size(resids));
        Bouts(:,vox) = B;
    end
else
    for vox = 1 : size(variables.lesion_dat,2) 
        check_for_interrupt(parameters)
        [B,~,resids] = regress(lesiondata(:,vox),modelspec); 
        lesiondataout(:,vox) = resids + repmat(B(1),size(resids));
        Bouts(:,vox) = B;
    end
end

variables.lesion_nuisance_model_betas = Bouts;
variables.lesion_dat2 = lesiondataout;    