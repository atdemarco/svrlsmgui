function variables = continuize_lesions(variables,parameters)
% 8/7/17 - added support for parallelization
% 8/7/17 - fixed bug - intercept term and adding its beta back into residuals 
% 11/16/17 - added check for user interrupt 
% 2/19/18 - added in-gui waitbar in parallelized and non-parallelized loop 
%           waitbar for the par loop involved replacing the pre-existing parfor code
%           with 'parfeval' which was introduced in MATLAB 2013b

lesiondataout = nan(size(variables.lesion_dat));
modelcols = variables.lesion_nuisance_model;
lesiondata = variables.lesion_dat;
const = ones(size(modelcols,1),1);
modelspec = [const modelcols];
Bouts = nan(size(modelspec,2),size(variables.lesion_dat,2));

if parameters.parallelize
    svrlsm_waitbar(parameters.waitbar,0,'Running lesion nuisance model (parallelized)...')
    batch_job_size = 2500; % this is going to be optimal for different systems/#cores
    nvox = size(lesiondata,2);
    njobs = ceil(nvox/batch_job_size); % gotta round up to capture all indices

    p = gcp(); % get current parallel pool
    for j = 1 : njobs
        this_job_start_index = ((j-1)*batch_job_size) + 1;
        this_job_end_index = min(this_job_start_index + batch_job_size-1,nvox); % need min so we don't go past valid indices
        job_indices = this_job_start_index:this_job_end_index;
        f(j) = parfeval(p, @parallel_lesion_batch_fcn, 2,lesiondata(:,job_indices),modelspec); % 2 outputs
    end
    
    Bouts = cell(1,njobs); %reserve space - note we want to accumulate in a row here
    lesiondataout = cell(1,njobs); %reserve space - note we want to accumulate in a row here

    for j = 1 : njobs
        check_for_interrupt(parameters) % allow user to interrupt
        [idx, value1,value2] = fetchNext(f);
        Bouts{idx} = value1; % combine these cells afterward
        lesiondataout{idx} = value2; % combine these cells afterward
        svrlsm_waitbar(parameters.waitbar,j/njobs) % update waitbar progress...
    end
    
    Bouts = cell2mat(Bouts); % combine afterward
    lesiondataout = cell2mat(lesiondataout); % combine afterward

else % not parallelized.
    
    svrlsm_waitbar(parameters.waitbar,0,'Running lesion nuisance model...')
    for vox = 1 : size(variables.lesion_dat,2) 
        if ~mod(vox,100), svrlsm_waitbar(parameters.waitbar,vox/size(variables.lesion_dat,2)); end 
        check_for_interrupt(parameters) % allow user to interrupt
        [B,~,resids] = regress(lesiondata(:,vox),modelspec); 
        lesiondataout(:,vox) = resids + repmat(B(1),size(resids));
        Bouts(:,vox) = B;
    end
end

svrlsm_waitbar(parameters.waitbar,0,'') % clear the waitbar

variables.lesion_nuisance_model_betas = Bouts;
variables.lesion_dat2 = lesiondataout;

% Parallel batch helper function for parallelizing the lesion nuisance model
function [Bouts,lesiondataout] = parallel_lesion_batch_fcn(lesiondata,modelspec)
    for vox = 1 : size(lesiondata,2) % run of the mill loop now...
        [B,~,resids] = regress(lesiondata(:,vox),modelspec);
        lesiondataout(:,vox) = resids + repmat(B(1),size(resids));
        Bouts(:,vox) = B;
    end
    
% old parallelized version without batch:
%     parfor vox = 1 : size(variables.lesion_dat,2)
%         check_for_interrupt(parameters)
%         [B,~,resids] = regress(lesiondata(:,vox),modelspec);
%         lesiondataout(:,vox) = resids + repmat(B(1),size(resids));
%         Bouts(:,vox) = B;
%     end
