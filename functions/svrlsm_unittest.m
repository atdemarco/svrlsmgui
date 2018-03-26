function svrlsm_unittest
% function to make sure all the basic functions are working right after
% changes... using a known dataset with stimulated scores...

% things to check:
% 1. pos and neg tails
% 2. regular cluster correct and cfwer correct
% 3. optim and non optim
% 4. libsvm svr, matlab svr, and mass univariate
% 5. parallel and nonparallel
% lesionvol correction?
% do use cache?

parameters = GetDefaultParameters;

parameters.analysis_out_path = '/home/crl/Documents/Projects/cfwer_svrlsm/svrlsmgui/output';

%parameters.analysis_name % < dynamic.

parameters.do_make_summary = false; % we'll test this in the future...
parameters.lesion_img_folder = '/home/crl/Documents/Projects/lesion_libraries/Binary_ANTS_to_OAT_X3refeed_unsmoothed_resampled_to_2pt5mm';
parameters.score_file = '/mnt/projects-current/Lesion_Tracing/lesion_libraries/rWABscores.csv';
parameters.score_name = 'AVComp';
parameters.cfwer_v_value = 1000; %mm3
parameters.PermNumVoxelwise = 500;
parameters.PermNumClusterwise = 500;
    
%for lesionvolcorr = {'Regress on Both','None'} % parameters.lesionvolcorrection
I=0;
for parallel = [false true]
    for optim = [false true]
        for cfwer = [false true]
            for tails = [1 0]
                for method = {'muvlsm','svr_matlab','svr_libsvm'}
                    parameters.parallelize = parallel;
                    parameters.optimization.do_optimize = optim;
                    parameters.do_CFWER = cfwer;
                    tailstring = myif(tails,'High scores are good','High scores are bad');
                    parameters.tails = tailstring;% tails{1};
                    
                    switch method{1}
                        case 'svr_matlab'
                            parameters.method.mass_univariate = 0;
                            parameters.useLibSVM = 0;
                        case 'svr_libsvm'
                            parameters.method.mass_univariate = 0;
                            parameters.useLibSVM = 1;
                        case 'muvlsm'
                            parameters.method.mass_univariate = 1;
                    end
                    
                    parameters.analysis_name = ['unit_par=' num2str(parallel) '_opt=' num2str(optim) '_cfwer=' num2str(cfwer) '_tail=' num2str(tails) '_met=' method{1}]; % update dynamically...
                    
                    if parameters.method.mass_univariate 
                        if parameters.optimization.do_optimize
                            continue % skip this, since we already have or will run this "non-optimized"
                        else % don't include optim in the name since that's not really a thing with the muvlsm
                            parameters.analysis_name = ['unit_par=' num2str(parallel) '_cfwer=' num2str(cfwer) '_tail=' num2str(tails) '_met=' method{1}]; % update dynamically...
                        end
                    end
                    I=I+1;
                    record(I).testnum = I;
                    record(I).parallel = parallel;
                    record(I).optim = optim;
                    record(I).cfwer = cfwer;
                    record(I).tails = tails; % tailstring;
                    record(I).method = method{1};
                    
                    disp(['Starting test ' num2str(I) '. -- ' parameters.analysis_name])
                    tic
                    RunAnalysisNoGUI(parameters);
                    toc
                    disp(['Completed test ' num2str(I) '. -- ' parameters.analysis_name])
                end
            end
        end
    end
end
disp('Tests completed:')
record = struct2table(record)

warning('< now check the test integrity vs...? standard analyses????')