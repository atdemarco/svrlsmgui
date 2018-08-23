function details = CheckIfNecessaryFilesAreInstalled(handles)
    % Try to add private functions to path....
    mypath = fileparts(which('svrlsmgui'));
    addpath(fullfile(mypath,'functions')) % when called 'private' this was unnecessary, but let's add it.
    % make sure nifti toolbox is on path.
    addpath(fullfile(mypath,'functions','nifti')) % nifti toolbox
    addpath(fullfile(mypath,'functions','libsvm-3.18')) % libsvm...
    addpath(fullfile(mypath,'functions','libsvm-3.18','matlab')) % libsvm's matlab directory...
    
    handles = UpdateProgress(handles,'Checking if necessary files are installed...',1);

    spm_found = which('spm','-all');
    if isempty(spm_found)
        handles = UpdateProgress(handles,'Warning: SPM12 is not installed and/or visible on MATLAB''s path.',1);
        details.spm = 0;
    elseif isempty(strfind(spm_found{1},'spm12'))
        handles = UpdateProgress(handles,'Warning: SPM found, but version 12 is not installed and/or at the top of MATLAB''s path.',1);
        details.spm = 0;
    else 
        handles = UpdateProgress(handles,'SPM12 is installed and visible on MATLAB''s path.',1);
        details.spm = 1;
        %topofpath = spm_found{1};
        %versions = {'spm5','spm8','spm12'};
        %found_version = find(~cellfun(@isempty,cellfun(@(x) strfind(topofpath,x), versions,'uni',false)));
        %details.spmversion = versions{found_version};
    end
    
    svmtrain_found = which('svmtrain','-all'); % nb: svmtrain is also the name of a statistics toolbox function.
    correct_svmtrains = cellfun(@(x) strfind(x,'libsvm'),svmtrain_found,'Uni',false);
    
    if ~isempty(correct_svmtrains{1})
        handles = UpdateProgress(handles,'libsvm is installed and visible on MATLAB''s path.',1);
        details.libsvm = 1;
    elseif ~all(isempty(correct_svmtrains)) % one of the function is right, but not at top of path
        handles = UpdateProgress(handles,'libsvm may be installed but it appears to be overloaded MATLAB''s path (svmtrain.m?).',1);
        details.libsvm = 0;
    else
        handles = UpdateProgress(handles,'libsvm is either not installed or not visible on MATLAB''s path.',1);
        details.libsvm = 0;
    end
    
    % Now for the MATLAB statistics svr functions
    matlab_stats_found = license('checkout','statistics_toolbox');
    if matlab_stats_found
        handles = UpdateProgress(handles,'MATLAB Statistics Toolbox license found.',1);
        details.stats_toolbox = 1;
    else
        handles = UpdateProgress(handles,'MATLAB Statistics Toolbox license not found.',1);
        details.stats_toolbox = 0;
    end
    
    % Can we parallelize?
    %if ~isempty(ver('distcomp')) && license('test','Distrib_Computing_Toolbox') && (feature('numcores') > 1) - Buxbaum group found ~isempty(ver('distcomp')) didn't work for them? 
    %if isstruct(ver('distcomp')) && license('test','Distrib_Computing_Toolbox') && (feature('numcores') > 1) % Josh found that isstruct(ver('distcomp')) fails.
    v = ver;
    if any(strcmp({v.Name},'Parallel Computing Toolbox')) && (feature('numcores') > 1) % 8/23/18 - check for toolbox with new first argument (the v.Name)
        handles = UpdateProgress(handles,'Parallelization available: Distributed Computing Toolbox installed, licensed, and > 1 core.',1);
        details.can_parallelize = true;
    else
        handles = UpdateProgress(handles,'Parallelization unavailable: Distributed Computing Toolbox not installed, or only 1 core.',1);
        details.can_parallelize = false;
    end    