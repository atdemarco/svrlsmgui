function results = nbytes_required_for_svrlsm_analysis(parameters,variables)
    %% Calculate how much space we think the analysis will take...
    results.bytes_precision = 4; % permutation data is single precision - 4 bytes per voxel.
    results.n_voxels_saved_after_each_permutation = numel(variables.m_idx); % these are only those voxels that exceed our minimum lesion cutoff.
    results.nbytes_required = results.bytes_precision * results.n_voxels_saved_after_each_permutation;
    
    results.nperms = 0;
    if parameters.DoPerformPermutationTesting
        results.nperms = parameters.PermNumVoxelwise;
        results.nbytes_required = results.nbytes_required * results.nperms;
    end
    
    results.also_store_pmap = parameters.do_CFWER;
    if results.also_store_pmap % then we double this value because we also need to store a file that mirrors the big beta binary final made up of p values...
        results.nbytes_required = results.nbytes_required * 2;
    end

    %% Check if we are likely to have enough space in our output directory for all the permutations, etc...
    results.outputdir = parameters.baseoutputdir; % all output will go into directories nested within here...
    FileObj = java.io.File(results.outputdir);
    results.free_bytes = FileObj.getFreeSpace;
    
    clear FileObj; % clean up before returning
    
    results.enough_space = results.free_bytes > results.nbytes_required;