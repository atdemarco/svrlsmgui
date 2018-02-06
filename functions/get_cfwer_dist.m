function [pval_null_dist,nextvaldist] = get_cfwer_dist(handles,parameters,variables)
    all_perm_pmaps = memmapfile(parameters.outfname_big_p,'Format','single');
    dataRef =  all_perm_pmaps.Data;

    pval_null_dist = nan(1,size(parameters.PermNumVoxelwise,1)); 
    nextvaldist = pval_null_dist;

    L = length(variables.m_idx); % for each voxel...
    for col = 1 : parameters.PermNumVoxelwise
        check_for_interrupt(parameters)
        cur_perm_image = dataRef(col:L:end); % index out each column using skips the length of the data...
        cur_perm_image_sorted_vals = sort(cur_perm_image);
        pval_null_dist(col) = cur_perm_image_sorted_vals(parameters.cfwer_v_value);
        nextvaldist(col) = sum(cur_perm_image <= pval_null_dist(col)) - parameters.cfwer_v_value; % difference between the desired and real value rank (resolution of the threshold)....
    end
