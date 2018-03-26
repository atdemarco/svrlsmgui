function variables = get_cfwer_dist(handles,parameters,variables)
    handles = UpdateProgress(handles,'Calculating CFWER thresholds...',1);
    
    % Create output dir for cfwer info
    success = CreateDirectory(variables.output_folder.cfwer); %#ok<NASGU>
    
    all_perm_pmaps = memmapfile(parameters.outfname_big_p,'Format','single');
    dataRef = all_perm_pmaps.Data;

    %% Convert the user submitted cutoff, which should be in mm3, to a voxel index (by volume)...
    % from the data spm gets from the image file, I think this is pixdim:
    cfwerinfo.pixdim = diag(variables.vo.mat(1:end-1,1:end-1)); % width, height, and depth of each voxel in this image in millimeters
    cfwerinfo.each_voxel_volume_mm3 = prod(cfwerinfo.pixdim); % each voxel's volume in cubic millimeters.
    cfwerinfo.v_value_in_mm3 = parameters.cfwer_v_value; % specified in the svrlsmgui interface...
    cfwerinfo.vox_limited_v_val = cfwerinfo.v_value_in_mm3 / cfwerinfo.each_voxel_volume_mm3; % converted here...
    cfwerinfo.vox_limited_v_val_floor = floor(cfwerinfo.vox_limited_v_val); % we floor this so you don't get free halfvoxels...
    cfwerinfo.v_value_in_mm3_actual_floor = cfwerinfo.each_voxel_volume_mm3 * cfwerinfo.vox_limited_v_val_floor; % the actual cubic millimeter of false positives we are admitting.
    cfwerinfo.cfwer_p_value = parameters.cfwer_p_value; % for ease of access in summary file.
    
    switch parameters.tailshort % tail
        case {'pos','neg'} % cfwerinfo.requested_v_index = nperms - cfwerinfo.vox_limited_v_val_floor;
            cfwerinfo.requested_v_index = cfwerinfo.vox_limited_v_val_floor;
%         case 'two'
%             error('not yet implemented')
%             cfwerinfo.requested_v_index.pos = floor(nperms - (cfwerinfo.actual_v_value_floor/2)); 
%             cfwerinfo.requested_v_index.neg = floor(cfwerinfo.actual_v_value_floor/2);
    end

    svrlsm_waitbar(parameters.waitbar,0,'Reconstructing null p-vols and calculating CFWER threshold.');
    
    %% First determine the whole-brain cutoff.
    %% This is done by reassembling each whole-volume of permutation p values. These are stored in the file with offsets 1 ... n_perms, skips of length(variables.m_idx), and number of elements of ...????
    L = length(variables.m_idx);
    nperms = parameters.PermNumVoxelwise;
    %disp(['numvoxels in each volume = ' num2str(L)])
    for col = 1 : nperms
        if ~mod(col,100) % update every 50 indices...
            svrlsm_waitbar(parameters.waitbar,col/nperms)
            check_for_interrupt(parameters)
        end
        
        cur_perm_image = dataRef(col:nperms:end);
        
     %  disp(['numvoxels in CUR volume = ' num2str(numel(cur_perm_image))])
            
        if parameters.SaveNullPMapsPreThresholding % need to update the inv option in this function V
            save_pre_thresholded_pmap_permutation(parameters,variables,cur_perm_image,col) % thresholded,
        end
        
        %% Sort the p-values in this null p-value volume, so small values are earlier in the list.
        cur_perm_image_sorted_vals = sort(cur_perm_image); %small to large, so earlier indices are smaller/more extreme p values.
        
        %% Figure out the "actual v" we will draw on for this null volume.
        switch parameters.tailshort % tail
            case {'pos','neg'}
%                 cfwerinfo.requested_v_index
%                 col
                cfwerinfo.p_of_requested_v(col) = cur_perm_image_sorted_vals(cfwerinfo.requested_v_index);
                cfwerinfo.n_found_p_of_requested_v(col) = numel(find(cur_perm_image_sorted_vals == cfwerinfo.p_of_requested_v(col))); % is there only 1 of this p value in our list (i.e. do we have adequate 'distributional resolution' to be sure of our rank here?)
                %disp('|')
%                disp(['col: ' num2str(col) ', the p of requested v (v=' num2str(cfwerinfo.requested_v_index) ') is ' num2str(cfwerinfo.p_of_requested_v(col)) ', nfound = ' num2str(cfwerinfo.n_found_p_of_requested_v(col)) '.'])
                cfwerinfo.adjust_updown(col) = 0; % neither up nor down
                if cfwerinfo.n_found_p_of_requested_v(col) > 1 % we have to find the next unique p value that's more or less extreme than the value..
                    if rand > .5 % then go to the left in the distribution
                        cfwerinfo.adjust_updown(col) = -1; % to show we went to the left.
                        more_extreme_voxels = cur_perm_image_sorted_vals < cfwerinfo.p_of_requested_v(col); % booleans of voxels with more extreme p values
                   %     disp(['num more extreme LESS THAN found:' num2str(sum(more_extreme_voxels)) '.'])
                        if ~any(more_extreme_voxels)
                            cfwerinfo.pval_null_dist.actual_v(col) = 1; % require to use the most extreme v.
                  %          disp(['since no more extreme, forced to use most extreme v = 1'])
                        else
                            cfwerinfo.pval_null_dist.actual_v(col) = find(more_extreme_voxels,1,'last'); % take the last one, and that's our v.
                 %           disp(['using the next most extreme v index = ' num2str(cfwerinfo.pval_null_dist.actual_v(col)) '.'])
                        end
                    else % go to the right in the distribution...
                        cfwerinfo.adjust_updown(col) = +1;% to show we went to the right.
                        more_extreme_voxels = cur_perm_image_sorted_vals > cfwerinfo.p_of_requested_v(col);
 %                       disp(['num more extreme GREATER THAN found:' num2str(sum(more_extreme_voxels)) '.'])
                        if ~any(more_extreme_voxels)
                             cfwerinfo.pval_null_dist.actual_v(col) = numel(cur_perm_image_sorted_vals); % require to use the most extreme v.
  %                           disp(['since no more extreme, forced to use most extreme v =' num2str(numel(cur_perm_image_sorted_vals))])
                        else
                            cfwerinfo.pval_null_dist.actual_v(col) = find(more_extreme_voxels,1,'first');
   %                         disp(['using the next most extreme v index = ' num2str(cfwerinfo.pval_null_dist.actual_v(col)) '.'])
                        end
                        
                    end
                else  % then we can get our exact requested v and p.
                    cfwerinfo.pval_null_dist.actual_v(col) = cfwerinfo.requested_v_index;
                end
                
                % draw the actual p value at that 'v' value.
                cfwerinfo.pval_null_dist.pval(col) = cur_perm_image_sorted_vals(cfwerinfo.pval_null_dist.actual_v(col));
    %            disp(['this col''s pval for the null dist =' num2str(cfwerinfo.pval_null_dist.pval(col)) '.'])
%             case 'two'
%                 error('this code has not been checked - will be updated soon')
%                 % negative side of the two-tail
%                 p_of_requested_v = cur_perm_image_sorted_vals(cfwerinfo.requested_v_index.neg);
%                 if numel(find(cur_perm_image_sorted_vals == p_of_requested_v)) > 1 % we have to find the next unique p value that's more extreme than the value..
%                     more_extreme_voxels = cur_perm_image_sorted_vals < p_of_requested_v; 
%                     cfwerinfo.pval_null_dist.actual_v.neg(col) = find(more_extreme_voxels,1,'last'); 
%                 else  % then we can get our exact requested v and p.
%                     cfwerinfo.pval_null_dist.actual_v.neg(col) = cfwerinfo.requested_v_index;
%                 end
%                 
%                 % draw the actual p value at that 'v' value.
%                 cfwerinfo.pval_null_dist.pval.neg(col) = cur_perm_image_sorted_vals(cfwerinfo.pval_null_dist.actual_v.neg(col));
% 
%                 % positive side of the two-tail
%                 p_of_requested_v = cur_perm_image_sorted_vals(cfwerinfo.requested_v_index.pos);
%                 if numel(find(cur_perm_image_sorted_vals == p_of_requested_v)) > 1 % we have to find the next unique p value that's more extreme than the value..
%                     more_extreme_voxels = cur_perm_image_sorted_vals > p_of_requested_v; % note we use greater than here.
%                     cfwerinfo.pval_null_dist.actual_v.pos(col) = find(more_extreme_voxels,1,'first'); % take the first next most extreme voxel value, and that index is our v.
%                 else  % then we can get our exact requested v and p.
%                     cfwerinfo.pval_null_dist.actual_v.pos(col) = cfwerinfo.requested_v_index;
%                 end
%                 
%                 % draw the actual p value at that 'v' value.
%                 cfwerinfo.pval_null_dist.pval.pos(col) = cur_perm_image_sorted_vals(cfwerinfo.pval_null_dist.actual_v.pos(col));
        end
    end
    
    %% Brainwide CFWER cutoffs...
    % Note: prctile() can produce values that do not strictly exist in the data. this seems not to be the end of the world...
    switch parameters.tailshort % tail
        case {'pos','neg'}
             cfwerinfo.cfwer_single_pval_answer = prctile(cfwerinfo.pval_null_dist.pval,100*parameters.cfwer_p_value); % this is the value that is used to threshold the data brain wide
             cfwerinfo.cfwer_single_pval_answer_unadjusted = prctile(cfwerinfo.p_of_requested_v,100*parameters.cfwer_p_value); % this is based on the *requested* v value no matter what. 
%         case 'two'
%             cfwerinfo.cfwer_single_pval_answer.pos = prctile(cfwerinfo.pval_null_dist.pos,100*(1-(parameters.cfwer_p_value/2))); 
%             cfwerinfo.cfwer_single_pval_answer.neg = prctile(cfwerinfo.pval_null_dist.neg,100*(parameters.cfwer_p_value/2)); 
    end
    
    %% Now let's go back through the null distributions and determine the distribution of # of suprathreshold clusters voxels of our null distributions
    svrlsm_waitbar(parameters.waitbar,0,'Counting suprathresholde CFWER voxels in permutation data.');
    for col = 1 : nperms % for each null volume...
        if ~mod(col,100) % update every X indices...
            check_for_interrupt(parameters)
            svrlsm_waitbar(parameters.waitbar,col/nperms);
        end
        
        cur_perm_image = dataRef(col:nperms:end);
        
        switch parameters.tailshort % tail
            case {'pos','neg'}
                include_mask = cur_perm_image <= cfwerinfo.cfwer_single_pval_answer; % highlight voxels the number of voxels that exceed our critical value.
                cfwerinfo.n_suprathreshold_vox_dist(col) = sum(include_mask(:)); 
                cfwerinfo.n_suprathreshold_vox_dist_unadjusted(col) = sum(cur_perm_image <= cfwerinfo.cfwer_single_pval_answer_unadjusted); % not adjusted for tie's in p-value rankings...
                
                % get the biggest null cluster under the adjusted whole-brain cutoff...
                cur_perm_image(~include_mask) = 0; % zero out voxels that didn't achieve our critical value...
                zerosmap = zeros(variables.vo.dim(1:3)); % make a zeros template....
                zerosmap(variables.m_idx) = cur_perm_image; % put it back into brain space so we can measure cluster extents...
                cfwerinfo.largest_cluster(col) = get_cur_largest_cluster_size_cfwer(parameters,variables,zerosmap);
                
%             case 'two'
%                 cfwerinfo.n_suprathreshold_vox_dist.pos(col) = sum(cur_perm_image >= cfwerinfo.cfwer_single_pval_answer.pos);
%                 cfwerinfo.n_suprathreshold_vox_dist.neg(col) = sum(cur_perm_image <= cfwerinfo.cfwer_single_pval_answer.neg);
        end
        
        if parameters.SaveNullPMapsPostThresholding
            save_post_thresholded_pmap_permutation(parameters,variables,cur_perm_image,include_mask,col);
        end
    end
    svrlsm_waitbar(parameters.waitbar,0,''); % reset
    
    %% save the cfwer results so later we can use it...
    fname = 'cfwerinfo.mat';
    fullcfwerout = fullfile(variables.output_folder.cfwer,fname);   
    save(fullcfwerout,'cfwerinfo')
    
    variables.files_created.cfwerinfo = fullcfwerout;
    variables.cfwerinfo = cfwerinfo;
    
    %% Save the resulting cluster lists
    variables.files_created.largest_clusters = fullfile(variables.output_folder.cfwer,'Largest null cluster list.mat');
    all_max_cluster_sizes = cfwerinfo.largest_cluster; % to keep it parallel with the regular non-cfwer clustering procedure.
    save(variables.files_created.largest_clusters,'all_max_cluster_sizes');