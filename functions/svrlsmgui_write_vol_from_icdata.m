function V = svrlsmgui_write_vol_from_icdata(parameters,variables,icdata)
% then redistribute the values back into the whole-brain volume...

%     assignin('base','variables',variables)
%     assignin('base','parameters',parameters)
%     assignin('base','icdata',icdata)
%     error('a')    

    V = variables.vo;
 
    [~,componentimg]=read_nifti(variables.files_created.ica_max_prob_ic_map_inds);

    % variables.l_idx_voxelwise = variables.l_idx;
    % variables.m_idx_voxelwise = variables.m_idx;

    
    tmp = zeros(variables.vo.dim(1:3));
    beta_map = tmp;

    for component = 1 : numel(icdata)
        curcomponentmask = componentimg==component;
        tmp(curcomponentmask) = icdata(component);
    end

    %tmp(variables.l_idx) = icdata; % return all lesion data to its l_idx indices.
    beta_map(variables.m_idx_voxelwise) = tmp(variables.m_idx_voxelwise); % m_idx -> m_idx

    % pass data on to SPM's real spm_write_vol()
    Y = beta_map;
    V = spm_write_vol(V,Y);