function can_use_cache_file = can_skip_generating_beta_perms(parameters,variables)
    % requirements:
    % 1. the file must exist already -- *this will ensure that the number of permutations is the same*
    % 2. the user allows us to draw on the cache (parameters.do_use_cache_when_available == true)
    %     [exist(parameters.outfname_big,'file') parameters.do_use_cache_when_available]

    can_use_cache_file=false; % default to false.
    if exist(parameters.outfname_big,'file') && parameters.do_use_cache_when_available
        % check the integrity (size) of the cache file vs what we are expecting
        if is_cache_file_size_ok(parameters,variables)
            % then we can use it...
            can_use_cache_file=true;

        else % if size fails to match, delete it.
            fclose all; % will this screw up something we have open?
            delete(parameters.outfname_big)
        end
    end


    function size_is_ok = is_cache_file_size_ok(parameters,variables)
    % partially written cache file will cause an error, so check size.
    info = dir(parameters.outfname_big);
    % expected size (single)
    bytes_in_a_single = 4;
    expected_size = parameters.PermNumVoxelwise * numel(variables.m_idx) * bytes_in_a_single; % this should be m_idx, not l_idx
    if info.bytes == expected_size
        size_is_ok = true;
    else
        size_is_ok = false;
    end