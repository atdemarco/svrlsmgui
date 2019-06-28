function [handles,parameters,variables] = cleanup(handles,parameters,variables)    
    handles = UpdateProgress(handles,'Cleaning up null data...',1);
    
    % Clean up as necessary
    if ~parameters.SavePermutationData
        fclose('all');
        delete(parameters.outfname_big); % try to delete the giant binary file with the permutation beta values
        if exist(parameters.outfname_big,'file') % if it still exists...
            fileInfo = dir(parameters.outfname_big); 
            warning(['Was not able to delete large binary file with raw null SVR betas in it. This file is ' sprintf('%.2f',fileInfo.bytes/2^30) ' GB so you may want to manually delete it and adjust permissions for the future.'])
        end
    end
        
    if parameters.do_CFWER && ~parameters.SavePermutationPData % then try to delete the binary file with our pvalues from the permutations
        delete(parameters.outfname_big_p); 
        if exist(parameters.outfname_big_p,'file') % if it still exists...
            fileInfo = dir(parameters.outfname_big_p); 
            warning(['Was not able to delete large binary file with raw null pvalues in it. This file is ' sprintf('%.2f',fileInfo.bytes/2^30) ' GB so you may want to manually delete it and adjust permissions for the future.'])
        end
    end
