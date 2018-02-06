function [handles,parameters,variables] = cleanup(handles,parameters,variables)    
    handles = UpdateProgress(handles,'Cleaning up null data...',1);
    % Clean up as necessary
    if ~parameters.SavePermutationData
        fclose('all');
        delete(parameters.outfname_big); % delete the monster bin file with raw permutation data in it.
        if exist(parameters.outfname_big,'file') % if it still exists...
            warning('Was not able to delete large binary file with raw permutation data in it. This file can be quite large, so you may want to manually clean up the file and adjust your permissions so that this is not a problem in the future.')
        end
    end
