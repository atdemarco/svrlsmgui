function OpenDirectoryInNativeWindow(dirpath)
    if isunix % this should at least work in Ubuntu...
        unix(['xdg-open ' dirpath ' &']);
    elseif ismac
        [~] = system(['open "' dirpath '"']);
    elseif ispc
        winopen(dirpath)
    else
        warndlg('Cannot open output directory because cannot determine OS in use.')
    end
