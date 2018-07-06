function openFileInSystemViewer(fpath)
    if ispc 
        winopen(fpath);
    elseif ismac
        system(['open ' fpath ' &']); % amper starts process in background
    elseif isunix % figure out what 'open' utility command is available on this sytem of a few potentials
        potentials={'xdg-open','gio open','gvfs-open','kde-open','gnome-open','dde-open', 'exo-open'};
        [successes,~]=cellfun(@(x) unix(x),potentials,'uni',false); % ignore second outputs.
        installed = find(cell2mat(successes)==1);  % returns 1 if found.
        if isempty(installed)
            msgbox('Unable to detect a known Unix file launcher utility')
        else
            command_to_use = potentials{installed(1)};
            unix([command_to_use ' ' fpath ' &']) % amper starts process in background
        end
    else
        msgbox('Couldn''t identify platform.')
    end
