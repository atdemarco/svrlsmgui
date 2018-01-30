function update_available = check_for_updates
    update_available = 0; % default to 0 in case of error(s)

    try
        recent_version_date = webread('https://raw.githubusercontent.com/atdemarco/svrlsmgui/master/version.txt');
        current_version_fname = fullfile(fileparts(which('svrlsmgui')),'version.txt');
        fileID = fopen(current_version_fname,'r');
        this_download_date = fscanf(fileID,'%s',10);
        fclose(fileID);
        if datenum(recent_version_date,'yyyy-mm-dd') < datenum(this_download_date,'yyyy-mm-dd')
            update_available = 1;
        end
    end

