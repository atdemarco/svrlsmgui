function parms = StartDocument(parms)
    % create new file
    parms.outhtmlfile = fullfile(parms.outdir,'overview.html');
    if exist(parms.outhtmlfile,'file')
        try %#ok<TRYNC>
            delete(parms.outhtmlfile); 
        end 
    end

    parms.fileID = fopen(parms.outhtmlfile,'w');

    % header 
    fprintf(parms.fileID,'<!DOCTYPE html>\n');
    fprintf(parms.fileID,'<html>');
    fprintf(parms.fileID,'<body>');
