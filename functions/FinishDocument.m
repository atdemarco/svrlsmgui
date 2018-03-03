function parms = FinishDocument(parms)
    % footer
    fprintf(parms.fileID,'</body>');
    fprintf(parms.fileID,'</html>');

    % clean up file
    fclose(parms.fileID);
    parms.fileID = [];
    
