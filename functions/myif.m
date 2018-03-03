function val = myif(teststatement,return_true,return_false)
    % convenience function
    if teststatement 
        val = return_true;
    else
        val = return_false;
    end
