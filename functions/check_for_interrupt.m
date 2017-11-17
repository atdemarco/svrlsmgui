function check_for_interrupt(parameters)
    % try to interrupt an analysis
    if parameters.runfromgui
        drawnow; % to allow a potentially queued cancel request to update the gcf
        if strcmp(get(gcf,'userdata'),'cancel') % check if a button press to cancel has actually been requested
            error('interrupted by user')
        end
    end