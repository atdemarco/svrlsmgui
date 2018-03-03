function stop = optim_outputfun(results,state)
    %% Output function to let the user know about progress.
    stop = false; % required. false - don't stop iterations.
    allobjs=findall(0);
    prog_rect_handle = allobjs(strcmp(arrayfun(@(x) x.Tag,allobjs,'uni',false),'progress_rectangle'));
    prog_text_handle = allobjs(strcmp(arrayfun(@(x) x.Tag,allobjs,'uni',false),'progress_text'));
    waitbar_handles = [prog_rect_handle prog_text_handle];
    switch state 
        %case 'initial'
        case 'iteration'
            svrlsm_waitbar(waitbar_handles,numel(results.IterationTimeTrace)/results.Options.MaxObjectiveEvaluations);
    end
    drawnow