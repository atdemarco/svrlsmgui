function svrlsm_waitbar(waitbar_handles,newval,message)
    if isempty(waitbar_handles)
        return
    end
%     if numel(waitbar_handles) ~= 2 % < caused by multiple svrlsmgui's open at once
%         error('Too many handles in input array')
%     end
    set_progress_percent(waitbar_handles(1),newval)
    if nargin > 2
        set_progress_text(waitbar_handles(2),message)
    end
    drawnow

% update percent progress
function set_progress_percent(rectangle_handle,newval)
    pos = get(rectangle_handle,'position');
    pos(3) = 100*newval; %update the width parameter....
    set(rectangle_handle,'position',pos)

% update message
function set_progress_text(txt_handle,message)
    set(txt_handle,'string',message)