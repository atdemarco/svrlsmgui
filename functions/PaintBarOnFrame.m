function im = PaintBarOnFrame(im,bar_xywh_percent,cmapname,colorminval,colormaxval,units,flip)
    % first put the bar on it - then we'll put text on the bar
    bar_wh_pix = bar_xywh_percent(3:4) .* [size(im,2) size(im,1)];
    imbar = zeros(round(bar_wh_pix))'; % blank correctly sized bar.
    cmap = eval([cmapname '(' num2str(size(imbar,2)) ')'])'; % match the width of the requested bar with cmap.

    if flip % for e.g., pmaps...
        cmap = fliplr(cmap);
    end
           
    imbar = repmat(imbar,1,1,3); % this has RGB layers now.
    for rgb = 1 : 3 
        imbar(:,:,rgb) = repmat(cmap(rgb,:),size(imbar,1),1);
    end

    % Outline the bar in white
    imbar(:,[1 end],:)=1; 
    imbar([1 end],:,:)=1;

    % normalize our color scale to the max of the im in which we will set it into 
    if max(im(:)) > 1
        imbar = double(imbar) .* double(max(im(:))); % so if we get in a 1-255 map we can still see our colors.
    end
    xyoffset = round(bar_xywh_percent(1:2) .* [size(im,2) size(im,1)]); % where to place the bar?
    im(xyoffset(2):xyoffset(2)+size(imbar,1)-1,xyoffset(1):xyoffset(1)+size(imbar,2)-1,:) = imbar;

    resize_amount = 2;

    % calculate offsets of label - center of colorbar in pixels.
    half_bar_height = round(size(imbar,1)/2);
    half_bar_width = round(size(imbar,2)/2);
    label_x_offset = resize_amount * (xyoffset(1) + half_bar_width); 
    label_y_offset = (resize_amount * (xyoffset(2) + half_bar_height)) - (.3*half_bar_height); % tiny adjustment at end.

    im = imresize(im, resize_amount,'bilinear','Antialiasing',true); % here we upsample so we get good looking font, I hope

    f = figure('visible','off');
    a = axes(f);
    imshow(im,'Parent',a);
    
    truesize(f); % one pixel per row/col
    
    % Draw the bar annotations
    text(label_x_offset,label_y_offset,units,'Color','k','FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle','FontSmoothing','off','Parent',a);
    text(label_x_offset - (1.1*resize_amount*half_bar_width),label_y_offset,num2str(colorminval),'Color','w','FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle','FontSmoothing','off','Parent',a); % min
    text(label_x_offset + (1.1*resize_amount*half_bar_width),label_y_offset,num2str(colormaxval),'Color','w','FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle','FontSmoothing','off','Parent',a); % max
    im = frame2im(getframe(a)); % do we want getframe a or f?
    close(f);