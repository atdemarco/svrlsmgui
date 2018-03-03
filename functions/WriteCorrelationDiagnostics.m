function WriteCorrelationDiagnostics(parms)
    if isfield(parms,'behavioralmodeldata') % earlier versions don't have this (added 0.08, 9/25/17)
        fprintf(parms.fileID,'<hr>');
        fprintf(parms.fileID,'<h2>Behavioral nuisance model diagnostics</h2>');
        if isempty(parms.behavioralmodeldata)
            fprintf(parms.fileID,'%s','No behavioral nuisance model was included, so no diagnostics to display.');
        else
            fprintf(parms.fileID,'%s','Correlation plot of variables included in behavioral nuisance model, including the primary behavioral predictor of interest:<br>');
            rawdata=table2array(parms.behavioralmodeldata);
            [r,p] = corrcoef(rawdata);

            corrplothandle = figure('visible','off');
            for subp = 1:numel(p)
                drawnow
                [i,j] = ind2sub(size(p),subp);
                curSP = subplot(size(p,1),size(p,2),subp,'parent',corrplothandle);
                currp = p(i,j);
                currr = r(i,j);
                if i == 1 % left column
                    ylabel(curSP,strrep(parms.behavioralmodeldata.Properties.VariableNames{j},'_',' ')); % so latex interp doesn't subscript things
                    hold on;
                end
                if j == size(p,2) % bottom row
                    xlabel(curSP,strrep(parms.behavioralmodeldata.Properties.VariableNames{i},'_',' ')); % so latex interp doesn't subscript things
                    hold on;
                end
                if i ~= j  % then we're off diagonal.
                  scatter(rawdata(:,i),rawdata(:,j),'parent',curSP);
                  hold on;
                  if currp < .05 % alpha for highlighting the r value red.
                      color = 'r';
                  else
                      color = 'k';
                  end
                   xoffset = max(get(curSP,'xlim')) - (.9*diff(get(curSP,'xlim'))); % 90% from the right
                   yoffset = max(get(curSP,'ylim')) - (.15*diff(get(curSP,'ylim'))); % 15% down from the top
                  text(xoffset,yoffset,sprintf('r=%0.2f',currr),'Color',color,'parent',curSP)
                else % we're on diagonal. show a histogram of this variable
                    histogram(rawdata(:,i),'parent',curSP); % i == j so we'll just use i.
                    hold on;
                end
            end        

            correl_im = getframe(corrplothandle); % capture whole figure.
            close(corrplothandle); % close the fig

            corrfname = 'behav_nuisance_correl_im.png';
            imwrite(correl_im.cdata,fullfile(parms.picturedir,corrfname));
            imstr = 'Correlation between variables in the behavioral nuisance model.';
            imtxt = ['<img src="images/' corrfname '" alt="' imstr '">'];
            fprintf(parms.fileID,'%s',imtxt);
        end
        fprintf(parms.fileID,'<br><br>'); % break before next section
    end