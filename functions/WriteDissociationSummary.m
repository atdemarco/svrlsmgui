function WriteDissociationSummary(parmsfile)
    warning('off', 'Images:initSize:adjustingMag');
    
    %% Read in the parameters file that saved our results and start the html output file
    parms = load(parmsfile);
    parms = parms.tosave;
    parms.outdir = parms.dissociation.output_folder.base;
    fields_to_copy = {'summary','do_CFWER','clusterwise_p','PermNumVoxelwise','DoPerformPermutationTesting','voxelwise_p'}; % inherit from first main analysis
    for f = 1 : numel(fields_to_copy)
        parms.(fields_to_copy{f}) = parms.dissociation.maineffects{1}.(fields_to_copy{f});
    end
    parms.method = parms.dissociation.maineffects{1}.method; % method for mass univariate vs not

    parms.targetSize = [79 95 79]; % this size gives us reasonable image output...
    parms = StartDocument(parms);

    %% Heading describing the analysis
    fprintf(parms.fileID,['<center><h1>SVR-LSMgui Dissociation Output Overview\n</h1></center>']);
    fprintf(parms.fileID,'<hr>');
    
    % Note if the analysis does not appear to be completed...
    % if ~parms.analysis_is_completed, fprintf(parms.fileID,'<p style="color:red;"><u>Note</u>: This analysis did not finish running.</p><br><br>\n'); end
    
    %% Write the top summary narrative...
    %if parms.summary.narrative, WriteSummaryNarrative(parms); end
    fprintf(parms.fileID,'<h2>Dissociation summary</h2>');
    behavA = parms.dissociation.maineffects{1}.double_dissociation_behaviors{1};
    behavB = parms.dissociation.maineffects{1}.double_dissociation_behaviors{2};
    descr = ['The dissociation analysis examined the spatial dissociation and colocalization of the lesion-symptom map for "' behavA '" and "' behavB '". The analysis looked for regions where ' behavA ' was more involved than ' behavB ' and vice versa. The analysis also looked for regions where both behaviors were involved.'];
    fprintf(parms.fileID,['<p>' descr '</p>']);
    
    %% For our pictures...
    parms.picturedir = fullfile(parms.outdir,'images');
    mkdir(parms.picturedir); % where we'll put our pictures.
    
    %% First, let's reslice an anatomical image from spm to the lesion space
    template = getReslicedTemplate(parms);
    files = getBasicFilesToShow(parms); % this includes many different files.

    %% Config what we'll show in our images
    slice_bounds_percent = [.30 .80]; % highest and lowest percent over which to draw slices.
    nslices = 10; % n of slices to show
    dimlen = size(template.img,3); % z axis
    slicebounds = floor([dimlen dimlen] .* slice_bounds_percent);
    slices_to_show = floor(linspace(slicebounds(1),slicebounds(2),nslices)); % nslices slices...
    bar_location = [.03 .84 .2 .1];
    image_widths = '100%'; 
    image_heights = '100%';
    cmapname='jet'; % normalize image to 255 to index out of color map -- scale using JET from -10 to +10
    cmap = eval([cmapname '(255)']);
    
    %% Show coverage map
    files_to_show={'A_valid_svrbs','B_valid_svrbs','lost_coverage'};
    filedescr = {['Valid voxels ' behavA],['Valid voxels ' behavB],'Lost coverage'};
    for f = 1 : numel(files_to_show)
        thisimg = files.(files_to_show{f}).img;
        raw_thisimg = thisimg~=0;
        beta_scale_max = 10; % so the cmap values will fall within -10 to 10.
        thisimg = thisimg + beta_scale_max;
        thisimg = ceil(255 * (thisimg ./ (2*beta_scale_max)));
        thisimg(thisimg==0) = 1; % this is a hack to avoid zeros because we can't index a zero out of the colormap...
        thisimg(thisimg>1) = 255;% so we get a binary map...

        imdata = [];
        for sl = 1 : numel(slices_to_show)
            curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl))); % Anatomical, not RGB
            curanatomicalslice=fliplr(curanatomicalslice);
            [R,G,B] = deal(curanatomicalslice); % tricky deal.
            % Make RGB by indexing out of colormap
            maptoshow = rot90(thisimg(:,:,slices_to_show(sl)));
            maptoshow = fliplr(maptoshow);

            relevant_pixels = fliplr(rot90(raw_thisimg(:,:,slices_to_show(sl))));
            if any(relevant_pixels(:)) % any suprathreshold voxels on this slice?
                R(relevant_pixels) = cmap(maptoshow(relevant_pixels),1);
                G(relevant_pixels) = cmap(maptoshow(relevant_pixels),2);
                B(relevant_pixels) = cmap(maptoshow(relevant_pixels),3);
            end

            curRGB = cat(3,R,G,B); % combine into RGB

            sep = ones(size(curanatomicalslice,1),1,3); % add vertical separator...
            imdata = [imdata sep curRGB]; % concat slice RGBs on horizontal axis
        end
        imdata = PaintBarOnFrame(imdata,bar_location,cmapname,0,1,'Coverage',0); % do not flip.
        imwrite(imdata,fullfile(parms.picturedir,[files_to_show{f} '.png']));
        fprintf(parms.fileID,'<hr>'); % break
        imstr = ['Coverage map (' filedescr{f} ')'];
        fprintf(parms.fileID,'<h2>%s</h2>',imstr);
        cur_alttext = imstr;
        imtxt = ['<img src="images/' files_to_show{f} '.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
        fprintf(parms.fileID,'%s\n',imtxt);
        fprintf(parms.fileID,'<br><br>\n');
    end
    
    
    %% Make slices for uncorrected beta threshold map - disjunction image, then conjunction image
    imdata = [];
    DISSOCIATION_TYPE = {'disjunction','conjunction'};
    FILE_NAMES = {'A_disjunct_B_svrb','A_conjunct_B_svrb'};
    for F = 1 : numel(FILE_NAMES)
        thisimg = files.(FILE_NAMES{F}).img;
        raw_thisimg = thisimg~=0;
        beta_scale_max = 10; % so the cmap values will fall within -10 to 10.
        thisimg = thisimg + beta_scale_max;
        thisimg = ceil(255 * (thisimg ./ (2*beta_scale_max)));
        thisimg(thisimg==0) = 1; % this is a hack to avoid zeros because we can't index a zero out of the colormap...
        imdata = [];

        for sl = 1 : numel(slices_to_show)
            curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl))); % Anatomical, not RGB
            curanatomicalslice=fliplr(curanatomicalslice);
            [R,G,B] = deal(curanatomicalslice); % tricky deal.
            % Make RGB by indexing out of colormap
            maptoshow = rot90(thisimg(:,:,slices_to_show(sl)));
            maptoshow = fliplr(maptoshow);

            relevant_pixels = fliplr(rot90(raw_thisimg(:,:,slices_to_show(sl))));
            if any(relevant_pixels(:)) % any suprathreshold voxels on this slice?
                R(relevant_pixels) = cmap(maptoshow(relevant_pixels),1);
                G(relevant_pixels) = cmap(maptoshow(relevant_pixels),2);
                B(relevant_pixels) = cmap(maptoshow(relevant_pixels),3);
            end

            curRGB = cat(3,R,G,B); % combine into RGB

            sep = ones(size(curanatomicalslice,1),1,3); % add vertical separator...
            imdata = [imdata sep curRGB]; % concat slice RGBs on horizontal axis
        end
        imdata = PaintBarOnFrame(imdata,bar_location,cmapname,-10,10,'svr-\beta (unthresholded)',0); % do not flip.
        imwrite(imdata,fullfile(parms.picturedir,[DISSOCIATION_TYPE{F} '_uncorr_beta_map.png']));
        fprintf(parms.fileID,'<hr>'); % break
        imstr = ['Unthresholded SVR-&beta; map (' DISSOCIATION_TYPE{F} ')'];
        fprintf(parms.fileID,'<h2>%s</h2>',imstr);
        cur_alttext = imstr;
        imtxt = ['<img src="images/' DISSOCIATION_TYPE{F} '_uncorr_beta_map.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
        fprintf(parms.fileID,'%s\n',imtxt);
        fprintf(parms.fileID,'<br><br>\n');
    end

    %% Now vox thresholded pval map... (NOW we will write z maps from -3 to +3)
    if parms.summary.voxelwise_thresholded
        if parms.DoPerformPermutationTesting && ~parms.method.mass_univariate
            voxelwisefiles = getVoxelwiseFilesToShow(parms); % this includes voxelwise thresholded files
        end
        field_names = {'disjunction_pos','disjunction_neg','conjunction'};
        field_names_descr = {'Disjunction (positive tail)','Disjunction (negative tail)','Conjunction'};
        for f = 1 : numel(field_names)
            fprintf(parms.fileID,'<hr>');
            fprintf(parms.fileID,['<h2>Voxelwise thresholded Z(P) Values ' field_names_descr{f} '</h2>']);
        
            if ~parms.DoPerformPermutationTesting % (voxelwisedir) % permutation testing wasn't conducted
                if ~parms.method.mass_univariate % then we are using svr-B maps...
                    fprintf(parms.fileID,'%s<br>','Permutation testing was not conducted so there is no threshold to apply to the uncorrected SVR-&beta; map.');
                else % we are using regular B maps...
                    fprintf(parms.fileID,'%s<br>','Permutation testing was not conducted so there is no threshold to apply to the uncorrected &beta; map.');
                end
            else % permutation was conducted so we have something to show >>
                curfieldname = field_names{f};
                curimg = voxelwisefiles.(curfieldname).img;
                raw_threshzmap_img = curimg~=0; % a mask.

                scalerange = 3;
                curimg(curimg > scalerange) = scalerange; % max scale at Z = 3.5..
                curimg(curimg < (-1*scalerange)) = -1*scalerange; % min scale at Z = -3.5..

                curimg = curimg + scalerange; % bring everything to zero.
                curimg = ceil(255 * (curimg ./ (2*scalerange)));
                curimg(curimg==0) = 1; % since we use these as indices into the cmap -- 0 gives error...
                imdata = [];
                for sl = 1 : numel(slices_to_show)
                    curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl))); % Anatomical, not RGB
                    curanatomicalslice=fliplr(curanatomicalslice);
                    [R,G,B] = deal(curanatomicalslice); % tricky deal.
                    corrpmap = fliplr(rot90(curimg(:,:,slices_to_show(sl)))); % Uncorrected p map, will make RGB by indexing out of colormap

                    relevant_pixels = find(fliplr(rot90(raw_threshzmap_img(:,:,slices_to_show(sl)))));
                    if any(relevant_pixels(:))
                        R(relevant_pixels) = cmap(corrpmap(relevant_pixels),1);
                        G(relevant_pixels) = cmap(corrpmap(relevant_pixels),2);
                        B(relevant_pixels) = cmap(corrpmap(relevant_pixels),3);
                    end
                    curRGB = cat(3,R,G,B); % combine into RGB
                    sep = ones(size(curanatomicalslice,1),1,3); % add vertical separator...
                    imdata = [imdata sep curRGB]; %#ok<*AGROW> % concat slice RGBs on horizontal axis
                end

                if parms.do_CFWER
                    imdata = PaintBarOnFrame(imdata,bar_location,cmapname,-1*scalerange,scalerange,['CFWER Z(p)-value (p < ' num2str(parms.voxelwise_p) ', ' num2str(parms.PermNumVoxelwise) ' perms)'],-1);
                    imstr = ['CFWER thresholded Z(p) values, p < ' num2str(parms.voxelwise_p) ' based on ' num2str(parms.PermNumVoxelwise) ' permutations (v = XXX, FWE = XXX)'];
                    imagename = 'cfwer_thresh_z_of_p_map.png';
                else % showing the voxelwise thresholded result from the regular svr-beta thresholding
                    imdata = PaintBarOnFrame(imdata,bar_location,cmapname,-1*scalerange,scalerange,['Z(p)-value (p < ' num2str(parms.voxelwise_p) ', ' num2str(parms.PermNumVoxelwise) ' perms)'],-1);
                    imstr = ['Voxelwise thresholded Z(P) values, p < ' num2str(parms.voxelwise_p) ' based on ' num2str(parms.PermNumVoxelwise) ' permutations.'];
                    imagename = 'vox_thresh_z_of_p_map.png';
                end

                imagename = [curfieldname '_' imagename]; % prefix with the type of dissociation
                
                imwrite(imdata,fullfile(parms.picturedir,imagename));
                fprintf(parms.fileID,'%s<br>',imstr);
                % write the img link html
                imtxt = ['<img src="images/' imagename '" alt="' imstr '" width="' image_widths '" height="' image_heights '">'];
                fprintf(parms.fileID,'%s',imtxt);
            end
        end
    end
    
    % Determine what output to summary depending on whether CFWER was chosen.
    if parms.summary.cfwer_diagnostics
        if parms.DoPerformPermutationTesting && parms.do_CFWER, WriteCFWERReport(parms); end
    end
    
    %dev1
    if parms.summary.clusterwise_thresholded
        if parms.DoPerformPermutationTesting && ~parms.method.mass_univariate
            clusterwisefiles = getClusterwiseFilesToShow(parms); % this includes clusterwise thresholded files
        end
        field_names = {'disjunction_pos','disjunction_neg','conjunction'};
        field_names_descr = {'Disjunction (positive tail)','Disjunction (negative tail)','Conjunction'};
        for f = 1 : numel(field_names)
            curfield = field_names{f};
            curdescr = field_names_descr{f};
            %% Now show cluster correction...
            fprintf(parms.fileID,'<hr>');
            fprintf(parms.fileID,['<h2>Clusterwise thresholded P(z)-map ' curdescr '</h2>']);
            if ~parms.DoPerformPermutationTesting
                imstr = 'Permutation testing was not conducted so there is no cluster correction data to display.';
                fprintf(parms.fileID,'%s<br>',imstr);
            elseif parms.DoPerformPermutationTesting && parms.do_CFWER
                imstr = 'Permutation testing was conducted, but because CFWER was chosen there is no cluster correction information to display.';
                fprintf(parms.fileID,'%s<br>',imstr);
            else % show clusters since we are just using regular null cluster distribution...
                clusterimg = clusterwisefiles.(curfield).img; %[clusteridx.hdr,clusteridx.img]=read_nifti(parms.files_created.significant_cluster_indices);
                clustertable = clusterwisefiles.(curfield).clustertable; %parms.files_created.clustertable;
                
                if exist(clustertable,'file')
                    cluster_table = readtable(clustertable);
                    nonsigclusters = cluster_table.clusterP > parms.clusterwise_p;
                    cluster_table(nonsigclusters,:) = []; % remove these rows...
                    last_significant_cluster = max(cluster_table.idx);
                    parms.nclusters = numel(parms.clusterwise_p);
                else
                    last_significant_cluster= 0;
                    parms.nclusters = 0;
                end

                if isempty(last_significant_cluster), last_significant_cluster=0; end

                nonsig_cluster_voxels = clusterimg > last_significant_cluster;

                %
                
                curimg = voxelwisefiles.(curfieldname).img;
                raw_threshzmap_img = curimg~=0; % a mask - we need to do this for the 3 different cluster outputs we're looping through...
                scalerange = 3;
                curimg(curimg > scalerange) = scalerange; % max scale at Z = 3.5..
                curimg(curimg < (-1*scalerange)) = -1*scalerange; % min scale at Z = -3.5..

                curimg = curimg + scalerange; % bring everything to zero.
                curimg = ceil(255 * (curimg ./ (2*scalerange)));
                curimg(curimg==0) = 1; % since we use these as indices into the cmap -- 0 gives error...

                %
                
                curimg(nonsig_cluster_voxels) = 0; % zero out voxels that aren't significant clusters...

                imdata = [];

                if last_significant_cluster == 1, plural=' was ';
                else, plural='s were ';
                end

                fprintf(parms.fileID,[num2str(last_significant_cluster) ' cluster' plural 'significant at p < ' strrep(num2str(parms.clusterwise_p),'0.','.') ' based on ' num2str(parms.PermNumVoxelwise) ' permutations.']);

                for sl = 1 : numel(slices_to_show)
                    curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl)));
                    curanatomicalslice = fliplr(curanatomicalslice);
                    [R,G,B] = deal(curanatomicalslice); % tricky deal.
                    clusteridxslice = fliplr(rot90(clusterimg(:,:,slices_to_show(sl)))); % get cluster idx slice to label clusters...
                    corrpmap = fliplr(rot90(curimg(:,:,slices_to_show(sl))));
                    raw_relevant_pixels = fliplr(rot90(raw_threshzmap_img(:,:,slices_to_show(sl))));
                    relevant_pixels = corrpmap~=0 & raw_relevant_pixels~=0; % mask by both so we don't get 0--> 128's...

                    if any(relevant_pixels(:))
                        R(relevant_pixels) = cmap(corrpmap(relevant_pixels),1);
                        G(relevant_pixels) = cmap(corrpmap(relevant_pixels),2);
                        B(relevant_pixels) = cmap(corrpmap(relevant_pixels),3);
                    end

                    curRGB = cat(3,R,G,B); % combine into RGB
                    rgbsize = size(curRGB);

                    f = figure('visible','off');  % new in v 2 to make text look better ...?
                    a=axes(f);
                    imshow(curRGB,'Parent',a); %truesize; % one pixel per row/col

                    % for each significant cluster, is there any relevant voxel on this frame?
                    for c = 1 : last_significant_cluster
                        if any(c==clusteridxslice(:)) % then this cluster appears on this slice, figure out where.
                            binarize_slice = clusteridxslice == c;
                            [~,x_max] = max(mean(binarize_slice,1)); % flipped for axes
                            [~,y_max] = max(mean(binarize_slice,2)); % flipped for axes
                            text(x_max,y_max,num2str(c),'Color','r','FontSize',10,'FontSmoothing','off');
                        end
                    end

                    [curRGB] = frame2im(getframe(a)); % pull frame
                    close(f)
                    curRGB = curRGB(1:rgbsize(1),1:rgbsize(2),:); % what pixels are we losing?
                    sep = 255*ones(size(curanatomicalslice,1),1,3); % add vertical separator...
                    imdata = [imdata sep curRGB]; % concat slice RGBs on horizontal axis
                end

                % we're still using jet for this plot...
                imdata = PaintBarOnFrame(imdata,bar_location,'jet',-1*scalerange,scalerange,['Z(p)-value (p < ' strrep(num2str(parms.clusterwise_p),'0.','.') ', ' num2str(parms.PermNumVoxelwise) ' perms)'],0);
                imwrite(imdata,fullfile(parms.picturedir,[curfield '_signif_cluster_slices.png']));

                % Label our output  in the html file and put in the html tags to make it show up.
                fprintf(parms.fileID,'<br><br>');

                imstr = ['Labeled significant clusters, p < ' num2str(parms.clusterwise_p) ' based on ' num2str(parms.PermNumVoxelwise) ' permutations.'];
                fprintf(parms.fileID,'%s<br>',imstr);
                cur_alttext = imstr;
                imtxt = ['<img src="images/' curfield '_signif_cluster_slices.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
                fprintf(parms.fileID,'%s',imtxt);
            end 
        end
        fprintf(parms.fileID,'<br><br>');
    end
    
    % Finish the html document
    FinishDocument(parms);
    htmloutpath = parms.outhtmlfile;
    
function template = getReslicedTemplate(parms)
    unthresholded_betamap_file = parms.variables.files_created.unthresholded_betamap; 
    spm_template_path = fullfile(fileparts(which('spm')),'canonical','single_subj_T1.nii');
    origtemplate_localpath = fullfile(fileparts(unthresholded_betamap_file),'template.nii');
    copyfile(spm_template_path,origtemplate_localpath); % copy the spm template to output directory.
    
    spm('defaults','fmri');
    spm_jobman('initcfg');
    
    matlabbatch{1}.spm.spatial.coreg.write.ref = {[unthresholded_betamap_file ',1']}; % lesion file is the target image space
    matlabbatch{1}.spm.spatial.coreg.write.source = {[origtemplate_localpath ',1']}; % SPM's single subject t1 template...
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r'; % will this overwrite?
    spm_jobman('run',matlabbatch);
    delete(origtemplate_localpath); % delete the local copy of the non-resliced template...
    [fpath,fname,ext]=fileparts(origtemplate_localpath);
    resliced_template_path = fullfile(fpath,['r' fname ext]);
    [template.hdr,template.img] = read_nifti(resliced_template_path); % read in resliced template...

    % Resize images for better looking output if appropriate and possible
    needToResize = ~all(size(template.img) == parms.targetSize);
    canImresize = ~isempty(which('imresize3'));
    if needToResize && canImresize, template.img = imresize3(template.img,parms.targetSize,'method','nearest'); end

    %% Configure underlay image data
    maxtemplatecolor = max(template.img(:)) * .75; % scale it down to enhance contrast.
    template.img = template.img / maxtemplatecolor;
    template.img = template.img - min(template.img);
    template.img = (template.img ./ max(template.img(:)));

function files = getBasicFilesToShow(parms)
    unthresholded_betamap_file = parms.variables.files_created.unthresholded_betamap; 
    base_dissoc_dir = fileparts(unthresholded_betamap_file);
    files_to_read = {'A_valid_svrbs','B_valid_svrbs','lost_coverage','A_disjunct_B_svrb','A_conjunct_B_svrb'};
    
    for f = 1 : numel(files_to_read)
        curfilename = files_to_read{f};
        fullfilepath = fullfile(base_dissoc_dir,[curfilename '.nii']);
        [files.(curfilename).hdr,files.(curfilename).img]=read_nifti(fullfilepath);
        needToResize = ~all(size(files.(curfilename).img) == parms.targetSize);
        canImresize = ~isempty(which('imresize3'));
        if needToResize && canImresize, files.(curfilename).img = imresize3(files.(curfilename).img,parms.targetSize,'method','nearest'); end
    end

function voxelwisefiles = getVoxelwiseFilesToShow(parms)
    unthresholded_betamap_file = parms.variables.files_created.unthresholded_betamap; 
    base_dissoc_dir = fileparts(unthresholded_betamap_file);
    files_to_read = {'disjunction\*_pos\Voxelwise thresholded Z map.nii','disjunction\*_neg\Voxelwise thresholded Z map.nii','conjunction\*\Voxelwise thresholded Z map.nii'};
    field_names = {'disjunction_pos','disjunction_neg','conjunction'};

    for f = 1 : numel(files_to_read)
        curfilename = files_to_read{f};
        curfilename = dir(fullfile(base_dissoc_dir,curfilename));
        curfilename = fullfile(curfilename(1).folder,curfilename(1).name);
        curfieldname = field_names{f};
        [voxelwisefiles.(curfieldname).hdr,voxelwisefiles.(curfieldname).img]=read_nifti(curfilename);
        needToResize = ~all(size(voxelwisefiles.(curfieldname).img) == parms.targetSize);
        canImresize = ~isempty(which('imresize3'));
        if needToResize && canImresize, voxelwisefiles.(curfieldname).img = imresize3(voxelwisefiles.(curfieldname).img,parms.targetSize,'method','nearest'); end
    end

function clusterwisefiles = getClusterwiseFilesToShow(parms)
    unthresholded_betamap_file = parms.variables.files_created.unthresholded_betamap; 
    base_dissoc_dir = fileparts(unthresholded_betamap_file);
    if parms.do_CFWER % only a single embedding
        error('finish the pathing for the cluster-corrected cfwer files')
        files_to_read = {'disjunction\*_pos\Voxelwise thresholded Z map.nii','disjunction\*_neg\Voxelwise thresholded Z map.nii','conjunction\*\Voxelwise thresholded Z map.nii'};
    else % double embedding
        files_to_read = {'disjunction\*_pos\Clust*\Significant clust indices.nii','disjunction\*_neg\Clust*\Significant clust indices.nii','conjunction\*\Clust*\Significant clust indices.nii'};
    end
    field_names = {'disjunction_pos','disjunction_neg','conjunction'};

    for f = 1 : numel(files_to_read)
        curfilename = files_to_read{f};
        curfilename = dir(fullfile(base_dissoc_dir,curfilename));
        curfilename = fullfile(curfilename(1).folder,curfilename(1).name);
        curfieldname = field_names{f};
        [clusterwisefiles.(curfieldname).hdr,clusterwisefiles.(curfieldname).img]=read_nifti(curfilename);

        clusterwisefiles.(curfieldname).clustertable = fullfile(fileparts(curfilename),'Table of clusters.txt');
        
        needToResize = ~all(size(clusterwisefiles.(curfieldname).img) == parms.targetSize);
        canImresize = ~isempty(which('imresize3'));
        if needToResize && canImresize, clusterwisefiles.(curfieldname).img = imresize3(clusterwisefiles.(curfieldname).img,parms.targetSize,'method','nearest'); end
    end
