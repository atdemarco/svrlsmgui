function htmloutpath = SummarizeAnalysis(parmsfile)
    warning('off', 'Images:initSize:adjustingMag');

    htmloutpath = [];
    parms=load(parmsfile);
    parms=parms.tosave;

    if isfield(parms,'do_make_summary') && ~parms.do_make_summary, return; end
    
    parms.outdir = fullfile(parms.baseoutputdir,[parms.score_name ', ' lower(parms.tails)]); %  force lower

    %% Start the html output file
    parms = StartDocument(parms);

%% Heading describing the analysis
fprintf(parms.fileID,['<center><h1>SVR-LSMgui (v' num2str(parms.gui_version) ') Output Overview\n</h1></center>']);
fprintf(parms.fileID,'<hr>'); %horizonal line and breaks

% Note if the analysis does not appear to be completed...
if ~parms.analysis_is_completed, fprintf(parms.fileID,'<p style="color:red;"><u>Note</u>: This analysis did not finish running.</p><br><br>\n'); end

%% Write the top summary narrative...
if parms.summary.narrative, WriteSummaryNarrative(parms); end

%% For our pictures...
parms.picturedir = fullfile(parms.outdir,'images');
mkdir(parms.picturedir); % where we'll put our pictures.

%% First, let's reslice an anatomical image from spm to the lesion space
lesionoverlapfile = parms.files_created.all_lesion_overlap;
spm_template_path = fullfile(fileparts(which('spm')),'canonical','single_subj_T1.nii');
origtemplate_localpath = fullfile(fileparts(lesionoverlapfile),'template.nii');
copyfile(spm_template_path,origtemplate_localpath); % copy the spm template to output directory.

spm('defaults','fmri');
spm_jobman('initcfg');

matlabbatch{1}.spm.spatial.coreg.write.ref = {[lesionoverlapfile ',1']}; % lesion file is the target image space
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

%% Resize images for better looking output if appropriate and possible
targetSize = [79 95 79]; % this size gives us reasonable image output...
needToResize = ~all(size(template.img) == targetSize);
canImresize = ~isempty(which('imresize3'));
if needToResize && canImresize
    disp('Resizing template for output...')
    template.img = imresize3(template.img,targetSize,'method','nearest');
end

% read lesion file
[lesionoverlapimg.hdr,lesionoverlapimg.img]=read_nifti(lesionoverlapfile);
nvoxels_any_lesion_val = nnz(lesionoverlapimg.img); % number of voxels with >0 lesions...
[minlesionmask.hdr,minlesionmask.img]=read_nifti(parms.files_created.min_lesion_mask);
nvox_meeting_lesion_minimum = nnz(minlesionmask.img(:)); % number of nonzero voxels in the mask...

if needToResize && canImresize
    disp('Resizing minimum lesion mask and lesion overlap mask for output...')
    minlesionmask.img = imresize3(minlesionmask.img,targetSize,'method','nearest');
    lesionoverlapimg.img = imresize3(lesionoverlapimg.img,targetSize,'method','nearest');
end

%% Config what we'll show in our images
slice_bounds_percent = [.30 .80]; % highest and lowest percent over which to draw slices.
nslices = 10; % n of slices to show

%% Configure underlay image data
maxtemplatecolor = max(template.img(:)) * .75; % scale it down to enhance contrast.
template.img = template.img / maxtemplatecolor;
template.img = template.img - min(template.img);
template.img = (template.img ./ max(template.img(:)));

%% Make slices for lesion overlap image
% normalize lesion overlap map to 255 to index out of color map
lesionoverlapimg.img = lesionoverlapimg.img - min(lesionoverlapimg.img(:));
lesionoverlapimg.img = ceil(255 * (lesionoverlapimg.img ./ max(lesionoverlapimg.img(:))));

dimlen = size(template.img,3); % z axis
slicebounds = floor([dimlen dimlen] .* slice_bounds_percent);
slices_to_show = floor(linspace(slicebounds(1),slicebounds(2),nslices)); % nslices slices...

imdata = [];

cmapname='hot';
cmap = eval([cmapname '(255)']); % colormap for lesion overlay.

for sl = 1 : numel(slices_to_show)
    % Anatomical, not RGB
    curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl)));
    curanatomicalslice = fliplr(curanatomicalslice);

    [R,G,B] = deal(curanatomicalslice); % tricky deal.

    % Lesion overlap, will make RGB by indexing out of colormap
    curlesionoverlapslice = rot90(lesionoverlapimg.img(:,:,slices_to_show(sl)));
    curlesionoverlapslice = fliplr(curlesionoverlapslice);
    
    relevant_pixels = find(curlesionoverlapslice>0);
    if any(relevant_pixels(:))
        R(relevant_pixels) = cmap(curlesionoverlapslice(relevant_pixels),1);
        G(relevant_pixels) = cmap(curlesionoverlapslice(relevant_pixels),2);
        B(relevant_pixels) = cmap(curlesionoverlapslice(relevant_pixels),3);
    end
    
    % Edge outline min lesion mask overlap slice in green now.
    curminslice=rot90(minlesionmask.img(:,:,slices_to_show(sl)));
    curminslice=fliplr(curminslice);
    BW = edge(curminslice,'Canny'); % edge detect.
    G(BW) = 1; % outline lesion min overlap region
    R(BW) = 0; % for contrast... 
    B(BW) = 0; % for contrast...
    curRGB = cat(3,R,G,B); % combine into RGB
    
    sep = ones(size(curanatomicalslice,1),1,3); % add vertical separator...
    imdata = [imdata sep curRGB]; % concat slice RGBs on horizontal axis
end

% We need to get the number of lesions before we can label our color bar...
[~,fname]=fileparts(lesionoverlapfile);
nlesionstotal = str2num(fname(strfind(fname,'=')+1:end));  %#ok<ST2NM>

bar_location = [.03 .84 .2 .1];

if parms.summary.lesion_overlap
    imdata = PaintBarOnFrame(imdata,bar_location,cmapname,1,nlesionstotal,'Overlap of lesions',0); % do not flip.
    imwrite(imdata,fullfile(parms.picturedir,'lesion_overlap.png'));

    %% Label our output in the html file and put in the html tags to make it show up.
    fprintf(parms.fileID,'<hr>'); % horizonal line
    fprintf(parms.fileID,'<h2>Lesion overlaps</h2>');
    imstr = ['Overlap of lesions in the patient sample (N = ' num2str(nlesionstotal) '), a green outline indicates regions meeting the minimum lesion overlap criterion for the analysis (N = '  num2str(parms.lesion_thresh) ') totaling ' num2str(nvox_meeting_lesion_minimum) ' or approximately ' num2str(round(100*(nvox_meeting_lesion_minimum/nvoxels_any_lesion_val))) '% of the ' num2str(nvoxels_any_lesion_val) ' voxels with any lesions present.<br>'];
    fprintf(parms.fileID,'%s',['<p>' imstr '</p>']);
    cur_alttext = imstr;
    image_widths = '100%';
    image_heights = '100%';
    imtxt = ['<img src="images/lesion_overlap.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
    fprintf(parms.fileID,'%s',imtxt);
    fprintf(parms.fileID,'<br><br>\n');
end

%% Now make slices for uncorrected beta threshold map
% normalize image to 255 to index out of color map -- scale using JET from -10 to +10
cmapname='jet';
cmap = eval([cmapname '(255)']);

[unthreshbetamap.hdr,unthreshbetamap.img]=read_nifti(parms.files_created.unthresholded_betamap);
if needToResize && canImresize
    disp('Resizing unthresholded beta map for output...')
    unthreshbetamap.img = imresize3(unthreshbetamap.img,targetSize,'method','nearest');
end

raw_unthreshbetamap_img = unthreshbetamap.img~=0;
beta_scale_max = 10; % so the cmap values will fall within -10 to 10.
unthreshbetamap.img = unthreshbetamap.img + beta_scale_max;
unthreshbetamap.img = ceil(255 * (unthreshbetamap.img ./ (2*beta_scale_max)));
unthreshbetamap.img(unthreshbetamap.img==0) = 1; % this is a hack to avoid zeros because we can't index a zero out of the colormap...
imdata = [];

if parms.summary.beta_map
    for sl = 1 : numel(slices_to_show)
        % Anatomical, not RGB
        curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl)));
        curanatomicalslice=fliplr(curanatomicalslice);

        [R,G,B] = deal(curanatomicalslice); % tricky deal.

        % Uncorrected beta map, will make RGB by indexing out of colormap
        uncorrbetamap = rot90(unthreshbetamap.img(:,:,slices_to_show(sl)));
        uncorrbetamap = fliplr(uncorrbetamap);

        relevant_pixels = fliplr(rot90(raw_unthreshbetamap_img(:,:,slices_to_show(sl))));
        if any(relevant_pixels(:)) % any suprathreshold voxels on this slice?
            R(relevant_pixels) = cmap(uncorrbetamap(relevant_pixels),1);
            G(relevant_pixels) = cmap(uncorrbetamap(relevant_pixels),2);
            B(relevant_pixels) = cmap(uncorrbetamap(relevant_pixels),3);
        end

        % Edge outline min lesion mask overlap slice in green now.
%         doOutline = false;
%         if doOutline
%             curminslice = rot90(minlesionmask.img(:,:,slices_to_show(sl))); %#ok<*UNRCH>
%             curminslice=fliplr(curminslice);
%             BW = edge(curminslice,'Canny'); % edge detect.
%             G(BW) = 1; % outline lesion min overlap region
%         end

        curRGB = cat(3,R,G,B); % combine into RGB

        sep = ones(size(curanatomicalslice,1),1,3); % add vertical separator...
        imdata = [imdata sep curRGB]; % concat slice RGBs on horizontal axis
    end
    imdata = PaintBarOnFrame(imdata,bar_location,cmapname,-10,10,'svr-\beta (unthresholded)',0); % do not flip.
    imwrite(imdata,fullfile(parms.picturedir,'uncorr_beta_map.png'));
    fprintf(parms.fileID,'<hr>'); % break
    imstr = 'Unthresholded SVR-&beta; map';
    fprintf(parms.fileID,'<h2>%s</h2>',imstr);
    cur_alttext = imstr;
    imtxt = ['<img src="images/uncorr_beta_map.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
    fprintf(parms.fileID,'%s\n',imtxt);
    fprintf(parms.fileID,'<br><br>\n');
end

%% Now vox thresholded pval map... (NOW we will write z maps from -3 to +3
cmapname='jet';% 'hot';%'jet';
cmap = eval([cmapname '(255)']);
    
if parms.summary.voxelwise_thresholded 
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Voxelwise thresholded Z(P) Values </h2>');
    
    if ~parms.DoPerformPermutationTesting % (voxelwisedir) % permutation testing wasn't conducted
        if ~parms.method.mass_univariate % then we are using svr-B maps...
            fprintf(parms.fileID,'%s<br>','Permutation testing was not conducted so there is no threshold to apply to the uncorrected SVR-&beta; map.');
        else % we are using regular B maps...
            fprintf(parms.fileID,'%s<br>','Permutation testing was not conducted so there is no threshold to apply to the uncorrected &beta; map.');
        end
    else % permutation was conducted so we have something to show >>
        %[threshpmap.hdr,threshpmap.img]=read_nifti(parms.files_created.thresholded_pmap);
        [threshzmap.hdr,threshzmap.img]=read_nifti(parms.files_created.thresholded_zmap);
        if needToResize && canImresize
            disp('Resizing thresholded z map image for output...')
            threshzmap.img = imresize3(threshzmap.img,targetSize,'method','nearest');
        end

        %raw_threshpmap_img = threshpmap.img~=0; % a mask.
        raw_threshzmap_img = threshzmap.img~=0; % a mask.
        
        scalerange = 3;
        threshzmap.img(threshzmap.img > scalerange) = scalerange; % max scale at Z = 3.5..
        threshzmap.img(threshzmap.img < (-1*scalerange)) = -1*scalerange; % min scale at Z = -3.5..
        
        threshzmap.img = threshzmap.img + scalerange; % bring everything to zero.
        threshzmap.img = ceil(255 * (threshzmap.img ./ (2*scalerange)));
        threshzmap.img(threshzmap.img==0) = 1; % since we use these as indices into the cmap -- 0 gives error...

        imdata = [];

        for sl = 1 : numel(slices_to_show)
            % Anatomical, not RGB
            curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl)));
            curanatomicalslice=fliplr(curanatomicalslice);

            [R,G,B] = deal(curanatomicalslice); % tricky deal.

            % Uncorrected p map, will make RGB by indexing out of colormap
            corrpmap = fliplr(rot90(threshzmap.img(:,:,slices_to_show(sl))));

            relevant_pixels = find(fliplr(rot90(raw_threshzmap_img(:,:,slices_to_show(sl)))));
            if any(relevant_pixels(:))
                R(relevant_pixels) = cmap(corrpmap(relevant_pixels),1);
                G(relevant_pixels) = cmap(corrpmap(relevant_pixels),2);
                B(relevant_pixels) = cmap(corrpmap(relevant_pixels),3);
            end

            % Edge outline min lesion mask overlap slice in green now.
            doOutline = true;
            if doOutline
                curminslice= rot90(minlesionmask.img(:,:,slices_to_show(sl)));
                curminslice=fliplr(curminslice);
                BW = edge(curminslice,'Canny'); % edge detect.
                G(BW) = 1; % outline lesion min overlap region
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

        imwrite(imdata,fullfile(parms.picturedir,imagename));
        fprintf(parms.fileID,'%s<br>',imstr);

        % write the img link html
        imtxt = ['<img src="images/' imagename '" alt="' imstr '" width="' image_widths '" height="' image_heights '">'];
        fprintf(parms.fileID,'%s',imtxt);
    end
end

% Determine what output to summary depending on whether CFWER was chosen.
if parms.summary.cfwer_diagnostics 
    if parms.DoPerformPermutationTesting && parms.do_CFWER, WriteCFWERReport(parms); end
end

if parms.summary.clusterwise_thresholded
    %% Now show cluster correction...
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Clusterwise thresholded P-map</h2>');

    if ~parms.DoPerformPermutationTesting
        imstr = 'Permutation testing was not conducted so there is no cluster correction data to display.';
        fprintf(parms.fileID,'%s<br>',imstr);
    elseif parms.DoPerformPermutationTesting && parms.do_CFWER
        imstr = 'Permutation testing was conducted, but because CFWE was chosen there is no cluster correction information to display.';
        fprintf(parms.fileID,'%s<br>',imstr);
    else % show clusters since we are just using regular null cluster distribution...
        [clusteridx.hdr,clusteridx.img]=read_nifti(parms.files_created.significant_cluster_indices);
        
        if needToResize && canImresize
            disp('Resizing cluster index image for output...')
            clusteridx.img = imresize3(clusteridx.img,targetSize,'method','nearest');
        end

        clustertable = parms.files_created.clustertable;
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

         nonsig_cluster_voxels = clusteridx.img > last_significant_cluster;
         %threshpmap.img(nonsig_cluster_voxels) = 0; % zero out voxels that aren't significant clusters...
         threshzmap.img(nonsig_cluster_voxels) = 0; % zero out voxels that aren't significant clusters...

        imdata = [];

        if last_significant_cluster == 1, plural=' was ';
        else, plural='s were ';
        end

        fprintf(parms.fileID,[num2str(last_significant_cluster) ' cluster' plural 'significant at p < ' strrep(num2str(parms.clusterwise_p),'0.','.') ' based on ' num2str(parms.PermNumVoxelwise) ' permutations.']);

        for sl = 1 : numel(slices_to_show)
            % Anatomical, not RGB
            curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl)));
            curanatomicalslice=fliplr(curanatomicalslice);

            [R,G,B] = deal(curanatomicalslice); % tricky deal.

            % get cluster idx slice to label clusters...
            clusteridxslice = fliplr(rot90(clusteridx.img(:,:,slices_to_show(sl))));

            % Uncorrected beta map, will make RGB by indexing out of colormap
            corrpmap = fliplr(rot90(threshzmap.img(:,:,slices_to_show(sl))));
            raw_relevant_pixels = fliplr(rot90(raw_threshzmap_img(:,:,slices_to_show(sl))));
            relevant_pixels = corrpmap~=0 & raw_relevant_pixels~=0; % mask by both so we don't get 0--> 128's...

            if any(relevant_pixels(:))
                R(relevant_pixels) = cmap(corrpmap(relevant_pixels),1);
                G(relevant_pixels) = cmap(corrpmap(relevant_pixels),2);
                B(relevant_pixels) = cmap(corrpmap(relevant_pixels),3);
            end

            % Outline edge of min lesion mask overlap slice in green now.
            curminslice= rot90(minlesionmask.img(:,:,slices_to_show(sl)));
            curminslice=fliplr(curminslice);
            BW = edge(curminslice,'Canny'); % edge detect.
            G(BW) = 1; % outline lesion min overlap region

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
        imwrite(imdata,fullfile(parms.picturedir,'signif_cluster_slices.png'));

        % Label our output  in the html file and put in the html tags to make it show up.
        fprintf(parms.fileID,'<br><br>');

        imstr = ['Labeled significant clusters, p < ' num2str(parms.clusterwise_p) ' based on ' num2str(parms.PermNumVoxelwise) ' permutations.'];
        fprintf(parms.fileID,'%s<br>',imstr);
        cur_alttext = imstr;
        imtxt = ['<img src="images/signif_cluster_slices.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
        fprintf(parms.fileID,'%s',imtxt);
    end
    fprintf(parms.fileID,'<br><br>');
end

%% Variable correlation diagnostics from behavioral nuisance model
if parms.summary.variable_diagnostics, WriteCorrelationDiagnostics(parms); end

%% Cluster correction stability plot
if parms.summary.cluster_stability, WriteClusterCorrectionStabilityPlot(parms); end

%% This is a report on the measure of the quality of the parameter, and uses optimalParameterReport()
if parms.summary.parameter_assessment, WriteOptimalParameterReport(parms); end

%% Write info about hyperparameter optimization record...
if parms.summary.hyperparameter_optimization_record , WriteHyperParamOptimReport(parms); end % this is a record of the hyperparameter optimization process

%% Write info about predicting behavior 
if parms.summary.predictions, WritePredictBehaviorReport(parms); end

% Finish the html document
FinishDocument(parms);
htmloutpath = parms.outhtmlfile;