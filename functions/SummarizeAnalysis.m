function htmloutpath = SummarizeAnalysis(parmsfile)

% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33% "
% To set the warning state, you must first know the message identifier for the one warning you want to enable. 
% Query the last warning to acquire the identifier.  For example: 
% warnStruct = warning('query', 'last');
% messageID = warnStruct.identifier
% messageID = MATLAB:concatenation:integerInteraction
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
if ~parms.analysis_is_completed
    fprintf(parms.fileID,'<p style="color:red;"><u>Note</u>: This analysis did not finish running.</p><br><br>\n');
end

% What software was used for SVR?
if parms.useLibSVM, svmtype = 'libSVM';
else svmtype = 'MATLAB'; end 

% Was the analysis parallelized?
if parms.parallelize,  parall = ''; 
else parall = 'not '; end 

% Was the analysis run from the gui?
if parms.runfromgui, runfromgui = '';
else runfromgui = ' not '; end 

% nsubs?
% anyone excluded?
nexcluded = numel(parms.excluded_subjects);
excluded_names = strjoin(parms.excluded_subjects,', ');
nsubs = parms.nsubjects; 

% is lesion volume and one_score correlated prior to correction?
[rho,pval] = corr(parms.one_score(:),parms.lesion_vol(:),'type','Pearson','tail','both');


%% Assemble narrative summary...

% Hypothesis direction and behavior variable name
descr = ['This analysis named ''' parms.analysis_name ''' tested the hypothesis (' lower(parms.tails) ') that there is a relationship between lesion status and the behavior score ''' parms.score_name '''.'];
descr = [ descr ' ' num2str(nsubs) ' subjects were listed for inclusion, and ' num2str(nexcluded) ' were excluded due to missing behavioral data.'];
if pval < 0.05
    descr = [ descr ' Prior to any correction, the behavior under investigation is significantly correlated with lesion volume across the patient group, (r = ' num2str(rho) ', p = ' num2str(pval) '), suggesting a lesion volume control may be appropriate for the analysis.'];
else
    descr = [ descr ' Prior to any correction, the behavior under investigation is not significantly correlated with lesion volume across the patient group, (r = ' num2str(rho) ', p = ' num2str(pval) '), suggesting a lesion volume control may not be appropriate for the analysis.'];
end

% Was it corrected for lesion volume?
if strcmp(parms.lesionvolcorrection,'None'), descr = [descr ' Data was not corrected for lesion volume.'];
else descr = [descr ' Data was corrected for lesion volume via ' parms.lesionvolcorrection '.']; 
end

% What were the covariates - were they used?
n_behav_covariates = numel(parms.control_variable_names);
if n_behav_covariates > 0
    if n_behav_covariates > 1
        waswere = 'covariates were';
        itthey = 'they';
    else
        waswere = 'covariate was';
        itthey = 'it';
    end
    
    % What text do we write for the covariate inclusion info...
    if ~parms.apply_covariates_to_behavior && ~parms.apply_covariates_to_lesion
        descr = [descr ' Although ' num2str(n_behav_covariates) ' behavioral ' waswere ' specified (' strjoin(parms.control_variable_names,', ') '), ' itthey ' were not included in a nuisance model for either the behavioral score or lesion data prior to SVR.'];
    elseif parms.apply_covariates_to_behavior && ~parms.apply_covariates_to_lesion
        descr = [descr ' ' num2str(n_behav_covariates) ' behavioral ' waswere ' specified (' strjoin(parms.control_variable_names,', ') ') and covariated out of the behavioral outcome variable ''' parms.score_name ''' in a nuisance model. No nuisance model was applied to the lesion data.'];
    elseif parms.apply_covariates_to_behavior && parms.apply_covariates_to_lesion
        descr = [descr ' ' num2str(n_behav_covariates) ' behavioral ' waswere ' specified (' strjoin(parms.control_variable_names,', ') ') and covariated out of both the behavioral outcome variable ''' parms.score_name ''' and the lesion data in nuisance models.'];
    elseif ~parms.apply_covariates_to_behavior && parms.apply_covariates_to_lesion
        descr = [descr ' ' num2str(n_behav_covariates) ' behavioral ' waswere ' specified (' strjoin(parms.control_variable_names,', ') ') and covariated out of the lesion data in a nuisance model. No nuisance model was applied to the behavioral outcome variable ''' parms.score_name '''.'];
    end
else % no covariates.
    descr = [descr ' No behavioral covariates were specified, so no nuisance model was applied to the behavioral score or the lesion data prior to SVR.'];
end
    
% Was permutation testing performed?
if parms.DoPerformPermutationTesting, descr = [descr ' The resulting SVR-&beta; values were thresholded at p < ' strrep(num2str(parms.voxelwise_p),'0.','.') ' and corrected for cluster size at p < ' strrep(num2str(parms.clusterwise_p),'0.','.') ', both based on ' num2str(parms.PermNumVoxelwise) ' permutations.'];
else descr = [ descr ' Permutation testing was not performed at either the voxel level or cluster level, meaning the results should be interpreted as tentative at very best.']; 
end

% What was the minimum lesion overlap set to?
descr = [descr ' The analysis was restricted to voxels with at least ' num2str(parms.lesion_thresh) ' overlapping lesions.'];

% Date run, what svr software was used, what gamma and cost was, 
descr = [descr [' The analysis was run on ' parms.datetime_run ' using ' svmtype '''s SVR procedures (' parall 'parallelized, ' runfromgui 'run from the GUI), with parameters gamma = ' num2str(parms.gamma) ' and cost = ' num2str(parms.cost) '.' ]];

% How long did it take to run?
one_hour = 60*60; % seconds.
if parms.time.runduration > one_hour % then report in hours and minutes...
    hours = parms.time.runduration/one_hour;
    minutes = round(rem(hours,1) * 60); % what percent of 60 minutes?
    if floor(hours) > 1, hplural = 's'; else hplural=''; end  %#ok<*SEPEX>
    if minutes > 1, mplural = 's'; else mplural=''; end 
    durationstring = [num2str(floor(hours)) ' hour' hplural ' and ' num2str(minutes) ' minute' mplural '.'];
else % report in minutes
    minutes = parms.time.runduration/60; % divide by 60 secs.
    seconds = round(rem(minutes,1) * 60); % what percent of 60 secs?
    if floor(minutes) > 1, mplural = 's'; else mplural=''; end
    if seconds > 1, splural = 's'; else splural=''; end 
    durationstring = [num2str(floor(minutes)) ' minute' mplural ' and ' num2str(seconds) ' second' splural '.'];
end
descr = [descr ' The analysis completed in ' durationstring];

% Was the output p file values inverted?
if parms.invert_p_map_flag, descr = [descr ' P map output values are inverted as p-1 for ease of thresholding during display (e.g., a value of .99 corresponds to p = .01).'];
else descr = [descr ' P map output values are not inverted (i.e. a voxel value of .01 corresponds to p = .01).'];
end

% Print the narrative summary to the file.
fprintf(parms.fileID,'<h2>Analysis summary</h2>');
fprintf(parms.fileID,['<p>' descr '</p>']);

%% For our pictures...
picturedir = fullfile(parms.outdir,'images');
mkdir(picturedir); % where we'll put our pictures.

%% First, let's reslice an anatomical image from spm to the lesion space
lesionoverlapfile = dir(fullfile(parms.outdir,'All lesion overlap n=*.nii'));
lesionoverlapfile = fullfile(lesionoverlapfile(1).folder,lesionoverlapfile(1).name);
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

%%
% read lesion file
[lesionoverlapimg.hdr,lesionoverlapimg.img]=read_nifti(lesionoverlapfile);
nvoxels_any_lesion_val = nnz(lesionoverlapimg.img); % number of voxels with >0 lesions...

minlesionmask = dir(fullfile(parms.outdir,'Minimum lesion mask n=*.nii')); % angle bracket removed 10/31/17 cause it creates bad filenames in Windows environment
[minlesionmask.hdr,minlesionmask.img]=read_nifti(fullfile(minlesionmask(1).folder,minlesionmask(1).name));

%% Config what we'll show in our images
slice_bounds_percent = [.30 .80]; % highest and lowest percent over which to draw slices.
nslices = 10; % n of slices to show

%% Configure underlay image data
maxtemplatecolor = max(template.img(:)) * .75; % scale it down to enhance brightness.
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
    curanatomicalslice=fliplr(curanatomicalslice);

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
    curminslice= rot90(minlesionmask.img(:,:,slices_to_show(sl)));
    curminslice=fliplr(curminslice);
    BW = edge(curminslice,'Canny'); % edge detect.
    G(BW) = 1; % outline lesion min overlap region
    
    curRGB = cat(3,R,G,B); % combine into RGB
    
    sep = ones(size(curanatomicalslice,1),1,3); % add vertical separator...
    imdata = [imdata sep curRGB]; % concat slice RGBs on horizontal axis
end

% We need to get the number of lesions before we can label our color bar...
[~,fname]=fileparts(lesionoverlapfile);
nlesionstotal = str2num(fname(strfind(fname,'=')+1:end));  %#ok<ST2NM>

bar_location = [.03 .84 .2 .1];
imdata = PaintBarOnFrame(imdata,bar_location,cmapname,1,nlesionstotal,'Overlap of lesions');
imwrite(imdata,fullfile(picturedir,'lesion_overlap.png'));

%% Label our output  in the html file and put in the html tags to make it show up.

fprintf(parms.fileID,'<hr>'); % horizonal line
fprintf(parms.fileID,'<h2>Lesion overlaps</h2>');

nvox_meeting_lesion_minimum = nnz(minlesionmask.img(:)); % number of nonzero voxels in the mask...

imstr = ['Overlap of lesions in the patient sample (N = ' num2str(nlesionstotal) '), a green outline indicates regions meeting the minimum lesion overlap criterion for the analysis (N = '  num2str(parms.lesion_thresh) ') totaling ' num2str(nvox_meeting_lesion_minimum) ' or approximately ' num2str(round(100*(nvox_meeting_lesion_minimum/nvoxels_any_lesion_val))) '% of the ' num2str(nvoxels_any_lesion_val) ' voxels with any lesions present.<br>'];
fprintf(parms.fileID,'%s',['<p>' imstr '</p>']);

cur_alttext = imstr;
image_widths = '100%'; %num2str(round(1.5*size(imdata,2)));
image_heights = '100%'; % num2str(round(1.5*size(imdata,1)));
imtxt = ['<img src="images/lesion_overlap.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
fprintf(parms.fileID,'%s',imtxt);
fprintf(parms.fileID,'<br><br>\n');

%% Now make slices for uncorrected beta threshold map
% normalize image to 255 to index out of color map -- scale using JET from -10 to +10
cmapname='jet';
cmap = eval([cmapname '(255)']);
[unthreshbetamap.hdr,unthreshbetamap.img]=read_nifti(fullfile(parms.outdir,'Beta map (unthresholded).nii'));
raw_unthreshbetamap_img = unthreshbetamap.img~=0;
beta_scale_max = 10; % so the cmap values will fall within -10 to 10.
unthreshbetamap.img = unthreshbetamap.img + beta_scale_max;
unthreshbetamap.img = ceil(255 * (unthreshbetamap.img ./ (2*beta_scale_max)));
unthreshbetamap.img(unthreshbetamap.img==0) = 1; % this is a hack to avoid zeros because we can't index a zero out of the colormap...
imdata = [];

for sl = 1 : numel(slices_to_show)
    % Anatomical, not RGB
    curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl)));
    curanatomicalslice=fliplr(curanatomicalslice);

    [R,G,B] = deal(curanatomicalslice); % tricky deal.

    % Uncorrected beta map, will make RGB by indexing out of colormap
    uncorrbetamap = rot90(unthreshbetamap.img(:,:,slices_to_show(sl)));
    uncorrbetamap = fliplr(uncorrbetamap);
    
    relevant_pixels = fliplr(rot90(raw_unthreshbetamap_img(:,:,slices_to_show(sl))));
    %relevant_pixels = uncorrbetamap>0;
    if any(relevant_pixels(:)) % any suprathreshold voxels on this slice?
        R(relevant_pixels) = cmap(uncorrbetamap(relevant_pixels),1);
        G(relevant_pixels) = cmap(uncorrbetamap(relevant_pixels),2);
        B(relevant_pixels) = cmap(uncorrbetamap(relevant_pixels),3);
    end
    
    % Edge outline min lesion mask overlap slice in green now.
    doOutline = false;
    
    if doOutline
        curminslice = rot90(minlesionmask.img(:,:,slices_to_show(sl))); %#ok<*UNRCH>
        curminslice=fliplr(curminslice);
        BW = edge(curminslice,'Canny'); % edge detect.
        G(BW) = 1; % outline lesion min overlap region
    end
    
    curRGB = cat(3,R,G,B); % combine into RGB
    
    sep = ones(size(curanatomicalslice,1),1,3); % add vertical separator...
    imdata = [imdata sep curRGB]; % concat slice RGBs on horizontal axis
end
imdata = PaintBarOnFrame(imdata,bar_location,cmapname,-10,10,'svr-\beta (unthresh)');
imwrite(imdata,fullfile(picturedir,'uncorr_beta_map.png'));

fprintf(parms.fileID,'<hr>'); % break

imstr = 'Unthresholded SVR-&beta; map';
fprintf(parms.fileID,'<h2>%s</h2>',imstr);
cur_alttext = imstr;
imtxt = ['<img src="images/uncorr_beta_map.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
fprintf(parms.fileID,'%s\n',imtxt);
fprintf(parms.fileID,'<br><br>\n');


%% Now thresholded svr-beta map -- again scale from -10 to +10.
% normalize image to 255 to index out of color map -- scale using JET from -10 to +10
cmapname='jet';
cmap = eval([cmapname '(255)']);
voxelwisedir = dir(fullfile(parms.outdir,'Voxelwise*'));

fprintf(parms.fileID,'<hr>');
fprintf(parms.fileID,'<h2>Voxelwise thresholded SVR-&beta; map</h2>');

if ~parms.DoPerformPermutationTesting % (voxelwisedir) % permutation testing wasn't conducted
    imstr = 'Permutation testing was not conducted so there is no threshold to apply to the uncorrected SVR-&beta; map.';
    fprintf(parms.fileID,'%s<br>',imstr);
else
    voxelwisedir = fullfile(voxelwisedir(1).folder,voxelwisedir(1).name);

    fullfile(parms.outdir,voxelwisedir,'Voxelwise thresholded beta map.nii')

    [threshbetamap.hdr,threshbetamap.img]=read_nifti(fullfile(voxelwisedir,'Voxelwise thresholded beta map.nii'));
    raw_threshbetamap_img = threshbetamap.img~=0;
    beta_scale_max = 10; % so the cmap values will fall within -10 to 10.
    threshbetamap.img = threshbetamap.img + beta_scale_max;
    threshbetamap.img = ceil(255 * (threshbetamap.img ./ (2*beta_scale_max)));
    threshbetamap.img(threshbetamap.img==0) = 1; % this is a hack to avoid zeros because we can't index a zero out of the colormap...

    imdata = [];

    for sl = 1 : numel(slices_to_show)
        % Anatomical, not RGB
        curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl)));
        curanatomicalslice=fliplr(curanatomicalslice);

        [R,G,B] = deal(curanatomicalslice); % tricky deal.

        % Uncorrected beta map, will make RGB by indexing out of colormap
        corrbetamap = fliplr(rot90(threshbetamap.img(:,:,slices_to_show(sl))));

        relevant_pixels = find(fliplr(rot90(raw_threshbetamap_img(:,:,slices_to_show(sl)))));

        if any(relevant_pixels(:))
            R(relevant_pixels) = cmap(corrbetamap(relevant_pixels),1);
            G(relevant_pixels) = cmap(corrbetamap(relevant_pixels),2);
            B(relevant_pixels) = cmap(corrbetamap(relevant_pixels),3);
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

    imdata = PaintBarOnFrame(imdata,bar_location,cmapname,-10,10,['svr-\beta (p < ' num2str(parms.voxelwise_p) ', ' num2str(parms.PermNumVoxelwise) ' perms)']);
    imwrite(imdata,fullfile(picturedir,'corr_beta_map.png'));

    imstr = ['Thresholded SVR-&beta; map, p < ' num2str(parms.voxelwise_p) ' based on ' num2str(parms.PermNumVoxelwise) ' permutations.'];
    fprintf(parms.fileID,'%s<br>',imstr);
    cur_alttext = imstr;

    % write the img link html
    imtxt = ['<img src="images/corr_beta_map.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
    fprintf(parms.fileID,'%s',imtxt);

end

%% Now show cluster correction...
fprintf(parms.fileID,'<hr>');
fprintf(parms.fileID,'<h2>Clusterwise thresholded SVR-&beta; map</h2>');

if ~parms.DoPerformPermutationTesting
    imstr = 'Permutation testing was not conducted so there is no cluster correction data to display.';
    fprintf(parms.fileID,'%s<br>',imstr);
else
    clusttablefile = fullfile(voxelwisedir,'Voxelwise thresholded beta map_clustidx.nii');
    %if ~exist(clusttablefile,'file'), return; end % cheap way to avoid an error.

    [clusteridx.hdr,clusteridx.img]=read_nifti(clusttablefile);

    clusterwisedir = dir(fullfile(voxelwisedir,'Clusterwise*'));
    clusterwisedir = fullfile(clusterwisedir(1).folder,clusterwisedir(1).name);
    cluster_table = readtable(fullfile(clusterwisedir,'Table of clusters.txt'));

    % [threshbetamap.hdr,threshbetamap.img]=read_nifti(fullfile(clusterwisedir,'Voxelwise thresholded beta map.nii'));

    nonsigclusters = cluster_table.clusterP > parms.clusterwise_p;
    cluster_table(nonsigclusters,:) = []; % remove these rows...

    last_significant_cluster = max(cluster_table.idx);

    if isempty(last_significant_cluster), last_significant_cluster=0; end

    nonsig_cluster_voxels = clusteridx.img > last_significant_cluster;
    threshbetamap.img(nonsig_cluster_voxels) = 0; % zero out voxels that aren't significant clusters...

    imdata = [];

    fprintf(parms.fileID,[num2str(last_significant_cluster) ' clusters were significant at p < ' strrep(num2str(parms.clusterwise_p),'0.','.') ' based on ' num2str(parms.PermNumVoxelwise) ' permutations.']);

    for sl = 1 : numel(slices_to_show)
        % Anatomical, not RGB
        curanatomicalslice = rot90(template.img(:,:,slices_to_show(sl)));
        curanatomicalslice=fliplr(curanatomicalslice);

        [R,G,B] = deal(curanatomicalslice); % tricky deal.

        % get cluster idx slice to label clusters...
        clusteridxslice = fliplr(rot90(clusteridx.img(:,:,slices_to_show(sl))));

        % Uncorrected beta map, will make RGB by indexing out of colormap
        corrbetamap = fliplr(rot90(threshbetamap.img(:,:,slices_to_show(sl))));

        raw_relevant_pixels = fliplr(rot90(raw_threshbetamap_img(:,:,slices_to_show(sl))));
        relevant_pixels = corrbetamap~=0 & raw_relevant_pixels~=0; % mask by both so we don't get 0--> 128's...

        if any(relevant_pixels(:))
            R(relevant_pixels) = cmap(corrbetamap(relevant_pixels),1);
            G(relevant_pixels) = cmap(corrbetamap(relevant_pixels),2);
            B(relevant_pixels) = cmap(corrbetamap(relevant_pixels),3);
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
        rgbsize = size(curRGB);

        f = figure('visible','off');
        imshow(curRGB)
        %truesize; % one pixel per row/col
        % for each significant cluster, is there any relevant voxel on this frame?
        for c = 1 : last_significant_cluster
            if any(c==clusteridxslice(:)) % then this cluster appears on this slice, figure out where.
                binarize_slice = clusteridxslice == c;
                [~,x_max] = max(mean(binarize_slice,1)); % flipped for axes
                [~,y_max] = max(mean(binarize_slice,2)); % flipped for axes
                text(x_max,y_max,num2str(c),'Color','r','FontSize',10,'FontSmoothing','off');
            end
        end

        [curRGB] = frame2im(getframe); % pull frame
        close(f)
        curRGB = curRGB(1:rgbsize(1),1:rgbsize(2),:); % what pixels are we losing?
        sep = 255*ones(size(curanatomicalslice,1),1,3); % add vertical separator...
        imdata = [imdata sep curRGB]; % concat slice RGBs on horizontal axis
    end

    % we're still using jet for this plot...
    imdata = PaintBarOnFrame(imdata,bar_location,'jet',-10,10,['svr-\beta (clust p < ' strrep(num2str(parms.clusterwise_p),'0.','.') ', ' num2str(parms.PermNumVoxelwise) ' perms)']);
    imwrite(imdata,fullfile(picturedir,'signif_cluster_slices.png'));

    % Label our output  in the html file and put in the html tags to make it show up.
    fprintf(parms.fileID,'<br><br>');

    imstr = ['Labeled significant clusters, p < ' num2str(parms.clusterwise_p) ' based on ' num2str(parms.PermNumVoxelwise) ' permutations.'];
    fprintf(parms.fileID,'%s<br>',imstr);
    cur_alttext = imstr;
    imtxt = ['<img src="images/signif_cluster_slices.png" alt="' cur_alttext '" width="' image_widths '" height="' image_heights '">'];
    fprintf(parms.fileID,'%s',imtxt);
end

fprintf(parms.fileID,'<br><br>');

parms.behavioralmodeldata

%% Variable correlation diagnostics from behavioral nuisance model - added v0.08 
if isfield(parms,'behavioralmodeldata') % earlier versions don't have this (added 0.08, 9/25/17)
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Behavioral nuisance model diagnostics</h2>');
    if isempty(parms.behavioralmodeldata)
        fprintf(parms.fileID,'%s','No behavioral nuisance model was included, so no diagnostics to display.');
    else
        fprintf(parms.fileID,'%s','Correlation plot of variables included in behavioral nuisance model, including the primary behavioral predictor of interest:');
        
        % close previous figs of this name because corrplot returns no handle.
        baseinfo = get(0); % we need to find the handle of the corr plot figure...
        children = baseinfo.Children; % it's one of the children of 0 
        close(children(strcmp({children.Tag},'corrPlotFigure'))); % same named figures will cause problem
        
        set(0,'DefaultFigureVisible','off'); % so we don't show the figure when we plot it
        corrplot(parms.behavioralmodeldata,'testR','on'); % do the plotting.
        set(0,'DefaultFigureVisible','on'); % back on;
        baseinfo = get(0); % we need to find the handle of the corr plot figure...
        children = baseinfo.Children; % it's one of the children of 0 
        figindex = find(strcmp({children.Tag},'corrPlotFigure')); % match the tag, which is hard coded
        corrplothandle = children(figindex); %#ok<FNDSB> % grab its handle;
        correl_im = getframe(corrplothandle); % capture whole figure.
        close(corrplothandle); % close the fig
        
        corrfname = 'behav_nuisance_correl_im.png';
        imwrite(correl_im.cdata,fullfile(picturedir,corrfname));
        imstr = 'Correlation between variables in the behavioral nuisance model.';
        fprintf(parms.fileID,'%s<br>',imstr);
        cur_alttext = imstr;
        imtxt = ['<img src="images/' corrfname '" alt="' cur_alttext '">'];
        fprintf(parms.fileID,'%s',imtxt);
        
    end
    fprintf(parms.fileID,'<br><br>'); % break before next section
end

%% Cluster correction stability plot
fprintf(parms.fileID,'<hr>');
fprintf(parms.fileID,'<h2>Cluster correction threshold stability</h2>');

if ~parms.DoPerformPermutationTesting
    fprintf(parms.fileID,'%s','Permutation testing was not conducted, so there is no cluster correction threshold stability to display.');
else
    assess_interval=100;
    clusters_to_show = 5;
    cluster_stability_im = plotClusterPermStability(clusterwisedir,assess_interval,clusters_to_show);

    imwrite(cluster_stability_im,fullfile(picturedir,'cluster_stabil.png'));

    imstr = ['Cluster size stability over ' num2str(parms.PermNumVoxelwise) ' permutations, assessed every ' num2str(assess_interval) ' permutations. Also plotted are up to ' num2str(clusters_to_show) ' clusters regardless of significance for comparison to critical threshold.'];
    fprintf(parms.fileID,'%s<br>',imstr);
    cur_alttext = imstr;
    imtxt = ['<img src="images/cluster_stabil.png" alt="' cur_alttext '">'];
    fprintf(parms.fileID,'%s',imtxt);
end

fprintf(parms.fileID,'<br><br>\n');

% Finish the html document
FinishDocument(parms)
htmloutpath = parms.outhtmlfile;

function parms = StartDocument(parms)
    % create new file
    parms.outhtmlfile = fullfile(parms.outdir,'overview.html');
    if exist(parms.outhtmlfile,'file')
        try delete(parms.outhtmlfile); end %#ok<TRYNC>
    end

    parms.fileID = fopen(parms.outhtmlfile,'w');

    % header 
    fprintf(parms.fileID,'<!DOCTYPE html>\n');
    fprintf(parms.fileID,'<html>');
    fprintf(parms.fileID,'<body>');




function parms = FinishDocument(parms)
    % footer
    fprintf(parms.fileID,'</body>');
    fprintf(parms.fileID,'</html>');

    % clean up file
    fclose(parms.fileID);
    parms.fileID = [];
    

function im = PaintBarOnFrame(im,bar_xywh_percent,cmapname,colorminval,colormaxval,units)
    % first put the bar on it - then we'll put text on the bar
    bar_wh_pix = bar_xywh_percent(3:4) .* [size(im,2) size(im,1)];
    imbar = zeros(round(bar_wh_pix))'; % blank correctly sized bar.
    cmap = eval([cmapname '(' num2str(size(imbar,2)) ')'])'; % match the width of the requested bar with cmap.

    imbar = repmat(imbar,1,1,3); % this has RGB layers now.
    for rgb = 1 : 3 
        imbar(:,:,rgb) = repmat(cmap(rgb,:),size(imbar,1),1);
    end

    % Outline the bar in white
    imbar(:,[1 end],:)=1; imbar([1 end],:,:)=1;

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
    imshow(im);
    
    truesize(f); % one pixel per row/col % added f...?
    
    % Draw the bar annotations
    text(label_x_offset,label_y_offset,units,'Color','k','FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle','FontSmoothing','off');
    text(label_x_offset - (1.1*resize_amount*half_bar_width),label_y_offset,num2str(colorminval),'Color','w','FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle','FontSmoothing','off'); % min
    text(label_x_offset + (1.1*resize_amount*half_bar_width),label_y_offset,num2str(colormaxval),'Color','w','FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle','FontSmoothing','off'); % max
    im = frame2im(getframe);
    close(f);