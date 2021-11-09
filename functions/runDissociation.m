function [success,handles] = runDissociation(hObject,eventdata,handles)
    success = 1; % by default?
    %% First, run the two main analyses we'll draw on
    handles.dissociation = [];
    handles.dissociation.starttime = now;
    handles = runMainAnalyses(hObject,eventdata,handles);
    handles = createBetaDeltaMapForDissociation(handles); % create new beta map (deltas)
    handles = createPvalsForDissociation(handles); % convert those betas into p values
    
function handles = runMainAnalyses(hObject,eventdata,handles)
    % handles.parameters.orig_name = handles.parameters.analysis_name
    for B = 1 : 2 % Run the full analysis for each main effect
        handles.parameters.SavePermutationData = true; % we'll need this data...
        curbehav = handles.parameters.double_dissociation_behaviors{B};
        handles = UpdateProgress(handles,['Dissocation analysis ' num2str(B) ': ' curbehav],1);
        handles.parameters.score_name = curbehav; % update score_name so we analyze the right behavior
        % dissoclabel = ['dissociation_mainpart_' num2str(B) 'of2'];
        % handles.parameters.analysis_name = dissoclabel; % dynamic directory field name - so we can reference from interaction module
        [success,handles] = RunSingleAnalysis(hObject,eventdata,handles);
        handles.dissociation.maineffects{B} = handles.parameters; % so we can figure out where we should read the data from for each main effect svrb map...
    end

function handles = createBetaDeltaMapForDissociation(handles)
    handles = UpdateProgress(handles,'Dissocation analysis - computing difference of beta maps...',1);

    dissociation = handles.dissociation; % for convenience.
    %% Where we'll put the results of the dissocation ...
    dissociation.output_folder = fullfile(dissociation.maineffects{1}.baseoutputdir,'dissocation');
    dissociation.delta_svrb_file = fullfile(dissociation.output_folder,'svrb_deltas_allperms.bin');
    
    %% What function will we use to derive the joint null distribution (and our real data)?
    dissociation.joint_distribution_fcn = @(x,y) x-y; % simple difference
    % dissociation.joint_distribution_fcn = @(x,y) min(x,y); %
    % dissociation.joint_distribution_fcn = @(x,y) sqrt(x.^2 .* y.^2); %

    %% Now create a svrb delta map form the two big base file beta map files
    % files = dir(fullfile(dissociation.maineffects{1}.baseoutputdir,'*\Voxwise*\Clustwise*\pmu_beta_maps_N_*.bin')); % these are the potentially large null svrb files
    for f = 1 : 2 % Grab a memmap handle to each
        permfile = dir(fullfile(fileparts(dissociation.maineffects{f}.parmsfile),'Voxwise*','Clustwise*','pmu_beta_maps*.bin'));
        cur_bigfname = fullfile(permfile.folder,permfile.name);
        all_perm_data{f} = memmapfile(cur_bigfname,'Format','single'); %#ok<AGROW>
    end

    idx_1_n = numel(dissociation.maineffects{1}.m_idx);
    idx_2_n = numel(dissociation.maineffects{2}.m_idx);

    [m_idx,ia,ib] = intersect(dissociation.maineffects{1}.m_idx,dissociation.maineffects{2}.m_idx);
    dissociation.m_idx = m_idx; % these are out .m_idx for the dissocation
    
    dissociation.orig_betas_A = dissociation.maineffects{1}.ori_beta_vals(ia); % need to use ia so we only get indices overlapping between the 2 analyses
    dissociation.orig_betas_B = dissociation.maineffects{2}.ori_beta_vals(ib); % need to use ib so we only get indices overlapping between the 2 analyses
    dissociation.orig_betas_delta = feval(dissociation.joint_distribution_fcn,dissociation.orig_betas_A,dissociation.orig_betas_B); % apply the function we use to calculate the joint distribution...
    
    nperms = dissociation.maineffects{1}.PermNumVoxelwise; % for each permutation, pull the full relevant data frame from each .bin file -- then extract the relevant indices that overlap.

    idx1_frame_inds = 1:idx_1_n; % we'll offset these by permutation count and reuse these indices.
    idx2_frame_inds = 1:idx_2_n; % we'll offset these by permutation count and reuse these indices.

    if ~exist(dissociation.output_folder,'dir'), mkdir(dissociation.output_folder); end

    assignin('base','dissociation',dissociation)

    %% Create and save the new permutation data
    fileID = fopen(dissociation.delta_svrb_file,'w');
    for p = 1 : nperms
        thisperm_offset_1 = (p-1)*idx_1_n; % so for the first permutation our offset is 0.
        thisperm_offset_2 =  (p-1)*idx_2_n; % so for the first permutation our offset is 0.

        result1_thisperm_frame = all_perm_data{1}.Data(thisperm_offset_1 + idx1_frame_inds);
        result2_thisperm_frame = all_perm_data{2}.Data(thisperm_offset_2 + idx2_frame_inds);
        result1_thisperm_relvox = result1_thisperm_frame(ia); % extract overlapping voxel index data
        result2_thisperm_relvox = result2_thisperm_frame(ib); % extract overlapping voxels index data
        
        % calculate/combine/subtract these. this makes our new beta values. -- or perform whatever function we have specified
        delta_thisperm_relvox = feval(dissociation.joint_distribution_fcn,result1_thisperm_relvox,result2_thisperm_relvox); 
        fwrite(fileID, delta_thisperm_relvox,'single'); % append to the big file
    end
    fclose(fileID); % close big file
    
    handles.dissociation = dissociation; % for convenience.
    
    %% Save the real delta results before we return
    tmp = dissociation.maineffects{1}.vo; % a template

    %% Base A information - make sure we have indices all correct
    tmp.fname = fullfile(dissociation.output_folder,'A_valid_svrbs.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.m_idx) = dissociation.orig_betas_A;
    svrlsmgui_write_vol(tmp, outimg);

    tmp.fname = fullfile(dissociation.output_folder,'A_all_svrbs.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.maineffects{1}.m_idx) = dissociation.maineffects{1}.ori_beta_vals; % all original analysis indices... - values shuld overlap with the intersect version
    svrlsmgui_write_vol(tmp, outimg);
    
    %% Base B information - make sure we have indices all correct
    tmp.fname = fullfile(dissociation.output_folder,'B_valid_svrbs.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.m_idx) = dissociation.orig_betas_B;
    svrlsmgui_write_vol(tmp, outimg);
    
    tmp.fname = fullfile(dissociation.output_folder,'B_all_svrbs.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.maineffects{2}.m_idx) = dissociation.maineffects{2}.ori_beta_vals; % all original analysis indices... - values shuld overlap with the intersect version
    svrlsmgui_write_vol(tmp, outimg);
    
    %% Write out the lost-coverage volume (places where both analyses don't have valid indices)
    lostinds = setxor(dissociation.maineffects{1}.m_idx,dissociation.maineffects{2}.m_idx);% indices that are only found in one map result (so not included in dissociation analysis)
    tmp.fname = fullfile(dissociation.output_folder,'lost_coverage.nii');
    outimg = zeros(tmp.dim);
    if any(lostinds), outimg(lostinds) = 1; end % note these lost indices...
    svrlsmgui_write_vol(tmp, outimg);
    
    %% Now the delta - make sure we have indices all correct
    tmp.fname = fullfile(dissociation.output_folder,'A_minus_B_svrb.nii');
    outimg = zeros(tmp.dim);
    outimg(dissociation.m_idx) = dissociation.orig_betas_delta;
    svrlsmgui_write_vol(tmp, outimg);
    
function handles = createPvalsForDissociation(handles)
     handles = UpdateProgress(handles,'Dissocation analysis - computing pvals of delta-beta maps...',1);
     parameters = handles.parameters; % for convenience - we can use the most recent waitbar, etc.
    
%      assignin('base','handles',handles)
%      assignin('base','parameters',parameters)
     
     if handles.parameters.runfromgui
         parameters.waitbar = [handles.progressaxes_rectangle handles.progressaxes_text];
     end
     
     dissociation = handles.dissociation;
    
    %% Ok now read back in our delta svrb file to convert to p values
    all_perm_data = memmapfile(dissociation.delta_svrb_file,'Format','single');

    dataRef = all_perm_data.Data; % will this eliminate some overhead
    L = length(dissociation.m_idx); % for each voxel feature  - this m_idx only refers to voxel features overlapping in both main effect analyses
     for col = 1 : L % For each voxel feature (each svrb-delta)
         if ~mod(col,50) % to reduce num of calls...
             check_for_interrupt(parameters)
             svrlsm_waitbar(parameters.waitbar,col/L)
         end
         
         curcol = dataRef(col:L:end); % index out each voxel column using skips the length of the data ---> this returns curcol, which has the length of npermutations
         observed_beta = dissociation.orig_betas_delta(col); % original observed beta value.
         
         if parameters.do_CFWER
             if col == 1 % open where we'll put our p-value converted volumetric data...
                 dissociation.outfname_big_p = fullfile(dissociation.output_folder,['pmu_p_maps_N_' num2str(length(variables.m_idx)) '.bin']);
                 fileID = fopen(dissociation.outfname_big_p,'w');
             end
             
             p_vec = betas2pvals(curcol,'pos'); %parameters.tailshort);
             fwrite(fileID, p_vec,'single'); % add our results from this voxel to the same giant file....
             
             if col == L, fclose(fileID); end % then it's the last permutation, and our file is written - close the cfwer permutation data output file
         end

        % Compute beta cutoff values and a pvalue map for the real observed betas
        voxp = dissociation.maineffects{1}.voxelwise_p;
        nperms = dissociation.maineffects{1}.PermNumVoxelwise;
        onetail_cutoff_index = median([1 round(voxp * nperms) nperms]);
        [one_tail_pos_alphas(col), pos_beta_map_cutoff(col)] = compare_real_beta(observed_beta,curcol,'pos',onetail_cutoff_index); % try to preclude p values of 0
     end
     
     dissociation.one_tail_pos_alphas = one_tail_pos_alphas;
     dissociation.pos_beta_map_cutoff = pos_beta_map_cutoff;
     
     handles.dissociation = dissociation; % for convenience
     
      assignin('base','handles',handles)
%      error('a')
% 
% 
%     %% now for each permutation compute the p values???
% 
% %         %% step 2: sort the betas in the huge data file data and create cutoff values
% %         % if desired, CFWER null p-maps are also created (nothing more) in this process
% %         [parameters,variables,thresholds] = step2(handles,parameters,variables,thresholds,all_perm_data);
% 
%     % svrlsm_waitbar(parameters.waitbar,0,''); % reset.
% 
%     % allhandles{1}.vo
%     % allhandles{1}.m_idx
% 
%     % files = dir('C:\Users\ad1470\Documents\GitHub\svrlsmgui\output\myanalysis\08-Nov-2021\*\Voxwise p005 (100 perms)\Clustwise p05 (100 perms)\pmu_beta_maps_N_100.bin')
% 
% 
%     % RUN DISSOCIATION FIRST
%     disp('run dissocation code...')
%     % Dissociation results will show review of main effects
%     % Overview of A results
%     % Overview of B results 
%     % Dissociation Results (A)
%     % a montage of slices where it's A but not B
% 
%     % the Null distribution is obtained by this function:
% 
%     % dissociation.null_dist_fct = @(x,y) min([x y]);
%     dissociation.null_dist_fct = @(x,y) x^2 + y^2;
% 
%     % a montage of slices where it's B but not A
% 
