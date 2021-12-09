function [parameters,variables] = svrlsm_prepare_ica(parameters,variables)
% We'll decompose the lesion data into ICs... and use percent damage to those ROIs as the lesion data values...

 % Reread in each subject in the analysis, the mask, and make concatenated ICA file for fsl
for ni= 1 : variables.SubNum % numel(variables.SubjectID)
    fname = [variables.SubjectID{ni}, '.nii'];
    fullfname = fullfile(parameters.lesion_img_folder, fname);
    svrlsm_waitbar(parameters.waitbar,ni / length(variables.SubjectID),sprintf('Rereading lesion file %s...',fname));
    %vo = spm_vol(fullfname); % The true voxel intensities of the jth image are given by: val*V.pinfo(1,j) + V.pinfo(2,j)
    %tmp = spm_read_vols(vo);
    [hdr,tmp]=read_nifti(fullfname); % cause i dunn how to output 4d files from the spm's header.
    tmp(isnan(tmp)) = 0; % Denan the image.
    tmp = tmp > 0;  % Binarize
    Ldat(:,:,:,ni) = uint8(tmp);
    check_for_interrupt(parameters)
end

% Make the ica output directory..
mkdir(variables.output_folder.ica)

% Write a 3D mask of ANY lesions!
fname = fullfile(variables.output_folder.ica,'3DAnyLesionMask.nii');
write_nifti(hdr,double(sum(Ldat,4)>0),fname); % double so it's not booleans.
variables.files_created.ica_3D_any_lesion_mask = fname; % so we know we made this file.

% Write out the 4D concatted lesions we'll need for FSL to run melodic...
hdr.dim(1) = 4;
hdr.dim(5) = variables.SubNum;
fname = fullfile(variables.output_folder.ica,'4DConcattedLesionMasks.nii');
write_nifti(hdr,Ldat,fname);
variables.files_created.ica_4d_all_lesions = fname; % so we know we made this file.

variables.output_folder.ica_icadata = fullfile(variables.output_folder.ica,'data.ica');

setup_ica % add the right stuff to the path...

wd=pwd;

cd(variables.output_folder.ica)

[~,cattedfroot]=fileparts(variables.files_created.ica_4d_all_lesions); % ica_concat_file); % -d or --dim should be the number of components, but adding it to the end of the call just causes it to fail and come back...
[~,outmaskfroot]=fileparts(variables.files_created.ica_3D_any_lesion_mask);

%% Compute the ICA decomposition...!
ncomps = 35; % variables.SubNum - 5;  % cause nsubs is too many, and nsubs-1 is too many too. will a different "method" allow more components to be deflated?
unix(sprintf(['melodic -i "' cattedfroot '" -v -d ' num2str(ncomps) ' -m "' outmaskfroot '" -a concat --Oall --nobet --report']));

% Write out the component overlaps in 3D...
[max_prob_ic_map_inds,max_prob_ic_map_vals] = write_max_ics(variables.output_folder.ica);
variables.files_created.ica_max_prob_ic_map_inds = max_prob_ic_map_inds;
variables.files_created.ica_max_prob_ic_map_vals = max_prob_ic_map_vals;

cd(wd); % go back to wd...

% now build our ica components - ica has already been run at this point.
%[~,cattedfroot]=fileparts(ica_concat_file);
%icadir = [cattedfroot '.ica'];

buildComponentType = 2;
disp(['Building ICA component type ' num2str(buildComponentType) '.'])
switch buildComponentType
    case 1
        %% read in the component loadings for each subject
%         reports = dir(fullfile(ica_path,icadir,'report','t*.txt'));
%         componentData=[];
%         for c = 1 : numel(reports)
%             tmp =readtable(fullfile(reports(c).folder,reports(c).name));
%             componentData(:,c) = tmp.Var1; % first column...
%         end
    case 2
        [~,alllesiondata] = read_nifti(variables.files_created.ica_4d_all_lesions); % ica_concat_file); % check for percent overlap with this...
        
        %% or write_max_ics() output ...
        %[~,componentimg]=read_nifti(fullfile(pwd,icadir,'stats','maxprob_ic_map_inds.nii'));
        %[~,componentimg]=read_nifti('maxprob_ic_map_inds.nii');
        
        [~,componentimg]=read_nifti(variables.files_created.ica_max_prob_ic_map_inds); % 'maxprob_ic_map_inds.nii');
        
        for s = 1:size(alllesiondata,4)
            curlesion = alllesiondata(:,:,:,s);
            for component = 1 : max(componentimg(:))
                curcomponentmask = componentimg==component;
                curcomponent_nvox = sum(curcomponentmask(:));
                curlesion_overlap_w_curcomp = curlesion&curcomponentmask;
                n_curlesion_overlap_w_curcomp = sum(curlesion_overlap_w_curcomp(:));
                tmp = n_curlesion_overlap_w_curcomp / curcomponent_nvox;
%                 if isnan(tmp)
%                     variables.SubjectID{s}
%                     component
%                     n_curlesion_overlap_w_curcomp 
%                     curcomponent_nvox
%                     disp('---')
%                 end
                componentData(s,component) = tmp; %#ok<AGROW>
            end
        end
    case 3 % multiply out the probabilities...
%         [~,alllesiondata] = read_nifti(ica_concat_file); % check for percent overlap with this...
%         componentData=[];
%         reports = dir(fullfile(ica_path,icadir,'report','t*.txt'));
%         ncomps = numel(reports);
%         for c = 1 : ncomps
%             [~,curcomponent] = read_nifti(fullfile(ica_path,icadir,'stats',['probmap_' num2str(c) '.nii']));
%             curcomponent_maxprob = sum(curcomponent(:)); % let the probabilities sum up without binarizing.
% 
%             for s = 1 : size(alllesiondata,4)
%                 curlesion = alllesiondata(:,:,:,s);
%                 multoutprob = curcomponent(curlesion>0); % just add up the probability values...
%                 componentData(s,c) = 100 * (sum(multoutprob(:)) / curcomponent_maxprob); % percent damage to this component prob map....
%             end
%         end
end

    % now that we've collected our component data, replace the lesion data matrix...
    variables.lesion_dat_voxelwise = variables.lesion_dat;
    
    variables.lesion_dat = componentData; % these are percent damage per component!
    variables.l_idx_voxelwise = variables.l_idx;
    variables.m_idx_voxelwise = variables.m_idx;
    
    variables.l_idx = variables.l_idx_voxelwise(1:ncomps); % We'll put these back later!
    variables.m_idx = variables.l_idx_voxelwise(1:ncomps); % We'll put these back later!
    
    %assignin('base','variables',variables)
    %assignin('base','parameters',parameters)

    % size(componentData)
    % tmp = readtable(lesion_text_file);
    % data.lesion_vol = tmp.lesion_vol;
    
function [max_prob_ic_map_inds,max_prob_ic_map_vals] = write_max_ics(ica_path)
    disp('Writing max IC maps')
    %wd = '/home/crl/Documents/Projects/svrica/concat_lesions_N=48.ica/stats';
    
    wd = fullfile(ica_path,'4DConcattedLesionMasks.ica','stats');
    cd(wd); % gotta go there cause the spaces and stuff in the folder names...
    
    files=dir('prob*.nii');
    
    if numel(files) < 10 % probably need to unzip them...
        unix('gunzip *.nii.gz')
        files=dir('prob*.nii');
    end
    
    for f = 1 : numel(files)
        curf = fullfile(files(f).folder,['probmap_' num2str(f) '.nii']);
        [hdr,img]=read_nifti(curf);
        allimgs(:,:,:,f) = img; %#ok<AGROW>
    end
    
    outimg = zeros(size(img));
    outimg2 = zeros(size(img));
    
    for x = 1 : size(allimgs,1)
        for y = 1 : size(allimgs,2)
            for z = 1 : size(allimgs,3)
                curvec = squeeze(allimgs(x,y,z,:));
                if any(curvec)
                    [val,ind] = max(curvec);
                    outimg(x,y,z) = ind;
                    outimg2(x,y,z) = val;
                end
            end
        end
    end
    max_prob_ic_map_inds = fullfile(ica_path,'maxprob_ic_map_inds.nii');
    write_nifti(hdr,outimg,max_prob_ic_map_inds)
    max_prob_ic_map_vals = fullfile(ica_path,'maxprob_ic_map_vals.nii');
    write_nifti(hdr,outimg2,max_prob_ic_map_vals)

function setup_ica
    setenv('PATH', [getenv('PATH') ':/usr/share/fsl/5.0/bin:/usr/lib/fsl/5.0']);
    setenv('LD_LIBRARY_PATH',[getenv('LD_LIBRARY_PATH') ':/usr/lib/fsl/5.0/'])
    setenv('FSLOUTPUTTYPE','NIFTI_GZ')
    setenv('FSLDIR','/usr/share/fsl/5.0')    