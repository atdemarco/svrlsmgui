function plotSvrlsmDisjunctionConfiguration(parms)
    
%     assignin('base','parms',parms)
% 
    % basedir = '~/Documents/svrlsmgui/output/yesno_prcnt_vs_fluency_prcnt_fullres_5k_moresubs/13-Jan-2022/';
    basedir = fileparts(parms.outdir);
    
    dissinfo = load(fullfile(basedir,'dissociation','Dissociation Parameters.mat'));

    score1 = dissinfo.tosave.variables.maineffects{1}.score_name;
    score2 = dissinfo.tosave.variables.maineffects{2}.score_name;

    %% Create a disjunction type map (DTM) for A>B and B>A

    %% Read in the mainA and mainB voxelwise p-maps
    tmp = dir(fullfile(basedir,[score1 '*/Voxwise*/Unthresholded P map.nii']));
    mainA_unthresh_voxpmap = niftiread(fullfile(tmp(1).folder,tmp(1).name));
    tmp = dir(fullfile(basedir,[score2 '*/Voxwise*/Unthresholded P map.nii']));
    mainB_unthresh_voxpmap = niftiread(fullfile(tmp(1).folder,tmp(1).name));

    mainA_signif = (mainA_unthresh_voxpmap > 0) & mainA_unthresh_voxpmap <= .05;
    mainA_subthresh = (mainA_unthresh_voxpmap > 0) & (mainA_unthresh_voxpmap < .1);
    mainA_antirelationship = mainA_unthresh_voxpmap >= .9;
    mainA_norelationship = (mainA_unthresh_voxpmap > .1) & (mainA_unthresh_voxpmap < .9);

    mainB_signif = (mainB_unthresh_voxpmap > 0) & mainB_unthresh_voxpmap <= .05;
    mainB_subthresh = (mainB_unthresh_voxpmap > 0) & (mainB_unthresh_voxpmap < .1);
    mainB_antirelationship = mainB_unthresh_voxpmap >= .9;
    mainB_norelationship = (mainB_unthresh_voxpmap > .1) & (mainB_unthresh_voxpmap < .9);

    %% For each disjunction tail, read in the 

    disjunction_tails = {'pos','neg'};
    for dt = 1 : numel(disjunction_tails)
        curtail = disjunction_tails{dt};
        tmp = dir(fullfile(basedir,'dissociation','disjunction',['Voxwise*_' curtail],'Clust*','Significant clust indices.nii'));
        clustfile = fullfile(tmp(1).folder,tmp(1).name);
        signif_cluster_map = niftiread(clustfile) > 0;
        hdr = niftiinfo(clustfile);   

        % Add configuration types
        switch curtail
            case 'neg'
                config1 = mainA_signif & mainB_norelationship;
                config2 = mainA_signif & mainB_signif;
                config3 = mainA_subthresh & mainB_antirelationship;
                config4 = mainA_norelationship & mainB_antirelationship;
            case 'pos'
                config1 = mainB_signif & mainA_norelationship;
                config2 = mainB_signif & mainA_signif;
                config3 = mainB_subthresh & mainA_antirelationship;
                config4 = mainB_norelationship & mainA_antirelationship;
        end

        disjunction_type_map_img = 5 .* zeros(size(signif_cluster_map)); % default all to 5 = unclassified...
        disjunction_type_map_img(config1) = 1;
        disjunction_type_map_img(config2) = 2;
        disjunction_type_map_img(config3) = 3;
        disjunction_type_map_img(config4) = 4;

        outdir = fileparts(clustfile);
        outimg = zeros(size(signif_cluster_map));
        outimg(signif_cluster_map)  = disjunction_type_map_img(signif_cluster_map); % only select voxels in significant clusters.
        niftiwrite(feval(hdr.Datatype,outimg),fullfile(outdir,'disjunction_type_map.nii'),hdr)
        %niftiwrite(feval(hdr.Datatype,disjunction_type_map_img),fullfile(outdir,'disjunction_type_map_allvox.nii'),hdr)
    end
