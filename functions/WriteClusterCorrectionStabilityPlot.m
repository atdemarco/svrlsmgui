function WriteClusterCorrectionStabilityPlot(parms)
    fprintf(parms.fileID,'<hr>');
    fprintf(parms.fileID,'<h2>Cluster correction threshold stability</h2>');

    if ~parms.DoPerformPermutationTesting
        fprintf(parms.fileID,'%s','Permutation testing was not conducted, so there is no cluster correction threshold stability to display.');
    elseif parms.DoPerformPermutationTesting && parms.do_CFWER
        imstr = 'Permutation testing was conducted, but CFWER was chosen. There is therefore no cluster stability information to display.';
        fprintf(parms.fileID,'%s<br>',imstr);
    else
        if parms.nclusters == 0
            fprintf(parms.fileID,'%s','Although permutation testing was conducted, there were zero voxels that passed voxelwise thresolding, and thus no clusters of any size to plot.');
        else
            assess_interval=100;
            clusters_to_show = 5;
            cluster_stability_im = plotClusterPermStability(parms,assess_interval,clusters_to_show);

            imwrite(cluster_stability_im,fullfile(parms.picturedir,'cluster_stabil.png'));

            imstr = ['Cluster size stability over ' num2str(parms.PermNumVoxelwise) ' permutations, assessed every ' num2str(assess_interval) ' permutations. Also plotted are up to ' num2str(clusters_to_show) ' clusters regardless of significance for comparison to critical threshold.'];
            fprintf(parms.fileID,'%s<br>',imstr);
            cur_alttext = imstr;
            imtxt = ['<img src="images/cluster_stabil.png" alt="' cur_alttext '">'];
            fprintf(parms.fileID,'%s',imtxt);
        end
    end

    fprintf(parms.fileID,'<br><br>\n');
