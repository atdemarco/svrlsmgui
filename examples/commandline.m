% run svrlsm on a few behaviors from a script

% make sure we can find files (svrlsmgui must be added to the path for this script)
which_svrlsmgui = which('svrlsmgui');
if isempty(which_svrlsmgui), error('SVRLSMgui is not added to the path, make sure to add it to the path to run this script.'); end
svrlsmguipath = fileparts(which_svrlsmgui);

% default files distributed with svrlsmgui
lesion_path = fullfile(svrlsmguipath,'default','lesion_imgs');
designfile_path = fullfile(svrlsmguipath,'default','PNT.csv');
output_path = fullfile(svrlsmguipath,'output');

% behaviors to run - just as a demo.
behaviors = {'randscore1','randscore2','randscore3'};

for B = behaviors
    
    current_behavior = B{1}; % this changes for each iteration
    tosave=[]; % clear it...
    tosave.analysis_root = output_path;
    tosave.score_file = designfile_path;
    tosave.score_name = current_behavior;
    tosave.lesion_img_folder = lesion_path;

    tosave.analysis_out_path = fullfile(output_path,'demo analysis');
    tosave.analysis_name = [current_behavior '_from_commandline']; % dynamic name as a demo
    
    tosave.lesion_thresh = 8;
    tosave.tails = 'One-tailed (negative)';
    tosave.lesionvolcorrection = 'Regress on Both'; 
    
    tosave.DoPerformPermutationTesting = true;
    tosave.voxelwise_p = 0.005;
    tosave.clusterwise_p = 0.05;

    tosave.parallelize = false;
    success = RunAnalysisNoGUI(tosave);
    
end