function scancheck_svrlsm_output(outputdir)

% to check...
analyses = dir(fullfile(outputdir,'*','*','*','Analysis Parameters.mat'));
if numel(analyses) == 0
    error(['No analyses found in the subdirectories of ' outputdir '...'])
else
    disp(['Checking ' num2str(numel(analyses)) ' analyses...'])
end

for f= 1 : numel(analyses)
    curanal_root = analyses(f).folder;
    f1 = dir(fullfile(curanal_root,'*.nii'));
    f2 = dir(fullfile(curanal_root,'*','*.nii'));
    f3 = dir(fullfile(curanal_root,'*','*','*.nii'));
    f4 = dir(fullfile(curanal_root,'*','*','*','*.nii'));
    allf = [f1;f2;f3;f4];
    for g= 1 : numel(allf)
        curf = fullfile(allf(g).folder,allf(g).name);
        [hdr,img]=read_nifti(curf);
        
        % check for nans
        img_ = img(:);
        if any(isnan(img_))
            fprintf('NANs: %s\n',curf)
        end
        
        if any(img_ < 0)
            fprintf('NEGs: %s\n',curf)
        end
        
        if any(isinf(img_))
            fprintf('INFs: %s\n',curf)
        end
        
        % nans
        % data mismatches...
        % empty images...
        % negatives...
    end
end