function V = svrlsmgui_write_vol(V,Y)
    % this is a wrapper for spm_write_vol
    % the main point of having it is to make sure data is written in a
    % format that behaves well in mricron or other programs.

    % float64 causes some weird behavior in some mricrons (?)
    % this is written out by default by SPM(12 at least) when the 
    % datatype is unknown. It is unknown in some cases of ITKSNAP output, too, maybe?
    
    % Adapted from spm_create_vol()
%     if ~isfield(V,'dt') % let's make sure it doesn't write out float64...
%         %V.dt = [spm_type('float64') spm_platform('bigend')];
%         %spm_type(V.dt(1));
%          V.dt = [spm_type('int8') spm_platform('bigend')]; % not float64.
%          Y=double(Y); % is this right?
%          warning('Writing out int8...')
%     end

    % pass data on to SPM's real spm_write_vol()
    V = spm_write_vol(V,Y);