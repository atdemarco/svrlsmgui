function copyimage(infname, outfname)

[inpathstr, inname, inext] = fileparts(infname);
[outpathstr, outname] = fileparts(outfname);

if isempty(outname)
    outname = inname;
end

switch inext
    case '.hdr'
    case '.img'
        copyfile(fullfile(inpathstr, [inname '.hdr']), fullfile(outpathstr, [outname '.hdr']));
        copyfile(fullfile(inpathstr, [inname '.img']), fullfile(outpathstr, [outname '.img']));
    case '.nii'
        copyfile(fullfile(inpathstr, [inname '.nii']), fullfile(outpathstr, [outname '.nii']));
    case ''
        if exist(fullfile(inpathstr, [inname '.nii'])) == 2
            copyfile(fullfile(inpathstr, [inname '.nii']), fullfile(outpathstr, [outname '.nii']));
        elseif exist(fullfile(inpathstr, [inname '.hdr'])) == 2
            copyfile(fullfile(inpathstr, [inname '.hdr']), fullfile(outpathstr, [outname '.hdr']));
            copyfile(fullfile(inpathstr, [inname '.img']), fullfile(outpathstr, [outname '.img']));
        else
            error('No file with valid extension found');
        end
    otherwise
        error('Unrecognized extension');
end
