function nii2pair(fname)

fname_base = fname(1:(end - 4));
hdr = read_nifti_hdr(fname);

fid = fopen(fname, 'r');
fseek(fid, hdr.vox_offset, 'bof');
img = fread(fid, inf, '*uint8');
fclose(fid);

hdr.vox_offset = 0;
hdr.magic = 'nii ';
hdr.magic(4) = 0;
write_nifti_hdr(hdr, [fname_base '.hdr']);

fid = fopen([fname_base '.img'], 'w');
fwrite(fid, img, 'uint8');
fclose(fid);

