function make_3d_brain
% X3DOM DEMOS ZOMG:
% https://www.x3dom.org/examples/ wow
% https://examples.x3dom.org/volren/volrenOpacityTestTF_aorta.xhtml
% https://doc.x3dom.org/tutorials/animationInteraction/transformations/example.html
% IMPORTANT: https://doc.x3dom.org/tutorials/animationInteraction/picking/example.html
% for thresholding?: https://doc.x3dom.org/tutorials/animationInteraction/onoutputchange/example.html

% less good:
% https://doc.x3dom.org/tutorials/animationInteraction/pickingBuffer/example.html

%cd(wd)

%fpath = 'C:\Users\andrew\Documents\mricron\templates';
fpath =  '/home/crl/Downloads/templates/ch2better.nii.gz';
img = niftiread(fpath);
% [hdr,img]=read_nifti(fullfile(fpath,'ch2bet.nii')); % niftiread
disp('1')
f=figure;
%colormap(map)
D=img;

amount=[5 5 5];

Ds = smooth3(D,'box',amount);
disp('2')
fv = patch(isosurface(Ds,5),'FaceColor',[1,.75,.65],'EdgeColor','none');
   
isonormals(Ds,fv)
disp('3')
p = patch(isocaps(D,5),... % hcap = patch(...
   'FaceColor','interp',...
   'EdgeColor','none');
disp('4')
set(p,'FaceColor','red','EdgeColor','none');

daspect([1,1,1]); view(3); axis tight; camlight; lighting gouraud;

title ='brain1';

% x3mesh(fv.faces, fv.vertices, 0.3, 'name', 'brain1', 'color', [0.5 0 1], ... 
%     'subheading', 'Click and drag to rotate. Scroll to zoom.')
disp('5')
x3mesh(fv.Faces, fv.Vertices, 0.5, 'name', 'brain1', 'color', [0.5 0 1], ... 
    'subheading', 'Click and drag to rotate. Scroll to zoom.','outdir',pwd)
disp('6')

% Required Inputs,
%   f : Faces of the input mesh
%   v : Vertices of the input mesh
% 
% Optional Inputs,
%   reduction   : factor to reduce size of mesh (default = 0.5) 
%                (useful because large meshes will render poorly and make large html files)
%   name        : File name and title of the html file (default = 'example')
%   subheading  : Additional text that can be added
%   color       : n x 3 vector specifying the RGB color of each vertex 
%                 If 1 x 3 vector then the whole mesh will
%                 have the same color. Values must be between 0 and 1
%                 e.g. [0.5 0 0.5]              
%   rotation    : set to:
%                        0, no rotation
%                        1, rotating mesh


