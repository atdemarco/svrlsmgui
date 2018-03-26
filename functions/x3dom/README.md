# matlab-to-html-plot

A simple matlab interface to plot meshes as interactive html plots

x3mesh is a small and simple script for converting a Matlab mesh consisting of faces and vertices into an x3dom object in an html file.

This allows 3d objects to be displayed and interacted with on the web (rotate and zoom).

Run demo1.m and demo2.m to try it out. The script produces a .html files which can be opened in a browser.

An example of the output can be found here:

http://www.birving.com/other/Example4.html

The script takes advantage of the x3dom web format. More details can be found here:

http://www.x3dom.org/

Not supported by internet explorer. Use a recent Firefox, Chrome or Safari browser.

Update:

1) Added better optional argument parser

2) Added option to set vertex color

3) Added option so set mesh to auto-rotate 
