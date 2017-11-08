# Start Here
This section is a non-technical description that will describe how to get the software running on your computer. If you encounter problems don't hesitate to report an issue on this page. I'll do my best to try to address problems as soon as possible.

1. Click the bright green button at the top right of this website that says "Clone or Download" and choose Download ZIP to download the master zip file. [This link](https://github.com/atdemarco/svrlsmgui/archive/master.zip) jumps right to that file.

2. Extract/unzip the file you downloaded, so that it now looks like a stand alone folder containing files.

3. Now that the master zip is downloaded and unzipped, make sure MATLAB is installed. If it's not installed, first [download it from MATLAB's website](https://www.mathworks.com/login?uri=https%3A%2F%2Fwww.mathworks.com%2Fdownloads%2Flatest_release) (may have to create a login). Your University may have a campus-wide "Total Academic Headcount" license, you can check on your University's website/IT services or use  MATLAB's [lookup tool](https://www.mathworks.com/academia/tah-support-program/eligibility/index.html).

4. Open MATLAB and [add to your MATLAB search path](https://www.mathworks.com/help/matlab/ref/addpath.html) the directory you extracted in step 2.

5. Add SPM to MATLAB's search path if it is not added.

6. Make sure you either have the MATLAB Machine Learning Toolbox or libSVM installed and visible to MATLAB's search path.

7. Open the svrlsmgui graphical interface by typing 'svrlsmgui' into the white MATLAB terminal and pressing enter or double-clicking the svrlsmgui.fig file contained in the directory you unzipped in step 2.

# What is SVR-LSM GUI?
A graphic multivariate lesion-symptom mapping toolbox.

Requirements:
*[The SPM12 toolbox](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
*MATLAB (note to self: add min version)
*MATLAB Parallel computing toolbox (for parallelization functionality)
*MATLAB Machine Learning Toolbox (for MATLAB's SVR functionality)
*libSVM (for libSVM SVR functionality, as in original Zhang et al., 2014)
