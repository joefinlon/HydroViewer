# HydroViewer
Python package for displaying particle images from optical array probe (OAP) data processed by the University of Illinois/Oklahoma OAP Processing Software (UIOOPS; https://github.com/joefinlon/UIOPS).

Please see the \*.sh files in the repository for examples on how to call the desired programs. Below is a brief description of the programs in this package. If you have a bug to report, send an email to jfinlon@uw.edu.

## Python Dependencies
-   Numpy (https://numpy.org/)
-   Xarray (http://xarray.pydata.org/en/stable/)
-   PIL (https://pillow.readthedocs.io/en/stable/)

## Supported Instruments
-   2DS/HVPS
-   Fast-2DC
-   CIP/PIP

## Features

### SelectHydroViewer
For help on running SelectHydroViewer in the background, follow the selectHydroViewer\_example.sh script.
-	Builds image buffers from particles meeting user-specified criteria (e.g., time, size, habit, and others)
-	Images outputted as pixel-for-pixel representation to preserve quality and reduce file size

### HydroViewer
For help on running HydroViewer in the background, follow the hydroViewer\_example.sh and hydroViewer\_example\_parallel.sh scripts.
-	Displays 6 image records at once
-	Automated saving of images for all image records in a file
-	Information on the frame number, time, and number of particles for each image record
-	Support for parallel plotting