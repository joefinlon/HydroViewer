# HydroViewer
Python package for displaying OAP images and particle properties.

HydroViewer is a Python package useful in displaying particle images from a variety of optical array probes (OAPs). The most stable version of the code will be on the GitHub repository, but the developer cannot guarantee that bugs do not exist. If you have a bug to report, send an email to finlon2@illinois.edu. The code will be updated over time with new user features, with such additions listed in the What’s New and Version History sections below.

Also included in this document is a brief overview of the features and how to run the software using files processed through the University of Illinois OAP Processing Software (UIOPS). The package includes 2 programs to display particle images: (1) HydroViewer for plotting to file 6 image records simultaneously in the background, and (2) HydroViewer2 for interactively plotting an image record and displaying particle properties for user-selected hydrometeors.

## What's New

### HydroViewer
Version 1.0

### HydroViewer2
- Added support for DMT probes (e.g., CIP, PIP)

## Features

### HydroViewer
-	Displays 6 image records at once
-	Automated saving of images for all image records in a file
-	Information on the frame number, time, and number of particles for each image record
-	Option to overlay an asterisk every 5 particles (helpful when associating images to the particle-by-particle data

### HydroViewer2
-	Displays one image record at a time with buttons to toggle through the file
-	Particle properties (e.g., maximum dimension, area, perimeter) shown after clicking on a particle within the image record
-	Ability to save images of user-selected particles

## HydroViewer Instructions

HydroViewer is designed to process 60,000 image records at a time. Processing more than this can have adverse effects such as excessive memory utilization. Because of this, it’s recommended to first check the number of frames in the image file and determine how many calls to the HydroViewer program are required for the flight.

The program has the following Python package dependencies: numpy, matplotlib, xarray, glob

HydroViewer is called using the following arguments:
-	campaign: string variable (e.g., ‘olympex’)
-	date: string variable in yyyymmdd format (e.g., ‘20151112’)
-	probeName: string variable [options – ‘2DS’, ‘HVPS’]
-	imageMode: 0 - all particles shaded black
-	annotate: 0 - no annotations ; 1 - a red asterisk is overlaid every 5 particles
- chunkNum: nth chunk (of 60,000 image records) to process within the file
-	inFile: file path to the decompressed image file (ouput data from the read_binary_*.m script)

The program is initiated with the following commands:
>> from hydroViewer import get_slice_endpoints, buffer_integrity, image_buffer, annotate_particle_incriments, initialize_inputs;
>> initialize_inputs(campaign, date, probeName, imageMode, annotate, chunkNum, inFile)

If all 60,000 frames are to be processed, it is recommended to run HydroViewer with 32 GB of memory allocated to one processor.

## HydroViewer2 Instructions

HydroViewer2 is designed to run within Jupyter notebook for the purpose of displaying particles from a single image record and allow the user to view particle properties for a user-selected particle.

The program has the following Python package dependencies: numpy, scipy, matplotlib, xarray, glob

HydroViewer2 is run by following the steps below:
1.	Run the first two cell blocks in the program (package import block and subroutines block)
2.	In the third cell block, modify the path to the particle-by-particle data file (partFile; output data from the imgProc_sm.m script) and modify the probe type (probeName). Run this cell.
3.	In the fourth cell block, modify the campaign (string variable; e.g., ‘olypex’), date (string variable; e.g., ‘20151112’), imageMode (0 - all particles shaded black), imageFile (file path to the decompressed image file; ouput data from the read_binary_*.m script), and frameStart (frame number to jump to within the image file) variables. Run this cell.
4.	Run the fifth cell block as well. Matplotlib should display an interactive plot which enables the capability to toggle between image records and register user clicks over a particle of interest. If desired: press the ‘Previous’ and ‘Next’ buttons to backtrack or advance one frame at a time through the data.
5.	With the desired frame displayed and the cell block highlighted, hover your cursor over a particle of interest and click within that region. The program records the position within the plot for later use. You may also notice that the cursor position in the lower right portion of the figure is displayed as you navigate around the image record.
6.	Run the sixth cell block. After 5-10 seconds, particle properties for the selected particle will be displayed immediately below the cell.
7.	In the seventh cell block, change the directory path as desired for the outFile variable. Separate images are saved to the specified directory depending on the probe used and the particle’s time, frame number, and particle number. Run this cell. The image (*.png) should appear in the specified folder. The pixel dimensions will be nDiodes x nSlices so that each pixel corresponds to a single pixel shadowed by the OAP.

## Known Bugs

### HydroViewer
-	Frame number in image filenames and frame numbers displayed above the image records do not directly correspond to the image record in the particle-by-particle data files (i.e., frame #1 corresponds to the second record number in the particle-by-particle data file)

### HydroViewer2
- When browsing through CIP/PIP image data using the previous/next buttons, the image record occasionally does not update --clicking the button again refreshes the display for an additional frame backward/forward

## Version History

### HydroViewer
Version 1.0
-	Initial code release

### HydroViewer2
Version 2.1
- Added support for DMT probes (e.g., CIP, PIP)

Version 2.0
- Initial code release
