#!/bin/bash

#SBATCH --job-name=hViewer
#SBATCH --time=1-0:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=48G
#SBATCH --output=%j-hViewer.out
#SBATCH --error=%j-hViewer.err

# Uncomment the line below and modify
#cd "/path/to/hydroViewer/script/"

#imgFile: Path to decompressed data generated from read_binary_* script in UIOOPS.
#plotDirectory: Path to save image strips to file.
#campaign: Name of project (e.g., 'olympex'). Allows for project-specific conditional statements to be added to plotting routines.
#probeName: Probe type ('Fast2DC','2DS', 'HVPS', 'CIP', 'PIP'). Used to determine decryption routines specific to the manufacturer.
#annotate: 0 - Do not annotate every 5 particles in record; 1 - Annotate
#chunkSize: Number of image records (frames) to plot in succession
#chunkNum: nth chunk of frames to process for current job (e.g., chunkNum=2 plots frames chunkSize+1 to 2*chunkSize)

# Uncomment the line below and modify
#python hydroViewer.py /path/to/UIOOPS/imageFile.nc /path/to/save/images/ socrates Fast2DC 0 60000 1