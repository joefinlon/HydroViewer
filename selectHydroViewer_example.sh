#!/bin/tcsh

#SBATCH --job-name=selHyVwr
#SBATCH --time=1-0:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=64G
#SBATCH --output=%j-selHyVwr.out
#SBATCH --error=%j-selHyVwr.err

set RunDir="/path/to/HydroViewer/directory/"
cd "/path/to/HydroViewer/directory/"

# EXAMPLE COMMAND TO PLOT ALL PARTICLES
#python selectHydroViewer.py image_filename particle_filename plot_directory campaign probe start_time end_time reject_criteria min_size max_size habit_criteria

# EXAMPLE COMMAND TO PLOT ACCEPTED DENDRITES BTWN 1830 & 1835 UTC FROM OLYMPEX 2DS W/ 0.75<D<1.00 MM & INT-ARR>1E-6 S
python selectHydroViewer.py /path/imagefile.nc /path/particlefile.nc /path/plot_dir/ olympex 2DS 183000 183500 0.75 1. 1e-6 100

# HABIT CODES:
#116 (tiny), 111 (oriented), 108 (linear), 97 (aggregate), 103 (graupel), 115 (sphere), 104 (hexagonal), 105 (irregular), 100 (dendrite)