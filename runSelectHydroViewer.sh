#!/bin/tcsh

#SBATCH --job-name=selHyVwr
#SBATCH --time=1-0:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=64G
#SBATCH --output=%j-selHyVwr.out
#SBATCH --error=%j-selHyVwr.err

set RunDir="/path/to/SelectHydroViewer/directory/"
cd "/path/to/SelectHydroViewer/directory/"

# EXAMPLE COMMAND TO PLOT ALL PARTICLES (PLEASE ALTER AS THIS WILL TAKE A LONG TIME)
#python -u -c 'from selectHydroViewer import probe_defaults,load_partData,get_partInds,get_imageData,get_slice_endpoints,image_buffer,hydro_viewer; hydro_viewer(imageFile="/path/to/image/data/fname.cdf",particleFile="/path/to/particle/data/fname.cdf",plotDirectory="/path/to/store/images/",campaign="",probeName="")'

# EXAMPLE COMMAND TO PLOT ACCEPTED DENDRITES BTWN 1830 & 1835 UTC FROM OLYMPEX 2DS W/ 0.75<D<1.00 MM & INT-ARR>1E-6 S
python -u -c 'from selectHydroViewer import probe_defaults,load_partData,get_partInds,get_imageData,get_slice_endpoints,image_buffer,hydro_viewer; hydro_viewer(imageFile="/path/to/image/data/fname.cdf",particleFile="/path/to/particle/data/fname.cdf",plotDirectory="/path/to/store/images/",campaign="olympex",probeName="2DS",startTime="183000",endTime="183500",rejStatus=[48, 104, 72, 117],minD=0.75,maxD=1.0,intArrThresh=1e-6,habitStatus=[100])'
