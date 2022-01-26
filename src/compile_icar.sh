#!/bin/bash
module load gnu
module load tmux
module load openmpi
module load opencoarrays
module load fftw
export OMP_NUM_THREADS=1
#make MODE=debugompslow
#make FFTW=$NCAR_ROOT_FFTW CAF_DIR=$NCAR_ROOT_OPENCOARRAYS NETCDF=$NCAR_ROOT_NETCDF F90=caf MODE=debugompslow
make FFTW=$NCAR_ROOT_FFTW CAF_DIR=$NCAR_ROOT_OPENCOARRAYS NETCDF=$NCAR_ROOT_NETCDF F90=caf 
