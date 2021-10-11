#!/bin/bash

module load tmux
module load openmpi
module load fftw
module load gnu
module load opencoarrays

export OPM_NUM_THREADS=1

make FFTW=$NCAR_ROOT_FFTW CAF_DIR=$NCAR_ROOT_OPENCOARRAYS NETCDF=$NCAR_ROOT_NETCDF F90=caf
