## Compiling ICAR

Edit the makefile to set the path to your compiled NetCDF and FFTW libraries

Also to set the compiler for your machine if necessary (defaults to gfortran)

    make clean
         # remove build by-products

    make
         # default (relatively high) optimization compile

    make install
        # compile if necessary, then install in the install directory [~/bin]

### Options
    MODE=fast           # more optimization, slower compile, WARNING:not safe optimizations
    MODE=profile        # set profiling options for gnu or intel compilers
    MODE=debug          # debug compile with optimizations
    MODE=debugslow      # debug compile w/o optimizations
    MODE=debugomp       # debug compile with optimizations and OpenMP parallelization
    MODE=debugompslow   # debug compile w/o optimizations but with OpenMP parallelization

    make doc
    # build doxygen documentation in docs/html

    make test
        # compiles various test programs (mpdata_test, fftshift_test, and calendar_test)

    add -jn to parallelize the compile over n processors

### Example
    make install MODE=debug -j4  # uses 4 processes to compile in debug mode


### Derecho System Specifics
Instructions for building ICAR on [Derecho](https://arc.ucar.edu/knowledge_base/74317833), a supercomputer at NCAR.
See [instructions](running/running.md#derecho-system-specifics) for running ICAR on Derecho.

#### Cray Compiler
``` bash
module load cce ncarcompilers cray-mpich netcdf fftw
make -j 4
```

#### GNU Compiler
``` bash
module load gcc ncarcompilers cray-mpich netcdf fftw caf/derecho-2.10.1
COMPILER=gnu make -j 4
```

#### Intel Compilers
  - __Note__: currently the classic Intel compiler is recommended for production runs but testing `ifx` and reporting issues would be useful.
  - __Note__: test `debugslow` mode is required for successful running and compilation, issue is being worked on.

##### Classic
``` bash
module load intel-classic ncarcompilers intel-mpi netcdf-mpi fftw-mpi
COMPILER=intel MODE=debugslow make -j 4
```
##### OneAPI
``` bash
module load intel-oneapi ncarcompilers intel-mpi netcdf-mpi fftw-mpi
COMPILER=intel MODE=debugslow F90=ifx FC=ifx make -j 4
```
