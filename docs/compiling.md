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
module load ncarenv cce cray-mpich netcdf fftw
make -j 4
```

#### GNU Compiler
``` bash
module load gcc ncarenv/23.09 cray-mpich netcdf fftw opencoarrays
COMPILER=gnu make -j 4
```


### Derecho Debugging
When running ICAR on larger domains there is a chance the program will run out of memory available for coarrays.
This may be difficult to diagnose because the program will stop without outputting any information but the following steps might help.
*NOTE*: The following debugging information is only for Cray compilers, GNU's OpenCoarrays is not setup to use the SHMEM library.

#### Set Helpful Debug Variables
Set the following variables, further detail and more variables can be found in the [Cray OpenSHMEMX](https://cray-openshmemx.readthedocs.io/en/latest/intro_shmem.html#cray-openshmemx-setup-and-running-specific-environment-variables) documentation.
* `SHMEM_MEMINFO_DISPLAY=1` to display information about the job's memory allocation during initialization
* `SHMEM_ABORT_ON_ERROR=1` and `MPICH_ABORT_ON_ERROR=1`, these are set when loading the ATP module

#### Increase Memory Available for Coarrays
The default symmetric heap size on Derecho is 64MB per process.
To increase the size set the variable `XT_SYMMETRIC_HEAP_SIZE` to an integer value with the suffix `M` for megabyte or `G` for gigabyte.
``` bash
export XT_SYMMETRIC_HEAP_SIZE=128M
```

#### Further Tools
* [Derecho debugging and profiling documentation](https://arc.ucar.edu/knowledge_base/149323810)
* Use [gdb](https://sourceware.org/gdb/current/onlinedocs/gdb) (with [cheatsheet](https://sourceware.org/gdb/current/onlinedocs/gdb)) or [gdb4hpc](https://cpe.ext.hpe.com/docs/debugging-tools/gdb4hpc.1.html)
* [Linaro Forge Tools](https://docs.linaroforge.com/23.0/html/forge/index.html) such as DTT or MAP.
* If the program returned `died from signal XYZ`, check the signal error against list shown by `$ kill -L` or `$ kill -l`

<!-- NOTE: removing Intel compiler information until tested more -->
<!-- #### Intel Compilers -->
<!--   - __Note__: currently the classic Intel compiler is recommended for production runs but testing `ifx` and reporting issues would be useful. -->
<!--   - __Note__: test `debugslow` mode is required for successful running and compilation, issue is being worked on. -->

<!-- ##### Classic -->
<!-- ``` bash -->
<!-- module load intel-classic ncarcompilers intel-mpi netcdf-mpi fftw-mpi -->
<!-- COMPILER=intel MODE=debugslow make -j 4 -->
<!-- ``` -->
<!-- ##### OneAPI -->
<!-- ``` bash -->
<!-- module load intel-oneapi ncarcompilers intel-mpi netcdf-mpi fftw-mpi -->
<!-- COMPILER=intel MODE=debugslow F90=ifx FC=ifx make -j 4 -->
<!-- ``` -->
