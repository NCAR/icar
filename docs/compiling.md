##Compiling ICAR
    
Edit the makefile to set the path to your compiled NetCDF and FFTW libraries
    
Also to set the compiler for your machine if necessary (defaults to gfortran)
    
    make clean
         # remove build by-products
        
    make
         # default (relatively high) optimization compile 
    
    make install
        # compile if necessary, then install in the install directory [~/bin]
    
###Options: 
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
    
###Example:
    make install MODE=debug -j4  # uses 4 processes to compile in debug mode
