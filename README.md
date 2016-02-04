#The Intermediate Complexity Atmospheric Research Model (ICAR)

ICAR is a simplified atmospheric model designed primarily for climate downscaling, atmospheric sensitivity tests, and hopefully educational uses. At this early stage, the model is still undergoing rapid development, and users are encouraged to get updates frequently. 

####Requirements
To run the model 3D time-varying atmospheric data are required, though an ideal test case can be generated for simple simulations as well.  See "Running the Model" below. There are some sample python scripts to help make input forcing files, but the WRF pre-processing system can also be used.  Low-resolution WRF output files can be used directly, various reanalysis and GCM output files can be used with minimal pre-processing (just get all the variables in the same netcdf file.)  In addition, a high-resolution netCDF topography file is required.  This will define the grid that ICAR will run on.  Finally and ICAR options file is used to specify various parameters for the model.  A sample options file is provided in the run/ directory. 

####Developing
For an outline of the basic code structure see the [ICAR code overview](https://github.com/NCAR/icar/blob/develop/docs/icar_code_overview.md)

For reference working with the model code and git, see the [ICAR and Git workflow](https://github.com/NCAR/icar/blob/develop/docs/howto/icar_and_git_howto.md). 

####Reference
Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016), *The Intermediate Complexity Atmospheric Research Model*, J. Hydrometeor, accepted.

---------------------------------------------------------------------------------
##Compiling the model: 
    
Edit the makefile to set the path to your compiled NetCDF and FFTW libraries
    
Also to set the compiler for your machine if necessary (defaults to gfortran)
    
    make clean
        remove build by-products
        
    make
        default (relatively high) optimization compile 
    
    make install
        compile if necessary, then install in the install directory [~/bin]
    
###Options: 
    MODE=fast           more optimization, slower compile, WARNING:not safe optimizations
    MODE=profile        set profiling options for gnu compiler in particular
    MODE=debug          debug compile with optimizations
    MODE=debugslow      debug compile w/o optimizations
    MODE=debugomp       debug compile with optimizations and OpenMP parallelization
    MODE=debugompslow   debug compile w/o optimizations but with OpenMP parallelization

    make test
        compiles various test programs (mpdata_test, fftshift_test, and calendar_test)
    
    add -jn to parallelize the compile over n processors 
    
###Example:
    make MODE=debug -j4  # uses 4 processes to compile in debug mode

---------------------------------------------------------------------------------
##Running the model:

###Example

$ icar options_file_name.nml

###Input files

####Primary settings file:
icar\_options.nml [other filename can be specified on the commandline e.g. icar your\_options\_file]
        This file specifies all other files and options.  See the example in run/icar\_options.nml

Other settings files can be defined in the primary settings file, or all settings can be included in one file
    
####Necessary netcdf files: 

1) Boundary / Initial conditions file(s) (from e.g. wrf, reanalysis, or GCM output).  These files are specified in the options file as boundary\_files. 

Must contain the following variables (optional variables are in square brackets) :
NOTE: any variable name can be used as specified in the options file. 

    QV      = Water Vapor mixing ratio                  (kg/kg)
    T       = Air Temperature                           (K with an optional offset)
    P       = Pressure                                  (Pa with an optional [PB] offset)
    U       = East-West wind                            (m/s)
    V       = North-South wind                          (m/s)
    HGT     = Terrain Height                            (m)
    Z       = 3D model level heights                    (m)
    LAT     = Latitude on mass (P/T/etc.) grid          (degrees)
    LONG    = Longitude on mass (P/T/etc.) grid         (degrees)
    [PB]    = base pressure to be added to P            (Pa)
    [QC]    = cloud water content mixing ratio          (kg/kg)
    [QI]    = cloud ice content mixing ratio            (kg/kg)
    [LATU]  = Latitude on the EW-wind grid              (degrees)
    [LONGU] = Longitude on the EW-wind grid             (degrees)
    [LATV]  = Latitude on the NS-wind grid              (degrees)
    [LONGV] = Longitude on the NS-wind grid             (degrees)
    [SST]   = Sea Surface Temperature                   (K)
    [SWD]   = Shortwave down at the surface             (W/m^2)
    [LWD]   = Longwave down at the surface              (W/m^2)
    [SH]    = Sensible heat flux from the surface       (W/m^2)
    [LH]    = Latent heat flux from the surface         (W/m^2)
    [PBLH]  = Specified height of PBL                   (m)

2) High-resolution file (all variables are on the high-resolution grid ICAR will run on).  This filename is specified in the options file as the (poorly named ) init\_conditions\_file. 

Must [optionally] contain:

    HGT     = Terrain Height                            (m)
    LAT     = Latitude on ICAR mass grid                (degrees)
    LONG    = Longitude on ICAR mass grid               (degrees)
    [LATU]  = Latitude on the ICAR EW-staggered wind grid         (degrees)
    [LONGU] = Longitude on the ICAR EW-staggered wind grid        (degrees)
    [LATV]  = Latitude on the ICAR NS-staggered wind grid         (degrees)
    [LONGV] = Longitude on the ICAR NS-staggered wind grid        (degrees)
    [LU]    = Land use cover classification for land surface model
    [SOIL]  = Soil type classification for land surface model

Must be in a (nearly) constant dx,dy projection (e.g. Lambert Conformal Conic)
            
**NB**: Prior to ICAR 1.0 (and we're not there yet), the namelist file in git repo is not always in sync with code! This *usually* just means that it may not include all current options, it will usually run. Develop branch more likely to be out of sync. 

---------------------------------------------------------------------------------
###Common errors:

1) **Segmentation fault**:
    Most likely due to your shell's stacksize limit (particularly with ifort). 
    
To *fix* it try the following: 

    in bash:
        ulimit -s unlimited
    in csh: 
        unlimit stacksize
        
**NB**: Some systems have a hard limit of 64MB (ulimit -s 65532), this *may* be enough depending on domain size. 
Reference: https://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors
