##Common errors

1) **Segmentation fault**:
    Most likely due to your shell's stacksize limit (particularly with ifort). 
    
To *fix* it try the following: 

    in bash:
        ulimit -s unlimited
    in csh: 
        unlimit stacksize

    **NB**: Some systems have a hard limit of 64MB (ulimit -s 65532), this *may* be enough depending on domain size. 
    Reference: https://software.intel.com/en-us/articles/determining-root-cause-of-sigsegv-or-sigbus-errors

2) **Common output file errors**
    If there are existing outputfiles with a different grid (e.g. number of levels, number of latitudes), ICAR will try to write to the existing file and will stop because the dimensions don't match. If you get an error such as : 
    
    NetCDF: Variable not found
    setup_varids: Searching for variable in existing file:nsq
    
Then ICAR is trying to write a variable (in this case nsq, the Brunt-Vaisala frequency) to the output file, but this variable does not exist is pre-existing output files.  If you get an error such as : 

    NetCDF: Start+count exceeds dimension bound
    output/icar_1990_01_01_00-00.nc:qv

Then the dimensions of the model domain have changed since the last run, and it is not possible to write the output to the existing files. 

In both cases, the solution is to delete (or move) the existing output files, or to modify your options file so that the output will be consistent with existing files.  

    
3) **dz_levels namelist error**
    "Fortran runtime error: Bad data for namelist object dz_levels"

    If compiled with gfortran, the namelist format for arrays has to be slightly different than is supplied in the run directory. An easy fix is to put the dz_levels array all onto one line. 
    
    
4) **other namelist error**
    If a newer namelist is used with an older version of code, you may get errors telling you that a given variable is not supported.  You can probably remove that line from the namelist and it will run correctly, though you should think about the variable that is being removed to decide what it means.  
