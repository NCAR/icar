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
    If there are existing outputfiles with a different grid (e.g. number of levels, number of latitudes), ICAR will try to write to the existing file and will stop because the dimensions don't match. 
    
3) **dz_levels namelist error**
    If compiled with gfortran, the namelist format for arrays has to be slightly different than is supplied in the run directory. An easy fix is to put the dz_levels array all onto one line. 
    
    
        
