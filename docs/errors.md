## Common errors

1) **Segmentation fault**:
    Possibly due to your shell's stacksize limit (particularly with ifort).

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

5) **Floating Point errors**
    Check your input data.  For example, if you are supplying Shortwave or longwave down at the surface, but those terms are 0, the LSM can cool off and the surface layer becomes too stable. This causes the surface fluxes to become numerically unstable, and eventually the system breaks.

6) **NetCDF related build error**
    If the NetCDF library supplied does not support the full NetCDF4 file format (based on HDF5) you will get the following error.  Note that using version 4 of the netcdf library does not imply full support for version 4 of the netcdf file format...

    `Error: Symbol ‘nf90_netcdf4’ at (1) has no IMPLICIT type`

    There are two options :
        1) The better option is to recompile the netcdf library, being sure that you first compile HDF5 and link the netcdf library to it.  You should be able to run
            `nc-config --has-nc4`
            and have it output “yes”

        2) You can edit the ICAR source code to not require netCDF4 support.  This might be easier, but it will mean that ICAR ca not save very large files (which is useful for the linear wind look up table).  To do this, simply change `or(NF90_CLOBBER,NF90_NETCDF4)` to `NF90_CLOBBER` in `io/lt_lut_io.f90`
