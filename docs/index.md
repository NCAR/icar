#The Intermediate Complexity Atmospheric Research Model (ICAR)

ICAR is a simplified atmospheric model designed primarily for climate downscaling, atmospheric sensitivity tests, and hopefully educational uses. ICAR combines the wind field from a lower resolution atmospheric model, such as a climate model, with an analytical solution for the effect of topography on that wind field to produce a high-resolution wind field. This wind field is then used to advect heat, moisture, and hydrometeor species (e.g. clouds) around the three dimensional domain, and combines this with microphysical calculations.  

Version 2 has been a major development and now requires that your fortran compiler supports coarrays (intel, gnu, cray).  Of note, cray may not work because it doesn't seem to properly support some uses of generic classes. Intel works on a single node, but is very slow across multiple nodes.  gfortran works well, but requires the [OpenCoarrays](https://github.com/sourceryinstitute/OpenCoarrays) library be installed.

##Requirements
To run the model 3D time-varying atmospheric data are required, though an ideal test case can be generated for simple simulations as well.  See "Running the Model" below. There are some sample python scripts to help make input forcing files, but the WRF pre-processing system can also be used.  Low-resolution WRF output files can be used directly, various reanalysis and GCM output files can be used with minimal pre-processing (just get all the variables in the same netcdf file.)  In addition, a high-resolution netCDF topography file is required.  This will define the grid that ICAR will run on.  Finally an ICAR options file is used to specify various parameters for the model.  A sample options file is provided in the run/ directory.

##Developing
For an outline of the basic code structure see the [ICAR code overview](icar_code_overview.md)

For reference working with the model code and git, see the [ICAR and Git workflow](howto/icar_and_git_howto.md).

You can see more documentation on the code by installing doxygen and running "make doc".  Documentation will be generated in docs/html

##References
Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016), *The Intermediate Complexity Atmospheric Research Model*, J. Hydrometeor, e-View http://dx.doi.org/10.1175/JHM-D-15-0155.1.

Rouson, D., Gutmann, E. D., Fanfarillo, A., Friesen, B. (2017) *Performance portability of an intermediate-complexity atmospheric research model in coarray Fortran*, Proceedings of the Second Annual PGAS Applications Workshop, 1-4
