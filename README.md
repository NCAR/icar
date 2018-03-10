# The Intermediate Complexity Atmospheric Research Model (ICAR)

[![Build Status](https://travis-ci.org/NCAR/icar.svg)](https://travis-ci.org/NCAR/icar)
[![Documentation Status](https://readthedocs.org/projects/icar/badge/)](http://icar.readthedocs.org/en/develop/)

ICAR is a simplified atmospheric model designed primarily for climate downscaling, atmospheric sensitivity tests, and hopefully educational uses. At this early stage, the model is still undergoing rapid development, and users are encouraged to get updates frequently. 

Documentation is (slowly) being built on [readthedocs](http://icar.readthedocs.org/en/develop/) and doxygen based documentation can be built now by running "make doc", and is available through [github-pages](http://NCAR.github.io/icar). 

#### Requirements
To run the model 3D time-varying atmospheric data are required, though an ideal test case can be generated for simple simulations as well.  See "Running the Model" below. There are some sample python scripts to help make input forcing files, but the WRF pre-processing system can also be used.  Low-resolution WRF output files can be used directly, various reanalysis and GCM output files can be used with minimal pre-processing (just get all the variables in the same netcdf file.)  In addition, a high-resolution netCDF topography file is required.  This will define the grid that ICAR will run on.  Finally and ICAR options file is used to specify various parameters for the model.  A sample options file is provided in the run/ directory. 

#### Developing
For an outline of the basic code structure see the [ICAR code overview](docs/icar_code_overview.md)

For reference working with the model code and git, see the [ICAR and Git workflow](docs/howto/icar_and_git_howto.md). 

#### Reference
Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016), *The Intermediate Complexity Atmospheric Research Model*, J. Hydrometeor, doi:[10.1175/JHM-D-15-0155.1](http://dx.doi.org/10.1175/JHM-D-15-0155.1).
