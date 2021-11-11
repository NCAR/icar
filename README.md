# The Intermediate Complexity Atmospheric Research Model (ICAR)

[![main](https://github.com/NCAR/icar/actions/workflows/icar-main-commit.yml/badge.svg)](https://github.com/NCAR/icar/actions/workflows/icar-main-commit.yml)
[![develop](https://github.com/NCAR/icar/actions/workflows/icar-develop-commit.yml/badge.svg)](https://github.com/NCAR/icar/actions/workflows/icar-develop-commit.yml)
[![Documentation Status](https://readthedocs.org/projects/icar/badge/)](http://icar.readthedocs.org/en/develop/)

ICAR is a simplified atmospheric model designed primarily for climate downscaling, atmospheric sensitivity tests, and hopefully educational uses. ICAR combines an analytical solution for flow over mountains (linear mountain wave theory) with the large scale flow for a driving model to predict the high resolution wind field. It then advects and heat and moisture through the domain while computing cloud microphysical effects. ICAR has includes a land surface model as well for land atmosphere interactions; ICAR can simulate open water fluxes, PBL mixing, surface radiation, and even parameterized convection.

In ICAR 2.0 (currently early alpha), ICAR supports parallelization across hundreds of computing nodes (the basic physics have been shown to scale up to nearly 100,000 processors) using coarray fortran. This version of the code has a significant overhaul of the original code base, and as a result not all functionality has been restored yet.

Documentation is (slowly) being built on [readthedocs](http://icar.readthedocs.org/en/develop/) and doxygen based documentation can be built now by running "make doc", and is available through [github-pages](http://NCAR.github.io/icar).

#### Requirements
To run the model 3D time-varying atmospheric data are required, though an ideal test case can be generated for simple simulations as well.  See "Running the Model" below. There are some sample python scripts to help make input forcing files, but the WRF pre-processing system can also be used.  Low-resolution WRF output files can be used directly, various reanalysis and GCM output files can be used with minimal pre-processing (just get all the variables in the same netcdf file.)  In addition, a high-resolution netCDF topography file is required.  This will define the grid that ICAR will run on.  Finally and ICAR options file is used to specify various parameters for the model.  A sample options file is provided in the run/ directory.

To run ICAR on more than one compute node requires a fortran compiler that supports the use of coarrays.  This includes ifort >= ~18, gfortran >= ~6.3 (with opencoarrays), and cray's fortran compiler. Note that ifort has often been extremely slow, cray's implementation is excellent but ICAR is not well tested with it, gfortran works very well, but some combinations of gfortran and opencoarrays may not work.

#### Developing
If you plan to make any major additions to ICAR, please get in touch, for minor changes feel free to just submit a pull request. The current workflow is to make changes and pull requests to the `develop` branch.

For an outline of the basic code structure see the [ICAR code overview](docs/icar_code_overview.md)

For reference working with the model code and git, see the [ICAR and Git workflow](docs/howto/icar_and_git_howto.md).

#### Reference
Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016), *The Intermediate Complexity Atmospheric Research Model*, J. Hydrometeor, doi:[10.1175/JHM-D-15-0155.1](http://dx.doi.org/10.1175/JHM-D-15-0155.1).
