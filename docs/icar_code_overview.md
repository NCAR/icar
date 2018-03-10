# ICAR Code Overview
**NOTE** Documentation is slowly being built up within the code using doxygen. After installing doxygen and dot, 
run "make doc" to create the API documentation in docs/html/

[Online ICAR docs](http://ncar.github.io/icar/ "ICAR docs")

This will be updated over time.  For now here is a cursory description to help point people to the right files. 

The main program resides in main/driver.f90

## Initialization
driver.f90 calls init\_model(), init\_physics(), and bc\_init()

### Reading parameters
see: main/init\_options.f90

### Domain Setup
see: main/init.f90 -> init\_model() -> init\_domain()

## Forcing code
see: main/boundary.f90 -> bc\_init(), bc\_update()

### Initial conditions
bc\_init() reads in the initial conditions for the atmosphere. 

### Boundary conditions
bc\_update() reads in the updated boundary conditions for the atmosphere, as well as the internal u,v, and p states. bc\_update() is VERY similar to bc\_init()

## Model time integration
see: main/time\_step.f90 -> step()

## Physics code
The various physics packages are stored in the physics/ directory.  This code is further broken down into microphysics (mp\_), advection (adv\_), planetary boundary layer (pbl\_), radiation (ra\_), land surface (lsm\_), and cumulus (cu\_) code, each of which has a *\_driver.f90 file that contains the routine called from time_step.f90.  In addition, the linear theory "dynamics" are in physics/linear\_wind.f90

## I/O and other utilities
### I/O
NetCDF input and output is handled by various helper routines in the io/ directory.  The io/output.f90 file contains the code to write the primary ICAR output files, while io\_routines.f90 contains various helper routines for both reading and writing files.  

### Utilities
Important utility functions are included for working with calendars (time.f90), geographic/horizontal interpolation (geo\_reader.f90), vertical interpolation (vinterp.f90), string manipulation (string.f90), and fft array manipulation (fftshift.f90).  

## Tests
A few rudimentary tests have been developed to help debug the code as it was being developed, and they are retained in the tests/ directory for future use.  
