# ICAR Code Overview
**NOTE** Documentation is slowly being built up within the code using doxygen. After installing doxygen and dot,

run `make doc` to create the API documentation in `docs/html/`

[Online ICAR docs](http://ncar.github.io/icar/ "ICAR docs")

This will be updated over time.  For now here is a cursory description to help point people to the right files.

The main program resides in `main/driver.f90`

## Initialization
`driver.f90` calls `init_model()`, `init_physics()`, and `bc_init()`

### Reading parameters
See: `main/init_options.f90`

### Domain Setup
See: `main/init.f90` -> `init_model()` -> `init_domain()`

## Forcing code
See: `main/boundary.f90` -> `bc_init()`, `bc_update()`

### Initial conditions
`bc_init()` reads in the initial conditions for the atmosphere.

### Boundary conditions
`bc_update()` reads in the updated boundary conditions for the atmosphere, as well as the internal u,v, and p states. `bc_update()` is VERY similar to `bc_init()`

## Model time integration
See: `main/time_step.f90` -> `step()`

## Physics code
The various physics packages are stored in the `physics/` directory.  This code is further broken down into microphysics (`mp_`), advection (`adv_`), planetary boundary layer (`pbl_`), radiation (`ra_`), land surface (`lsm_`), and cumulus (`cu_`) code, each of which has a `*_driver.f90` file that contains the routine called from `time_step.f90`.  In addition, the linear theory "dynamics" are in `physics/linear\_wind.f90`

## I/O and other utilities
### I/O
NetCDF input and output is handled by various helper routines in the `io/` directory.  The `io/output.f90` file contains the code to write the primary ICAR output files, while `io_routines.f90` contains various helper routines for both reading and writing files.

### Utilities
Important utility functions are included for working with calendars (`time.f90`), geographic/horizontal interpolation (`geo_reader.f90`), vertical interpolation (`vinterp.f90`), string manipulation (`string.f90`), and fft array manipulation (`fftshift.f90`).

## Tests
A few rudimentary tests have been developed to help debug the code as it was being developed, and they are retained in the `tests/` directory for future use.
