!>-----------------------------------------
!! Main Program
!!
!! Initialize options and memory in init_model
!! Read initial conditions in bc_init (from a restart file if requested)
!! initialize physics packages in init_physics (e.g. tiedke and thompson if used)
!! If this run is a restart run, then set start to the restart timestep
!!      in otherwords, ntimesteps is the number of BC updates from the beginning of the entire model
!!      run, not just from the begining of this restart run
!! calculate model time in seconds based on the time between BC updates (in_dt)
!! Calculate the next model output time from current model time + output time delta (out_dt)
!!
!! Finally, loop until ntimesteps are reached updating boundary conditions and stepping the model forward
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!-----------------------------------------
program icar
    use time_obj,           only : Time_type
    use options_obj,        only : options_t
    use domain_obj,         only : domain_t
    use forcing_obj,        only : forcing_t
    use output_obj,         only : options_t

    implicit none

    type(options_t)     :: options
    type(domain_t)      :: domain
    type(forcing_t)     :: boundary
    type(output_t)      :: output
    type(output_t)      :: restart_state
    type(Time_type)     :: next_output

    call options%read_options()
    !-----------------------------------------
    !  Model Initialization
    !
    !   initialize model including options, terrain, lat, lon data.
    write(*,*) "Initializing Domain"
    call domain%init(options)

    j=

    ! read initial conditions from the boundary file
    write(*,*) "Initializing Boundary conditions"
    call boundary%init(options, domain)

    next_output= domain%model_time + options%output_dt

    write(*,*) "Initializing Physics packages"
    call domain%init_physics(options, domain)

    ! initialize the output module
    ! this can't be called until after bc_init, so that rain accumulations can be initialized in a restart
    call output%init(domain, options)

    ! write the initial state of the model (primarily useful for debugging)
    if (.not.options%restart) then
        call output%write(options, domain%model_time)
    endif

    !-----------------------------------------
    !-----------------------------------------
    !  Time Loop
    !
    !   note that a timestep here is a forcing input timestep O(1-3hr), not a physics timestep O(20-100s)
    do while (domain%model_time < options%end_time)
        write(*,*) ""
        write(*,*) " ----------------------------------------------------------------------"
        write(*,*) "  Model time = ", trim(domain%model_time%as_string())
        write(*,*) "   End  time = ", trim(options%end_time%as_string())

        ! update boundary conditions (dXdt variables) so we can integrate to the next step
        call boundary%update(domain, options)
        write(*,*) "  Next input = ", trim(boundary%next_domain%model_time%as_string())

        ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
        call domain%integrate(domain, options, boundary, next_output)

    end do
!
!-----------------------------------------

end program icar

! This is the Doxygen mainpage documentation.  This should be moved to another file at some point.

!>------------------------------------------
!!  @mainpage
!!
!!  @section Introduction
!!  ICAR is a simplified atmospheric model designed primarily for climate downscaling, atmospheric sensitivity tests,
!!  and hopefully educational uses. At this early stage, the model is still undergoing rapid development, and users
!!  are encouraged to get updates frequently.
!!
!!  @section Running_ICAR
!!  To run the model 3D time-varying atmospheric data are required, though an ideal test case can be generated for
!!  simple simulations as well. There are some sample python scripts to help make input forcing files, but the WRF
!!  pre-processing system can also be used. Low-resolution WRF output files can be used directly, various reanalysis
!!  and GCM output files can be used with minimal pre-processing (just get all the variables in the same netcdf file.)
!!  In addition, a high-resolution netCDF topography file is required. This will define the grid that ICAR will run on.
!!  Finally an ICAR options file is used to specify various parameters for the model. A sample options file is provided
!!  in the run/ directory.
!!
!!  @section Developing
!!  This document provides the primary API and code structure documentation. The code is based on github.com/NCAR/icar
!!  Developers are encouraged to fork the main git repository and maintain their own git repository from which to
!!  issue pull requests.
!!
!!  @section Reference
!!  Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016),
!!  The Intermediate Complexity Atmospheric Research Model, J. Hydrometeor, doi:<a href="http://dx.doi.org/10.1175/JHM-D-15-0155.1">10.1175/JHM-D-15-0155.1</a>.
!!
!!------------------------------------------
