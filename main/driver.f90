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
    use time,               only : calendar_date, date_to_mjd         ! Convert between date and modified Julian Day
    use init,               only : init_model, init_physics           ! Initialize model (not initial conditions)
    use boundary_conditions,only : bc_init,bc_update,bc_find_step     ! Boundary and initial conditions
    use data_structures          ! *_type datatypes                   ! Data-types and physical "constants"
    use output,             only : write_domain, output_init          ! Used to output initial model state
    use time_step,          only : step                               ! Advance the model forward in time
    use string,             only : str                                ! Convert real,integer,double to string
    
    implicit none
    type(options_type) :: options
    type(domain_type)  :: domain
    type(bc_type)      :: boundary
    integer            :: i,nx,ny,start_point
    integer            :: year, month, day, hour, minute, second
    double precision   :: model_time,next_output
        
!-----------------------------------------
!  Model Initialization
!
!   initialize model including options, terrain, lat, lon data. 
    call init_model(options,domain,boundary)
    
    ! set up the timeing for the model
    if (options%restart) then
        start_point=options%restart_step-1
    else
        start_point=bc_find_step(options)
    endif
    
    model_time=start_point * DBLE(options%in_dt) + options%time_zero
    domain%model_time=model_time
    next_output=model_time+options%out_dt
    
    call calendar_date(model_time/86400.0D+0 + 50000, year, month, day, hour, minute, second)
    domain%current_month=month
    
    ! read initial conditions from the boundary file
    write(*,*) "Initializing Boundary conditions"
    call bc_init(domain, boundary, options)
    
    write(*,*) "Initializing Physics packages"
    call init_physics(options,domain)
    
    ! initialize the output module
    ! this can't be called until after bc_init, so that rain accumulations can be initialized in a restart
    call output_init(domain, options)
    
    ! write the initial state of the model (primarily useful for debugging)
    if (.not.options%restart) then
        call write_domain(domain,options,nint((model_time-options%time_zero)/options%out_dt))
    endif
!
!-----------------------------------------
!-----------------------------------------
!  Time Loop
!
!   note that a timestep here is a forcing input timestep O(1-3hr), not a physics timestep O(20-100s)
    do i=start_point,options%ntimesteps
        write(*,*) ""
        write(*,*) " ----------------------------------------------------------------------"
        write(*,*) "Timestep:", i-options%time_step_zero, "  of ", options%ntimesteps-options%time_step_zero
        write(*,*) "  Model time=", trim(str((model_time-options%time_zero)/3600.0,fmt="(F10.2)")) ,"hrs"
        call calendar_date(model_time/86400.0D+0 + 50000, year, month, day, hour, minute, second)
        domain%current_month=month
        write(*,'(A,i4,"/",i2.2"/"i2.2" "i2.2":"i2.2":"i2.2)') "  Date = ",year,month,day,hour,minute,second
        
        ! update boundary conditions (dXdt variables) so we can integrate to the next step
        call bc_update(domain,boundary,options)
        
        ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
        call step(domain,options,boundary,model_time,next_output)

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
