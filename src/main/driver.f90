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

    use iso_fortran_env,    only : output_unit
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use output_interface,   only : output_t
    use time_step,          only : step                               ! Advance the model forward in time
    use initialization,     only : init_model, init_physics
    use timer_interface,    only : timer_t
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use wind,               only : update_winds
    use restart_interface,  only : restart_model
    use icar_constants,     only : kVARS

    use land_surface,               only : lsm_init

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    type(boundary_t):: boundary, add_cond
    ! type(output_t)  :: dataset
    ! type(output_t)  :: surface_dataset
    type(output_t)  :: restart_dataset
    type(output_t)  :: output_dataset
    type(timer_t)   :: initialization_timer, total_timer, input_timer, output_timer, physics_timer
    type(Time_type) :: next_output
    type(time_delta_t) :: small_time_delta

    character(len=1024) :: file_name, restart_file_name
    character(len=49)   :: file_date_format = '(I4,"-",I0.2,"-",I0.2,"_",I0.2,"-",I0.2,"-",I0.2)'
    integer :: i, restart_counter
    integer :: output_vars(18)


    call small_time_delta%set(1)

    call total_timer%start()
    call initialization_timer%start()
    !-----------------------------------------
    !  Model Initialization
    !
    ! Reads config options and initializes domain and boundary conditions
    call init_model(options, domain, boundary, add_cond)  ! added boundary structure for external files (additional conditions)

    if (this_image()==1) write(*,*) "Setting up output files"
    ! should be combined into a single setup_output call
    
    if (this_image()==1 .and. options%parameters%frames_per_outfile<2) then
        print*,"  frames per output file should be 2 or more. Currently: ", options%parameters%frames_per_outfile
    else
        if (this_image()==1) print*,"  frames per output file= ", options%parameters%frames_per_outfile
    end if    
    
    call restart_dataset%set_domain(domain)
    call restart_dataset%add_variables(options%vars_for_restart, domain)

    call output_dataset%set_domain(domain)
    call output_dataset%add_variables(options%output_options%vars_for_output, domain)

    if (options%parameters%restart) then
        if (this_image()==1) write(*,*) "Reading restart data"
        call restart_model(domain, restart_dataset, options)
    endif

    ! physics drivers need to be initialized after restart data are potentially read in.
    call init_physics(options, domain)

    !-----------------------------------------
    !-----------------------------------------
    !  Time Loop
    !
    !   note that a timestep here is a forcing input timestep O(1-3hr), not a physics timestep O(20-100s)
    write(file_name, '(A,I6.6,"_",A,".nc")')      &
            trim(options%output_options%output_file), &
            this_image(),                         &
            trim(domain%model_time%as_string(file_date_format))

    write(restart_file_name, '(A,I6.6,"_",A,".nc")')      &
            trim(options%output_options%restart_file), &
            this_image(),                         &
            trim(domain%model_time%as_string(file_date_format))


    call initialization_timer%stop()

    i=1
    call output_timer%start()
    call restart_dataset%save_file(trim(restart_file_name), i, domain%model_time)
    restart_counter = 0
    call output_dataset%save_file(trim(file_name), i, domain%model_time)

    next_output = domain%model_time + options%output_options%output_dt
    call output_timer%stop()
    i = i + 1

    if (this_image()==1) write(*,*) "Initialization complete, beginning physics integration."

    do while (domain%model_time < options%parameters%end_time)
        call flush(output_unit)

        ! -----------------------------------------------------
        !
        !  Read input data if necessary
        !
        ! -----------------------------------------------------
        call input_timer%start()
        if (boundary%current_time <= (domain%model_time + small_time_delta) ) then
            if (this_image()==1) write(*,*) ""
            if (this_image()==1) write(*,*) " ----------------------------------------------------------------------"
            if (this_image()==1) write(*,*) "Updating Boundary conditions"
            call boundary%update_forcing(options)
            call domain%interpolate_forcing(boundary, update=.True.)
            call update_winds(domain, options)

            ! Make the boundary condition dXdt values into units of [X]/s
            call domain%update_delta_fields(boundary%current_time - domain%model_time)
        endif
        call input_timer%stop()



        ! -----------------------------------------------------
        !
        !  Integrate physics forward in time
        !
        ! -----------------------------------------------------
        if (this_image()==1) write(*,*) "Running Physics"
        if (this_image()==1) write(*,*) "  Model time = ", trim(domain%model_time%as_string())
        if (this_image()==1) write(*,*) "   End  time = ", trim(options%parameters%end_time%as_string())
        if (this_image()==1) write(*,*) "  Next Input = ", trim(boundary%current_time%as_string())
        if (this_image()==1) write(*,*) "  Next Output= ", trim(next_output%as_string())

        ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
        call physics_timer%start()
        call step(domain, step_end(boundary%current_time, next_output), options)
        call physics_timer%stop()


        ! -----------------------------------------------------
        !
        !  Write output data if it is time
        !
        ! -----------------------------------------------------
        ! This is a bit of a hack until the output object is set up better to handle files with specified number of steps per file (or months or...)
        ! ideally this will just become timer_start, save_file()...
        ! the output object needs a pointer to model_time and will have to know how to create output file names
        call output_timer%start()
        if ((domain%model_time + small_time_delta) >= next_output) then
            if (this_image()==1) write(*,*) "Writing output file"
            ! if (i>24) then
            if (i>options%parameters%frames_per_outfile) then
                write(file_name, '(A,I6.6,"_",A,".nc")')    &
                    trim(options%output_options%output_file),   &
                    this_image(),                           &
                    trim(domain%model_time%as_string(file_date_format))
                i = 1
            endif
            call output_dataset%save_file(trim(file_name), i, next_output)

            restart_counter = restart_counter + 1
            if (restart_counter == options%output_options%restart_count) then
                restart_counter = 0

                write(restart_file_name, '(A,I6.6,"_",A,".nc")')      &
                        trim(options%output_options%restart_file), &
                        this_image(),                         &
                        trim(domain%model_time%as_string(file_date_format))

                call restart_dataset%save_file(trim(restart_file_name), 1, next_output)
            endif

            next_output = next_output + options%output_options%output_dt
            i = i + 1
        endif

        call output_timer%stop()

    end do
    !
    !-----------------------------------------
    call total_timer%stop()

    if (this_image()==1) then
        write(*,*) ""
        write(*,*) "Model run from : ",trim(options%parameters%start_time%as_string())
        write(*,*) "           to  : ",trim(options%parameters%end_time%as_string())
        write(*,*) "Domain : ",trim(options%parameters%init_conditions_file)
        write(*,*) "Number of images:",num_images()
        write(*,*) ""
        write(*,*) "First image timing:"
        write(*,*) "total   : ", trim(total_timer%as_string())
        write(*,*) "init    : ", trim(initialization_timer%as_string())
        write(*,*) "input   : ", trim(input_timer%as_string())
        write(*,*) "output  : ", trim(output_timer%as_string())
        write(*,*) "physics : ", trim(physics_timer%as_string())
    endif

contains

    function step_end(time1, time2) result(min_time)
        implicit none
        type(Time_type), intent(in) :: time1
        type(Time_type), intent(in) :: time2
        type(Time_type) :: min_time

        if (time1 <= time2 ) then
            min_time = time1
        else
            min_time = time2
        endif

    end function

end program

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
