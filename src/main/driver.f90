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

    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use output_interface,   only : output_t
    use time_step,          only : step                               ! Advance the model forward in time
    use initialization,     only : init_model
    use timer_interface,    only : timer_t

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    type(boundary_t):: boundary
    type(output_t)  :: dataset
    type(timer_t)   :: initialization_timer, total_timer, input_timer, output_timer, physics_timer

    character(len=1024) :: file_name
    integer :: i

    call total_timer%start()
    call initialization_timer%start()
    !-----------------------------------------
    !  Model Initialization
    !
    ! Reads config options and initializes domain and boundary conditions
    call init_model(options, domain, boundary)

    if (this_image()==1) print*, "Setting up output files"
    ! should be combined into a single setup_output call
    call dataset%set_domain(domain)
    call dataset%add_variables(options%vars_for_restart, domain)

    !-----------------------------------------
    !-----------------------------------------
    !  Time Loop
    !
    !   note that a timestep here is a forcing input timestep O(1-3hr), not a physics timestep O(20-100s)
    write(file_name, '(I4.4,"_",A,".nc")') this_image(), trim(domain%model_time%as_string())

    do i=1,len_trim(file_name)
        if (file_name(i:i)==" ") file_name = file_name(:i-1)//"_"//file_name(i+1:)
        if (file_name(i:i)=="/") file_name = file_name(:i-1)//"-"//file_name(i+1:)
        if (file_name(i:i)==":") file_name = file_name(:i-1)//"-"//file_name(i+1:)
    end do
    file_name = trim(options%parameters%output_file)//trim(file_name)
    i=1

    call initialization_timer%stop()

    call output_timer%start()
    call dataset%save_file(trim(file_name), i, domain%model_time)
    call output_timer%stop()
    i = i + 1

    do while (domain%model_time < options%parameters%end_time)

        if (this_image()==1) write(*,*) ""
        if (this_image()==1) write(*,*) " ----------------------------------------------------------------------"
        if (this_image()==1) print*, "Updating Boundary conditions"
        call input_timer%start()
        call boundary%update_forcing(options)
        call domain%interpolate_forcing(boundary, update=.True.)
        call input_timer%stop()

        if (this_image()==1) print*, "Running Physics"
        if (this_image()==1) write(*,*) "  Model time = ", trim(domain%model_time%as_string())
        if (this_image()==1) write(*,*) "   End  time = ", trim(options%parameters%end_time%as_string())
        if (this_image()==1) write(*,*) "  Input time = ", trim(boundary%current_time%as_string())


        ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
        call physics_timer%start()
        call step(domain, boundary%current_time, options)
        call physics_timer%stop()

        ! This is an ugly hack until the output object is set up better to handle multiple time steps per file
        ! (that may just need "unlimited" specified in variables?)
        if (this_image()==1) print*, "Writing output file"
        if (i>24) then
            write(file_name, '(I4.4,"_",A,".nc")') this_image(), trim(domain%model_time%as_string())

            do i=1,len_trim(file_name)
                if (file_name(i:i)==" ") file_name = file_name(:i-1)//"_"//file_name(i+1:)
                if (file_name(i:i)=="/") file_name = file_name(:i-1)//"-"//file_name(i+1:)
                if (file_name(i:i)==":") file_name = file_name(:i-1)//"-"//file_name(i+1:)
            end do
            file_name = trim(options%parameters%output_file)//trim(file_name)
            i = 1
        endif

        if (this_image()==1) print*, trim(file_name)
        call output_timer%start()
        call dataset%save_file(trim(file_name), i, domain%model_time)
        call output_timer%stop()
        i = i + 1

    end do
    !
    !-----------------------------------------
    call total_timer%stop()

    if (this_image()==1) then
        print*, ""
        print*, "Model run from : ",trim(options%parameters%start_time%as_string())
        print*, "           to  : ",trim(options%parameters%end_time%as_string())
        print*, "Domain : ",trim(options%parameters%init_conditions_file)
        print*, "Number of images:",num_images()
        print*, ""
        print*, "First image timing:"
        print*, "total   : ", trim(total_timer%as_string())
        print*, "init    : ", trim(initialization_timer%as_string())
        print*, "input   : ", trim(input_timer%as_string())
        print*, "output  : ", trim(output_timer%as_string())
        print*, "physics : ", trim(physics_timer%as_string())
    endif

    sync all

    if (this_image()==num_images()) then
        print*, ""
        print*, "Last image timing:"
        print*, "total   : ", trim(total_timer%as_string())
        print*, "init    : ", trim(initialization_timer%as_string())
        print*, "input   : ", trim(input_timer%as_string())
        print*, "output  : ", trim(output_timer%as_string())
        print*, "physics : ", trim(physics_timer%as_string())
    endif

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
