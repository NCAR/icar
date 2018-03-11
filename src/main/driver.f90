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
    use wind,               only : update_winds

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    type(boundary_t):: boundary
    type(output_t)  :: dataset

    character(len=1024) :: file_name
    integer :: i
    integer :: omp_get_max_threads

    !-----------------------------------------
    !  Model Initialization
    !
    if (this_image()==1) then
        print*, "Number of coarray image:",num_images()
        print*, "Max number of OpenMP Threads:",omp_get_max_threads()
        print*, "Reading options"
    endif
    call options%init()

    if (this_image()==1) print*, "Initializing Domain"
    call domain%init(options)

    if (this_image()==1) print*, "Initializing boundary condition data structure"
    call boundary%init(options)

    if (this_image()==1) print*, "Reading Initial conditions from boundary dataset"
    call domain%get_initial_conditions(boundary, options)

    call update_winds(domain, options)

    if (this_image()==1) print*, "Setting up output files"
    ! should be combined into a single setup_output call
    call dataset%set_domain(domain)
    call dataset%add_variables(options%vars_for_restart, domain)

    !-----------------------------------------
    !-----------------------------------------
    !  Time Loop
    !
    !   note that a timestep here is a forcing input timestep O(1-3hr), not a physics timestep O(20-100s)
    do while (domain%model_time < options%parameters%end_time)

        call boundary%update_forcing(options)
        call domain%interpolate_forcing(boundary, update=.True.)

        if (this_image()==1) write(*,*) ""
        if (this_image()==1) write(*,*) " ----------------------------------------------------------------------"
        if (this_image()==1) write(*,*) "  Model time = ", trim(domain%model_time%as_string())
        if (this_image()==1) write(*,*) "   End  time = ", trim(options%parameters%end_time%as_string())
        if (this_image()==1) write(*,*) "  Input time = ", trim(boundary%current_time%as_string())


        ! this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
        call step(domain, boundary%current_time, options)


        ! This is an ugly hack until the output object is set up better to handle multiple time steps per file
        ! (that may just need "unlimited" specified in variables?)
        if (this_image()==1) print*, "Writing output file"
        write(file_name, '("icar_restart_output_",I3.3,"_",A,".nc")') this_image(), trim(domain%model_time%as_string())

        do i=1,len_trim(file_name)
            if (file_name(i:i)==" ") file_name = file_name(:i-1)//"_"//file_name(i+1:)
            if (file_name(i:i)=="/") file_name = file_name(:i-1)//"-"//file_name(i+1:)
            if (file_name(i:i)==":") file_name = file_name(:i-1)//"-"//file_name(i+1:)
        end do

        call dataset%save_file("output/"//trim(file_name))

    end do
    !
    !-----------------------------------------

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
