submodule(options_interface) options_implementation

    use icar_constants,             only : kMAINTAIN_LON, MAXFILELENGTH, MAXVARLENGTH, MAX_NUMBER_FILES, MAXLEVELS, &
                                           kNO_STOCHASTIC, kVERSION_STRING, kMAX_FILE_LENGTH, kMAX_NAME_LENGTH, pi, &
                                           kWATER_LAKE, &
                                           kWIND_LINEAR, kLINEAR_ITERATIVE_WINDS, kITERATIVE_WINDS, kCONSERVE_MASS
    use io_routines,                only : io_newunit, check_writeable_path
    use time_io,                    only : find_timestep_in_file
    use time_delta_object,          only : time_delta_t
    use time_object,                only : Time_type
    use string,                     only : str
    use model_tracking,             only : print_model_diffs

    use convection,                 only : cu_var_request
    use land_surface,               only : lsm_var_request
    use radiation,                  only : ra_var_request
    use microphysics,               only : mp_var_request
    use advection,                  only : adv_var_request
    use wind,                       only : wind_var_request
    use planetary_boundary_layer,   only : pbl_var_request

    use output_metadata,            only : get_varname

    implicit none


contains


    !> ----------------------------------------------------------------------------
    !!  Read all namelists from the options file specified on the command line
    !!
    !!  Reads the commandline (or uses default icar_options.nml filename)
    !!  Checks that the version of the options file matches the version of the code
    !!  Reads each namelist successively, all options are stored in supplied options object
    !!
    !! ----------------------------------------------------------------------------

    !> -------------------------------
    !! Initialize an options object
    !!
    !! Allocated coarray options types, reads namelists on image 1, and distributes data to all images
    !!
    !! -------------------------------
    module subroutine init(this)
        implicit none
        class(options_t),   intent(inout)  :: this

!       reads a series of options from a namelist file and stores them in the
!       options data structure
        character(len=MAXFILELENGTH) :: options_filename
        integer :: i

        options_filename = get_options_file()
        if (this_image()==1) write(*,*) "Using options file = ", trim(options_filename)

        call version_check(         options_filename,   this%parameters)
        call physics_namelist(      options_filename,   this)
        call var_namelist(          options_filename,   this%parameters)
        call parameters_namelist(   options_filename,   this%parameters)
        call output_namelist(       options_filename,   this%output_options)
        call model_levels_namelist( options_filename,   this%parameters)

        call lt_parameters_namelist(    this%parameters%lt_options_filename,    this)
        call block_parameters_namelist( this%parameters%block_options_filename, this)
        call mp_parameters_namelist(    this%parameters%mp_options_filename,    this)
        call adv_parameters_namelist(   this%parameters%adv_options_filename,   this)
        call lsm_parameters_namelist(   this%parameters%lsm_options_filename,   this)
        call cu_parameters_namelist(    this%parameters%cu_options_filename,    this)
        call bias_parameters_namelist(  this%parameters%bias_options_filename,  this)
        call rad_parameters_namelist(   this%parameters%rad_options_filename,   this)

        if (this%parameters%restart) then
            ! if (this_image()==1) write(*,*) "  (opt) Restart = ", this%parameters%restart
            call init_restart_options(options_filename, this%parameters)
            this%parameters%start_time = this%parameters%restart_time
        endif

        call filename_namelist(options_filename, this%parameters)

        ! check for any inconsistencies in the options requested
        call options_check(this)

        call collect_physics_requests(this)

    end subroutine init

    !> -------------------------------
    !! Call all physics driver var_request routines
    !!
    !! var_request routines allow physics modules to requested
    !! which variables they need to have allocated, advected, and written in restart files
    !!
    !! -------------------------------
    subroutine collect_physics_requests(options)
        type(options_t) :: options

        call ra_var_request(options)
        call lsm_var_request(options)
        call pbl_var_request(options)
        call cu_var_request(options)
        call mp_var_request(options)
        call adv_var_request(options)
        call wind_var_request(options)
        call pbl_var_request(options)

    end subroutine

    !> -------------------------------
    !! Add list of new variables to a list of variables
    !!
    !! Adds one to the associated index of the list and returns an error
    !! Sets Error/=0 if any of the variables suggested are outside the bounds of the list
    !!
    !! -------------------------------
    subroutine add_to_varlist(varlist, varids, error)
        implicit none
        integer, intent(inout)  :: varlist(:)
        integer, intent(in)     :: varids(:)
        integer, intent(out), optional  :: error

        integer :: i, ierr

        ierr=0
        do i=1,size(varids)
            if (varids(i) <= size(varlist)) then
                varlist( varids(i) ) = varlist( varids(i) ) + 1
            else
                if (this_image()==1) write(*,*) "WARNING: trying to add var outside of permitted list:",varids(i), size(varlist)
                ierr=1
            endif
        enddo

        if (present(error)) error=ierr

    end subroutine add_to_varlist


    !> -------------------------------
    !! Add a set of variable(s) to the internal list of variables to be allocated
    !!
    !! Sets error /= 0 if an error occurs in add_to_varlist
    !!
    !! -------------------------------
    module subroutine alloc_vars(this, input_vars, var_idx, error)
        class(options_t),  intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error

        integer :: ierr

        ierr=0
        if (present(var_idx)) then
            call add_to_varlist(this%vars_to_allocate,[var_idx], ierr)
        endif

        if (present(input_vars)) then
            call add_to_varlist(this%vars_to_allocate,input_vars, ierr)
        endif

        if (present(error)) error=ierr

    end subroutine alloc_vars


    !> -------------------------------
    !! Add a set of variable(s) to the internal list of variables to be output in a restart file
    !!
    !! Sets error /= 0 if an error occurs in add_to_varlist
    !!
    !! -------------------------------
    module subroutine restart_vars(this, input_vars, var_idx, error)
        class(options_t),  intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error

        integer :: ierr

        ierr=0
        if (present(var_idx)) then
            call add_to_varlist(this%vars_for_restart,[var_idx], ierr)
        endif

        if (present(input_vars)) then
            call add_to_varlist(this%vars_for_restart,input_vars, ierr)
        endif

        if (present(error)) error=ierr

    end subroutine


    !> -------------------------------
    !! Add a set of variable(s) to the internal list of variables to be advected
    !!
    !! Sets error /= 0 if an error occurs in add_to_varlist
    !!
    !! -------------------------------
    module subroutine advect_vars(this, input_vars, var_idx, error)
        class(options_t),  intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error

        integer :: ierr

        ierr=0
        if (present(var_idx)) then
            call add_to_varlist(this%vars_to_advect,[var_idx], ierr)
        endif

        if (present(input_vars)) then
            call add_to_varlist(this%vars_to_advect,input_vars, ierr)
        endif

        if (present(error)) error=ierr

    end subroutine



    !> ----------------------------------------------------------------------------
    !!  Read in the name of the options files from the command line
    !!
    !!  @retval     options_file    The name of the options parameter file to use.
    !!                              Default = icar_options.nml
    !!
    !! ----------------------------------------------------------------------------
    function get_options_file() result(options_file)
        implicit none
        character(len=MAXFILELENGTH) :: options_file

        ! Internal variables
        integer :: error
        logical :: file_exists
        ! default options filename
        character(len=*), parameter :: default_options_file = "icar_options.nml"

        ! if a commandline argument was supplied, read the options filename from there
        if (command_argument_count()>0) then
            ! read the commandline argument
            call get_command_argument(1, options_file, status=error)
            ! if there was any problem revert to the default filename
            if (error > 0) then
                options_file = default_options_file

            ! error -1 means the filename supplied was too long
            elseif (error == -1) then
                if (this_image()==1) write(*,*) "Options filename = ", trim(options_file), " ...<cutoff>"
                if (this_image()==1) write(*,*) "Maximum filename length = ", MAXFILELENGTH
                stop "ERROR: options filename too long"
            endif

        ! If not arguments were supplied use the default filename
        else
            options_file = default_options_file
        endif

        ! Check that the options file actually exists
        INQUIRE(file=trim(options_file), exist=file_exists)

        ! if options file does not exist, print an error and quit
        if (.not.file_exists) then
            if (this_image()==1) write(*,*) "Using options file = ", trim(options_file)
            stop "Options file does not exist. "
        endif
    end function



    !> -------------------------------
    !! Check the version number in the namelist file and compare to the current model version
    !!
    !! If the namelist version doesn't match, print the differences between that version and this
    !! and STOP execution
    !!
    !! -------------------------------
    subroutine version_check(filename,options)
        character(len=*),            intent(in)     :: filename
        type(parameter_options_type),intent(inout)  :: options

        character(len=MAXVARLENGTH) :: version,comment
        integer                     :: name_unit

        namelist /model_version/ version,comment


        !default comment:
        comment="Model testing"

        ! read namelists
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=model_version)
        close(name_unit)

        if (version.ne.kVERSION_STRING) then
            if (this_image()==1) write(*,*) "Model version does not match namelist version"
            if (this_image()==1) write(*,*) "  Model version: ",kVERSION_STRING
            if (this_image()==1) write(*,*) "  Namelist version: ",trim(version)
            call print_model_diffs(version)
            stop
        endif
        options%version = version
        options%comment = comment

        if (this_image()==1) write(*,*) "  Model version: ",trim(version)

    end subroutine version_check

    !> -------------------------------
    !! Checks options in the options data structure for consistency
    !!
    !! Stops or prints a large warning depending on warning level requested and error found
    !!
    !! -------------------------------
    subroutine options_check(options)
        ! Minimal error checking on option settings
        implicit none
        type(options_t), intent(inout)::options

        if (options%parameters%t_offset.eq.(-9999)) then
            ! if (options%parameters%warning_level>0) then
            !     if (this_image()==1) write(*,*) "WARNING, WARNING, WARNING"
            !     if (this_image()==1) write(*,*) "WARNING, Using default t_offset=0"
            !     if (this_image()==1) write(*,*) "WARNING, WARNING, WARNING"
            ! endif
            options%parameters%t_offset = 0
        endif


        ! convection can modify wind field, and ideal doesn't rebalance winds every timestep
        if ((options%physics%convection.ne.0).and.(options%parameters%ideal)) then
            if (options%parameters%warning_level>3) then
                if (this_image()==1) write(*,*) ""
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                if (this_image()==1) write(*,*) "WARNING, Running convection in ideal mode may be bad..."
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%parameters%warning_level==10) then
                if (this_image()==1) write(*,*) "Set warning_level<10 to continue"
                stop
            endif
        endif

        ! wind calculations almost require fixed_dz_advection settings
        if ((options%physics%windtype.eq.kITERATIVE_WINDS).and.(.not.options%parameters%fixed_dz_advection)) then
            if (options%parameters%warning_level>3) then
                if (this_image()==1) write(*,*) ""
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                if (this_image()==1) write(*,*) "WARNING, wind=3 setting is best used with fixed_dz_advection=.True."
                if (this_image()==1) write(*,*) "WARNING, setting fixed_dz_advection=.True."
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                options%parameters%fixed_dz_advection = .True.
            endif
            if (options%parameters%warning_level==10) then
                if (this_image()==1) write(*,*) "Set warning_level<10 to continue"
                stop
            endif
        endif
        if (((options%physics%windtype.eq.kWIND_LINEAR).or.(options%physics%windtype.eq.kLINEAR_ITERATIVE_WINDS)).and.(options%parameters%fixed_dz_advection)) then
            if (options%parameters%warning_level>3) then
                if (this_image()==1) write(*,*) ""
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                if (this_image()==1) write(*,*) "WARNING, wind=1 or 5 setting is best used with fixed_dz_advection=.False."
                if (this_image()==1) write(*,*) "WARNING, setting fixed_dz_advection=.False."
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                if (this_image()==1) write(*,*) ""
                options%parameters%fixed_dz_advection = .False.
            endif
            if (options%parameters%warning_level==10) then
                if (this_image()==1) write(*,*) "Set warning_level<10 to continue"
                stop
            endif
        endif
        if ((options%physics%windtype.eq.0).and.(options%parameters%fixed_dz_advection)) then
            if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            if (this_image()==1) write(*,*) "WARNING setting fixed_dz_advection=False for wind=0"
            if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            options%parameters%fixed_dz_advection = .False.
        endif
        if ((options%physics%windtype.eq.kCONSERVE_MASS).and.(.not.options%parameters%fixed_dz_advection)) then
            if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            if (this_image()==1) write(*,*) "WARNING setting fixed_dz_advection=True for wind=2"
            if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            options%parameters%fixed_dz_advection = .True.
        endif

        ! if using a real LSM, feedback will probably keep hot-air from getting even hotter, so not likely a problem
        if ((options%physics%landsurface>1).and.(options%physics%boundarylayer==0)) then
            if (options%parameters%warning_level>2) then
                if (this_image()==1) write(*,*) ""
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                if (this_image()==1) write(*,*) "WARNING, Running LSM without PBL may overheat the surface and CRASH the model. "
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%parameters%warning_level>=7) then
                if (this_image()==1) write(*,*) "Set warning_level<7 to continue"
                stop
            endif
        endif
        ! if using perscribed LSM fluxes, no feedbacks are present, so the surface layer is likely to overheat.
        if ((options%physics%landsurface==1).and.(options%physics%boundarylayer==0)) then
            if (options%parameters%warning_level>0) then
                if (this_image()==1) write(*,*) ""
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                if (this_image()==1) write(*,*) "WARNING, Prescribed LSM fluxes without a PBL may overheat the surface and CRASH. "
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%parameters%warning_level>=5) then
                if (this_image()==1) write(*,*) "Set warning_level<5 to continue"
                stop
            endif
        endif

        ! prior to v 0.9.3 this was assumed, so throw a warning now just in case.
        if ((options%parameters%z_is_geopotential .eqv. .False.).and. &
            (options%parameters%zvar=="PH")) then
            if (options%parameters%warning_level>1) then
                if (this_image()==1) write(*,*) ""
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                if (this_image()==1) write(*,*) "WARNING z variable is not assumed to be geopotential height when it is 'PH'."
                if (this_image()==1) write(*,*) "WARNING If z is geopotential, set z_is_geopotential=True in the namelist."
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%parameters%warning_level>=7) then
                if (this_image()==1) write(*,*) "Set warning_level<7 to continue"
                stop
            endif
        endif
        if ((options%parameters%z_is_geopotential .eqv. .True.).and. &
            (options%parameters%z_is_on_interface .eqv. .False.)) then
            if (options%parameters%warning_level>1) then
                if (this_image()==1) write(*,*) ""
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
                if (this_image()==1) write(*,*) "WARNING geopotential height is no longer assumed to be on model interface levels."
                if (this_image()==1) write(*,*) "WARNING To interpolate geopotential, set z_is_on_interface=True in the namelist. "
                if (this_image()==1) write(*,*) "WARNING WARNING WARNING"
            endif
        endif

    end subroutine options_check


    !> ------------------
    !!  Determine the filename to be used for this particular image/process based on the restart filename, restart time, and image number
    !!
    !!  Uses the same calculation that is used to get the output filename when writing.
    !! -------------------
    function get_image_filename(image_number, initial_filename, restart_time) result(file_name)
        implicit none
        integer,            intent(in) :: image_number
        character(len=*),   intent(in) :: initial_filename
        type(Time_type),    intent(in) :: restart_time

        character(len=kMAX_FILE_LENGTH) :: file_name
        integer :: n, i

        character(len=49)   :: file_date_format = '(I4,"-",I0.2,"-",I0.2,"_",I0.2,"-",I0.2,"-",I0.2)'

        write(file_name, '(A,I6.6,"_",A,".nc")') trim(initial_filename), image_number, trim(restart_time%as_string(file_date_format))

        n = len(trim(file_name))
        file_name(n-10:n-9) = "00"

    end function get_image_filename


    !> -------------------------------
    !! Initialize the restart options
    !!
    !! Reads the restart namelist if this is a restart run
    !!
    !! -------------------------------
    subroutine init_restart_options(filename, options)
        ! initialize the restart specifications
        ! read in the namelist, and calculate the restart_step if appropriate
        character(len=*),               intent(in)    :: filename    ! name of the file containing the restart namelist
        type(parameter_options_type),   intent(inout) :: options     ! options data structure to store output

        ! Local variables
        character(len=MAXFILELENGTH) :: restart_file    ! file name to read restart data from
        integer :: restart_step                         ! time step relative to the start of the restart file
        integer :: restart_date(6)                      ! date to restart
        integer :: name_unit                            ! logical unit number for namelist
        type(Time_type) :: restart_time, time_at_step   ! restart date as a modified julian day
        real :: input_steps_per_day                     ! number of input time steps per day (for calculating restart_time)

        namelist /restart_info/ restart_step, restart_file, restart_date

        restart_date=[-999,-999,-999,-999,-999,-999]
        restart_step=-999

        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=restart_info)
        close(name_unit)

        if (minval(restart_date)<0) then
            if (this_image()==1) write(*,*) "------ Invalid restart_date ERROR ------"
            stop "restart_date must be specified in the namelist"
        endif

        ! calculate the modified julian day for th restart date given
        call restart_time%init(options%calendar)
        call restart_time%set(restart_date(1), restart_date(2), restart_date(3), &
                              restart_date(4), restart_date(5), restart_date(6))

        ! find the time step that most closely matches the requested restart time (<=)
        restart_file = get_image_filename(this_image(), restart_file, restart_time)
        restart_step = find_timestep_in_file(restart_file, 'time', restart_time, time_at_step)

        ! check to see if we actually updated the restart date and print if in a more verbose mode
        if (options%debug) then
            if (restart_time /= time_at_step) then
                if (this_image()==1) write(*,*) " updated restart date: ", trim(time_at_step%as_string())
            endif
        endif

        restart_time = time_at_step

        if (options%debug) then
            if (this_image()==1) write(*,*) " ------------------ "
            if (this_image()==1) write(*,*) "RESTART INFORMATION"
            if (this_image()==1) write(*,*) "mjd",         restart_time%mjd()
            if (this_image()==1) write(*,*) "date:",       trim(restart_time%as_string())
            if (this_image()==1) write(*,*) "date",        restart_date
            if (this_image()==1) write(*,*) "file",   trim(restart_file)
            if (this_image()==1) write(*,*) "forcing step",restart_step
            if (this_image()==1) write(*,*) " ------------------ "
        endif

        ! save the parameters in the master options structure
        options%restart_step_in_file = restart_step
        options%restart_file         = restart_file
        options%restart_date         = restart_date
        options%restart_time         = restart_time

        if (options%debug) then
            if (this_image()==1) write(*,*) " step in restart file",options%restart_step_in_file
        endif

    end subroutine init_restart_options


    !> -------------------------------
    !! Read physics options to use from a namelist file
    !!
    !! -------------------------------
    subroutine physics_namelist(filename,options)
        implicit none
        character(len=*),intent(in)     :: filename
        type(options_t), intent(inout)  :: options

        integer :: name_unit
!       variables to be used in the namelist
        integer :: pbl, lsm, water, mp, rad, conv, adv, wind

!       define the namelist
        namelist /physics/ pbl, lsm, water, mp, rad, conv, adv, wind

!       default values for physics options (advection+linear winds+simple_microphysics)
        pbl = 0 ! 0 = no PBL,
                ! 1 = Stupid PBL (only in LSM=1 module),
                ! 2 = local PBL diffusion after Louis 1979
                ! 3 = YSU PBL (not complete)

        lsm = 0 ! 0 = no LSM,
                ! 1 = Fluxes from GCM,
                ! 2 = simple LSM, (not complete)
                ! 3 = Noah LSM

        water =0! 0 = no open water fluxes,
                ! 1 = Fluxes from GCM, (needs lsm=1)
                ! 2 = Simple fluxes (needs SST in forcing data)
                ! 3 = WRF's lake model (needs lake depth in hi-res data)

        mp  = 1 ! 0 = no MP,
                ! 1 = Thompson et al (2008),
                ! 2 = "Linear" microphysics
                ! 3 = Morrison
                ! 4 = WSM5

        rad = 0 ! 0 = no RAD,
                ! 1 = Surface fluxes from GCM, (radiative cooling ~1K/day in LSM=1 module),
                ! 2 = cloud fraction based radiation + radiative cooling
                ! 3 = RRTMG

        conv= 0 ! 0 = no CONV,
                ! 1 = Tiedtke scheme
                ! 2 = Kain-Fritsch scheme

        adv = 1 ! 0 = no ADV,
                ! 1 = upwind advection scheme
                ! 2 = MPDATA

        wind= 1 ! 0 = no LT,
                ! 1 = linear theory wind perturbations
                ! 2 = terrain induced horizontal accelleration

!       read the namelist
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=physics)
        close(name_unit)

!       store options
        options%physics%boundarylayer = pbl
        options%physics%convection    = conv
        options%physics%advection     = adv
        options%physics%landsurface   = lsm
        options%physics%watersurface  = water
        options%physics%microphysics  = mp
        options%physics%radiation     = rad
        options%physics%windtype      = wind

    end subroutine physics_namelist


    !> -------------------------------
    !! Check that a required input variable is present
    !!
    !! If not present, halt the program
    !!
    !! -------------------------------
    subroutine require_var(inputvar, var_name)
        implicit none
        character(len=*), intent(in) :: inputvar
        character(len=*), intent(in) :: var_name

        if (trim(inputvar)=="") then
            if (this_image()==1) write(*,*) "Variable: ",trim(var_name), " is required."
            stop
        endif

    end subroutine require_var

    !> -------------------------------
    !! Initialize the variable names to be written to standard output
    !!
    !! Reads the output_list namelist
    !!
    !! -------------------------------
    subroutine output_namelist(filename, options)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(output_options_type), intent(inout) :: options

        integer :: name_unit, i, j, status
        real    :: outputinterval, restartinterval

        character(len=kMAX_FILE_LENGTH) :: output_file, restart_file
        character(len=kMAX_NAME_LENGTH) :: names(kMAX_STORAGE_VARS)
        character (len=MAXFILELENGTH) :: output_file_frequency

        namelist /output_list/ names, outputinterval, restartinterval, &
                               output_file, restart_file, output_file_frequency

        output_file         = "icar_out_"
        restart_file        = "icar_rst_"
        names(:)            = ""
        outputinterval      =  3600
        restartinterval     =  24 ! in units of outputintervals

        open(io_newunit(name_unit), file=filename)
        read(name_unit, nml=output_list)
        close(name_unit)

        do j=1, kMAX_STORAGE_VARS
            if (trim(names(j)) /= "") then
                do i=1, kMAX_STORAGE_VARS
                    if (trim(get_varname(i)) == trim(names(j))) then
                        call add_to_varlist(options%vars_for_output, [i])
                    endif
                enddo
            endif
        enddo

        options%out_dt     = outputinterval
        call options%output_dt%set(seconds=outputinterval)

        if (trim(output_file_frequency) /= "") then
            options%output_file_frequency = output_file_frequency
        else
            ! if outputing at half-day or longer intervals, create monthly files
            if (outputinterval>=43200) then
                options%output_file_frequency="monthly"
            ! if outputing at half-hour or longer intervals, create daily files
            else if (outputinterval>=1800) then
                options%output_file_frequency="daily"
            ! otherwise create a new output file every timestep
            else
                options%output_file_frequency="every step"
            endif
        endif

        options%rst_dt = outputinterval * restartinterval
        call options%restart_dt%set(seconds=options%rst_dt)

        if (restartinterval<0) then
            options%restart_count = restartinterval
        else
            options%restart_count = max(24, nint(restartinterval))
        endif

        ! check if directory paths to output and restart file strings exist
        ! stop and report error if they do not rather than failing at write step
        call check_writeable_path(output_file)
        call check_writeable_path(restart_file)
        options%output_file = output_file
        options%restart_file = restart_file

    end subroutine output_namelist


    !> -------------------------------
    !! Initialize the variable names to be read
    !!
    !! Reads the var_list namelist
    !!
    !! -------------------------------
    subroutine var_namelist(filename,options)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(parameter_options_type), intent(inout) :: options
        integer :: name_unit, i, j
        character(len=MAXVARLENGTH) :: landvar,lakedepthvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar,zbvar,  &
                                        hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,           &
                                        pvar,pbvar,tvar,qvvar,qcvar,qivar,qrvar,qgvar,qsvar,hgtvar,shvar,lhvar,pblhvar,   &
                                        psvar, pslvar, snowh_var, &
                                        soiltype_var, soil_t_var,soil_vwc_var,swe_var, soil_deept_var,           &
                                        vegtype_var,vegfrac_var, vegfracmax_var, albedo_var, lai_var, canwat_var, linear_mask_var, nsq_calibration_var,  &
                                        swdown_var, lwdown_var, sst_var, rain_var, time_var, sinalpha_var, cosalpha_var, &
                                        lat_ext, lon_ext, swe_ext, hsnow_ext, rho_snow_ext, tss_ext, tsoil2D_ext, tsoil3D_ext, z_ext, time_ext

        namelist /var_list/ pvar,pbvar,tvar,qvvar,qcvar,qivar,qrvar,qgvar,qsvar,hgtvar,shvar,lhvar,pblhvar,   &
                            landvar,lakedepthvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar,zbvar, &
                            psvar, pslvar, snowh_var, &
                            hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,           &
                            soiltype_var, soil_t_var,soil_vwc_var,swe_var,soil_deept_var,           &
                            vegtype_var,vegfrac_var, vegfracmax_var, albedo_var, lai_var, canwat_var, linear_mask_var, nsq_calibration_var,  &
                            swdown_var, lwdown_var, sst_var, rain_var, time_var, sinalpha_var, cosalpha_var, &
                            lat_ext, lon_ext, swe_ext, hsnow_ext, rho_snow_ext, tss_ext, tsoil2D_ext, tsoil3D_ext,  z_ext, time_ext

        ! no default values supplied for variable names
        hgtvar=""
        latvar=""
        lonvar=""
        time_var=""
        uvar=""
        ulat=""
        ulon=""
        vvar=""
        vlat=""
        vlon=""
        pslvar=""
        psvar=""
        pvar=""
        pbvar=""
        tvar=""
        qvvar=""
        qcvar=""
        qivar=""
        qrvar=""
        qsvar=""
        qgvar=""
        zvar=""
        zbvar=""
        shvar=""
        lhvar=""
        swdown_var=""
        lwdown_var=""
        sst_var=""
        pblhvar=""
        hgt_hi=""
        landvar=""
        lakedepthvar=""
        lat_hi=""
        lon_hi=""
        ulat_hi=""
        ulon_hi=""
        vlat_hi=""
        vlon_hi=""
        soiltype_var=""
        soil_t_var=""
        soil_vwc_var=""
        swe_var=""
        snowh_var=""
        soil_deept_var=""
        vegtype_var=""
        vegfrac_var=""
        vegfracmax_var=""
        albedo_var=""
        lai_var=""
        canwat_var=""
        linear_mask_var=""
        nsq_calibration_var=""
        rain_var=""
        sinalpha_var=""
        cosalpha_var=""
        lat_ext=""
        lon_ext=""
        swe_ext=""
        hsnow_ext=""
        rho_snow_ext=""
        tss_ext=""
        tsoil2D_ext=""
        tsoil3D_ext=""
        z_ext = ""
        time_ext = ""

        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=var_list)
        close(name_unit)

        call require_var(lonvar, "Longitude")
        call require_var(latvar, "Latitude")
        call require_var(lat_hi, "High-res Lat")
        call require_var(lon_hi, "High-res Lon")
        call require_var(hgt_hi, "High-res HGT")
        call require_var(time_var, "Time")

        options%compute_p = .False.
        if ((pvar=="") .and. ((pslvar/="") .or. (psvar/=""))) options%compute_p = .True.
        if (options%compute_p) then
            if ((pslvar == "").and.(hgtvar == "")) then
                write(*,*) "ERROR: if surface pressure is used to compute air pressure, then surface height must be specified"
                error stop
            endif
        endif

        options%compute_z = .False.
        if ((zvar=="") .and. ((pslvar/="") .or. (psvar/=""))) options%compute_z = .True.
        if (options%compute_z) then
            if (pvar=="") then
                if (this_image()==1) write(*,*) "ERROR: either pressure (pvar) or atmospheric level height (zvar) must be specified"
                error stop
            endif
        endif


        options%vars_to_read(:) = ""
        i=1
        ! 2D geometry variable names (for coarse model)
        options%hgtvar      = hgtvar    ; options%vars_to_read(i) = hgtvar;     options%dim_list(i) = 2;    i = i + 1
        options%latvar      = latvar    ! these variables are read explicitly so not added to vars_to_read list
        options%lonvar      = lonvar    ! these variables are read explicitly so not added to vars_to_read list
        options%time_var    = time_var  ! these variables are read explicitly so not added to vars_to_read list

        ! U varname and associated lat/lon var names
        options%uvar        = uvar      ; options%vars_to_read(i) = uvar;       options%dim_list(i) = 3;    i = i + 1
        if (ulat=="") ulat=latvar
        if (ulon=="") ulon=lonvar
        options%ulat        = ulat      !; options%vars_to_read(i) = ulat;       options%dim_list(i) = 2;    i = i + 1
        options%ulon        = ulon      !; options%vars_to_read(i) = ulon;       options%dim_list(i) = 2;    i = i + 1

        ! V varname and associated lat/lon var names
        options%vvar        = vvar      ; options%vars_to_read(i) = vvar;       options%dim_list(i) = 3;    i = i + 1
        if (vlat=="") vlat=latvar
        if (vlon=="") vlon=lonvar
        options%vlat        = vlat      !; options%vars_to_read(i) = vlat;       options%dim_list(i) = 2;    i = i + 1
        options%vlon        = vlon      !; options%vars_to_read(i) = vlon;       options%dim_list(i) = 2;    i = i + 1

        ! Primary model variable names
        options%pbvar       = pbvar     ; options%vars_to_read(i) = pbvar;      options%dim_list(i) = 3;    i = i + 1
        if (options%compute_p) then
            pvar = "air_pressure_computed"
            options%pvar        = pvar  ; options%vars_to_read(i) = pvar;       options%dim_list(i) = -3;   i = i + 1
        else
            options%pvar        = pvar  ; options%vars_to_read(i) = pvar;       options%dim_list(i) = 3;    i = i + 1
        endif
        options%psvar       = psvar     ; options%vars_to_read(i) = psvar;      options%dim_list(i) = 2;    i = i + 1
        options%pslvar      = pslvar    ; options%vars_to_read(i) = pslvar;     options%dim_list(i) = 2;    i = i + 1
        options%tvar        = tvar      ; options%vars_to_read(i) = tvar;       options%dim_list(i) = 3;    i = i + 1
        options%qvvar       = qvvar     ; options%vars_to_read(i) = qvvar;      options%dim_list(i) = 3;    i = i + 1
        options%qcvar       = qcvar     ; options%vars_to_read(i) = qcvar;      options%dim_list(i) = 3;    i = i + 1
        options%qivar       = qivar     ; options%vars_to_read(i) = qivar;      options%dim_list(i) = 3;    i = i + 1
        options%qrvar       = qrvar     ; options%vars_to_read(i) = qrvar;      options%dim_list(i) = 3;    i = i + 1
        options%qsvar       = qsvar     ; options%vars_to_read(i) = qsvar;      options%dim_list(i) = 3;    i = i + 1
        options%qgvar       = qgvar     ; options%vars_to_read(i) = qgvar;      options%dim_list(i) = 3;    i = i + 1

        ! vertical coordinate
        ! if (options%time_varying_z) then
        if (options%compute_z) then
            zvar = "height_computed"
            options%vars_to_read(i) = zvar;      options%dim_list(i) = -3;    i = i + 1
        else
            options%vars_to_read(i) = zvar;      options%dim_list(i) = 3;    i = i + 1
        endif
        ! endif
        options%zvar        = zvar      ! this could get reassigned from "" to "height_computed" above
        options%zbvar       = zbvar     !; options%vars_to_read(i) = zbvar;      options%dim_list(i) = 3;    i = i + 1

        ! 2D model variables (e.g. Land surface and PBL height)
        options%shvar       = shvar     ; options%vars_to_read(i) = shvar;      options%dim_list(i) = 2;    i = i + 1
        options%lhvar       = lhvar     ; options%vars_to_read(i) = lhvar;      options%dim_list(i) = 2;    i = i + 1
        options%pblhvar     = pblhvar   ; options%vars_to_read(i) = pblhvar;    options%dim_list(i) = 2;    i = i + 1

        ! Shortwave and longwave down at the surface
        options%swdown_var  = swdown_var; options%vars_to_read(i) = swdown_var; options%dim_list(i) = 2;    i = i + 1
        options%lwdown_var  = lwdown_var; options%vars_to_read(i) = lwdown_var; options%dim_list(i) = 2;    i = i + 1

        ! Sea surface temperature
        options%sst_var     = sst_var   ; options%vars_to_read(i) = sst_var;    options%dim_list(i) = 2;    i = i + 1
        options%rain_var    = rain_var  ; options%vars_to_read(i) = rain_var;   options%dim_list(i) = 2;    i = i + 1

        !------------------------------------------------------
        ! these variables are read explicitly so not added to vars_to_read list
        !------------------------------------------------------
        ! variable names for the high resolution domain
        options%hgt_hi          = hgt_hi
        options%landvar         = landvar
        options%lakedepthvar    = lakedepthvar
        options%lat_hi          = lat_hi
        options%lon_hi          = lon_hi
        options%ulat_hi         = ulat_hi
        options%ulon_hi         = ulon_hi
        options%vlat_hi         = vlat_hi
        options%vlon_hi         = vlon_hi

        options%sinalpha_var    = sinalpha_var
        options%cosalpha_var    = cosalpha_var

        ! soil and vegetation parameters
        options%soiltype_var       = soiltype_var
        options%soil_t_var         = soil_t_var
        options%soil_vwc_var       = soil_vwc_var
        options%swe_var            = swe_var
        options%snowh_var          = snowh_var
        options%soil_deept_var     = soil_deept_var
        options%vegtype_var        = vegtype_var
        options%vegfrac_var        = vegfrac_var
        options%vegfracmax_var     = vegfracmax_var
        options%albedo_var         = albedo_var
        options%lai_var            = lai_var
        options%canwat_var         = canwat_var

        ! optional calibration variables for linear wind solution
        options%linear_mask_var     = linear_mask_var
        options%nsq_calibration_var = nsq_calibration_var

        !------------------------------------------------------
        ! external variables for initialization  -  mainly snow-related, can be extended later on. swe only for now?
        !------------------------------------------------------
        options%ext_var_list(:) = ""
        j=1
        options%lat_ext         = lat_ext
        options%lon_ext         = lon_ext
        options%swe_ext         = swe_ext      ; options%ext_var_list(j) = swe_ext;       options%ext_dim_list(j) = 2;    j = j + 1
        options%rho_snow_ext    = rho_snow_ext ; options%ext_var_list(j) = rho_snow_ext;  options%ext_dim_list(j) = 2;    j = j + 1
        options%hsnow_ext       = hsnow_ext    ; options%ext_var_list(j) = hsnow_ext;     options%ext_dim_list(j) = 2;    j = j + 1
        options%tss_ext         = tss_ext      ; options%ext_var_list(j) = tss_ext;       options%ext_dim_list(j) = 2;    j = j + 1
        options%tsoil2D_ext     = tsoil2D_ext    ; options%ext_var_list(j) = tsoil2D_ext;     options%ext_dim_list(j) = 2;    j = j + 1
        options%tsoil3D_ext     = tsoil3D_ext    ; options%ext_var_list(j) = tsoil3D_ext;     options%ext_dim_list(j) = 3;    j = j + 1
        ! options%z_ext      = z_ext   ; options%ext_var_list(j) = z_ext;       options%ext_dim_list(j) = 3;    j = j + 1
        options%time_ext        = time_ext    ; options%ext_var_list(j) = time_ext;      options%ext_dim_list(j) = 1;    j = j + 1




    end subroutine var_namelist


    !> -------------------------------
    !! Initialize the main parameter options
    !!
    !! Reads parameters for the ICAR simulation
    !! These include setting flags that request other namelists be read
    !!
    !! -------------------------------
    subroutine parameters_namelist(filename,options)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(parameter_options_type), intent(inout) :: options

        integer :: name_unit
        type(time_delta_t) :: dt
        ! parameters to read

        real    :: dx, dxlow, outputinterval, restartinterval, inputinterval, t_offset, smooth_wind_distance, frames_per_outfile, agl_cap
        real    :: cfl_reduction_factor, rh_limit, sst_min_limit, cp_limit
        integer :: ntimesteps, wind_iterations
        integer :: longitude_system
        integer :: nz, n_ext_winds,buffer, warning_level, cfl_strictness
        logical :: ideal, readz, readdz, interactive, debug, external_winds, surface_io_only, &
                   mean_winds, mean_fields, restart, advect_density, z_is_geopotential, z_is_on_interface,&
                   high_res_soil_state, use_agl_height, time_varying_z, t_is_potential, qv_is_spec_humidity, &
                   qv_is_relative_humidity, &
                   use_mp_options, use_lt_options, use_adv_options, use_lsm_options, use_bias_correction, &
                   use_block_options, use_cu_options, use_rad_options

        character(len=MAXFILELENGTH) :: date, calendar, start_date, forcing_start_date, end_date
        integer :: year, month, day, hour, minute, second
        character(len=MAXFILELENGTH) :: mp_options_filename, lt_options_filename, &
                                        adv_options_filename, lsm_options_filename, &
                                        bias_options_filename, block_options_filename, &
                                        cu_options_filename, rad_options_filename


        namelist /parameters/ ntimesteps, wind_iterations, outputinterval, frames_per_outfile, inputinterval, surface_io_only,                &
                              dx, dxlow, ideal, readz, readdz, nz, t_offset, rh_limit, sst_min_limit, cp_limit, &
                              debug, warning_level, interactive, restart,                                &
                              external_winds, buffer, n_ext_winds, advect_density, smooth_wind_distance, &
                              mean_winds, mean_fields, z_is_geopotential, z_is_on_interface,             &
                              date, calendar, high_res_soil_state, t_is_potential,                       &
                              qv_is_relative_humidity, qv_is_spec_humidity,                              &
                              use_agl_height, agl_cap, start_date, forcing_start_date, end_date,         &
                              time_varying_z,  longitude_system,            &
                              cfl_reduction_factor,     cfl_strictness,     &
                              mp_options_filename,      use_mp_options,     &
                              block_options_filename,   use_block_options,  &
                              lt_options_filename,      use_lt_options,     &
                              lsm_options_filename,     use_lsm_options,    &
                              adv_options_filename,     use_adv_options,    &
                              bias_options_filename,    use_bias_correction,&
                              cu_options_filename,      use_cu_options,     &
                              rad_options_filename,     use_rad_options

!       default parameters
        surface_io_only     = .False.
        mean_fields         = .False.
        mean_winds          = .False.
        external_winds      = .False.
        n_ext_winds         = 1
        t_offset            = (-9999)
        rh_limit            = -1
        sst_min_limit       = 273.15
        cp_limit            = 500
        buffer              = 0
        advect_density      = .False.
        t_is_potential      = .True.
        qv_is_spec_humidity = .False.
        qv_is_relative_humidity=.False.
        z_is_geopotential   = .False.
        z_is_on_interface   = .False.
        wind_iterations     = 100
        dxlow               = 100000
        restart             = .False.
        ideal               = .False.
        debug               = .False.
        interactive         = .False.
        warning_level       = -9999
        readz               = .False.
        readdz              = .True.
        nz                  =  MAXLEVELS
        smooth_wind_distance= -9999
        calendar            = "gregorian"
        inputinterval       = 3600
        high_res_soil_state = .False.
        use_agl_height      = .False.
        agl_cap             = 300
        date                = ""
        start_date          = ""
        forcing_start_date  = ""
        end_date            = ""
        time_varying_z      = .False.
        cfl_reduction_factor=  0.9
        cfl_strictness      =  3
        inputinterval       =  3600
        outputinterval      =  3600
        frames_per_outfile  =  24
        longitude_system    = kMAINTAIN_LON

        ! flag set to read specific parameterization options
        use_mp_options=.False.
        mp_options_filename = filename

        use_lt_options=.False.
        lt_options_filename = filename

        use_adv_options=.False.
        adv_options_filename = filename

        use_cu_options=.False.
        cu_options_filename = filename

        use_lsm_options=.False.
        lsm_options_filename = filename

        use_rad_options=.False.
        rad_options_filename = filename

        use_bias_correction=.False.
        bias_options_filename = filename

        use_block_options=.False.
        block_options_filename = filename

        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=parameters)
        close(name_unit)

        if (ideal) then
            if (this_image()==1) write(*,*) " Running Idealized simulation "
        endif

        if ((trim(date)=="").and.(trim(start_date)/="")) date = start_date

        if (trim(forcing_start_date)=="") then
            if (trim(date)/="") then
                forcing_start_date=date
            endif
        else
            date=forcing_start_date
        endif
        if (warning_level==-9999) then
            if (debug) then
                warning_level=2
            else
                warning_level=5
            endif
        endif


        if (smooth_wind_distance.eq.(-9999)) then
            smooth_wind_distance=dxlow*2
            if (this_image()==1) write(*,*) " Default smoothing distance = lowdx*2 = ", smooth_wind_distance
        endif

        options%t_offset=t_offset

        if (rh_limit > 5) then
            write(*,*) "RH limit >5 specified, assuming it is in %"
            rh_limit = rh_limit/100
        endif
        options%rh_limit=rh_limit

        options%sst_min_limit = sst_min_limit

        options%cp_limit = cp_limit / 3600.0 ! convert from mm/hr to mm/s

        if (smooth_wind_distance<0) then
            write(*,*) " Wind smoothing must be a positive number"
            write(*,*) " smooth_wind_distance = ",smooth_wind_distance
            if (warning_level>4) then
                stop
            else
                smooth_wind_distance = dxlow*2
            endif
        endif
        options%smooth_wind_distance = smooth_wind_distance
        options%longitude_system = longitude_system

        options%in_dt      = inputinterval
        call options%input_dt%set(seconds=inputinterval)

        options%out_dt     = outputinterval
        call options%output_dt%set(seconds=outputinterval)
        ! if outputing at half-day or longer intervals, create monthly files
        ! if (outputinterval>=43200) then
        !     options%output_file_frequency="monthly"
        ! ! if outputing at half-hour or longer intervals, create daily files
        ! else if (outputinterval>=1800) then
        !     options%output_file_frequency="daily"
        ! ! otherwise create a new output file every timestep
        ! else
        !     options%output_file_frequency="every step"
        ! endif

        ! options%paramters%frames_per_outfile : this may cause trouble with the above, but a nicer way
        options%frames_per_outfile = frames_per_outfile

        ! options%surface_io_only = surface_io_only


        options%calendar=calendar

        call options%initial_time%init(calendar)
        call options%initial_time%set(date)
        if (start_date=="") then
            options%start_time = options%initial_time
        else
            call options%start_time%init(calendar)
            call options%start_time%set(start_date)
        endif
        if (trim(end_date)/="") then
            call options%end_time%init(calendar)
            call options%end_time%set(end_date)
        else
            call dt%set(seconds=ntimesteps * inputinterval)
            options%end_time = options%start_time + dt
        endif

        options%dx               = dx
        options%dxlow            = dxlow
        options%ideal            = ideal
        options%readz            = readz
        options%readdz           = readdz
        options%buffer           = buffer
        options%mean_winds       = mean_winds
        options%mean_fields      = mean_fields
        options%advect_density   = advect_density
        options%debug            = debug
        options%interactive      = interactive
        options%warning_level    = warning_level
        options%use_agl_height   = use_agl_height
        options%agl_cap          = agl_cap

        options%qv_is_relative_humidity = qv_is_relative_humidity
        options%qv_is_spec_humidity= qv_is_spec_humidity
        options%t_is_potential   = t_is_potential
        options%z_is_geopotential= z_is_geopotential
        options%z_is_on_interface= z_is_on_interface

        options%external_winds = external_winds
        options%ext_winds_nfiles = n_ext_winds
        options%restart = restart

        options%nz = nz
        options%wind_iterations = wind_iterations
        options%high_res_soil_state = high_res_soil_state
        options%time_varying_z = time_varying_z

        options%cfl_reduction_factor = cfl_reduction_factor
        options%cfl_strictness = cfl_strictness

        options%use_mp_options      = use_mp_options
        options%mp_options_filename = mp_options_filename

        options%use_lt_options      = use_lt_options
        options%lt_options_filename = lt_options_filename

        options%use_cu_options      = use_cu_options
        options%cu_options_filename = cu_options_filename

        options%use_adv_options     = use_adv_options
        options%adv_options_filename= adv_options_filename

        options%use_lsm_options     = use_lsm_options
        options%lsm_options_filename= lsm_options_filename

        options%use_rad_options     = use_rad_options
        options%rad_options_filename= rad_options_filename

        options%use_bias_correction  = use_bias_correction
        options%bias_options_filename= bias_options_filename

        options%use_block_options     = use_block_options
        options%block_options_filename= block_options_filename

        ! options are updated when complete
    end subroutine parameters_namelist

    !> -------------------------------
    !! Initialize the microphysics options
    !!
    !! Reads the mp_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine mp_parameters_namelist(mp_filename,options)
        implicit none
        character(len=*),   intent(in)    :: mp_filename
        type(options_t),    intent(inout) :: options
        integer :: name_unit

        real    :: Nt_c, TNO, am_s, rho_g, av_s, bv_s, fv_s, av_g, bv_g, av_i
        real    :: Ef_si, Ef_rs, Ef_rg, Ef_ri
        real    :: C_cubes, C_sqrd, mu_r, t_adjust
        logical :: Ef_rw_l, EF_sw_l
        integer :: top_mp_level
        real    :: local_precip_fraction
        integer :: update_interval

        namelist /mp_parameters/ Nt_c, TNO, am_s, rho_g, av_s, bv_s, fv_s, av_g, bv_g, av_i,    &   ! thompson microphysics parameters
                                Ef_si, Ef_rs, Ef_rg, Ef_ri,                                     &   ! thompson microphysics parameters
                                C_cubes, C_sqrd, mu_r, Ef_rw_l, Ef_sw_l, t_adjust,              &   ! thompson microphysics parameters
                                top_mp_level, local_precip_fraction, update_interval

        ! because mp_options could be in a separate file (should probably set all namelists up to have this option)
        if (options%parameters%use_mp_options) then
            if (trim(mp_filename)/=trim(get_options_file())) then
                call version_check(mp_filename, options%parameters)
            endif
        endif

        ! set default parameters
        Nt_c  = 100.e6      !  50, 100,500,1000
        TNO   = 5.0         !  0.5, 5, 50
        am_s  = 0.069       ! 0.052 (Heymsfield), 0.02 (Mitchell), 0.01.
                            ! Note that these values are converted to mks units. Was given as cgs units in Morrison p3 code
        rho_g = 500.0       ! 800, 500, 200
        av_s  = 40.0        ! 11.72 (Locatelli and Hobbs)
        bv_s  = 0.55        ! 0.41
        fv_s  = 100.0       ! 0
        av_g  = 442.0       ! 19.3   from "Cloud-Resolving Modelling of Convective Processes, by Gao and Li,
        bv_g  = 0.89        ! 0.37
        av_i  = 1847.5      ! 700 (Ikawa and Saito)
        Ef_si = 0.05
        Ef_rs = 0.95        ! 1
        Ef_rg = 0.75        ! 1
        Ef_ri = 0.95        ! 1
        C_cubes = 0.5       ! 0.25 Based on Thesis paper "Validation and Improvements of Simulated
                            !      Cloud Microphysics and Orographic Precipitation over the Pacific Northwest"
        C_sqrd  = 0.3
        mu_r    = 0.        ! 1, 2, 5
        t_adjust= 0.0       ! -5, 10, 15
        Ef_rw_l = .False.   ! True sets ef_rw = 1, insted of max 0.95
        Ef_sw_l = .False.   ! True sets ef_rw = 1, insted of max 0.95

        update_interval       = 0 ! update every time step
        top_mp_level          = 0 ! if <=0 just use the actual model top
        local_precip_fraction = 1 ! if <1: the remaining fraction (e.g. 1-x) of precip gets distributed to the surrounding grid cells

        ! read in the namelist
        if (options%parameters%use_mp_options) then
            open(io_newunit(name_unit), file=mp_filename)
            read(name_unit, nml=mp_parameters)
            close(name_unit)
        endif

        ! store the data back into the mp_options datastructure
        options%mp_options%Nt_c     = Nt_c
        options%mp_options%TNO      = TNO
        options%mp_options%am_s     = am_s
        options%mp_options%rho_g    = rho_g
        options%mp_options%av_s     = av_s
        options%mp_options%bv_s     = bv_s
        options%mp_options%fv_s     = fv_s
        options%mp_options%av_g     = av_g
        options%mp_options%bv_g     = bv_g
        options%mp_options%av_i     = av_i
        options%mp_options%Ef_si    = Ef_si
        options%mp_options%Ef_rs    = Ef_rs
        options%mp_options%Ef_rg    = Ef_rg
        options%mp_options%Ef_ri    = Ef_ri
        options%mp_options%mu_r     = mu_r
        options%mp_options%t_adjust = t_adjust
        options%mp_options%C_cubes  = C_cubes
        options%mp_options%C_sqrd   = C_sqrd
        options%mp_options%Ef_rw_l  = Ef_rw_l
        options%mp_options%Ef_sw_l  = Ef_sw_l

        if (top_mp_level < 0) top_mp_level = options%parameters%nz + top_mp_level

        options%mp_options%update_interval      = update_interval
        options%mp_options%top_mp_level         = top_mp_level
        options%mp_options%local_precip_fraction= local_precip_fraction

    end subroutine mp_parameters_namelist

    !> -------------------------------
    !! Initialize the blocking options
    !!
    !! Reads the block_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine block_parameters_namelist(filename, options)
        implicit none

        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options

        integer :: name_unit

        real    :: blocking_contribution  ! fractional contribution of flow blocking perturbation that is added [0-1]
        real    :: smooth_froude_distance ! distance (m) over which Froude number is smoothed
        integer :: n_smoothing_passes     ! number of times the smoothing window is applied
        real    :: block_fr_max           ! max froude no at which flow is only partially blocked above, no blocking
        real    :: block_fr_min           ! min froude no at which flow is only partially blocked below, full blocking
        logical :: block_flow             ! use a blocking parameterization


        ! define the namelist
        namelist /block_parameters/ smooth_froude_distance, &
                                    n_smoothing_passes,     &
                                    block_fr_max,           &
                                    block_fr_min,           &
                                    blocking_contribution,  &
                                    block_flow

        ! because block_options could be in a separate file
        if (options%parameters%use_block_options) then
            if (trim(filename)/=trim(get_options_file())) then
                call version_check(filename,options%parameters)
            endif
        endif

        ! set default values
        smooth_froude_distance  = 6000
        n_smoothing_passes      = 3
        block_fr_max            = 0.75
        block_fr_min            = 0.5
        blocking_contribution   = 0.5
        block_flow              = .False.

        ! read the namelist options
        if (options%parameters%use_block_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=block_parameters)
            close(name_unit)
        endif

        ! copy values into data type
        associate( opt => options%block_options )
            opt%smooth_froude_distance = smooth_froude_distance
            opt%n_smoothing_passes     = n_smoothing_passes
            opt%block_fr_max           = block_fr_max
            opt%block_fr_min           = block_fr_min
            opt%blocking_contribution  = blocking_contribution
            opt%block_flow             = block_flow
        end associate

    end subroutine block_parameters_namelist


    !> -------------------------------
    !! Initialize the Linear Theory options
    !!
    !! Reads the lt_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine lt_parameters_namelist(filename, options)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options

        integer :: name_unit

        integer :: vert_smooth
        logical :: variable_N           ! Compute the Brunt Vaisala Frequency (N^2) every time step
        logical :: smooth_nsq               ! Smooth the Calculated N^2 over vert_smooth vertical levels
        integer :: buffer                   ! number of grid cells to buffer around the domain MUST be >=1
        integer :: stability_window_size    ! window to average nsq over
        real    :: max_stability            ! limits on the calculated Brunt Vaisala Frequency
        real    :: min_stability            ! these may need to be a little narrower.
        real    :: linear_contribution      ! multiplier on uhat,vhat before adding to u,v
        real    :: linear_update_fraction   ! controls the rate at which the linearfield updates (should be calculated as f(in_dt))

        real    :: N_squared                ! static Brunt Vaisala Frequency (N^2) to use
        logical :: remove_lowres_linear     ! attempt to remove the linear mountain wave from the forcing low res model
        real    :: rm_N_squared             ! static Brunt Vaisala Frequency (N^2) to use in removing linear wind field
        real    :: rm_linear_contribution   ! fractional contribution of linear perturbation to wind field to remove from the low-res field

        logical :: spatial_linear_fields    ! use a spatially varying linear wind perturbation
        logical :: linear_mask              ! use a spatial mask for the linear wind field
        logical :: nsq_calibration          ! use a spatial mask to calibrate the nsquared (brunt vaisala frequency) field

        ! Look up table generation parameters
        real    :: dirmax, dirmin
        real    :: spdmax, spdmin
        real    :: nsqmax, nsqmin
        integer :: n_dir_values, n_nsq_values, n_spd_values
        real    :: minimum_layer_size       ! Minimum vertical step to permit when computing LUT.
                                            ! If model layers are thicker, substepping will be used.

        ! parameters to control reading from or writing an LUT file
        logical :: read_LUT, write_LUT
        character(len=MAXFILELENGTH) :: u_LUT_Filename, v_LUT_Filename, LUT_Filename
        logical :: overwrite_lt_lut

        ! define the namelist
        namelist /lt_parameters/ variable_N, smooth_nsq, buffer, stability_window_size, max_stability, min_stability, &
                                 linear_contribution, linear_update_fraction, N_squared, vert_smooth, &
                                 remove_lowres_linear, rm_N_squared, rm_linear_contribution, &
                                 spatial_linear_fields, linear_mask, nsq_calibration, minimum_layer_size, &
                                 dirmax, dirmin, spdmax, spdmin, nsqmax, nsqmin, n_dir_values, n_nsq_values, n_spd_values, &
                                 read_LUT, write_LUT, u_LUT_Filename, v_LUT_Filename, overwrite_lt_lut, LUT_Filename

         ! because lt_options could be in a separate file
         if (options%parameters%use_lt_options) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options%parameters)
             endif
         endif


        ! set default values
        variable_N = .True.
        smooth_nsq = .True.
        buffer = 50                    ! number of grid cells to buffer around the domain MUST be >=1
        stability_window_size = 10     ! window to average nsq over
        vert_smooth = 10
        max_stability = 6e-4           ! limits on the calculated Brunt Vaisala Frequency
        min_stability = 1e-7           ! these may need to be a little narrower.

        N_squared = 3e-5               ! static Brunt Vaisala Frequency (N^2) to use
        linear_contribution = 1        ! fractional contribution of linear perturbation to wind field (e.g. u_hat multiplied by this)
        remove_lowres_linear = .False. ! attempt to remove the linear mountain wave from the forcing low res model
        rm_N_squared = 3e-5            ! static Brunt Vaisala Frequency (N^2) to use in removing linear wind field
        rm_linear_contribution = 1     ! fractional contribution of linear perturbation to wind field to remove from the low-res field

        linear_update_fraction = 0.2   ! fraction of linear perturbation to add each time step
        spatial_linear_fields = .True. ! use a spatially varying linear wind perturbation
        linear_mask = .False.          ! use a spatial mask for the linear wind field
        nsq_calibration = .False.      ! use a spatial mask to calibrate the nsquared (brunt vaisala frequency) field

        ! look up table generation parameters
        dirmax = 2*pi
        dirmin = 0
        spdmax = 30
        spdmin = 0
        nsqmax = log(max_stability)
        nsqmin = log(min_stability)
        n_dir_values = 24
        n_nsq_values = 5
        n_spd_values = 6
        minimum_layer_size = 100

        read_LUT = .False.
        write_LUT = .True.
        u_LUT_Filename = "MISSING"
        v_LUT_Filename = "MISSING"
        LUT_Filename   = "MISSING"
        overwrite_lt_lut = .True.

        ! read the namelist options
        if (options%parameters%use_lt_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=lt_parameters)
            close(name_unit)
        endif

        ! store everything in the lt_options structure
        associate(opt => options%lt_options )
            opt%buffer = buffer
            opt%stability_window_size = stability_window_size
            opt%max_stability = max_stability
            opt%min_stability = min_stability
            opt%variable_N = variable_N
            opt%smooth_nsq = smooth_nsq

            if (vert_smooth<0) then
                write(*,*) " Vertical smoothing must be a positive integer"
                write(*,*) " vert_smooth = ",vert_smooth
                stop
            endif
            opt%vert_smooth=vert_smooth

            opt%N_squared = N_squared
            opt%linear_contribution = linear_contribution
            opt%remove_lowres_linear = remove_lowres_linear
            opt%rm_N_squared = rm_N_squared
            opt%rm_linear_contribution = rm_linear_contribution
            opt%linear_update_fraction = linear_update_fraction
            opt%spatial_linear_fields = spatial_linear_fields
            opt%linear_mask = linear_mask
            opt%nsq_calibration = nsq_calibration
            opt%dirmax = dirmax
            opt%dirmin = dirmin
            opt%spdmax = spdmax
            opt%spdmin = spdmin
            opt%nsqmax = nsqmax
            opt%nsqmin = nsqmin
            opt%n_dir_values = n_dir_values
            opt%n_nsq_values = n_nsq_values
            opt%n_spd_values = n_spd_values
            opt%minimum_layer_size = minimum_layer_size
            opt%read_LUT = read_LUT
            opt%write_LUT = write_LUT
            if (trim(u_LUT_Filename)=="MISSING") then
                if (trim(LUT_Filename)=="MISSING") then
                    u_LUT_Filename="Linear_Theory_LUT.nc"
                else
                    u_LUT_Filename=LUT_Filename
                endif
            endif

            ! check if directory paths to LUT file strings exist, if not fails before write step
            call check_writeable_path(u_LUT_Filename)
            call check_writeable_path(v_LUT_Filename)
            opt%u_LUT_Filename = u_LUT_Filename
            opt%v_LUT_Filename = v_LUT_Filename
            opt%overwrite_lt_lut = overwrite_lt_lut

        end associate

    end subroutine lt_parameters_namelist


    !> -------------------------------
    !! Initialize the advection options
    !!
    !! Reads the adv_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine adv_parameters_namelist(filename, options)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options

        type(adv_options_type)::adv_options
        integer :: name_unit

        logical :: boundary_buffer          ! apply some smoothing to the x and y boundaries in MPDATA
        logical :: flux_corrected_transport ! use the flux corrected transport option in MPDATA
        integer :: mpdata_order             ! MPDATA order of correction (e.g. 1st=upwind, 2nd=classic, 3rd=better)

        ! define the namelist
        namelist /adv_parameters/ boundary_buffer, flux_corrected_transport, mpdata_order

         ! because adv_options could be in a separate file
         if (options%parameters%use_adv_options) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options%parameters)
             endif
         endif


        ! set default values
        boundary_buffer = .False.
        flux_corrected_transport = .True.
        mpdata_order = 2

        ! read the namelist options
        if (options%parameters%use_adv_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=adv_parameters)
            close(name_unit)
        endif

        ! store everything in the adv_options structure
        adv_options%boundary_buffer = boundary_buffer
        adv_options%flux_corrected_transport = flux_corrected_transport
        adv_options%mpdata_order = mpdata_order

        ! copy the data back into the global options data structure
        options%adv_options = adv_options
    end subroutine adv_parameters_namelist


    !> -------------------------------
    !! Initialize the convection scheme options
    !!
    !! Reads the cu_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine cu_parameters_namelist(filename, options)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options

        type(cu_options_type) :: cu_options
        integer :: name_unit

        real :: tendency_fraction, tend_qv_fraction, tend_qc_fraction, tend_th_fraction, tend_qi_fraction
        real :: stochastic_cu


        ! define the namelist
        namelist /cu_parameters/ tendency_fraction, tend_qv_fraction, tend_qc_fraction, tend_th_fraction, tend_qi_fraction, &
                                 stochastic_cu

        ! because adv_options could be in a separate file
        if (options%parameters%use_cu_options) then
            if (trim(filename)/=trim(get_options_file())) then
                call version_check(filename,options%parameters)
            endif
        endif

        ! set default values
        stochastic_cu       = kNO_STOCHASTIC

        tendency_fraction   = 1.0
        tend_qv_fraction    = -1.0
        tend_qc_fraction    = -1.0
        tend_th_fraction    = -1.0
        tend_qi_fraction    = -1.0

        ! read the namelist options
        if (options%parameters%use_cu_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=cu_parameters)
            close(name_unit)
        endif

        ! if not set separately, default to the global tendency setting
        if (tend_qv_fraction < 0) tend_qv_fraction = tendency_fraction
        if (tend_qc_fraction < 0) tend_qc_fraction = tendency_fraction
        if (tend_th_fraction < 0) tend_th_fraction = tendency_fraction
        if (tend_qi_fraction < 0) tend_qi_fraction = tendency_fraction

        ! store everything in the cu_options structure
        cu_options%stochastic_cu        = stochastic_cu
        cu_options%tendency_fraction    = tendency_fraction
        cu_options%tend_qv_fraction     = tend_qv_fraction
        cu_options%tend_qc_fraction     = tend_qc_fraction
        cu_options%tend_th_fraction     = tend_th_fraction
        cu_options%tend_qi_fraction     = tend_qi_fraction

        ! copy the data back into the global options data structure
        options%cu_options = cu_options
    end subroutine cu_parameters_namelist



    !> -------------------------------
    !! Sets the default value for each of three land use categories depending on the LU_Categories input
    !!
    !! -------------------------------
    subroutine set_default_LU_categories(options, urban_category, ice_category, water_category, LU_Categories, lake_category)
        ! if various LU categories were not defined in the namelist (i.e. they == -1) then attempt
        ! to define default values for them based on the LU_Categories variable supplied.
        implicit none
        type(options_t),    intent(inout)::options
        integer, intent(inout) :: urban_category, ice_category, water_category, lake_category
        character(len=MAXVARLENGTH), intent(in) :: LU_Categories

        if (trim(LU_Categories)=="MODIFIED_IGBP_MODIS_NOAH") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15
            if (water_category==-1) water_category = 17
            if (lake_category==-1) lake_category = 21

        elseif (trim(LU_Categories)=="USGS") then
            if (urban_category==-1) urban_category = 1
            if (ice_category==-1)   ice_category = -1
            if (water_category==-1) water_category = 16
            ! if (lake_category==-1) lake_category = 16  ! No separate lake category!
            if((options%physics%watersurface==kWATER_LAKE) .AND. (this_image()==1)) then
                write(*,*) "WARNING: Lake model selected (water=3), but USGS LU-categories has no lake category"
            endif

        elseif (trim(LU_Categories)=="USGS-RUC") then
            if (urban_category==-1) urban_category = 1
            if (ice_category==-1)   ice_category = 24
            if (water_category==-1) water_category = 16
            if (lake_category==-1) lake_category = 28
            ! also note, lakes_category = 28
            ! write(*,*) "WARNING: not handling lake category (28)"

        elseif (trim(LU_Categories)=="MODI-RUC") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15
            if (water_category==-1) water_category = 17
            if (lake_category==-1) lake_category = 21
            ! also note, lakes_category = 21
            ! write(*,*) "WARNING: not handling lake category (21)"

        elseif (trim(LU_Categories)=="NLCD40") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15 ! and 22?
            ! if (water_category==-1) water_category = 17 ! and 21 'Open Water'
            if(options%physics%watersurface==kWATER_LAKE) write(*,*) "WARNING: Lake model selected (water=3), but NLCD40 LU-categories has no lake category"
            write(*,*) "WARNING: not handling all varients of categories (e.g. permanent_snow=15 is, but permanent_snow_ice=22 is not)"
        endif

    end subroutine set_default_LU_categories


    !> -------------------------------
    !! Initialize the bias correction options
    !!
    !! Reads the bias_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine bias_parameters_namelist(filename, options)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options

        type(bias_options_type)     ::  bias_options
        integer :: name_unit

        character(len=MAXFILELENGTH):: bias_correction_filename ! file containing bias correction data
        character(len=MAXVARLENGTH) :: rain_fraction_var        ! name of variable containing the fraction to multiply rain by

        ! define the namelist
        namelist /bias_parameters/ bias_correction_filename, rain_fraction_var

         ! because adv_options could be in a separate file
         if (options%parameters%use_bias_correction) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options%parameters)
             endif
         endif


        ! set default values
        bias_correction_filename = options%parameters%init_conditions_file
        rain_fraction_var        = "rain_fraction"

        ! read the namelist options
        if (options%parameters%use_bias_correction) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=bias_parameters)
            close(name_unit)
        endif

        ! store everything in the bias_options structure
        bias_options%filename           = bias_correction_filename
        bias_options%rain_fraction_var  = rain_fraction_var

        ! copy the data back into the global options data structure
        options%bias_options = bias_options
    end subroutine bias_parameters_namelist


    !> -------------------------------
    !! Initialize the land surface model options
    !!
    !! Reads the lsm_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine lsm_parameters_namelist(filename, options)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout)::options

        type(lsm_options_type) :: lsm_options
        integer :: name_unit

        character(len=MAXVARLENGTH) :: LU_Categories ! Category definitions (e.g. USGS, MODIFIED_IGBP_MODIS_NOAH)
        real    :: lh_feedback_fraction
        real    :: sh_feedback_fraction
        real    :: sfc_layer_thickness
        real    :: dz_lsm_modification
        real    :: wind_enhancement
        real    :: max_swe
        logical :: monthly_vegfrac                   ! read in 12 months of vegfrac data
        logical :: monthly_albedo                    ! same for albedo (requires vegfrac be monthly)
        integer :: update_interval                   ! minimum number of seconds between LSM updates
        integer :: urban_category                    ! index that defines the urban category in LU_Categories
        integer :: ice_category                      ! index that defines the ice category in LU_Categories
        integer :: water_category                    ! index that defines the water category in LU_Categories
        integer :: lake_category                    ! index that defines the lake category in (some) LU_Categories

        ! define the namelist
        namelist /lsm_parameters/ LU_Categories, lh_feedback_fraction, sh_feedback_fraction, update_interval, &
                                  urban_category, ice_category, water_category, lake_category,                &
                                  monthly_vegfrac, monthly_albedo, sfc_layer_thickness, dz_lsm_modification,  &
                                  wind_enhancement, max_swe

         ! because adv_options could be in a separate file
         if (options%parameters%use_lsm_options) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options%parameters)
             endif
         endif


        ! set default values
        LU_Categories   = "MODIFIED_IGBP_MODIS_NOAH"
        update_interval = 300 ! 5 minutes
        monthly_vegfrac = .False.
        monthly_albedo = .False.

        ! default values for these will be set after reading LU_Categories
        urban_category  = -1
        ice_category    = -1
        water_category  = -1
        lake_category   = -1

        lh_feedback_fraction = 1.0
        sh_feedback_fraction = 0.625
        sfc_layer_thickness  = 400.0
        dz_lsm_modification  = 0.5
        wind_enhancement = 1.5

        max_swe = 1e10

        ! read the namelist options
        if (options%parameters%use_lsm_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=lsm_parameters)
            close(name_unit)
        endif

        call set_default_LU_categories(options, urban_category, ice_category, water_category, LU_Categories, lake_category)

        ! store everything in the lsm_options structure
        lsm_options%LU_Categories        = LU_Categories
        lsm_options%monthly_vegfrac      = monthly_vegfrac
        lsm_options%monthly_albedo       = monthly_albedo
        lsm_options%update_interval      = update_interval
        lsm_options%urban_category       = urban_category
        lsm_options%ice_category         = ice_category
        lsm_options%water_category       = water_category
        lsm_options%lake_category        = lake_category
        lsm_options%lh_feedback_fraction = lh_feedback_fraction
        lsm_options%sh_feedback_fraction = sh_feedback_fraction
        lsm_options%sfc_layer_thickness  = sfc_layer_thickness
        lsm_options%dz_lsm_modification  = dz_lsm_modification
        lsm_options%wind_enhancement     = wind_enhancement
        lsm_options%max_swe              = max_swe

        ! copy the data back into the global options data structure
        options%lsm_options = lsm_options
    end subroutine lsm_parameters_namelist


    !> -------------------------------
    !! Initialize the radiation model options
    !!
    !! Reads the rad_parameters namelist or sets default values
    !! -------------------------------
    subroutine rad_parameters_namelist(filename, options)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout)::options

        type(rad_options_type) :: rad_options
        integer :: name_unit

        integer :: update_interval_rrtmg             ! minimum number of seconds between RRTMG updates
        integer :: icloud                            ! how RRTMG interacts with clouds
        logical :: read_ghg
        logical :: use_simple_sw
        ! define the namelist
        namelist /rad_parameters/ update_interval_rrtmg, icloud, read_ghg, use_simple_sw


         ! because adv_options could be in a separate file
         if (options%parameters%use_rad_options) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options%parameters)
             endif
         endif


        ! set default values
        update_interval_rrtmg = 1800 ! 30 minutes
        icloud          = 3    ! effective radius from microphysics scheme
        read_ghg        = .false.
        use_simple_sw   = .false.

        ! read the namelist options
        if (options%parameters%use_rad_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=rad_parameters)
            close(name_unit)
        endif

        ! store everything in the radiation_options structure
        rad_options%update_interval_rrtmg = update_interval_rrtmg
        rad_options%icloud                = icloud
        rad_options%read_ghg              = read_ghg
        rad_options%use_simple_sw         = use_simple_sw

        ! copy the data back into the global options data structure
        options%rad_options = rad_options
    end subroutine rad_parameters_namelist



    !> -------------------------------
    !! Set up model levels, either read from a namelist, or from a default set of values
    !!
    !! Reads the z_info namelist or sets default values
    !!
    !! -------------------------------
    subroutine model_levels_namelist(filename,options)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(parameter_options_type), intent(inout) :: options

        integer :: name_unit, this_level
        real, allocatable, dimension(:) :: dz_levels
        real, dimension(45) :: fulldz
        logical :: space_varying, fixed_dz_advection, dz_modifies_wind, sleve, use_terrain_difference

        real :: flat_z_height, decay_rate_L_topo, decay_rate_S_topo, sleve_n
        integer :: terrain_smooth_windowsize, terrain_smooth_cycles

        namelist /z_info/ dz_levels, space_varying, dz_modifies_wind, flat_z_height, fixed_dz_advection, sleve, terrain_smooth_windowsize, terrain_smooth_cycles, decay_rate_L_topo, decay_rate_S_topo, sleve_n, use_terrain_difference

        this_level=1
        space_varying = .False.
        fixed_dz_advection = .False.
        dz_modifies_wind = .False.
        flat_z_height = -1
        sleve = .False.
        terrain_smooth_windowsize = 3
        terrain_smooth_cycles = 5
        decay_rate_L_topo = 2.
        decay_rate_S_topo = 6.
        sleve_n = 1.2
        use_terrain_difference = .False.

        ! read the z_info namelist if requested
        if (options%readdz) then
            allocate(dz_levels(MAXLEVELS))
            dz_levels=-999

            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=z_info)
            close(name_unit)

            ! if nz wasn't specified in the namelist, we assume a HUGE number of levels
            ! so now we have to figure out what the actual number of levels read was
            if (options%nz==MAXLEVELS) then
                do this_level=1,MAXLEVELS-1
                    if (dz_levels(this_level+1)<=0) then
                        options%nz=this_level
                        exit
                    endif
                end do
                options%nz=this_level
            endif

            ! allocate(options%dz_levels(options%nz))
            if (minval(dz_levels(1:options%nz)) < 0) then
                if (this_image()==1) write(*,*) "WARNING: dz seems to be less than 0, this is not physical and is probably an error in the namelist"
                if (this_image()==1) write(*,*) "Note that the gfortran compiler will not read dz_levels spread across multiple lines. "
                error stop
            endif

            options%dz_levels = -1
            options%dz_levels(1:options%nz)=dz_levels(1:options%nz)
            deallocate(dz_levels)
        else
        ! if we are not reading dz from the namelist, use default values from a WRF run
        ! default mean layer thicknesses from a 36km WRF run over the "CO-headwaters" domain
            fulldz=[36.,   51.,   58.,   73.,   74.,  111.,  113.,  152.,  155.,  157.,  160.,  245., &
                   251.,  258.,  265.,  365.,  379.,  395.,  413.,  432.,  453.,  476.,  503.,  533., &
                   422.,  443.,  467.,  326.,  339.,  353.,  369.,  386.,  405.,  426.,  450.,  477., &
                   455.,  429.,  396.,  357.,  311.,  325.,  340.,  356.,  356.]
           if (options%nz>45) then
               options%nz=45
           endif
            ! allocate(options%dz_levels(options%nz))
            options%dz_levels = -1
            options%dz_levels(1:options%nz) = fulldz(1:options%nz)
        endif

        options%dz_modifies_wind = dz_modifies_wind
        options%space_varying_dz = space_varying
        options%flat_z_height = flat_z_height
        options%fixed_dz_advection = fixed_dz_advection
        options%sleve = sleve
        options%terrain_smooth_windowsize = terrain_smooth_windowsize
        options%terrain_smooth_cycles = terrain_smooth_cycles
        options%decay_rate_L_topo = decay_rate_L_topo  ! decay_rate_large_scale_topography
        options%decay_rate_S_topo = decay_rate_S_topo ! decay_rate_small_scale_topography !
        options%sleve_n = sleve_n
        options%use_terrain_difference = use_terrain_difference


        if (dz_modifies_wind) then
            write(*,*) "WARNING: setting dz_modifies_wind to true is not recommended, use wind = 2 instead"
            write(*,*) "if you want to continue and enable this, you will need to change this code in the options_obj"
            error stop
        endif

    end subroutine model_levels_namelist

    !> ----------------------------------------------------------------------------
    !!  Read in the name of the boundary condition files from a text file
    !!
    !!  @param      filename        The name of the text file to read
    !!  @param[out] forcing_files   An array to store the filenames in
    !!  @retval     nfiles          The number of files read.
    !!
    !! ----------------------------------------------------------------------------
    function read_forcing_file_names(filename, forcing_files) result(nfiles)
        implicit none
        character(len=*) :: filename
        character(len=MAXFILELENGTH), dimension(MAX_NUMBER_FILES) :: forcing_files
        integer :: nfiles
        integer :: file_unit
        integer :: i, error
        character(len=MAXFILELENGTH) :: temporary_file

        open(unit=io_newunit(file_unit), file=filename)
        i=0
        error=0
        do while (error==0)
            read(file_unit, *, iostat=error) temporary_file
            if (error==0) then
                i=i+1
                forcing_files(i) = temporary_file
            endif
        enddo
        close(file_unit)
        nfiles = i
        ! print out a summary
        if (this_image()==1) write(*,*) "  Boundary conditions files to be used:"
        if (nfiles>10) then
            if (this_image()==1) write(*,*) "    nfiles=", trim(str(nfiles)), ", too many to print."
            if (this_image()==1) write(*,*) "    First file:", trim(forcing_files(1))
            if (this_image()==1) write(*,*) "    Last file: ", trim(forcing_files(nfiles))
        else
            do i=1,nfiles
                if (this_image()==1) write(*,*) "      ",trim(forcing_files(i))
            enddo
        endif

    end function read_forcing_file_names

    !> -------------------------------
    !! Initialize the list of input files
    !!
    !! Reads the file_list namelist or sets default values
    !!
    !! -------------------------------
    subroutine filename_namelist(filename, options)
        ! read in filenames from up to two namelists
        ! init_conditions_file = high res grid data
        ! boundary_files       = list of files for boundary conditions (low-res)
        ! ext_wind_files       = list of files to read wind data on the high-res grid (optional)
        ! linear_mask_file     = file to read a high-res fractional mask for the linear perturbation
        implicit none
        character(len=*),             intent(in)    :: filename
        type(parameter_options_type), intent(inout) :: options

        character(len=MAXFILELENGTH) :: init_conditions_file, output_file, forcing_file_list, &
                                        linear_mask_file, nsq_calibration_file, external_files

        character(len=MAXFILELENGTH), allocatable :: boundary_files(:), ext_wind_files(:)
        integer :: name_unit, nfiles, i

        ! set up namelist structures

        namelist /files_list/ init_conditions_file, output_file, boundary_files, forcing_file_list, &
                              linear_mask_file, nsq_calibration_file, external_files

        namelist /ext_winds_info/ ext_wind_files

        linear_mask_file="MISSING"
        nsq_calibration_file="MISSING"
        allocate(boundary_files(MAX_NUMBER_FILES))
        boundary_files="MISSING"
        forcing_file_list="MISSING"
        external_files="MISSING"

        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=files_list)
        close(name_unit)

        options%external_files = external_files

        options%init_conditions_file=init_conditions_file

        i=1
        do while (trim(boundary_files(i))/=trim("MISSING"))
            i=i+1
        end do
        nfiles=i-1

        if ((nfiles==0).and.(trim(forcing_file_list)/="MISSING")) then
            nfiles = read_forcing_file_names(forcing_file_list, boundary_files)
        else
            if (nfiles==0) then
                stop "No boundary conditions files specified."
            endif
        endif

        allocate(options%boundary_files(nfiles))
        options%boundary_files(1:nfiles) = boundary_files(1:nfiles)
        deallocate(boundary_files)

        if (trim(linear_mask_file)=="MISSING") then
            linear_mask_file = options%init_conditions_file
        endif
        options%linear_mask_file=linear_mask_file
        if (trim(nsq_calibration_file)=="MISSING") then
            nsq_calibration_file = options%init_conditions_file
        endif
        options%nsq_calibration_file=nsq_calibration_file

        if (options%external_winds) then
            allocate(ext_wind_files(options%ext_winds_nfiles))

            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=ext_winds_info)
            close(name_unit)

            allocate(options%ext_wind_files(options%ext_winds_nfiles))
            options%ext_wind_files(1:options%ext_winds_nfiles) = ext_wind_files(1:options%ext_winds_nfiles)
            deallocate(ext_wind_files)
        endif
    end subroutine filename_namelist


end submodule
