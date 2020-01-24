!> ----------------------------------------------------------------------------
!!  Read model options from namelist structures
!!  Also checks commandline arguments to an options filename
!!
!!  Entry point is init_options, everything else follows
!!
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!  Trude Eidhammer added the ability to read microphysics parameters
!!
!! ----------------------------------------------------------------------------
module initialize_options
    use data_structures
    use io_routines,                only : io_nearest_time_step, io_newunit
    use model_tracking,             only : print_model_diffs
    use time,                       only : date_to_mjd, parse_date, time_init
    use time_io,                    only : find_timestep_in_file
    use time_object,                only : Time_type
    use time_delta_object,          only : time_delta_t
    use string,                     only : str

    implicit none

    private
    public :: init_options
contains


    subroutine init_options(options)
!       reads a series of options from a namelist file and stores them in the
!       options data structure
        implicit none
        type(options_type), intent(inout) :: options
        character(len=MAXFILELENGTH) :: options_filename

        options_filename=get_options_file()
        write(*,*) "Using options file = ", trim(options_filename)

        call version_check(         options_filename,   options)
        call physics_namelist(      options_filename,   options)
        call var_namelist(          options_filename,   options)
        call parameters_namelist(   options_filename,   options)
        call model_levels_namelist( options_filename,   options)

        call lt_parameters_namelist(    options%lt_options_filename,    options)
        call block_parameters_namelist( options%block_options_filename, options)
        call mp_parameters_namelist(    options%mp_options_filename,    options)
        call adv_parameters_namelist(   options%adv_options_filename,   options)
        call lsm_parameters_namelist(   options%lsm_options_filename,   options)
        call bias_parameters_namelist(  options%bias_options_filename,  options)

        if (options%restart) then
            call init_restart_options(options_filename,options)
        endif

        call filename_namelist(options_filename, options)
        ! check for any inconsistencies in the options requested
        call options_check(options)
    end subroutine init_options


    function get_options_file() result(options_file)
        implicit none
        character(len=MAXFILELENGTH) ::options_file
        integer :: error
        logical :: file_exists

        if (command_argument_count()>0) then
            call get_command_argument(1,options_file, status=error)
            if (error>0) then
                options_file="icar_options.nml"
            elseif (error==-1) then
                write(*,*) "Options filename = ", trim(options_file), " ...<cutoff>"
                write(*,*) "Maximum filename length = ", MAXFILELENGTH
                stop "ERROR: options filename too long"
            endif
        else
            options_file="icar_options.nml"
        endif
        INQUIRE(file=trim(options_file), exist=file_exists)
        if (.not.file_exists) then
            write(*,*) "Using options file = ", trim(options_file)
            stop "Options file does not exist. "
        endif
    end function



! check the version number in the namelist file and compare to the current model version
! if the namelist version doesn't match, print the differences between that version and this
! and STOP execution
    subroutine version_check(filename,options)
        character(len=*),intent(in) :: filename
        type(options_type),intent(inout)::options
        character(len=MAXVARLENGTH) :: version,comment
        integer:: name_unit

        namelist /model_version/ version,comment
        !default comment:
        comment="Model testing"
        ! read namelists
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=model_version)
        close(name_unit)
        if (version.ne.kVERSION_STRING) then
            write(*,*) "Model version does not match namelist version"
            write(*,*) "  Model version: ",kVERSION_STRING
            write(*,*) "  Namelist version: ",trim(version)
            call print_model_diffs(version)
            stop
        endif
        options%version=version
        options%comment=comment
        write(*,*) "  Model version: ",trim(version)
    end subroutine version_check

    subroutine options_check(options)
        ! Minimal error checking on option settings
        implicit none
        type(options_type), intent(inout)::options

        if (options%t_offset.eq.(-9999)) then
            if (options%warning_level>0) then
                write(*,*) "WARNING, WARNING, WARNING"
                write(*,*) "WARNING, Using default t_offset=0"
                write(*,*) "WARNING, WARNING, WARNING"
            endif
            options%t_offset = 0
        endif


        ! convection can modify wind field, and ideal doesn't rebalance winds every timestep
        if ((options%physics%convection.ne.0).and.(options%ideal)) then
            if (options%warning_level>3) then
                write(*,*) ""
                write(*,*) "WARNING WARNING WARNING"
                write(*,*) "WARNING, Running convection in ideal mode may be bad..."
                write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%warning_level==10) then
                write(*,*) "Set warning_level<10 to continue"
                stop
            endif
        endif
        ! if using a real LSM, feedback will probably keep hot-air from getting even hotter, so not likely a problem
        if ((options%physics%landsurface>1).and.(options%physics%boundarylayer==0)) then
            if (options%warning_level>2) then
                write(*,*) ""
                write(*,*) "WARNING WARNING WARNING"
                write(*,*) "WARNING, Running LSM without PBL may overheat the surface and CRASH the model. "
                write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%warning_level>=7) then
                write(*,*) "Set warning_level<7 to continue"
                stop
            endif
        endif
        ! if using perscribed LSM fluxes, no feedbacks are present, so the surface layer is likely to overheat.
        if ((options%physics%landsurface==1).and.(options%physics%boundarylayer==0)) then
            if (options%warning_level>0) then
                write(*,*) ""
                write(*,*) "WARNING WARNING WARNING"
                write(*,*) "WARNING, Prescribed LSM fluxes without a PBL may overheat the surface and CRASH. "
                write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%warning_level>=5) then
                write(*,*) "Set warning_level<5 to continue"
                stop
            endif
        endif

        ! prior to v 0.9.3 this was assumed, so throw a warning now just in case.
        if ((options%z_is_geopotential .eqv. .False.).and.(options%zvar=="PH")) then
            if (options%warning_level>1) then
                write(*,*) ""
                write(*,*) "WARNING WARNING WARNING"
                write(*,*) "WARNING z variable is not assumed to be geopotential height when it is 'PH'."
                write(*,*) "WARNING If z is geopotential, set z_is_geopotential=True in the namelist."
                write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%warning_level>=7) then
                write(*,*) "Set warning_level<7 to continue"
                stop
            endif
        endif
        if ((options%z_is_geopotential .eqv. .True.).and.(options%z_is_on_interface .eqv. .False.)) then
            if (options%warning_level>1) then
                write(*,*) ""
                write(*,*) "WARNING WARNING WARNING"
                write(*,*) "WARNING geopotential height is no longer assumed to be on model interface levels."
                write(*,*) "WARNING To interpolate geopotential, set z_is_on_interface=True in the namelist. "
                write(*,*) "WARNING WARNING WARNING"
            endif
        endif

    end subroutine options_check

    subroutine init_restart_options(filename, options)
        ! initialize the restart specifications
        ! read in the namelist, and calculate the restart_step if appropriate
        character(len=*),  intent(in)    :: filename    ! name of the file containing the restart namelist
        type(options_type),intent(inout) :: options     ! options data structure to store output

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
            write(*,*) "------ Invalid restart_date ERROR ------"
            stop "restart_date must be specified in the namelist"
        endif

        ! calculate the modified julian day for th restart date given
        call restart_time%init(options%calendar)
        call restart_time%set(restart_date(1), restart_date(2), restart_date(3), &
                              restart_date(4), restart_date(5), restart_date(6))

        ! find the time step that most closely matches the requested restart time (<=)
        restart_step = find_timestep_in_file(restart_file, 'time', restart_time, time_at_step)

        ! check to see if we actually udpated the restart date and print if in a more verbose mode
        if (options%debug) then
            if (restart_time /= time_at_step) then
                write(*,*) " updated restart date: ", trim(time_at_step%as_string())
            endif
        endif

        restart_time = time_at_step

        if (options%debug) then
            write(*,*) " ------------------ "
            write(*,*) "RESTART INFORMATION"
            write(*,*) "mjd",         restart_time%mjd()
            write(*,*) "date:",       trim(restart_time%as_string())
            write(*,*) "date",        restart_date
            write(*,*) "file",   trim(restart_file)
            write(*,*) "forcing step",restart_step
            write(*,*) " ------------------ "
        endif

        ! save the parameters in the master options structure
        options%restart_step_in_file = restart_step
        options%restart_file         = restart_file
        options%restart_date         = restart_date
        options%restart_time         = restart_time

        if (options%debug) write(*,*) " step in restart file",options%restart_step_in_file

    end subroutine init_restart_options


!   read physics options to use from a namelist file
    subroutine physics_namelist(filename,options)
        implicit none
        character(len=*),intent(in) :: filename
        type(options_type), intent(inout) :: options
        integer :: name_unit

!       variables to be used in the namelist
        integer::pbl,lsm,water,mp,rad,conv,adv,wind
!       define the namelist
        namelist /physics/ pbl,lsm,water,mp,rad,conv,adv,wind

!       default values for physics options (advection+linear winds+simple_microphysics)
        pbl = 0 ! 0 = no PBL,
                ! 1 = Stupid PBL (only in LSM=1 module),
                ! 2 = local PBL diffusion after Louis 1979
                ! 3 = YSU PBL (not complete)

        lsm = 0 ! 0 = no LSM,
                ! 1 = Fluxes from GCM,
                ! 2 = simple LSM, (not complete)
                ! 3 = Noah LSM

        water = 0 ! 0 = no open water fluxes,
                ! 1 = Fluxes from GCM, (needs lsm=1)
                ! 2 = Simple fluxes (needs SST in forcing data)

        mp  = 1 ! 0 = no MP,
                ! 1 = Thompson et al (2008),
                ! 2 = "Linear" microphysics

        rad = 0 ! 0 = no RAD,
                ! 1 = Surface fluxes from GCM, (radiative cooling ~1K/day in LSM=1 module),
                ! 2 = cloud fraction based radiation + radiative cooling

        conv= 0 ! 0 = no CONV,
                ! 1 = Tiedke scheme
                ! 2 = Kain-Fritsch scheme

        adv = 1 ! 0 = no ADV,
                ! 1 = upwind advection scheme
                ! 2 = MPDATA

        wind= 1 ! 0 = no LT,
                ! 1 = linear theory wind perturbations

!       read the namelist
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=physics)
        close(name_unit)

!       store options
        options%physics%boundarylayer=pbl
        options%physics%convection=conv
        options%physics%advection=adv
        options%physics%landsurface=lsm
        options%physics%watersurface=water
        options%physics%microphysics=mp
        options%physics%radiation=rad
        options%physics%windtype=wind

    end subroutine physics_namelist

    subroutine var_namelist(filename,options)
        implicit none
        character(len=*),intent(in) :: filename
        type(options_type), intent(inout) :: options
        integer :: name_unit
        character(len=MAXVARLENGTH) :: landvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar,zbvar,  &
                                        hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,           &
                                        pvar,pbvar,tvar,qvvar,qcvar,qivar,hgtvar,shvar,lhvar,pblhvar,   &
                                        soiltype_var, soil_t_var,soil_vwc_var,soil_deept_var,           &
                                        vegtype_var,vegfrac_var, linear_mask_var, nsq_calibration_var,  &
                                        swdown_var, lwdown_var, sst_var, rain_var, time_var,            &
                                        qrvar, qsvar, qgvar, qnivar, qnrvar, nvar

        namelist /var_list/ pvar,pbvar,tvar,qvvar,qcvar,qivar,hgtvar,shvar,lhvar,pblhvar,   &
                            landvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar,zbvar, &
                            hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,           &
                            soiltype_var, soil_t_var,soil_vwc_var,soil_deept_var,           &
                            vegtype_var,vegfrac_var, linear_mask_var, nsq_calibration_var,  &
                            swdown_var, lwdown_var, sst_var, rain_var, time_var,            &
                            qrvar, qsvar, qgvar, qnivar, qnrvar, nvar

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
        pvar=""
        pbvar=""
        tvar=""
        qvvar=""
        qcvar=""
        qivar=""
        qrvar=""         ! jh - added as optional field
        qsvar=""         ! jh - added as optional field
        qgvar=""         ! jh - added as optional field
        qnivar=""        ! jh - added as optional field
        qnrvar=""        ! jh - added as optional field
        nvar=""          ! jh - added as optional field
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
        lat_hi=""
        lon_hi=""
        ulat_hi=""
        ulon_hi=""
        vlat_hi=""
        vlon_hi=""
        soiltype_var=""
        soil_t_var=""
        soil_vwc_var=""
        soil_deept_var=""
        vegtype_var=""
        vegfrac_var=""
        linear_mask_var=""
        nsq_calibration_var=""
        rain_var=""

        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=var_list)
        close(name_unit)

        ! 2D geometry variable names (for coarse model)
        options%hgtvar      = hgtvar
        options%latvar      = latvar
        options%lonvar      = lonvar
        options%time_var    = time_var

        ! U varname and associated lat/lon var names
        options%uvar        = uvar
        if (ulat=="") ulat=latvar
        if (ulon=="") ulon=lonvar
        options%ulat        = ulat
        options%ulon        = ulon

        ! V varname and associated lat/lon var names
        options%vvar        = vvar
        if (vlat=="") vlat=latvar
        if (vlon=="") vlon=lonvar
        options%vlat        = vlat
        options%vlon        = vlon

        ! Primary model variable names
        options%pbvar       = pbvar
        options%pvar        = pvar
        options%tvar        = tvar
        options%qvvar       = qvvar
        options%qcvar       = qcvar
        options%qivar       = qivar
        
        options%qrvar       = qrvar    ! jh - added as optional field
        options%qsvar       = qsvar    ! jh - added as optional field
        options%qgvar       = qgvar    ! jh - added as optional field
        options%qnivar      = qnivar   ! jh - added as optional field
        options%qnrvar      = qnrvar   ! jh - added as optional field
        
        options%nvar        = nvar     ! jh - added as optional field
        
        ! vertical coordinate
        options%zvar        = zvar
        options%zbvar       = zbvar

        ! 2D model variables (e.g. Land surface and PBL height)
        options%shvar       = shvar
        options%lhvar       = lhvar
        options%pblhvar     = pblhvar

        ! Shortwave and longwave down at the surface
        options%swdown_var  = swdown_var
        options%lwdown_var  = lwdown_var

        ! Sea surface temperature
        options%sst_var     = sst_var
        options%rain_var    = rain_var

        !------------------------------------------------------
        ! variable names for the high resolution domain
        options%hgt_hi          = hgt_hi
        options%landvar         = landvar
        options%lat_hi          = lat_hi
        options%lon_hi          = lon_hi
        options%ulat_hi         = ulat_hi
        options%ulon_hi         = ulon_hi
        options%vlat_hi         = vlat_hi
        options%vlon_hi         = vlon_hi

        ! soil and vegetation parameters
        options%soiltype_var    = soiltype_var
        options%soil_t_var      = soil_t_var
        options%soil_vwc_var    = soil_vwc_var
        options%soil_deept_var  = soil_deept_var
        options%vegtype_var     = vegtype_var
        options%vegfrac_var     = vegfrac_var

        ! optional calibration variables for linear wind solution
        options%linear_mask_var     = linear_mask_var
        options%nsq_calibration_var = nsq_calibration_var
    end subroutine var_namelist

    subroutine parameters_namelist(filename,options)
        implicit none
        character(len=*),intent(in) :: filename
        type(options_type), intent(inout) :: options
        integer :: name_unit

        real    :: dx, dxlow, outputinterval, inputinterval, t_offset, smooth_wind_distance
        real    :: cfl_reduction_factor
        integer :: ntimesteps
        type(time_delta_t) :: dt
        integer :: nz, n_ext_winds,buffer, warning_level, cfl_strictness
        logical :: ideal, readz, readdz, interactive, debug, external_winds, surface_io_only, &
                   mean_winds, mean_fields, restart, advect_density, z_is_geopotential, z_is_on_interface,&
                   high_res_soil_state, use_agl_height, time_varying_z, t_is_potential, qv_is_spec_humidity, &
                   qv_is_relative_humidity, &
                   use_mp_options, use_lt_options, use_adv_options, use_lsm_options, use_bias_correction, &
                   use_block_options

        character(len=MAXFILELENGTH) :: date, calendar, start_date, forcing_start_date, end_date
        integer :: year, month, day, hour, minute, second
        character(len=MAXFILELENGTH) :: mp_options_filename, lt_options_filename, &
                                        adv_options_filename, lsm_options_filename, &
                                        bias_options_filename, block_options_filename

        namelist /parameters/ ntimesteps,outputinterval,inputinterval, surface_io_only, &
                              dx,dxlow,ideal,readz,readdz,nz,t_offset,debug, interactive, &
                              external_winds,buffer,n_ext_winds,advect_density,smooth_wind_distance, &
                              mean_winds,mean_fields,restart, z_is_geopotential, z_is_on_interface,&
                              date, calendar, high_res_soil_state,warning_level, t_is_potential,  &
                              qv_is_relative_humidity, qv_is_spec_humidity,  &
                              use_agl_height, start_date, forcing_start_date, end_date, time_varying_z, &
                              cfl_reduction_factor, cfl_strictness,         &
                              mp_options_filename,      use_mp_options,     &
                              block_options_filename,   use_block_options,  &
                              lt_options_filename,      use_lt_options,     &
                              lsm_options_filename,     use_lsm_options,    &
                              adv_options_filename,     use_adv_options,    &
                              bias_options_filename,    use_bias_correction

!       default parameters
        surface_io_only=.False.
        mean_fields=.False.
        mean_winds=.False.
        external_winds=.False.
        n_ext_winds=1
        t_offset=(-9999)
        buffer=0
        advect_density=.False.
        t_is_potential=.True.
        qv_is_spec_humidity=.False.
        qv_is_relative_humidity=.False.
        z_is_geopotential=.False.
        z_is_on_interface=.False.
        dxlow=100000
        restart=.False.
        ideal=.False.
        debug=.False.
        interactive=.False.
        warning_level=-9999
        readz=.False.
        readdz=.True.
        nz  = MAXLEVELS
        smooth_wind_distance=-9999
        calendar="gregorian"
        high_res_soil_state=.False.
        use_agl_height=.False.
        start_date=""
        forcing_start_date=""
        end_date=""
        time_varying_z=.False.
        cfl_reduction_factor = 0.9
        cfl_strictness = 3
        inputinterval = 3600
        outputinterval = 3600

        ! flag set to read specific parameterization options
        use_mp_options=.False.
        mp_options_filename = filename

        use_lt_options=.False.
        lt_options_filename = filename

        use_adv_options=.False.
        adv_options_filename = filename

        use_lsm_options=.False.
        lsm_options_filename = filename

        use_bias_correction=.False.
        bias_options_filename = filename

        use_block_options=.False.
        block_options_filename = filename

        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=parameters)
        close(name_unit)

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
            write(*,*) " Default smoothing distance = lowdx*2 = ", smooth_wind_distance
        endif

        options%t_offset=t_offset
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

        options%in_dt      = inputinterval
        call options%input_dt%set(seconds=inputinterval)
        options%out_dt     = outputinterval
        call options%output_dt%set(seconds=outputinterval)
        ! if outputing at half-day or longer intervals, create monthly files
        if (outputinterval>=43200) then
            options%output_file_frequency="monthly"
        ! if outputing at half-hour or longer intervals, create daily files
        else if (outputinterval>=300) then
            options%output_file_frequency="daily"
        ! otherwise create a new output file every timestep
        else
            options%output_file_frequency="every step"
        endif
        options%surface_io_only = surface_io_only


        options%calendar=calendar
        call time_init(calendar)
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
        options%dx = dx
        options%dxlow = dxlow
        options%ideal = ideal
        if (ideal) then
            write(*,*) " Running Idealized simulation (time step does not advance)"
        endif
        options%readz = readz
        options%readdz = readdz
        options%buffer = buffer
        options%mean_winds = mean_winds
        options%mean_fields = mean_fields
        options%advect_density = advect_density
        options%debug = debug
        options%interactive = interactive
        options%warning_level = warning_level
        options%use_agl_height = use_agl_height
        options%t_is_potential = t_is_potential
        options%qv_is_relative_humidity = qv_is_relative_humidity
        options%qv_is_spec_humidity = qv_is_spec_humidity
        options%z_is_geopotential = z_is_geopotential
        options%z_is_on_interface = z_is_on_interface

        options%external_winds = external_winds
        options%ext_winds_nfiles = n_ext_winds
        options%restart = restart

        options%nz = nz

        options%high_res_soil_state = high_res_soil_state
        options%time_varying_z = time_varying_z
        options%cfl_reduction_factor = cfl_reduction_factor
        options%cfl_strictness = cfl_strictness


        options%use_mp_options = use_mp_options
        options%mp_options_filename  = mp_options_filename   ! trude added

        options%use_lt_options = use_lt_options
        options%lt_options_filename  = lt_options_filename

        options%use_adv_options = use_adv_options
        options%adv_options_filename  = adv_options_filename

        options%use_lsm_options = use_lsm_options
        options%lsm_options_filename  = lsm_options_filename

        options%use_bias_correction = use_bias_correction
        options%bias_options_filename  = bias_options_filename

        options%use_block_options = use_block_options
        options%block_options_filename  = block_options_filename

        ! options are updated when complete
    end subroutine parameters_namelist

    subroutine mp_parameters_namelist(mp_filename,options)
        implicit none
        character(len=*),intent(in) :: mp_filename
        type(options_type), intent(inout) :: options
        integer :: name_unit

        real :: Nt_c, TNO, am_s,rho_g,av_s,bv_s,fv_s,av_g,bv_g,av_i,Ef_si,Ef_rs,Ef_rg,Ef_ri
        real :: C_cubes, C_sqrd, mu_r, t_adjust
        logical :: Ef_rw_l, EF_sw_l
        integer :: top_mp_level
        real :: local_precip_fraction
        integer :: update_interval

        namelist /mp_parameters/ Nt_c,TNO, am_s, rho_g, av_s,bv_s,fv_s,av_g,bv_g,av_i,Ef_si,Ef_rs,Ef_rg,Ef_ri,&     ! trude added Nt_c, TNO
                              C_cubes,C_sqrd, mu_r, Ef_rw_l, Ef_sw_l, t_adjust, top_mp_level, local_precip_fraction, update_interval

        ! because mp_options could be in a separate file (shoudl probably set all namelists up to have this option)
        if (options%use_mp_options) then
            if (trim(mp_filename)/=trim(get_options_file())) then
                call version_check(mp_filename,options)
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

        update_interval = 0 ! update every time step
        top_mp_level = 0    ! if <=0 just use the actual model top
        local_precip_fraction = 1 ! if <1: the remaining fraction (e.g. 1-x) of precip gets distributed to the surrounding grid cells
        ! read in the namelist
        if (options%use_mp_options) then
            open(io_newunit(name_unit), file=mp_filename)
            read(name_unit,nml=mp_parameters)
            close(name_unit)
        endif

        ! store the data back into the mp_options datastructure
        options%mp_options%Nt_c = Nt_c
        options%mp_options%TNO = TNO
        options%mp_options%am_s = am_s
        options%mp_options%rho_g = rho_g
        options%mp_options%av_s = av_s
        options%mp_options%bv_s = bv_s
        options%mp_options%fv_s = fv_s
        options%mp_options%av_g = av_g
        options%mp_options%bv_g = bv_g
        options%mp_options%av_i = av_i
        options%mp_options%Ef_si = Ef_si
        options%mp_options%Ef_rs = Ef_rs
        options%mp_options%Ef_rg = Ef_rg
        options%mp_options%Ef_ri = Ef_ri
        options%mp_options%mu_r = mu_r
        options%mp_options%t_adjust = t_adjust
        options%mp_options%C_cubes = C_cubes
        options%mp_options%C_sqrd = C_sqrd
        options%mp_options%Ef_rw_l = Ef_rw_l
        options%mp_options%Ef_sw_l = Ef_sw_l

        options%mp_options%update_interval = update_interval
        if (top_mp_level < 0) top_mp_level = options%nz + top_mp_level
        options%mp_options%top_mp_level = top_mp_level
        options%mp_options%local_precip_fraction = local_precip_fraction
    end subroutine mp_parameters_namelist

    subroutine block_parameters_namelist(filename, options)
        implicit none

        character(len=*),   intent(in)   :: filename
        type(options_type), intent(inout):: options

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
        if (options%use_block_options) then
            if (trim(filename)/=trim(get_options_file())) then
                call version_check(filename,options)
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
        if (options%use_block_options) then
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


    subroutine lt_parameters_namelist(filename, options)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_type), intent(inout):: options

        integer :: name_unit

        integer :: vert_smooth
        logical :: variable_N               ! Compute the Brunt Vaisala Frequency (N^2) every time step
        logical :: N_from_forcing           ! If true N is calculated from the forcing data at every forcing timestep instead of using the atmosperic fields of ICAR (standard: False)
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
        namelist /lt_parameters/ variable_N, N_from_forcing, smooth_nsq, buffer, stability_window_size, max_stability, min_stability, &
                                 linear_contribution, linear_update_fraction, N_squared, vert_smooth, &
                                 remove_lowres_linear, rm_N_squared, rm_linear_contribution, &
                                 spatial_linear_fields, linear_mask, nsq_calibration, minimum_layer_size, &
                                 dirmax, dirmin, spdmax, spdmin, nsqmax, nsqmin, n_dir_values, n_nsq_values, n_spd_values, &
                                 read_LUT, write_LUT, u_LUT_Filename, v_LUT_Filename, overwrite_lt_lut, LUT_Filename

         ! because lt_options could be in a separate file
         if (options%use_lt_options) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options)
             endif
         endif


        ! set default values
        variable_N = .True.
        N_from_forcing = .False.       ! standard is set to false to preserve original behaviour
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
        if (options%use_lt_options) then
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
            opt%N_from_forcing = N_from_forcing
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

            opt%u_LUT_Filename = u_LUT_Filename
            opt%v_LUT_Filename = v_LUT_Filename
            opt%overwrite_lt_lut = overwrite_lt_lut

        end associate
        
        if (options%debug) then
                if (options%lt_options%N_from_forcing) then
                        write(*,*) "DEBUG"
                        write(*,*) "DEBUG experimental linear theory option activated"
                        write(*,*) "DEBUG N_from_forcing is set to true. Calculating N from forcing."
                        if (.NOT. options%lt_options%variable_N) then
                                write(*,*) "DEBUG however, setting is useless since variable_N is set to false!"
                        endif
                        write(*,*) "DEBUG"
                endif
        endif

    end subroutine lt_parameters_namelist


    subroutine adv_parameters_namelist(filename, options)
        implicit none
        character(len=*), intent(in) :: filename
        type(options_type), intent(inout)::options
        type(adv_options_type)::adv_options

        integer :: name_unit

        logical :: boundary_buffer          ! apply some smoothing to the x and y boundaries in MPDATA
        logical :: flux_corrected_transport ! use the flux corrected transport option in MPDATA
        integer :: mpdata_order             ! MPDATA order of correction (e.g. 1st=upwind, 2nd=classic, 3rd=better)
        ! jhorak: options to choose the boundary condition used at the model top added.
        ! option values and the associated boundary condition are:
        ! 0 ... zero gradient (default behaviour if not set)
        ! 1 ... constant gradient BC (or adjusted correspondingly if this would lead to a negative value)
        ! 2 ... zero value 
        ! 3 ... zero value but on the flux divergence in the topmost vertical layer for downdrafts
        integer :: bc_top                   ! type of boundary condition to apply to microphysics fields (qv,qc,qi,qr,qs,qg) at the model top
        integer :: bc_th_top                ! type of boundary condition to apply to potential temperature at the model




        ! define the namelist
        namelist /adv_parameters/ boundary_buffer, flux_corrected_transport, mpdata_order, bc_top, bc_th_top

         ! because adv_options could be in a separate file
         if (options%use_adv_options) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options)
             endif
         endif


        ! set default values
        boundary_buffer = .False.
        flux_corrected_transport = .True.
        mpdata_order = 2
        bc_top = 0      
        bc_th_top = 0

        ! read the namelist options
        if (options%use_adv_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=adv_parameters)
            close(name_unit)
        endif

        ! store everything in the adv_options structure
        adv_options%boundary_buffer = boundary_buffer
        adv_options%flux_corrected_transport = flux_corrected_transport
        adv_options%mpdata_order = mpdata_order
        adv_options%bc_top = bc_top
        adv_options%bc_th_top = bc_th_top

        if (options%debug) then
                if (adv_options%bc_top > 0) then
                        write(*,*) "DEBUG"
                        write(*,*) "DEBUG experimental boundary condition applied to microphysics fields at model"
                        if (adv_options%bc_top == 1) then
                                write(*,*) "DEBUG constant gradient boundary condition for microphysics fields"
                        endif
                        if (adv_options%bc_top == 2) then
                                write(*,*) "DEBUG zero value on inflow, zero-gradient on outflow for microphysics fields"
                        endif
                        write(*,*) "DEBUG"
                endif
                if (adv_options%bc_th_top > 0) then
                        write(*,*) "DEBUG"
                        write(*,*) "DEBUG experimental boundary condition applied at model"
                        if (adv_options%bc_th_top == 1) then
                                write(*,*) "DEBUG constant gradient boundary condition for potential temperature"
                        endif
                        write(*,*) "DEBUG"
                endif
        endif

        ! copy the data back into the global options data structure
        options%adv_options = adv_options
    end subroutine adv_parameters_namelist


    subroutine set_default_LU_categories(urban_category, ice_category, water_category, LU_Categories)
        ! if various LU categories were not defined in the namelist (i.e. they == -1) then attempt
        ! to define default values for them based on the LU_Categories variable supplied.
        implicit none
        integer, intent(inout) :: urban_category, ice_category, water_category
        character(len=MAXVARLENGTH), intent(in) :: LU_Categories

        if (trim(LU_Categories)=="MODIFIED_IGBP_MODIS_NOAH") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15
            if (water_category==-1) water_category = 17
        elseif (trim(LU_Categories)=="USGS") then
            if (urban_category==-1) urban_category = 1
            if (ice_category==-1)   ice_category = -1
            if (water_category==-1) water_category = 16
        elseif (trim(LU_Categories)=="USGS-RUC") then
            if (urban_category==-1) urban_category = 1
            if (ice_category==-1)   ice_category = 24
            if (water_category==-1) water_category = 16
            ! also note, lakes_category = 28
            write(*,*) "WARNING: not handling lake category (28)"
        elseif (trim(LU_Categories)=="MODI-RUC") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15
            if (water_category==-1) water_category = 17
            ! also note, lakes_category = 21
            write(*,*) "WARNING: not handling lake category (21)"
        elseif (trim(LU_Categories)=="NLCD40") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15 ! and 22?
            if (water_category==-1) water_category = 17 ! and 21
            write(*,*) "WARNING: not handling all varients of categories (e.g. permanent_snow=15 is, but permanent_snow_ice=22 is not)"
        endif

    end subroutine set_default_LU_categories


    subroutine bias_parameters_namelist(filename, options)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_type), intent(inout):: options

        type(bias_options_type)     ::  bias_options

        integer :: name_unit

        character(len=MAXFILELENGTH):: bias_correction_filename ! file containing bias correction data
        character(len=MAXVARLENGTH) :: rain_fraction_var        ! name of variable containing the fraction to multiply rain by

        ! define the namelist
        namelist /bias_parameters/ bias_correction_filename, rain_fraction_var

         ! because adv_options could be in a separate file
         if (options%use_bias_correction) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options)
             endif
         endif


        ! set default values
        bias_correction_filename = options%init_conditions_file
        rain_fraction_var        = "rain_fraction"

        ! read the namelist options
        if (options%use_bias_correction) then
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


    subroutine lsm_parameters_namelist(filename, options)
        implicit none
        character(len=*), intent(in) :: filename
        type(options_type), intent(inout)::options
        type(lsm_options_type)::lsm_options

        integer :: name_unit

        character(len=MAXVARLENGTH) :: LU_Categories ! Category definitions (e.g. USGS, MODIFIED_IGBP_MODIS_NOAH)
        logical :: monthly_vegfrac                   ! read in 12 months of vegfrac data
        integer :: update_interval                   ! minimum number of seconds between LSM updates
        integer :: urban_category                    ! index that defines the urban category in LU_Categories
        integer :: ice_category                      ! index that defines the ice category in LU_Categories
        integer :: water_category                    ! index that defines the water category in LU_Categories

        ! define the namelist
        namelist /lsm_parameters/ LU_Categories, update_interval, monthly_vegfrac, &
                                  urban_category, ice_category, water_category

         ! because adv_options could be in a separate file
         if (options%use_lsm_options) then
             if (trim(filename)/=trim(get_options_file())) then
                 call version_check(filename,options)
             endif
         endif


        ! set default values
        LU_Categories   = "MODIFIED_IGBP_MODIS_NOAH"
        update_interval = 300 ! 5 minutes
        monthly_vegfrac = .False.

        ! default values for these will be set after reading LU_Categories
        urban_category  = -1
        ice_category    = -1
        water_category  = -1

        ! read the namelist options
        if (options%use_lsm_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,nml=lsm_parameters)
            close(name_unit)
        endif

        call set_default_LU_categories(urban_category, ice_category, water_category, LU_Categories)

        ! store everything in the lsm_options structure
        lsm_options%LU_Categories   = LU_Categories
        lsm_options%monthly_vegfrac = monthly_vegfrac
        lsm_options%update_interval = update_interval
        lsm_options%urban_category  = urban_category
        lsm_options%ice_category    = ice_category
        lsm_options%water_category  = water_category

        ! copy the data back into the global options data structure
        options%lsm_options = lsm_options
    end subroutine lsm_parameters_namelist


    ! set up model levels, either read from a namelist, or from a default set of values
    subroutine model_levels_namelist(filename,options)
        implicit none
        character(len=*), intent(in) :: filename
        type(options_type), intent(inout) :: options

        integer :: name_unit, this_level
        real, allocatable, dimension(:) :: dz_levels
        real, dimension(45) :: fulldz

        namelist /z_info/ dz_levels
        this_level=1

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
            allocate(options%dz_levels(options%nz))

            options%dz_levels(1:options%nz)=dz_levels(1:options%nz)

            if (minval(options%dz_levels)<1) then
                print*, "NB: gfortran doesn't read namelist arrays on multiple lines (check dz_levels)"
                stop "ERROR: model levels must be > 1m vertical"
            endif

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
            allocate(options%dz_levels(options%nz))
            options%dz_levels = fulldz(1:options%nz)
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
        write(*,*) "  Boundary conditions files to be used:"
        if (nfiles>10) then
            write(*,*) "    nfiles=", trim(str(nfiles)), ", too many to print."
            write(*,*) "    First file:", trim(forcing_files(1))
            write(*,*) "    Last file: ", trim(forcing_files(nfiles))
        else
            do i=1,nfiles
                write(*,*) "      ",trim(forcing_files(i))
            enddo
        endif

    end function read_forcing_file_names

    subroutine filename_namelist(filename, options)
        ! read in filenames from up to two namelists
        ! init_conditions_file = high res grid data
        ! boundary_files       = list of files for boundary conditions (low-res)
        ! ext_wind_files       = list of files to read wind data on the high-res grid (optional)
        ! linear_mask_file     = file to read a high-res fractional mask for the linear perturbation
        implicit none
        character(len=*), intent(in) :: filename
        type(options_type), intent(inout) :: options

        character(len=MAXFILELENGTH) :: init_conditions_file, output_file, forcing_file_list, &
                                        linear_mask_file, nsq_calibration_file
        character(len=MAXFILELENGTH), allocatable :: boundary_files(:), ext_wind_files(:)
        integer :: name_unit, nfiles, i

        ! set up namelist structures
        namelist /files_list/ init_conditions_file, output_file, boundary_files, forcing_file_list, &
                              linear_mask_file, nsq_calibration_file
        namelist /ext_winds_info/ ext_wind_files

        linear_mask_file="MISSING"
        nsq_calibration_file="MISSING"
        allocate(boundary_files(MAX_NUMBER_FILES))
        boundary_files="MISSING"
        forcing_file_list="MISSING"

        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=files_list)
        close(name_unit)

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
        options%boundary_files=boundary_files(1:nfiles)
        deallocate(boundary_files)

        options%output_file=output_file
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
            options%ext_wind_files=ext_wind_files
            deallocate(ext_wind_files)
        endif
    end subroutine filename_namelist

end module initialize_options
