module initialize_options
    use data_structures
    use io_routines,                only : io_nearest_time_step, io_newunit
    use model_tracking,             only : print_model_diffs
    use time,                       only : date_to_mjd, parse_date, time_init

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
        
        call version_check(options_filename,options)
        call physics_namelist(options_filename,options)
        call var_namelist(options_filename,options)
        call parameters_namelist(options_filename,options)
        call model_levels_namelist(options_filename, options)
        
        ! ++ trude
        !       read in mp_options file
        write(*,*) "Initializing mp_Options"
        call init_mp_options(options%mp_options_filename,options)
        ! -- trude
        
        if (options%restart) then
            call init_restart_options(options_filename,options)
        endif

        call filename_namelist(options_filename, options)
        ! check for any inconsistencies in the options requested
        call options_check(options)
    end subroutine init_options


    function get_options_file()
        implicit none
        character(len=MAXFILELENGTH) ::get_options_file
        integer :: error
        logical :: file_exists
    
        if (command_argument_count()>0) then
            call get_command_argument(1,get_options_file, status=error)
            if (error>0) then
                get_options_file="icar_options.nml"
            elseif (error==-1) then
                write(*,*) "Options filename = ", trim(get_options_file), " ...<cutoff>"
                write(*,*) "Maximum filename length = ", MAXFILELENGTH
                stop("ERROR: options filename too long")
            endif
        else
            get_options_file="icar_options.nml"
        endif
        write(*,*) "Options filename = ", trim(get_options_file)
        INQUIRE(file=trim(get_options_file), exist=file_exists)
        if (.not.file_exists) then
            stop("Options file does not exist. ")
        endif
    end function


! ++ trude  
    subroutine init_mp_options(mp_options_filename,options)
!       reads a series of options from a namelist file and stores them in the 
!       options data structure
        implicit none
        character(len=*), intent(in) :: mp_options_filename
        type(options_type), intent(inout) :: options
        
        if (options%use_mp_options) then
            call version_check(mp_options_filename,options)
        endif
        call parameters_mp_namelist(mp_options_filename,options)
! NBNBNB trude, check if I can comment these 3 lines out        
!       if (options%restart) then
!           call init_restart(options_filename,options)
!       endif
    end subroutine init_mp_options
! -- trude


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
        if (version.ne."0.8.1") then
            write(*,*) "Model version does not match namelist version"
            write(*,*) "  Model version: 0.8.1"
            write(*,*) "  Namelist version: ",trim(version)
            call print_model_diffs(version)
            stop
        endif
        options%version=version
        options%comment=comment
        write(*,*) "Model version: ",trim(version)
    end subroutine version_check

    subroutine options_check(options)
        ! Minimal error checking on option settings
        implicit none
        type(options_type), intent(inout)::options
    
        if (options%t_offset.eq.(-9999)) then
            if (options%warning_level>0) then
                write(*,*), "WARNING, WARNING, WARNING"
                write(*,*), "WARNING, Using default t_offset=300"
                write(*,*), "WARNING, WARNING, WARNING"
            endif
            options%t_offset=300
        endif

        ! convection can modify wind field, and ideal doesn't rebalance winds every timestep
        if ((options%physics%convection.ne.0).and.(options%ideal)) then
            if (options%warning_level>1) then
                write(*,*) "WARNING WARNING WARNING"
                write(*,*) "WARNING, Running convection in ideal mode may be bad..."
                write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%warning_level==10) then
                stop
            endif
        endif
        if ((options%physics%landsurface>0).and.(options%physics%boundarylayer==0)) then
            if (options%warning_level>0) then
                write(*,*) "WARNING WARNING WARNING"
                write(*,*) "WARNING, Running LSM without PBL may overheat the surface and CRASH the model. "
                write(*,*) "WARNING WARNING WARNING"
            endif
            if (options%warning_level>=5) then
                write(*,*) "Set warning_level<5 to continue"
                stop
            endif
        endif
    
    end subroutine options_check
    
    subroutine init_restart_options(filename, options)
        character(len=*), intent(in) :: filename
        type(options_type),intent(inout) :: options
        
        character(len=MAXFILELENGTH) :: init_conditions_file, output_file,restart_file
        integer :: restart_step, restart_date(6), name_unit
        real*8 :: restart_mjd
        
        namelist /restart_info/ restart_step, restart_file, restart_date
        
        restart_date=[-999,-999,-999,-999,-999,-999]
        restart_step=-999
        
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=restart_info)
        close(name_unit)
        
        options%restart_step=restart_step
        options%restart_file=restart_file
        options%restart_date=restart_date
        
        
        restart_mjd = date_to_mjd(restart_date(1), restart_date(2), restart_date(3), &
                                  restart_date(4), restart_date(5), restart_date(6))
        write(*,*) "mjd",restart_mjd
        write(*,*) "date",restart_date
        write(*,*) "file",trim(restart_file)
        write(*,*) "forcing step",restart_step
        
        if (restart_step==-999) then
            restart_step=FLOOR((restart_mjd - options%initial_mjd + 1e-6) / (options%in_dt/86400.0)) + 1
            options%restart_step=restart_step
            write(*,*) "updated forcing step",restart_step
        endif
        ! in case the supplied restart date doesn't line up with an input forcing step, recalculate
        ! the restart date (mjd) based off the nearest input step
        restart_mjd = options%initial_mjd + (restart_step-1)*(options%in_dt/86400.0)
        write(*,*) "updated mjd",restart_mjd
        
        ! now find the closest previous output step to the current restart date
        options%restart_step_in_file=io_nearest_time_step(restart_file, restart_mjd)
        write(*,*) "step in restart file",options%restart_step_in_file
        
    end subroutine init_restart_options
    
    
!   read physics options to use from a namelist file
    subroutine physics_namelist(filename,options)
        implicit none
        character(len=*),intent(in) :: filename
        type(options_type), intent(inout) :: options
        integer :: name_unit
        
!       variables to be used in the namelist
        integer::pbl,lsm,mp,rad,conv,adv,wind
!       define the namelist
        namelist /physics/ pbl,lsm,mp,rad,conv,adv,wind
        
!       default values for physics options (advection+linear winds+simple_microphysics)
        pbl=2   ! 0 = no PBL, 1 = STUPID PBL (in LSM module), 2 = local PBL diffusion
        lsm=3   ! 0 = no LSM, 1 = Fluxes from GCM, 2 = simple LSM, 3 = Noah LSM
        mp=2    ! 0 = no MP,  1 = Thompson et al (2008), 2 = "Linear" microphysics
        rad=2   ! 0 = no RAD, 1 = radiative cooling 1K/day (in LSM), 2 = cloud fraction based radiation
        conv=0  ! 0 = no CONV,1 = Tiedke scheme
        adv=1   ! 0 = no ADV, 1 = upwind advection scheme
        wind=1  ! 0 = no LT,  1 = linear theory wind perturbations
        
!       read the namelist
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=physics)
        close(name_unit)

!       store options
        options%physics%boundarylayer=pbl
        options%physics%convection=conv
        options%physics%advection=adv
        options%physics%landsurface=lsm
        options%physics%microphysics=mp
        options%physics%radiation=rad
        options%physics%windtype=wind
        
    end subroutine physics_namelist
    
    subroutine var_namelist(filename,options)
        implicit none
        character(len=*),intent(in) :: filename
        type(options_type), intent(inout) :: options
        integer :: name_unit
        character(len=MAXVARLENGTH) :: landvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar,&
                                        hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,     &
                                        pvar,pbvar,tvar,qvvar,qcvar,qivar,hgtvar,shvar,lhvar,pblhvar,&
                                        soiltype_var, soil_t_var,soil_vwc_var,soil_deept_var, &
                                        vegtype_var,vegfrac_var, linear_mask_var
                                        
        namelist /var_list/ pvar,pbvar,tvar,qvvar,qcvar,qivar,hgtvar,shvar,lhvar,pblhvar,&
                            landvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar, &
                            hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi, &
                            soiltype_var, soil_t_var,soil_vwc_var,soil_deept_var, &
                            vegtype_var,vegfrac_var, linear_mask_var
        
        hgtvar="HGT"
        latvar="XLAT"
        lonvar="XLONG"
        uvar="U"
        ulat="XLAT_U"
        ulon="XLONG_U"
        vvar="V"
        vlat="XLAT_V"
        vlon="XLONG_V"
        pvar="P"
        pbvar="PB"
        tvar="T"
        qvvar="QVAPOR"
        qcvar="QCLOUD"
        qivar="QICE"
        zvar="Z"
        shvar="HFLX"
        lhvar="LVFLX"
        pblhvar="PBLH"
        hgt_hi="HGT"
        landvar="XLAND"
        lat_hi="XLAT"
        lon_hi="XLONG"
        ulat_hi="XLAT_U"
        ulon_hi="XLONG_U"
        vlat_hi="XLAT_V"
        vlon_hi="XLONG_V"
        soiltype_var="" !"SOILTYPE"
        soil_t_var="" !"TSOIL"
        soil_vwc_var="" !"SOILSMC"
        soil_deept_var="" !"DEEPT"
        vegtype_var="" !"VEGTYPE"
        vegfrac_var="" !"VEGFRAC"
        linear_mask_var="data"
        
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=var_list)
        close(name_unit)

!       2D geometry variable names (for coarse model)
        options%hgtvar=hgtvar
        options%latvar=latvar
        options%lonvar=lonvar
!       U varname and associated lat/lon var names
        options%uvar=uvar
        options%ulat=ulat
        options%ulon=ulon
!       V varname and associated lat/lon var names
        options%vvar=vvar
        options%vlat=vlat
        options%vlon=vlon
!       Primary model variable names
        options%pbvar=pbvar
        options%pvar=pvar
        options%tvar=tvar
        options%qvvar=qvvar
        options%qcvar=qcvar
        options%qivar=qivar
!       vertical coordinate
        options%zvar=zvar
!       2D model variables (e.g. Land surface and PBL height)       
        options%shvar=shvar
        options%lhvar=lhvar
        options%pblhvar=pblhvar
        
!       separate variable names for the high resolution domain
        options%hgt_hi=hgt_hi
        options%landvar=landvar
        options%lat_hi=lat_hi
        options%lon_hi=lon_hi
        options%ulat_hi=ulat_hi
        options%ulon_hi=ulon_hi
        options%vlat_hi=vlat_hi
        options%vlon_hi=vlon_hi
        
!       soil and vegetation parameters
        options%soiltype_var=soiltype_var
        options%soil_t_var=soil_t_var
        options%soil_vwc_var=soil_vwc_var
        options%soil_deept_var=soil_deept_var
        options%vegtype_var=vegtype_var
        options%vegfrac_var=vegfrac_var
        
        options%linear_mask_var=linear_mask_var
    end subroutine var_namelist
    
    subroutine parameters_namelist(filename,options)
        implicit none
        character(len=*),intent(in) :: filename
        type(options_type), intent(inout) :: options
        integer :: name_unit
        
        real    :: dx, dxlow, outputinterval, inputinterval, t_offset, smooth_wind_distance
        real    :: rotation_scale_height, N_squared, rm_N_squared,linear_contribution, rm_linear_contribution, linear_update_fraction
        integer :: ntimesteps, nfiles, xmin, xmax, ymin, ymax, vert_smooth
        integer :: nz, n_ext_winds,buffer, warning_level
        logical :: ideal, readz, readdz, debug, external_winds, remove_lowres_linear, variable_N, &
                   mean_winds, mean_fields, restart, advect_density, high_res_soil_state, &
                   use_agl_height, spatial_linear_fields, time_varying_z, linear_mask, use_mp_options
        character(len=MAXFILELENGTH) :: date, calendar, start_date
        integer :: year, month, day, hour, minute, second
! ++ trude
        character(len=MAXFILELENGTH) :: mp_options_filename
! -- trude

        namelist /parameters/ ntimesteps,outputinterval,inputinterval,dx,dxlow,ideal,readz,readdz,nz,t_offset,debug,nfiles, &
                              external_winds,buffer,n_ext_winds,advect_density,smooth_wind_distance, &
                              remove_lowres_linear,mean_winds,mean_fields,restart,xmin,xmax,ymin,ymax,vert_smooth, &
                              date, calendar, high_res_soil_state,rotation_scale_height,warning_level, variable_N, &
                              N_squared,rm_N_squared,linear_contribution,rm_linear_contribution, use_agl_height,  &
                              spatial_linear_fields, linear_update_fraction, start_date, time_varying_z, linear_mask, &
                              mp_options_filename, use_mp_options    ! trude added
        
!       default parameters
        mean_fields=.False.
        mean_winds=.False.
        external_winds=.False.
        n_ext_winds=1
        t_offset=(-9999)
        buffer=0
        remove_lowres_linear=.False.
        advect_density=.True.
        restart=.False.
        ideal=.False.
        debug=.False.
        warning_level=-9999
        readz=.False.
        readdz=.True.
        xmin=   1
        ymin=   1
        xmax= (-1)
        ymax= (-1)
        nz  = MAXLEVELS
        smooth_wind_distance=-9999
        vert_smooth=2
        calendar="gregorian"
        high_res_soil_state=.True.
        rotation_scale_height=2000.0
        N_squared=6.e-5
        rm_N_squared=6.e-5
        variable_N=.False.
        linear_contribution=1.0
        rm_linear_contribution=1.0
        linear_update_fraction=0.25
        use_agl_height=.True.
        spatial_linear_fields=.False.
        start_date=""
        time_varying_z=.False.
        linear_mask=.False.
        use_mp_options=.False.

        mp_options_filename = 'mp_options.nml'    ! trude added
        
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=parameters)
        close(name_unit)
        
        if (warning_level==-9999) then
            if (debug) then
                warning_level=2
            else
                warning_level=5
            endif
        endif
            
        
        if (smooth_wind_distance.eq.(-9999)) then
            smooth_wind_distance=dxlow*3
            write(*,*), "Default smoothing distance = lowdx*3 = ", smooth_wind_distance
        endif
        
        options%t_offset=t_offset
        if (smooth_wind_distance<0) then 
            write(*,*) "Wind smoothing must be a positive number"
            write(*,*) "smooth_wind_distance = ",smooth_wind_distance
            stop
        endif
        options%smooth_wind_distance=smooth_wind_distance
        if (vert_smooth<0) then 
            write(*,*) "Vertical smoothing must be a positive integer"
            write(*,*) "vert_smooth = ",vert_smooth
            stop
        endif
        options%vert_smooth=vert_smooth
        
        options%nfiles=nfiles
        options%ntimesteps=ntimesteps
        options%in_dt=inputinterval
        options%out_dt=outputinterval
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
        call parse_date(date, year, month, day, hour, minute, second)
        call time_init(calendar)
        options%initial_mjd=date_to_mjd(year, month, day, hour, minute, second)
        if (start_date=="") then
            options%start_mjd=options%initial_mjd
        else
            call parse_date(start_date, year, month, day, hour, minute, second)
            options%start_mjd=date_to_mjd(year, month, day, hour, minute, second)
        endif
            
        options%time_zero=((options%initial_mjd-50000) * 86400.0)
        options%dx=dx
        options%dxlow=dxlow
        options%ideal=ideal
        if (ideal) then
            write(*,*) "Running Idealized simulation (time step does not advance)"
        endif
        options%readz=readz
        options%readdz=readdz
        options%buffer=buffer
        options%remove_lowres_linear=remove_lowres_linear
        options%mean_winds=mean_winds
        options%mean_fields=mean_fields
        options%advect_density=advect_density
        options%debug=debug
        options%warning_level=warning_level
        options%rotation_scale_height=rotation_scale_height
        options%use_agl_height=use_agl_height
        options%spatial_linear_fields=spatial_linear_fields
        
        options%external_winds=external_winds
        options%ext_winds_nfiles=n_ext_winds
        options%restart=restart
        
        options%nz=nz
        options%xmin=xmin
        options%xmax=xmax
        options%ymin=ymin
        options%ymax=ymax
        
        options%N_squared=N_squared
        options%rm_N_squared=rm_N_squared
        options%variable_N=variable_N
        options%linear_contribution=linear_contribution
        options%rm_linear_contribution=rm_linear_contribution
        options%linear_mask = linear_mask
        options%linear_update_fraction = linear_update_fraction

        options%high_res_soil_state=high_res_soil_state
        options%time_varying_z=time_varying_z
        
        options%use_mp_options = use_mp_options
        options%mp_options_filename  = mp_options_filename   ! trude added

    end subroutine parameters_namelist

! ++ trude
    subroutine parameters_mp_namelist(mp_filename,options)
        implicit none
        character(len=*),intent(in) :: mp_filename
        type(options_type), intent(inout) :: options
        integer :: name_unit

        real :: Nt_c, TNO, am_s,rho_g,av_s,bv_s,fv_s,av_g,bv_g,av_i,Ef_si,Ef_rs,Ef_rg,Ef_ri
        real :: C_cubes, C_sqrd, mu_r, t_adjust
        logical :: Ef_rw_l, EF_sw_l

        namelist /parameters/ Nt_c,TNO, am_s, rho_g, av_s,bv_s,fv_s,av_g,bv_g,av_i,Ef_si,Ef_rs,Ef_rg,Ef_ri,&     ! trude added Nt_c, TNO
                              C_cubes,C_sqrd, mu_r, Ef_rw_l, Ef_sw_l, t_adjust
!       default parameters
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
        
        if (options%use_mp_options) then
            open(io_newunit(name_unit), file=mp_filename)
            read(name_unit,nml=parameters)
            close(name_unit)
        endif
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
    end subroutine parameters_mp_namelist
! -- trude
    
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
            allocate(dz_levels(options%nz))
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
            options%dz_levels=fulldz(1:options%nz)
        endif
        
    end subroutine model_levels_namelist
    
    
    subroutine filename_namelist(filename, options)
        ! read in filenames from up to two namelists
        ! init_conditions_file = high res grid data
        ! boundary_files       = list of files for boundary conditions (low-res)
        ! ext_wind_files       = list of files to read wind data on the high-res grid (optional)
        ! linear_mask_file     = file to read a high-res fractional mask for the linear perturbation
        implicit none
        character(len=*), intent(in) :: filename
        type(options_type), intent(inout) :: options

        character(len=MAXFILELENGTH) :: init_conditions_file, output_file, linear_mask_file
        character(len=MAXFILELENGTH), allocatable :: boundary_files(:), ext_wind_files(:)
        integer :: name_unit
        
        ! set up namelist structures
        namelist /files_list/ init_conditions_file, output_file, boundary_files, linear_mask_file
        namelist /ext_winds_info/ ext_wind_files
        
        linear_mask_file="MISSING"
        allocate(boundary_files(options%nfiles))
        
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=files_list)
        close(name_unit)
        
        options%init_conditions_file=init_conditions_file
        
        allocate(options%boundary_files(options%nfiles))
        options%boundary_files=boundary_files
        deallocate(boundary_files)
        
        options%output_file=output_file
        if (trim(linear_mask_file)=="MISSING") then
            linear_mask_file = options%init_conditions_file
        endif
        options%linear_mask_file=linear_mask_file
        
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
