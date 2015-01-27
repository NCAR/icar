!> ----------------------------------------------------------------------------
!! 	Model Initialization includes allocating memory for boundary and domain 
!! 		data structures.  It reads all of the options from the namelist
!! 		file (or files).  It also reads in Lat/Lon and Terrain data.  This module
!!		also sets up geographic (and vertical) look uptables for the forcing data
!! 		Finally, there is a driver routine to initialize all model physics packages
!! 
!! 	This module was initially written to read WRF output files as input.
!! 		This should serve as a basis for any additional file types, and can be 
!! 		readily modified.  At some point this could be modified to check the 
!! 		type if input file being specified, and call an appropriate routine. 
!! 		e.g. if options%inputtype=="WRF" then call wrf_init/update()
!! 
!!   The module has been updated to allow arbitrary named variables
!!       this allows the use of e.g. ERAi, but still is not as flexible as it could be
!!       e.g. if a forcing doesn't contain one or more variable types it will not work. 
!!       Need to let varname=""
!!
!!   The use of various python wrapper scripts in helpers/ makes it easy to add new
!!       datasets, and make them conform to the expectations of the current system. 
!!		For now there are no plans to near term plans to substantially modify this. 
!!
!!	Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module init
	use data_structures
	use io_routines,				only : io_read2d, io_read2di, io_read3d, io_newunit, &
		 								   io_write3d,io_write3di, io_nearest_time_step
	use geo,						only : geo_LUT, geo_interp, geo_interp2d
	use vertical_interpolation,     only : vLUT
	use microphysics,				only : mp_init
	use convection,				    only : init_convection
	use planetary_boundary_layer,   only : pbl_init
	use radiation,					only : radiation_init
	use land_surface,				only : lsm_init
	use model_tracking,			    only : print_model_diffs
	use wind,						only : init_winds
	use time,						only : date_to_mjd, parse_date, time_init
	
	implicit none
	private
	public::init_model,init_physics
	
contains
	subroutine init_model(options,domain,boundary)
		implicit none
		type(options_type), intent(inout) :: options
		type(domain_type), intent(inout):: domain
		type(bc_type), intent(inout):: boundary
		character(len=MAXFILELENGTH) :: options_filename
		
		options_filename=get_options_file()
! 		read in options file
		write(*,*) "Initializing Options"
		call init_options(options_filename,options)
! 		allocate and initialize the domain
		write(*,*) "Initializing Domain"
		call init_domain(options,domain)
! 		allocate and initialize the boundary conditions structure (includes 3D grids too...)
!		this might be more apropriately though of as a forcing data structure (for low res model)
		write(*,*) "Initializing Boundaries"
		call init_bc(options,domain,boundary)
		write(*,*) "Finished basic initialization"
		
	end subroutine init_model
	
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
	
	subroutine init_physics(options,domain)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(inout) :: domain

		! initialize microphysics code (e.g. compute look up tables in Thompson et al)
		call mp_init(options) !this could easily be moved to init_model...
		
		call init_convection(domain,options)
		
		call pbl_init(domain,options)
		
		call radiation_init(domain,options)
		
		call lsm_init(domain,options)
		
	end subroutine init_physics
	
! 	read physics options to use from a namelist file
	subroutine physics_namelist(filename,options)
		implicit none
		character(len=*),intent(in) :: filename
		type(options_type), intent(inout) :: options
		integer :: name_unit
		
! 		variables to be used in the namelist
		integer::pbl,lsm,mp,rad,conv,adv,wind
! 		define the namelist
		namelist /physics/ pbl,lsm,mp,rad,conv,adv,wind
		
! 		default values for physics options (advection+linear winds+simple_microphysics)
		pbl=2   ! 0 = no PBL, 1 = STUPID PBL (in LSM module), 2 = local PBL diffusion
		lsm=3   ! 0 = no LSM, 1 = Fluxes from GCM, 2 = simple LSM, 3 = Noah LSM
		mp=2	! 0 = no MP,  1 = Thompson et al (2008), 2 = "Linear" microphysics
		rad=2   ! 0 = no RAD, 1 = radiative cooling 1K/day (in LSM), 2 = cloud fraction based radiation
		conv=0  ! 0 = no CONV,1 = Tiedke scheme
		adv=1   ! 0 = no ADV, 1 = upwind advection scheme
		wind=1  ! 0 = no LT,  1 = linear theory wind perturbations
		
! 		read the namelist
		open(io_newunit(name_unit), file=filename)
		read(name_unit,nml=physics)
		close(name_unit)

! 		store options
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
		                                vegtype_var,vegfrac_var
										
		namelist /var_list/ pvar,pbvar,tvar,qvvar,qcvar,qivar,hgtvar,shvar,lhvar,pblhvar,&
							landvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar, &
							hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi, &
							soiltype_var, soil_t_var,soil_vwc_var,soil_deept_var, &
		                    vegtype_var,vegfrac_var
		
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
		
		open(io_newunit(name_unit), file=filename)
		read(name_unit,nml=var_list)
		close(name_unit)

! 		2D geometry variable names (for coarse model)
		options%hgtvar=hgtvar
		options%latvar=latvar
		options%lonvar=lonvar
! 		U varname and associated lat/lon var names
		options%uvar=uvar
		options%ulat=ulat
		options%ulon=ulon
! 		V varname and associated lat/lon var names
		options%vvar=vvar
		options%vlat=vlat
		options%vlon=vlon
! 		Primary model variable names
		options%pbvar=pbvar
		options%pvar=pvar
		options%tvar=tvar
		options%qvvar=qvvar
		options%qcvar=qcvar
		options%qivar=qivar
! 		vertical coordinate
		options%zvar=zvar
! 		2D model variables (e.g. Land surface and PBL height)		
		options%shvar=shvar
		options%lhvar=lhvar
		options%pblhvar=pblhvar
		
! 		separate variable names for the high resolution domain
		options%hgt_hi=hgt_hi
		options%landvar=landvar
		options%lat_hi=lat_hi
		options%lon_hi=lon_hi
		options%ulat_hi=ulat_hi
		options%ulon_hi=ulon_hi
		options%vlat_hi=vlat_hi
		options%vlon_hi=vlon_hi
		
! 		soil and vegetation parameters
		options%soiltype_var=soiltype_var
		options%soil_t_var=soil_t_var
		options%soil_vwc_var=soil_vwc_var
		options%soil_deept_var=soil_deept_var
		options%vegtype_var=vegtype_var
		options%vegfrac_var=vegfrac_var
	end subroutine var_namelist
	
	subroutine parameters_namelist(filename,options)
		implicit none
		character(len=*),intent(in) :: filename
		type(options_type), intent(inout) :: options
		integer :: name_unit
		
		real    :: dx, dxlow, outputinterval, inputinterval, t_offset, smooth_wind_distance
		real    :: rotation_scale_height, N_squared,linear_contribution
		integer :: ntimesteps, nfiles, xmin, xmax, ymin, ymax, vert_smooth
		integer :: nz, n_ext_winds,buffer, warning_level
		logical :: ideal, readz, readdz, debug, external_winds, remove_lowres_linear, variable_N, &
		           mean_winds, mean_fields, restart, add_low_topo, advect_density, high_res_soil_state
		character(len=MAXFILELENGTH) :: date, calendar
		integer :: year, month, day, hour, minute, second
		
		
		namelist /parameters/ ntimesteps,outputinterval,inputinterval,dx,dxlow,ideal,readz,readdz,nz,t_offset,debug,nfiles, &
							  external_winds,buffer,n_ext_winds,add_low_topo,advect_density,smooth_wind_distance, &
							  remove_lowres_linear,mean_winds,mean_fields,restart,xmin,xmax,ymin,ymax,vert_smooth, &
							  date, calendar, high_res_soil_state,rotation_scale_height,warning_level, variable_N, &
							  N_squared,linear_contribution
		
! 		default parameters
		mean_fields=.False.
		mean_winds=.False.
		external_winds=.False.
		n_ext_winds=1
		t_offset=(-9999)
		buffer=0
		remove_lowres_linear=.False.
		add_low_topo=.False.
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
		smooth_wind_distance=-9999
		vert_smooth=2
		calendar="gregorian"
		high_res_soil_state=.True.
		rotation_scale_height=2000.0
		N_squared=6.37e-5
		variable_N=.False.
		linear_contribution=1.0
		
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
			
		
		if (t_offset.eq.(-9999)) then
			write(*,*), "WARNING, WARNING, WARNING"
			write(*,*), "WARNING, WARNING, WARNING"
			write(*,*), ""
			write(*,*), "	Using default t_offset=300"
			write(*,*), ""
			write(*,*), "WARNING, WARNING, WARNING"
			write(*,*), "WARNING, WARNING, WARNING"
			t_offset=300
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
		options%add_low_topo=add_low_topo
		options%mean_winds=mean_winds
		options%mean_fields=mean_fields
		options%advect_density=advect_density
		options%debug=debug
		options%warning_level=warning_level
		options%rotation_scale_height=rotation_scale_height
		
		options%external_winds=external_winds
		options%ext_winds_nfiles=n_ext_winds
		options%restart=restart
		
		options%nz=nz
		options%xmin=xmin
		options%xmax=xmax
		options%ymin=ymin
		options%ymax=ymax
		
		options%N_squared=N_squared
		options%variable_N=variable_N
		options%linear_contribution=linear_contribution

		options%high_res_soil_state=high_res_soil_state
		
		
	end subroutine parameters_namelist
	
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
		type(options_type), intent(in)::options
		
		! convection can modify wind field, and ideal doesn't rebalance winds every timestep
		if ((options%physics%convection.ne.0).and.(options%ideal)) then
			if (options%warning_level>1) then
				write(*,*) "WARNING WARNING WARNING"
				write(*,*) "WARNING, running convection in ideal mode may be bad..."
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
	
	subroutine init_restart(options_filename, options)
		character(len=*), intent(in) :: options_filename
		type(options_type),intent(inout) :: options
		
		character(len=MAXFILELENGTH) :: init_conditions_file, output_file,restart_file
		integer :: restart_step, restart_date(6), name_unit
		real*8 :: restart_mjd
		
		namelist /restart_info/ restart_step, restart_file, restart_date
		
		restart_date=[-999,-999,-999,-999,-999,-999]
		restart_step=-999
		
		open(io_newunit(name_unit), file=options_filename)
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
		
	end subroutine init_restart
	
	subroutine init_options(options_filename,options)
! 		reads a series of options from a namelist file and stores them in the 
! 		options data structure
		implicit none
		character(len=*), intent(in) :: options_filename
		type(options_type), intent(inout) :: options
		
		character(len=MAXFILELENGTH) :: init_conditions_file, output_file
		character(len=MAXFILELENGTH),allocatable:: boundary_files(:),ext_wind_files(:)
		real, allocatable, dimension(:) :: dz_levels
   		real,dimension(45)::fulldz
		integer :: name_unit
		
! 		set up namelist structures
		namelist /z_info/ dz_levels
		namelist /files_list/ init_conditions_file,output_file,boundary_files
		namelist /ext_winds_info/ ext_wind_files
		
		call version_check(options_filename,options)
		call physics_namelist(options_filename,options)
		call var_namelist(options_filename,options)
		call parameters_namelist(options_filename,options)
		
		! read the z_info namelist if requested
		if (options%readdz) then
			allocate(dz_levels(options%nz),options%dz_levels(options%nz))
			
			open(io_newunit(name_unit), file=options_filename)
			read(name_unit,nml=z_info)
			close(name_unit)
			
			options%dz_levels=dz_levels
			deallocate(dz_levels)
		else
		! if we are not reading dz from the namelist, use default values from a WRF run
		! default mean layer thicknesses from a 36km WRF run over the "CO-headwaters" domain
			fulldz=[36.,   51.,   58.,   73.,   74.,  111.,  113.,  152.,  155.,  157.,  160.,  245., &
				   251.,  258.,  265.,  365.,  379.,  395.,  413.,  432.,  453.,  476.,  503.,  533., &
				   422.,  443.,  467.,  326.,  339.,  353.,  369.,  386.,  405.,  426.,  450.,  477., &
				   455.,  429.,  396.,  357.,  311.,  325.,  340.,  356.,  356.]
 			allocate(options%dz_levels(options%nz))
			options%dz_levels=fulldz(1:options%nz)
		endif
		
		if (options%restart) then
			call init_restart(options_filename,options)
		endif

		
		allocate(boundary_files(options%nfiles))
		open(io_newunit(name_unit), file=options_filename)
		read(name_unit,nml=files_list)
		close(name_unit)
		
		options%init_conditions_file=init_conditions_file
		
		allocate(options%boundary_files(options%nfiles))
		options%boundary_files=boundary_files
		deallocate(boundary_files)
		
		options%output_file=output_file
		
		if (options%external_winds) then
			allocate(ext_wind_files(options%ext_winds_nfiles))
			
			open(io_newunit(name_unit), file=options_filename)
			read(name_unit,nml=ext_winds_info)
			close(name_unit)
			
			allocate(options%ext_wind_files(options%ext_winds_nfiles))
			options%ext_wind_files=ext_wind_files
			deallocate(ext_wind_files)
		endif
		
		! check for any inconsistencies in the options requested
		call options_check(options)
	end subroutine init_options
	
! 	Allow running over a sub-domain, by removing the outer N grid cells from all sides of the domain (lat,lon,terrain)
	subroutine remove_edges(domain,edgesize)
		implicit none
		type(domain_type), intent(inout) :: domain
		integer, intent(in)::edgesize
		integer::nx1,ny1,nx2,ny2,nz
		real,allocatable,dimension(:,:)::temp_data
		
		nx1=size(domain%lat,1)
		ny1=size(domain%lat,2)
		nx2=nx1-(edgesize*2)
		ny2=ny1-(edgesize*2)
		
		allocate(temp_data(nx1,ny1))
		
		temp_data=domain%lat
		deallocate(domain%lat)
		allocate(domain%lat(nx2,ny2))
		domain%lat=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
		
		temp_data=domain%lon
		deallocate(domain%lon)
		allocate(domain%lon(nx2,ny2))

		domain%lon=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		temp_data=domain%terrain
		deallocate(domain%terrain)
		allocate(domain%terrain(nx2,ny2))
		domain%terrain=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		temp_data=domain%landmask
		deallocate(domain%landmask)
		allocate(domain%landmask(nx2,ny2))
		domain%landmask=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		
		deallocate(temp_data)
		
		nx1=size(domain%u_geo%lat,1)
		ny1=size(domain%u_geo%lat,2)
		nx2=nx1-(edgesize*2)
		ny2=ny1-(edgesize*2)
		
		allocate(temp_data(nx1,ny1))
		temp_data=domain%u_geo%lat
		deallocate(domain%u_geo%lat)
		allocate(domain%u_geo%lat(nx2,ny2))
		domain%u_geo%lat=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
		
		temp_data=domain%u_geo%lon
		deallocate(domain%u_geo%lon)
		allocate(domain%u_geo%lon(nx2,ny2))
		domain%u_geo%lon=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		deallocate(temp_data)
		
		nx1=size(domain%v_geo%lat,1)
		ny1=size(domain%v_geo%lat,2)
		nx2=nx1-(edgesize*2)
		ny2=ny1-(edgesize*2)
		
		allocate(temp_data(nx1,ny1))
		temp_data=domain%v_geo%lat
		deallocate(domain%v_geo%lat)
		allocate(domain%v_geo%lat(nx2,ny2))
		domain%v_geo%lat=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
		
		temp_data=domain%v_geo%lon
		deallocate(domain%v_geo%lon)
		allocate(domain%v_geo%lon(nx2,ny2))
		domain%v_geo%lon=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
		
		deallocate(temp_data)
	end subroutine remove_edges
	
! 	allocate all arrays in domain
	subroutine domain_allocation(domain,nx,nz,ny,nsoil)
		implicit none
		type(domain_type), intent(inout) :: domain
		integer,intent(in)::nx,nz,ny
		integer,intent(in),optional :: nsoil
		integer :: ns
		
		ns=4
		if (present(nsoil)) then
			ns=nsoil
		endif
		
		! atmosphere allocation
		allocate(domain%p(nx,nz,ny))		! air pressure [Pa]
		domain%p=100000
		allocate(domain%u(nx+1,nz,ny))		! eastward wind [m/s]
		domain%u=0
		allocate(domain%v(nx,nz,ny+1))		! northward wind [m/s]
		domain%v=0
		allocate(domain%w(nx,nz,ny))		! vertical wind [grid/s]
		domain%w=0
		allocate(domain%ur(nx+1,nz,ny))		! eastward wind * density [m/s kg/m^3]
		domain%ur=0
		allocate(domain%vr(nx,nz,ny+1))		! northward wind * density[m/s kg/m^3]
		domain%vr=0
		allocate(domain%wr(nx,nz,ny))		! vertical wind * density [grid/s kg/m^3]
		domain%wr=0
		allocate(domain%th(nx,nz,ny))		! potential temperature [K]
		domain%th=280
		allocate(domain%qv(nx,nz,ny))		! water vapor [kg/kg]
		domain%qv=0.0002
		allocate(domain%cloud(nx,nz,ny))	! liquid cloud water content mixing ratio [kg/kg]
		domain%cloud=0
		allocate(domain%ice(nx,nz,ny))		! frozen cloud water content mixing ratio [kg/kg]
		domain%ice=0
		allocate(domain%nice(nx,nz,ny))		! cloud ice number concentration [cm-3]
		domain%nice=0
		allocate(domain%qrain(nx,nz,ny))	! rain mixing ratio [kg/kg]
		domain%qrain=0
		allocate(domain%nrain(nx,nz,ny))	! rain drop number concentration [cm-3]
		domain%nrain=0
		allocate(domain%qsnow(nx,nz,ny))	! snow  mixing ratio [kg/kg]
		domain%qsnow=0
		allocate(domain%qgrau(nx,nz,ny))	! graupel mixing ratio [kg/kg]
		domain%qgrau=0
		allocate(domain%pii(nx,nz,ny))		! exner function
		domain%pii=1
		allocate(domain%qv_adv_tendency(nx,nz,ny)) ! advective qv tendency [qv/s]
		domain%qv_adv_tendency=0
		allocate(domain%qv_pbl_tendency(nx,nz,ny)) ! qv tendency from PBL scheme [qv/s]
		domain%qv_pbl_tendency=0
		allocate(domain%rho(nx,nz,ny))		! air density [kg/m^3]
		domain%rho=1
		allocate(domain%cloudfrac(nx,ny))	! cloud fraction
		domain%cloudfrac=0
		
		! land-atm flux allocation
		allocate(domain%rain(nx,ny))		! accumulated total rainfall [kg/m^2]
		domain%rain=0
		allocate(domain%crain(nx,ny))		! accumulated convective rainfall
		domain%crain=0
		allocate(domain%snow(nx,ny))		! accumulated snow fall
		domain%snow=0
		allocate(domain%graupel(nx,ny))		! accumulated graupel fall
		domain%graupel=0
		allocate(domain%current_rain(nx,ny))! rain fall in current time step
		domain%current_rain=0
		allocate(domain%current_snow(nx,ny))! snow fall in current time step
		domain%current_snow=0
		allocate(domain%swdown(nx,ny))		! shortwave down at surface
		domain%swdown=0
		allocate(domain%lwdown(nx,ny))		! longwave down at surface
		domain%lwdown=0
		allocate(domain%lwup(nx,ny))		! longwave up from surface
		domain%lwup=0

		allocate(domain%sensible_heat(nx,ny)) ! sensible heat flux from surface
		domain%sensible_heat=0
		allocate(domain%latent_heat(nx,ny))	! latent heat flux from surface
		domain%latent_heat=0
		allocate(domain%ground_heat(nx,ny))	! ground heat flux into ground
		domain%ground_heat=0
		allocate(domain%pbl_height(nx,ny))	! planetary boundary layer height (not always used)
		domain%pbl_height=0
		
		! land surface allocation
		allocate(domain%soil_t(nx,ns,ny))	! 3D soil temperature
		domain%soil_t=280
		allocate(domain%soil_vwc(nx,ns,ny))	! 3D soil volumetric water content
		domain%soil_vwc=0.25
		
		allocate(domain%soil_tdeep(nx,ny))		! deep soil temperature
		domain%soil_tdeep=280
		allocate(domain%soil_totalmoisture(nx,ny))	! soil column total moisture content
		domain%soil_totalmoisture=500 ! =2000mm * 0.25 (vwc)
		allocate(domain%skin_t(nx,ny))			! skin temperature
		domain%skin_t=280
		allocate(domain%snow_swe(nx,ny))		! snow water equivalent
		domain%snow_swe=0
		allocate(domain%vegfrac(nx,ny))			! vegetation cover fraction (%)
		domain%vegfrac=50 !% veg cover
		allocate(domain%canopy_water(nx,ny))	! canopy water content
		domain%canopy_water=0
		allocate(domain%soil_type(nx,ny))		! USGS soil type
		domain%soil_type=6 ! Loam
		allocate(domain%veg_type(nx,ny))		! Vegetation type
		domain%veg_type=7  ! grassland

	end subroutine domain_allocation
	
	subroutine copy_z(input,output,interpolate_dim)
		implicit none
		class(interpolable_type), intent(in) :: input
		class(interpolable_type), intent(inout) :: output
		integer,intent(in)::interpolate_dim
		
		integer::nxi,nyi,nzi,nxo,nyo,nzo
		
		! dimensions of the input data
		nxi=size(input%z,1)
		nzi=size(input%z,2)
		nyi=size(input%z,3)
		! dimensions of the output data
		nxo=size(output%lat,1)
		nzo=nzi
		nyo=size(output%lat,2)
		
		allocate(output%z(nxo,nzo,nyo))
		if (interpolate_dim==1) then
			output%z(2:nxo-1,:,:)=(input%z(1:nxi-1,:,:) + input%z(2:nxi,:,:))/2
			output%z(1,:,:)=input%z(1,:,:)
			output%z(nxo,:,:)=input%z(nxi,:,:)
		else if (interpolate_dim==3) then
			output%z(:,:,2:nyo-1)=(input%z(:,:,1:nyi-1) + input%z(:,:,2:nyi))/2
			output%z(:,:,1)=input%z(:,:,1)
			output%z(:,:,nyo)=input%z(:,:,nyi)
		else
			write(*,*) "Can not interpolate z data over z dimension"
		endif
		
	end subroutine copy_z

	subroutine init_domain_land(domain,options)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(inout):: domain
		integer,dimension(:,:),allocatable :: temp_idata
		real,dimension(:,:,:),allocatable :: temp_rdata
		real,dimension(:,:),allocatable :: temp_rdata_2d
		integer:: buffer,nx,ny,nz
		buffer=options%buffer
		nx=size(domain%veg_type,1)+buffer*2
		ny=size(domain%veg_type,2)+buffer*2
		nz=size(domain%soil_t,2)
		
		
		! Veg cover fraction = 2D real
		if (options%vegfrac_var.ne."") then
			call io_read2d(options%init_conditions_file,options%vegfrac_var,temp_rdata_2d,1)
			domain%vegfrac=temp_rdata_2d(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
			deallocate(temp_rdata_2d)
		endif

		! Veg TYPE = 2D integer
		if (options%vegtype_var.ne."") then
			call io_read2di(options%init_conditions_file,options%vegtype_var,temp_idata,1)
			domain%veg_type=temp_idata(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
			deallocate(temp_idata)
		endif
		
		! Soil TYPE = 2D integer
		if (options%soiltype_var.ne."") then
			call io_read2di(options%init_conditions_file,options%soiltype_var,temp_idata,1)
			domain%soil_type=temp_idata(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
			deallocate(temp_idata)
		endif
		
		! Soil Volumetric Water Content = 3D real
		if (options%soil_vwc_var.ne."") then
			call io_read3d(options%init_conditions_file,options%soil_vwc_var,temp_rdata,1)  
			domain%soil_vwc=reshape(temp_rdata(1+buffer:nx-buffer,1+buffer:ny-buffer,1:nz) &	! subset the data by buffer grid cells
									,[nx-buffer*2,nz,ny-buffer*2],order=[1,3,2])				! and reshape to move the z axis to the middle
			deallocate(temp_rdata)
		endif
		
		! Soil Temperature = 3D real
		if (options%soil_t_var.ne."") then
			call io_read3d(options%init_conditions_file,options%soil_t_var,temp_rdata,1)
			domain%soil_t=reshape(temp_rdata(1+buffer:nx-buffer,1+buffer:ny-buffer,1:nz) &  ! subset the data by buffer grid cells
									,[nx-buffer*2,nz,ny-buffer*2],order=[1,3,2])			! and reshape to move the z axis to the middle
			deallocate(temp_rdata)
		endif
		
		! Deep Soil Temperature = 2D real
		if (options%soil_deept_var.ne."") then
			call io_read2d(options%init_conditions_file,options%soil_deept_var,temp_rdata_2d,1)
			domain%soil_tdeep=temp_rdata_2d(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
			deallocate(temp_rdata_2d)
		endif
		
		
	end subroutine init_domain_land
		

! 	initialize the domain e.g. lat,lon,terrain, 3D z coordinate
	subroutine init_domain(options, domain)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(inout):: domain
		real,dimension(:,:,:),allocatable::temporary_z
		integer:: ny,nz,nx,i,buf
		
! 		these are the only required variables on a high-res grid, lat, lon, and terrain elevation
		call io_read2d(options%init_conditions_file,options%hgt_hi,domain%terrain,1)
		call io_read2d(options%init_conditions_file,options%lat_hi,domain%lat,1)
		call io_read2d(options%init_conditions_file,options%lon_hi,domain%lon,1)
		call io_read2d(options%init_conditions_file,options%ulat_hi,domain%u_geo%lat,1)
		call io_read2d(options%init_conditions_file,options%ulon_hi,domain%u_geo%lon,1)
		call io_read2d(options%init_conditions_file,options%vlat_hi,domain%v_geo%lat,1)
		call io_read2d(options%init_conditions_file,options%vlon_hi,domain%v_geo%lon,1)
		if (options%landvar/="") then
			call io_read2d(options%init_conditions_file,options%landvar,domain%landmask,1)
		else
			nx=size(domain%lat,1)
			ny=size(domain%lat,2)
			allocate(domain%landmask(nx,ny))
			domain%landmask=1 !if we weren't supplied a landmask field, assume all is land (what we care about anyway)
		endif
		
		if(options%buffer>0) then
			call remove_edges(domain,options%buffer)
		endif
! 		use the lat variable to define the x and y dimensions for all other variables
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
! 		assumes nz is defined in the options
		nz=options%nz
		
! 		if a 3d grid was also specified, then read those data in
		if ((options%readz).and.(options%ideal).and.(options%zvar.ne."")) then
			
			call io_read3d(options%init_conditions_file,options%zvar, domain%z)
! 			dz also has to be calculated from the 3d z file
			buf=options%buffer
			allocate(domain%dz(nx,nz,ny))
			allocate(temporary_z(nx,nz,ny))
			do i=1,nz-1
				domain%dz(:,i,:)=domain%z(buf+1:nx+buf,buf+1:ny+buf,i+1)-domain%z(buf+1:nx+buf,buf+1:ny+buf,i)
				temporary_z(:,i,:)=domain%z(buf+1:nx+buf,buf+1:ny+buf,i)
			enddo
			temporary_z(:,nz,:)=domain%z(buf+1:nx+buf,buf+1:ny+buf,nz)
			domain%dz(:,nz,:)=domain%dz(:,nz-1,:)
			deallocate(domain%z)
			allocate(domain%z(nx,nz,ny))
			domain%z=temporary_z
			deallocate(temporary_z)
		else
! 			otherwise, set up the z grid to be evenly spaced in z using the terrain +dz/2 for the base
! 			and z[1]+i*dz for the res
			allocate(domain%z(nx,nz,ny))
			allocate(domain%dz(nx,nz,ny))
			domain%dz(:,1,:)=options%dz_levels(1)
			domain%z(:,1,:)=domain%terrain+options%dz_levels(1)/2
			do i=2,nz
				domain%z(:,i,:)=domain%z(:,i-1,:)+(options%dz_levels(i)+options%dz_levels(i-1))/2
				domain%dz(:,i,:)=options%dz_levels(i)
			enddo
			
		endif
		call copy_z(domain,domain%u_geo,interpolate_dim=1)
		call copy_z(domain,domain%v_geo,interpolate_dim=3)
		
! 		all other variables should be allocated and initialized to 0
		call domain_allocation(domain,nx,nz,ny)
		
! 		initializing land
		call init_domain_land(domain,options)
		
! 		store dx in domain as well as options, read as an option, but it is more appropriate in domain
		domain%dx=options%dx
		
		call init_winds(domain,options)
	end subroutine init_domain
	
! 	allocate arrays in boundary condition data structure
	subroutine boundary_allocate(boundary,nx,nz,ny)
		implicit none
		type(bc_type), intent(inout) :: boundary
		integer,intent(in)::nx,nz,ny
		
		allocate(boundary%du_dt(nx+1,nz,ny))
		boundary%du_dt=0
		allocate(boundary%dv_dt(nx,nz,ny+1))
		boundary%dv_dt=0
		allocate(boundary%dw_dt(nx,nz,ny))
		boundary%dw_dt=0
		allocate(boundary%dp_dt(nx,nz,ny))
		boundary%dp_dt=0
		allocate(boundary%drho_dt(nx,nz,ny))
		boundary%drho_dt=0
		allocate(boundary%dth_dt(nz,max(nx,ny),4))
		boundary%dth_dt=0
		allocate(boundary%dqv_dt(nz,max(nx,ny),4))
		boundary%dqv_dt=0
		allocate(boundary%dqc_dt(nz,max(nx,ny),4))
		boundary%dqc_dt=0
		allocate(boundary%dlh_dt(nx,ny))
		boundary%dlh_dt=0
		allocate(boundary%dsh_dt(nx,ny))
		boundary%dsh_dt=0
		allocate(boundary%dpblh_dt(nx,ny))
		boundary%dpblh_dt=0
	end subroutine boundary_allocate
	
! 	initialize the boundary condition data structure e.g. lat,lon,terrain,3D Z coord
	subroutine init_bc_data(options,boundary,domain)
		implicit none
		type(options_type), intent(in) :: options
		type(bc_type), intent(inout):: boundary
		type(domain_type), intent(in):: domain
		real,dimension(45)::fulldz
		real,dimension(:,:,:),allocatable::zbase
		integer::nx,ny,nz,i
		
! 		these variables are required for any boundary/forcing file type
		call io_read2d(options%boundary_files(1),options%latvar,boundary%lat)
		call io_read2d(options%boundary_files(1),options%lonvar,boundary%lon)
		call io_read2d(options%boundary_files(1),options%ulat,boundary%u_geo%lat)
		call io_read2d(options%boundary_files(1),options%ulon,boundary%u_geo%lon)
		call io_read2d(options%boundary_files(1),options%vlat,boundary%v_geo%lat)
		call io_read2d(options%boundary_files(1),options%vlon,boundary%v_geo%lon)
		call io_read2d(options%boundary_files(1),options%hgtvar,boundary%terrain)
		
! 		read in the vertical coordinate
		call io_read3d(options%boundary_files(1),options%zvar,zbase)
		nx=size(zbase,1)
		ny=size(zbase,2)
		nz=size(zbase,3)
		allocate(boundary%lowres_z(nx,nz,ny))
! 		write(*,*) "WARNING: HACK in init.f90 boundary init : z=z*9.8"
! 		boundary%lowres_z=reshape(zbase*9.8,[nx,nz,ny],order=[1,3,2])
		boundary%lowres_z=reshape(zbase,[nx,nz,ny],order=[1,3,2])
		deallocate(zbase)
		if (options%zvar=="PH") then
			call io_read3d(options%boundary_files(1),"PHB",zbase)
			
			boundary%lowres_z=(boundary%lowres_z+reshape(zbase,[nx,nz,ny],order=[1,3,2])) / gravity
			deallocate(zbase)
		endif
		
! 		nx=size(boundary%lat,1)
! 		nz=options%nz
! 		ny=size(boundary%lat,2)
! 		allocate(boundary%z(nx,nz,ny))
! 		allocate(boundary%dz(nx,nz,ny))
! 		! assumes that vertical level thicknesses are identical in the forcing and model domain!
! 		fulldz(:nz)=domain%dz(1,:nz,1)
! 		! set up the bc_data dz and z variables (vertical levels)
! 		boundary%dz(:,1,:)=fulldz(1)
! 		! z is calculated from dz, dz is thickness of layers.
! 		boundary%z(:,1,:)=boundary%terrain+fulldz(1)/2
! 		do i=2,nz
! 			boundary%z(:,i,:)=boundary%z(:,i-1,:)+(fulldz(i)+fulldz(i-1))/2
! 			boundary%dz(:,i,:)=fulldz(i)
! 		enddo
		
		! all other structures must be allocated and initialized, but will be set on a high-res grid
		! u/v are seperate so we can read them on the low res grid and adjust/rm-linearwinds before interpolating
		! this also makes it easier to change how these variables are read from various forcing model file structures
		nx=size(boundary%u_geo%lon,1)
		ny=size(boundary%u_geo%lon,2)
		allocate(boundary%u(nx,nz,ny))
		boundary%u=0
		nx=size(boundary%v_geo%lon,1)
		ny=size(boundary%v_geo%lon,2)
		allocate(boundary%v(nx,nz,ny))
		boundary%v=0
		
		nz=options%nz
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
		
		call boundary_allocate(boundary,nx,nz,ny)
	end subroutine init_bc_data
	
	
! 	sets up the data used to rotate the wind field when using an external
! 	high-res wind field from e.g. a high res WRF run
    subroutine setup_extwinds(domain)
        implicit none
        type(wind_type),intent(inout)::domain
        integer::nx,ny
        
        nx=size(domain%terrain,1)
        ny=size(domain%terrain,2)
		
		! dzdx/y used in rotating windfield back to terrain following grid in a simple fashion
		allocate(domain%dzdx(nx-1,ny))
		allocate(domain%dzdy(nx,ny-1))
		domain%dzdx=sqrt((domain%terrain(2:nx,:)-domain%terrain(1:nx-1,:))**2+domain%dx**2)/domain%dx
		domain%dzdy=sqrt((domain%terrain(:,2:ny)-domain%terrain(:,1:ny-1))**2+domain%dx**2)/domain%dx
    end subroutine
	
	
! 	initialize the external wind system (GEOLUT)
	subroutine init_ext_winds(options,bc)
		implicit none
		type(options_type), intent(in) :: options
		type(bc_type),intent(inout) :: bc
			
		real, allocatable, dimension(:,:,:) :: u,v
		real, allocatable, dimension(:,:) :: lat,lon
		real, allocatable, dimension(:,:) :: temporary_terrain
		
		call io_read2d(options%ext_wind_files(1),options%hgt_hi,bc%ext_winds%terrain,1)
		call io_read2d(options%ext_wind_files(1),options%latvar,bc%ext_winds%lat)
		call io_read2d(options%ext_wind_files(1),options%lonvar,bc%ext_winds%lon)
		call io_read2d(options%ext_wind_files(1),options%ulat_hi,bc%ext_winds%u_geo%lat)
		call io_read2d(options%ext_wind_files(1),options%ulon_hi,bc%ext_winds%u_geo%lon)
		call io_read2d(options%ext_wind_files(1),options%vlat_hi,bc%ext_winds%v_geo%lat)
		call io_read2d(options%ext_wind_files(1),options%vlon_hi,bc%ext_winds%v_geo%lon)
		write(*,*) "Setting up ext wind geoLUTs"
		call geo_LUT(bc%next_domain%u_geo, bc%ext_winds%u_geo)
		call geo_LUT(bc%next_domain%v_geo, bc%ext_winds%v_geo)
		
		! if the external wind file has a different shape in either dimension, compute the GEOLUT and interpolate terrain
		! this was particularly problematic for ideal WRF simluation comparisons
		if ((size(bc%ext_winds%terrain,1)/=size(bc%next_domain%terrain,1)).or. &
		    (size(bc%ext_winds%terrain,2)/=size(bc%next_domain%terrain,2))) then
			call geo_LUT(bc%next_domain, bc%ext_winds)
			allocate(temporary_terrain(size(bc%ext_winds%terrain,1),size(bc%ext_winds%terrain,2)))
			temporary_terrain=bc%ext_winds%terrain
			deallocate(bc%ext_winds%terrain)
			allocate(bc%ext_winds%terrain(size(bc%next_domain%terrain,1),size(bc%next_domain%terrain,2)))
			call geo_interp2d(bc%ext_winds%terrain,temporary_terrain,bc%ext_winds%geolut)
			deallocate(temporary_terrain)
		endif
		
! 		force all weight to be on the first x,y pair...
! 		this assumes the "external winds" file is on the exact same grid as the high res model grid
		bc%ext_winds%u_geo%geolut%w(2:,:,:)=0
		bc%ext_winds%u_geo%geolut%w(1,:,:)=1
		bc%ext_winds%v_geo%geolut%w(2:,:,:)=0
		bc%ext_winds%v_geo%geolut%w(1,:,:)=1
		bc%ext_winds%dx=bc%next_domain%dx
		call setup_extwinds(bc%ext_winds)
	end subroutine init_ext_winds
	
	subroutine swap_z(bc)
		type(bc_type), intent(inout) :: bc
		real,allocatable,dimension(:,:,:) :: tempz
		integer::nx,nz,ny
		nx=size(bc%lowres_z,1)
		nz=size(bc%lowres_z,2)
		ny=size(bc%lowres_z,3)
		allocate(tempz(nx,nz,ny))
		tempz=bc%lowres_z
		
		nx=size(bc%z,1)
		nz=size(bc%z,2)
		ny=size(bc%z,3)
		deallocate(bc%lowres_z)
		allocate(bc%lowres_z(nx,nz,ny))
		bc%lowres_z=bc%z

		nx=size(tempz,1)
		nz=size(tempz,2)
		ny=size(tempz,3)
		deallocate(bc%z)
		allocate(bc%z(nx,nz,ny))
		bc%z=tempz
		deallocate(tempz)
		
	end subroutine swap_z

	subroutine move_lut(inputgeo,outputgeo)
		! move the contents of one geographic look up table into another
		! moving implies deletion of the initial data, so inputgeo is destroyed
		type(geo_look_up_table), intent(inout) :: inputgeo,outputgeo
		integer::nx,ny,nz
		nx=size(inputgeo%x,1)
		ny=size(inputgeo%x,2)
		nz=size(inputgeo%x,3)
		
		allocate(outputgeo%x(nx,ny,nz))
		allocate(outputgeo%y(nx,ny,nz))
		allocate(outputgeo%w(nx,ny,nz))
		
		outputgeo%x=inputgeo%x
		outputgeo%y=inputgeo%y
		outputgeo%w=inputgeo%w
		
		call destroy_lut(inputgeo)
	end subroutine move_lut
	
	subroutine destroy_lut(geolut)
		! deallocate all memory associated with a geographic look up table
		type(geo_look_up_table), intent(inout) :: geolut
		deallocate(geolut%x)
		deallocate(geolut%y)
		deallocate(geolut%w)
	end subroutine destroy_lut
	
! 	initialize the boundary condiditions (init data structures and GEOLUT)
!	initializes external winds if necessary, adds low-res terrain to the high-res domain if desired
	subroutine init_bc(options,domain,boundary)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(inout):: domain
		type(bc_type), intent(inout):: boundary
		type(geo_look_up_table) :: u_temp_geo,v_temp_geo
		integer::i,nx,ny,nz
			
		boundary%dx=options%dxlow
! 		set up base data
		call init_bc_data(options,boundary,domain)
		call init_domain(options,boundary%next_domain) !set up a domain to hold the forcing for the next time step
		
! 		create the geographic look up table used to calculate boundary forcing data
		write(*,*) "Setting up domain geographic Look Up Tables"
		! set up a look up table from low-res grid center to high-res u-offset coordinates
		call geo_LUT(domain%u_geo,boundary)
		call move_lut(boundary%geolut,u_temp_geo)
		! set up a look up table from low-res grid center to high-res v-offset coordinates
		call geo_LUT(domain%v_geo,boundary)
		call move_lut(boundary%geolut,v_temp_geo)
		call geo_LUT(domain,boundary)
		call geo_LUT(domain%u_geo,boundary%u_geo)
		call geo_LUT(domain%v_geo,boundary%v_geo)
		
		if (options%external_winds) then
			call init_ext_winds(options,boundary)
		endif
		
		! interpolate the low-res terrain to the high-res grid for pressure adjustments. 
		! the correct way would probably be to adjust all low-res pressures to Sea level before interpolating
		! then pressure adjustments all occur from SLP. 
		! This should be done on a separate lowres terrain grid so the embedded high res terrain grid 
		! can also be used in pressure adjustments on each time step...
		nx=size(domain%terrain,1)
		ny=size(domain%terrain,2)
		allocate(boundary%lowres_terrain(nx,ny))
		call geo_interp2d(boundary%lowres_terrain,boundary%terrain,boundary%geolut)
		
		nz=size(boundary%lowres_z,2)
		allocate(boundary%z(nx,nz,ny))
		call geo_interp(boundary%z,boundary%lowres_z,boundary%geolut,.false.)
		
		nx=size(domain%u_geo%lat,1)
		ny=size(domain%u_geo%lat,2)
		allocate(boundary%u_geo%z(nx,nz,ny))
		call geo_interp(boundary%u_geo%z,boundary%lowres_z,u_temp_geo,.false.)
		
		nx=size(domain%v_geo%lat,1)
		ny=size(domain%v_geo%lat,2)
		allocate(boundary%v_geo%z(nx,nz,ny))
		call geo_interp(boundary%v_geo%z,boundary%lowres_z,v_temp_geo,.false.)
		
		call destroy_lut(v_temp_geo)
		call destroy_lut(u_temp_geo)
		
		if (options%add_low_topo) then
			domain%terrain=domain%terrain+(boundary%lowres_terrain-sum(boundary%lowres_terrain) &
											 /size(boundary%lowres_terrain))/2.0
			do i=1,size(domain%z,2)
				domain%z(:,i,:)=domain%z(:,i,:)+(boundary%lowres_terrain-sum(boundary%lowres_terrain) &
												 /size(boundary%lowres_terrain))/2.0
			enddo
		endif
		write(*,*) "Setting up vertical interpolation Look Up Tables"
		call vLUT(domain,boundary)
		call vLUT(domain%u_geo,boundary%u_geo)
		call vLUT(domain%v_geo,boundary%v_geo)
! 		call io_write3di("vlutz1.nc","z",boundary%vert_lut%z(1,:,:,:))
! 		call io_write3di("vlutz2.nc","z",boundary%vert_lut%z(2,:,:,:))
! 		call io_write3d("vlutw1.nc","w",boundary%vert_lut%w(1,:,:,:))
! 		call io_write3d("vlutw2.nc","w",boundary%vert_lut%w(2,:,:,:))
!
! 		call io_write3di("u_vlutz.nc","z",boundary%u_geo%vert_lut%z(1,:,:,:))
! 		call io_write3d("u_vlutw.nc","w",boundary%u_geo%vert_lut%w(1,:,:,:))
! 		call io_write3di("v_vlutz.nc","z",boundary%v_geo%vert_lut%z(1,:,:,:))
! 		call io_write3d("v_vlutw.nc","w",boundary%v_geo%vert_lut%w(1,:,:,:))
!
! 		call io_write3d("bc_hires_z_u.nc","data",boundary%u_geo%z)
! 		call io_write3d("bc_hires_z_v.nc","data",boundary%v_geo%z)
! 		call io_write3d("domain_z_u.nc","data",domain%u_geo%z)
! 		call io_write3d("domain_z_v.nc","data",domain%v_geo%z)
! 		call io_write3d("bc_lowresz.nc","data",boundary%lowres_z)
! 		call io_write3d("bc_hiresz.nc","data",boundary%z)
! 		call io_write3d("domain_z.nc","data",domain%z)

		call swap_z(boundary)
		
	end subroutine init_bc
end module
