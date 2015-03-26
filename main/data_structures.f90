!>------------------------------------------------
!! 
!! Contains type definitions for a variety of model data strucutres
!! Also defines model constants (e.g. gravity, and MAXFILELENGTH)
!!
!! General Field Definitions
!!
!! ---- 3D fields ---- NX x NZ x NY
!! u     = wind in east direction        					[m/s]
!! v     = wind in north direction       					[m/s]
!! w     = wind in vertical direction    					[m/s] (possibly scaled by dx/dz)
!! 
!! p     = pressure                      					[pa]
!! th    = potential temperature         					[K]
!!
!! qv    = water vapor (mixing ratio)    					[kg/kg]
!! cloud = cloud water                   					[kg/kg]
!! ice   = cloud ice                     					[kg/kg]
!! qrain = rain mixing ratio             					[kg/kg]
!! qsnow = snow mixing ratio             					[kg/kg]
!! qgrau = graupel mixing ratio          					[kg/kg]
!! nice  = ice number concentration      					[1/cm^3]
!! nrain = rain number concentration     					[1/cm^3]
!!
!! ---- 2D fields ---- NX x NY
!! 		---- moisture fluxes ----
!! rain  = rain+crain+snow+graupel       					[mm]
!! crain = convective rain at surface    					[mm]
!! snow  = snow at surface               					[mm]
!! graupel = graupel at surface          					[mm]
!!
!! 		---- energy fluxes ----
!! sensible_heat = Sensible heat flux from surface			[W/m^2]
!! latent_heat   = Latent heat flux from surface				[W/m^2]
!! pbl_height    = Height of the planetary boundary layer	[m]
!!
!! 		---- Radiation variables ----
!! cloudfrac		= Cloud fraction 							[0-1]
!! swdown		= Shortwave down at land surface			[W/m^2]
!! lwdown		= Longwave down at land surface				[W/m^2]
!! lwup			= Lonwave up from the land surface			[W/m^2]
!!
!! ---- Land Surface variables ---- 
!! 	3D fields ---- NX x NZ x NY
!! soil_t 		= 3D Soil temperature						[K]
!! soil_vwc		= 3D Soil volumetric water content			[m^3/m^3]
!!
!!   2D fields ---- NX x NY
!! skin_t 		= Land surface skin temperature				[K]
!! soil_tdeep	= Temperature at the soil column bottom		[K]
!! vegfrac		= vegetation cover fraction 				[%]
!! snow_swe 		= Snow water equivalent on the land surface	[mm]
!! soil_totalmoisture = Soil column total water content 		[mm]
!! soil_type 	= Soil type (index for USGS classification in SOILPARM.TBL)	[1-19]
!! veg_type 		= Vegetation type (index for VEGPARM.TBL)					[1-27]
!! landmask      = Map of Land vs Water grid cells			[0,1,2]
!!
!! ---- NOTE ----
!! dX_dt variables are the increment in boundary conditions between internal model time steps
!! some of these are 2d, some are 3d
!! 
!! ---- model structure ----
!! terrain = surface elevation           [m]
!! z = model layer height (at mid point) [m]
!! dz = layer thickness                  [m]
!!
!! sintheta = sine of the angle between grid and geographic coords   []
!! costheta = cosine of the angle between grid and geographic coords []
!! fzs      = buffered FFT(terrain) for linear wind calculations   
!!
!!	Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module data_structures
	use, intrinsic :: iso_c_binding ! needed for fftw compatible complex types
	implicit none

!------------------------------------------------
! Model constants (string lengths)
!------------------------------------------------
	integer,parameter::MAXFILELENGTH=100 ! maximum file name length
	integer,parameter::MAXVARLENGTH=100  ! maximum variable name length
	
!------------------------------------------------
! Physical Constants
!------------------------------------------------
	real, parameter :: LH_vaporization=2260000.0 ! J/kg
	! should be calculated as 2.5E6 + (2106.0 - 4218.0)*temp_degC ?
    real, parameter :: R   = 287.058 ! J/(kg K) specific gas constant for air
    real, parameter :: cp  = 1012.0  ! J/kg/K   specific heat capacity of moist STP air? 
	real, parameter :: gravity= 9.81    ! m/s^2    gravity
	real, parameter :: pi  = 3.1415927 ! pi
	real, parameter :: stefan_boltzmann = 5.67e-8 ! the Stefan-Boltzmann constant
	
!------------------------------------------------
! 	various data structures for use in geographic interpolation routines
!------------------------------------------------
	! contains the location of a specific grid point
	type position
		integer::x,y
	end type position
	! contains location of surrounding 4 grid cells
	type fourpos
		integer::x(4),y(4)
	end type fourpos
	
	! a geographic look up table for spatial interpolation, from x,y with weight w
	type geo_look_up_table
		! x,y index positions, [n by m by 4] where there are 4 surrounding low-res points 
		! for every high resolution point grid point to interpolate to
		integer,allocatable,dimension(:,:,:)::x,y
		! weights to use for each of the 4 surrounding gridpoints.  Sum(over axis 3) must be 1.0
		real,allocatable,dimension(:,:,:)::w
	end type geo_look_up_table

	!------------------------------------------------
	! A look up table for vertical interpolation. from z with weight w
	!------------------------------------------------
	type vert_look_up_table
		! z index positions for all x,y,z points (x 2 for above and below z levels)
		integer,allocatable,dimension(:,:,:,:)::z
		! weights to use for each of the two surrounding points.  Sum (over axis 1) must be 1.0
		real,allocatable,dimension(:,:,:,:)::w
	end type vert_look_up_table

	!------------------------------------------------
	! generic interpolable type so geo interpolation routines will work on winds, domain, or boundary conditions. 
	!------------------------------------------------
	type interpolable_type
		real, allocatable, dimension(:,:) :: lat,lon
		real, allocatable, dimension(:,:,:) :: z
		type(vert_look_up_table)::vert_lut
		type(geo_look_up_table)::geolut
		logical :: dx_errors_printed=.False.
		logical :: dy_errors_printed=.False.
	end type interpolable_type


	!------------------------------------------------
	! type to contain external wind fields, only real addition is nfiles... maybe this could be folded in elsewhere?
	!------------------------------------------------
	type, extends(interpolable_type) :: wind_type
		real, allocatable, dimension(:,:,:) :: u,v
		type(interpolable_type)				:: u_geo,v_geo
		real, allocatable, dimension(:,:) 	:: terrain,dzdx,dzdy
		real :: dx
		integer :: nfiles
	end type wind_type

	!------------------------------------------------
	! generic linearizable type so we can add linear wind field to domain or remove it from low-res (BC) U/V
	!------------------------------------------------
	type, extends(interpolable_type) :: linearizable_type
		! linear theory computes u,v at z.  Trying rho to mitigate boussinesq approx... 
		real, allocatable, dimension(:,:,:)	:: u,v,dz,rho,th
		type(interpolable_type)				:: u_geo,v_geo
		real, allocatable, dimension(:,:)	:: terrain,dzdx,dzdy
		complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: fzs !FFT(terrain)
		real::dx
	end type linearizable_type
	
	!------------------------------------------------
	! All fields needed in the domain defined in detail above
	!------------------------------------------------
	type, extends(linearizable_type) :: domain_type
		! 3D atmospheric fields
		real, allocatable, dimension(:,:,:) :: p,w,pii,ur,vr,wr
		real, allocatable, dimension(:,:,:) :: qv,cloud,ice,nice,qrain,nrain,qsnow,qgrau
		! 3D atmospheric field tendencies
		real, allocatable, dimension(:,:,:) :: qv_adv_tendency,qv_pbl_tendency
		! 3D soil field
		real, allocatable, dimension(:,:,:) :: soil_t, soil_vwc
		! 2D fields, primarily fluxes to/from the land surface
		! precip fluxes
		real, allocatable, dimension(:,:)   :: rain,crain,snow,graupel
		real, allocatable, dimension(:,:)   :: current_rain, current_snow
		! radiative fluxes (and cloud fraction)
		real, allocatable, dimension(:,:)   :: swdown, lwdown, cloudfrac, lwup
		! turbulent fluxes (and ground heat flux)
		real, allocatable, dimension(:,:)   :: sensible_heat,latent_heat,ground_heat
		! domain parameters (and PBL height)
		real, allocatable, dimension(:,:)   :: pbl_height,landmask ! store PBL height (if available) and the land-sea mask
		real, allocatable, dimension(:,:)   :: sintheta, costheta !rotations about the E-W, N-S grid
		! land surface state and parameters
		real, allocatable, dimension(:,:)   :: soil_tdeep, skin_t, soil_totalmoisture, snow_swe
		real, allocatable, dimension(:,:)   :: vegfrac,canopy_water
		integer, allocatable, dimension(:,:):: soil_type,veg_type
		! current model time step length (should this be somewhere else?)
		real::dt
		! current model time (seconds from options%time_zero)
		real*8 :: model_time
	end type domain_type

	!------------------------------------------------
	! boundary conditions type, must be linearizable so we can remove low res linear wind field
	!------------------------------------------------
	type, extends(linearizable_type) :: bc_type
		! not sure these are used anymore...
		real, allocatable, dimension(:,:,:) :: p,qv
		! dX_dt variables are the change in variable X between two forcing time steps
		! wind and pressure dX_dt fields applied to full 3d grid, others applied only to boundaries
		real, allocatable, dimension(:,:,:) :: du_dt,dv_dt,dw_dt,dp_dt,drho_dt,dth_dt,dqv_dt,dqc_dt
		! sh, lh, and pblh fields are only 2d. These are only used with LSM option 1 and are derived from forcing file
		real, allocatable, dimension(:,:) :: dsh_dt,dlh_dt,dpblh_dt
		! store the low resolution versionf of terrain and atmospheric elevations
		real,allocatable,dimension(:,:)::lowres_terrain
		real,allocatable,dimension(:,:,:)::lowres_z
		! store the full high-res 3D grid for the next time step to compute dXdt fields
		! includes high res versions of low res terrain and z
		type(domain_type)::next_domain
		! if we are using an external wind field, store them here temporarily... 
		! does this need to be separate from next_domain other than the nfiles attribute?
		type(wind_type)::ext_winds
	end type bc_type

	!------------------------------------------------
	! type to store integer options for each physics package
	!------------------------------------------------
	type physics_type
		integer::microphysics
		integer::advection
		integer::boundarylayer
		integer::landsurface
		integer::radiation
		integer::convection
		integer::windtype
	end type physics_type

        !! ++ trude
        type mp_options_type
                 real :: Nt_c
                 real :: TNO
                 real :: am_s
                 real :: rho_g
                 real :: av_s, bv_s, fv_s, av_i
                 real :: av_g, bv_g
                 real :: Ef_si, Ef_rs, Ef_rg, Ef_ri
                 real :: C_cubes, C_sqrd
                 real :: mu_r
                 real :: t_adjust
                 logical :: Ef_rw_l, EF_sw_l
       end type mp_options_type
        !! -- trude
	!------------------------------------------------
	! store all model options
	!------------------------------------------------
	type options_type
                character (len=MAXVARLENGTH) :: version,comment

		! file names
		character (len=MAXFILELENGTH) :: init_conditions_file
		character (len=MAXFILELENGTH), dimension(:), allocatable::boundary_files,ext_wind_files
		character (len=MAXFILELENGTH) :: output_file,restart_file,output_file_frequency

		! variable names from init/BC/wind/... files
		character (len=MAXVARLENGTH) :: landvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon, &
										hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi, &
										pvar,pbvar,tvar,qvvar,qcvar,qivar,qrvar,qsvar,qgvar,hgtvar, &
										shvar,lhvar,pblhvar,zvar, &
										soiltype_var, soil_t_var,soil_vwc_var,soil_deept_var, &
										vegtype_var,vegfrac_var
!! ++ trude
 		character(len=MAXFILELENGTH) :: mp_options_filename
!! -- trude

		! various boolean options
		logical :: ideal 				! this is an ideal simulation, forcing will be held constant
		logical :: readz 				! read atmospheric grid elevations from file
		logical :: readdz				! read atm model layer thicknesses from namelist
		logical :: debug 				! outputs a little more information at runtime (not much at present)
		logical :: external_winds 		! read a high res 3d wind field from an external file (e.g. a high res WRF run)
		logical :: remove_lowres_linear ! attempt to remove the linear mountain wave from the forcing low res model
		logical :: mean_winds 			! use only a mean wind field across the entire model domain
		logical :: mean_fields			! use only a mean forcing field across the model boundaries 
		logical :: restart 				! this is a restart run, read model conditions from a restart file
		logical :: add_low_topo 		! option to add low resolution topography back to the high res model to 
										! mitigate low-res wind field downdrafts over mountains (don't use)
		logical :: advect_density 		! properly incorporate density into the advection calculations. 
										! Doesn't play nice with linear winds
		logical :: high_res_soil_state  ! read the soil state from the high res input file not the low res file
		logical :: variable_N			! Compute the Brunt Vaisala Frequency (N^2) every time step

		integer :: buffer				! buffer to remove from all sides of the high res grid supplied
		integer :: ymin,ymax,xmin,xmax 	! never implemented : would permit buffers of different distances on all sides
		integer :: vert_smooth 			! number of model levels to smooth winds over in the vertical
		! various integer parameters/options
		integer :: ntimesteps 			! total number of time steps to be simulated
		integer :: nz 					! number of model vertical levels
		integer :: nfiles 				! number of forcing files to read from namelist
		integer :: ext_winds_nfiles 	! number of extrenal wind filenames to read from namelist
		integer :: restart_step 		! step in forcing data to begin running
		integer :: restart_date(6) 		! date to initialize from (y,m,d, h,m,s)
		integer :: restart_step_in_file ! step in restart file to initialize from

		! various real parameters/options
		real :: dx 						! grid cell width [m]
		real :: dxlow 					! forcing model grid cell width [m]
		real :: in_dt 					! time step between forcing inputs [s]
		real :: out_dt					! time step between output [s]
		real :: outputinterval 			! time steps per output
		real :: inputinterval  			! time steps per input
		real :: smooth_wind_distance 	! distance over which to smooth the forcing wind field (m)
		real :: N_squared				! static Brunt Vaisala Frequency (N^2) to use
		real :: linear_contribution     ! fractional contribution of linear perturbation to wind field (e.g. u_hat multiplied by this)
		logical :: spatial_linear_fields! use a spatially varying linear wind perturbation
		
		! date/time parameters
		double precision :: initial_mjd ! Modified Julian Day of the first model time step [days]
		double precision :: time_zero   ! Starting model initial time step (mjd-50000)*3600 [s]

		real :: t_offset				! offset to temperature because WRF outputs potential temperature-300
		real, allocatable, dimension(:)::dz_levels ! model layer thicknesses to be read from namelist
		real :: rotation_scale_height   ! height to minimize wind rotation into the terrain following grid below [m]
		logical :: use_agl_height       ! interpolate from forcing to model layers using Z above ground level, not sea level
		
		! defines which physics package to be used. 
		type(physics_type)::physics
!! ++ trude
		! parameterization options
                type(mp_options_type)::mp_options
!! -- trude
		integer :: warning_level        ! level of warnings to issue when checking options settings 0-10.  
										! 0  = Don't print anything
										! 1  = print serious warnings
		! (DEFAULT if debug=True)		! 2  = print all warnings
										! 3-4 ... nothing specified
		! (DEFAULT if debug=False)		! 5  = Stop for options that are likely to break the model (print all warnings) 
										! 6-8... nothing specified
										! 9  = stop on serious warnings 
										! 10 = stop on all warnings
	end type options_type
end module data_structures	
