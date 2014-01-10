module data_structures
	use, intrinsic :: iso_c_binding

	integer,parameter::MAXFILELENGTH=100
	integer,parameter::MAXVARLENGTH=100

	real, parameter :: LH_vaporization=2260000.0 ! J/kg
    real, parameter :: R=287.058   ! J/(kg K) specific gas constant for air
    real, parameter :: cp = 1012.0 ! J/kg/K   specific heat capacity of moist STP air? 
	real, parameter :: g = 9.81    ! m/s^2    gravity
	
! 	various data structures for use in geographic interpolation routines
! 	contains the location of a specific grid point
	type position
		integer::x,y
	end type position
! 	contains location of surrounding 4 grid cells
	type fourpos
		integer::x(4),y(4)
	end type fourpos
	
! 	a geographic look up table for spatial interpolation, from x,y with weight w
	type geo_look_up_table
		integer,allocatable,dimension(:,:,:)::x,y
		real,allocatable,dimension(:,:,:)::w
	end type geo_look_up_table
	
!   generic interpolable type so geo interpolation routines will work on winds, domain, or boundary conditions. 	
	type interpolable_type
		real, allocatable, dimension(:,:) :: lat,lon
		type(geo_look_up_table)::geolut
	end type interpolable_type

	!------------------------------------------------
! 
! General Field Definitions
!
! ---- 3D fields ----
! u     = wind in east direction        [m/s]
! v     = wind in north direction       [m/s]
! w     = wind in vertical direction    [m/s] (scaled by dx/dz)
! 
! p     = pressure                      [pa]
! th    = potential temperature         [K]
!
! qv    = vapor pressure (mixing ratio) [kg/kg]
! cloud = cloud water                   [kg/kg]
! ice   = cloud ice                     [kg/kg]
! qrain = rain mixing ratio             [kg/kg]
! qsnow = snow mixing ratio             [kg/kg]
! qgrau = graupel mixing ratio          [kg/kg]
! nice  = ice number concentration      [1/cm^3]
! nrain = rain number concentration     [1/cm^3]
!
! ---- 2D fields ----
! rain  = rain+snow+graupel at surface  [mm]
! snow  = snow at surface               [mm]
! graupel = graupel at surface          [mm]
!
! ---- model structure ----
! terrain = surface elevation           [m]
! z = model layer height (at mid point) [m]
! dz = layer thickness                  [m]
!------------------------------------------------

! 	type to contain external wind fields, only real addition is nfiles... maybe this could be folded in elsewhere?
	type, extends(interpolable_type) :: wind_type
		real, allocatable, dimension(:,:,:) :: u,v
		type(interpolable_type)::u_geo,v_geo
		real, allocatable, dimension(:,:) :: terrain,dzdx,dzdy
		real :: dx
		integer :: nfiles
	end type wind_type

! 	generic linearizable type so we can add linear wind field to domain or remove it from low-res (BC) U/V
	type, extends(interpolable_type) :: linearizable_type
! 		linear theory computes u,v at z
		real, allocatable, dimension(:,:,:):: u,v,dz,z
		type(interpolable_type)::u_geo,v_geo
		real, allocatable, dimension(:,:) :: terrain,dzdx,dzdy
		complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: fzs !FFT(terrain)
		real::dx
	end type linearizable_type
	
! 	All fields needed in the domain
	type, extends(linearizable_type) :: domain_type
		real, allocatable, dimension(:,:,:) :: p,th,w
		real, allocatable, dimension(:,:,:) :: qv,cloud,ice,nice,qrain,nrain,qsnow,qgrau
		real, allocatable, dimension(:,:) :: rain,crain,snow,graupel,sensible_heat,latent_heat,pbl_height
		real::dt
	end type domain_type

! 	boundary conditions type, must be linearizable so we can remove low res linear wind field
	type, extends(linearizable_type) :: bc_type
! 		not sure these are used anymore...
		real, allocatable, dimension(:,:,:) :: p,th,qv
! 		dXdt variables are the change in variable X between two forcing time steps
! 		wind and pressure dXdt fields applied to full 3d grid, others applied only to boundaries
		real, allocatable, dimension(:,:,:) :: dudt,dvdt,dwdt,dpdt,dthdt,dqvdt,dqcdt
! 		sh, lh, and pblh fields are only 2d
		real, allocatable, dimension(:,:) :: dshdt,dlhdt,dpblhdt
		real,allocatable,dimension(:,:)::lowres_terrain
		real,allocatable,dimension(:,:,:)::lowres_z
! 		store the full 3D grid for the next time step to compute dXdt fields
		type(domain_type)::next_domain
! 		if we are using external winds, store them here temporarily... does this need to be separate from next_domain other than nfiles?
		type(wind_type)::ext_winds
	end type bc_type

! 	type to store integer options for each physics package (not all used at present)
	type physics_type
		integer::microphysics
		integer::advection
		integer::boundarylayer
		integer::landsurface
		integer::radiation
		integer::convection
		integer::windtype
	end type physics_type
	
! 	store all model options
	type options_type
! 		file names
		character (len=MAXFILELENGTH) :: init_conditions_file
		character (len=MAXFILELENGTH), allocatable::boundary_files(:),ext_wind_files(:)
		character (len=MAXFILELENGTH) :: output_file,restart_file
! 		variable names from init/BC/wind/... files
		character (len=MAXVARLENGTH) :: latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon, &
										hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi, &
										pvar,tvar,qvvar,qcvar,qivar,qrvar,qsvar,qgvar,hgtvar, &
										shvar,lhvar,pblhvar
! 		various boolean options
		logical :: ideal,readz, debug, external_winds,remove_lowres_linear,&
				   mean_winds,mean_fields,restart,add_low_topo
! 		buffer to remove from all sides of the high res grid supplied
		integer :: buffer=0
! 		various integer parameters/options
		integer :: ntimesteps,nz,nfiles,ext_winds_nfiles,restart_step
! 		various real parameters/options
		real :: dx,dxlow,in_dt,out_dt,outputinterval,inputinterval
! 		defines which physics package to be used. 
		type(physics_type)::physics
	end type options_type
end module data_structures	