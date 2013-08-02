module data_structures
	use, intrinsic :: iso_c_binding

	integer,parameter::MAXFILELENGTH=100
	integer,parameter::MAXVARLENGTH=100
	
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
	
! 	type to contain external wind fields, only real addition is nfiles... maybe this could be folded in elsewhere?
	type, extends(interpolable_type) :: wind_type
		real, allocatable, dimension(:,:,:) :: u,v
		integer :: nfiles
	end type wind_type

! 	generic linearizable type so we can add linear wind field to domain or remove it from low-res (BC) U/V
	type, extends(interpolable_type) :: linearizable_type
! 		linear theory computes u,v at z
		real, allocatable, dimension(:,:,:):: u,v,dz,z
		real, allocatable, dimension(:,:) :: terrain,dzdx,dzdy
		complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: fzs
		real::dx
	end type linearizable_type
	
! 	All fields needed in the domain
	type, extends(linearizable_type) :: domain_type
		real, allocatable, dimension(:,:,:) :: p,th,w
		real, allocatable, dimension(:,:,:) :: qv,cloud,ice,nice,qrain,nrain,qsnow,qgrau
		real, allocatable, dimension(:,:) :: rain,snow,graupel
		real::dt
	end type domain_type

! 	boundary conditions type, must be linearizable so we can remove low res linear wind field
	type, extends(linearizable_type) :: bc_type
! 		not sure these are used anymore...
		real, allocatable, dimension(:,:,:) :: p,th,qv
! 		wind and pressure dXdt fields applied to full 3d grid, others applied only to boundaries
		real, allocatable, dimension(:,:,:) :: dudt,dvdt,dwdt,dpdt,dthdt,dqvdt,dqcdt
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
		character (len=MAXVARLENGTH) :: latvar,lonvar,uvar,vvar,pvar,thvar,qvvar,qcvar,qivar,qrvar,qsvar,qgvar
! 		various boolean options
		logical :: readz, debug, external_winds,remove_lowres_linear,mean_winds,mean_fields,restart
! 		buffer to remove from all sides of the high res grid supplied
		integer :: buffer=0
! 		various integer parameters/options
		integer :: ntimesteps,nz,nfiles,ext_winds_nfiles,restart_step
! 		various real parameters/options
		real :: dx,io_dt,outputinterval,dz
! 		defines which physics package to be used. 
		type(physics_type)::physics
	end type options_type
end module data_structures	