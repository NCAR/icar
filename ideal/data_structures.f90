module data_structures
	use, intrinsic :: iso_c_binding

	integer,parameter::MAXFILELENGTH=100
	integer,parameter::MAXVARLENGTH=100
	
	type position
		integer::x,y
	end type position
	type fourpos
		integer::x(4),y(4)
	end type fourpos
	
	type geo_look_up_table
		integer,allocatable,dimension(:,:,:)::x,y
		real,allocatable,dimension(:,:,:)::w
	end type geo_look_up_table
	
	
	type interpolable_type
		real, allocatable, dimension(:,:) :: lat,lon
		type(geo_look_up_table)::geolut
	end type interpolable_type
		
	type, extends(interpolable_type) :: wind_type
		real, allocatable, dimension(:,:,:) :: u,v
		integer :: nfiles
	end type wind_type
	
	type, extends(interpolable_type) :: linearizable_type
		real, allocatable, dimension(:,:,:):: u,v,dz,z
		real, allocatable, dimension(:,:) :: terrain,dzdx,dzdy
		complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: fzs
		real::dx
	end type linearizable_type
	
	type, extends(linearizable_type) :: domain_type
		real, allocatable, dimension(:,:,:) :: p,th,w
		real, allocatable, dimension(:,:,:) :: qv,cloud,ice,nice,qrain,nrain,qsnow,qgrau
		real, allocatable, dimension(:,:) :: rain,snow,graupel
		real::dt
	end type domain_type

	type, extends(linearizable_type) :: bc_type
		real, allocatable, dimension(:,:,:) :: p,th,qv
		real, allocatable, dimension(:,:,:) :: dudt,dvdt,dwdt,dpdt,dthdt,dqvdt,dqcdt
		type(domain_type)::next_domain
		type(wind_type)::ext_winds
	end type bc_type

	
	type physics_type
		integer::microphysics
		integer::advection
		integer::boundarylayer
		integer::landsurface
		integer::radiation
		integer::convection
		integer::windtype
	end type physics_type
		
	type options_type
		character (len=MAXFILELENGTH) :: init_conditions_file
		character (len=MAXFILELENGTH), allocatable::boundary_files(:),ext_wind_files(:)
		character (len=MAXFILELENGTH) :: output_file,restart_file
		character (len=MAXVARLENGTH) :: latvar,lonvar
		logical :: readz, debug, external_winds,remove_lowres_linear,mean_winds,mean_fields,restart
		integer :: buffer=0
		integer :: ntimesteps,nz,nfiles,ext_winds_nfiles,restart_step
		real :: dx,io_dt,outputinterval,dz
		type(physics_type)::physics
	end type options_type
end module data_structures	