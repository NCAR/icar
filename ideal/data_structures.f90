module data_structures
	use, intrinsic :: iso_c_binding
	
	type position
		integer::x,y
	end type position
	type four_pos
		integer::x(4),y(4)
	end type four_pos
	
	type geo_look_up_table
		integer,allocatable,dimension(:,:,:)::x,y
		real,allocatable,dimension(:,:,:)::w
	end type geo_look_up_table
	
	type domain_type
		real, allocatable, dimension(:,:,:) :: u,v,w,p,th,dz,z
		real, allocatable, dimension(:,:,:) :: qv,cloud,ice,nice,qrain,nrain,qsnow,qgrau
		real, allocatable, dimension(:,:) :: terrain,rain,snow,graupel,dzdx,dzdy
		real, allocatable, dimension(:,:) :: lat,lon
		complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: fzs
		real::dx,dt
	end type domain_type
	
	type bc_type
		real, allocatable, dimension(:,:) :: lat,lon,terrain
		type(geo_look_up_table)::geolut
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
		character (len=100) :: init_conditions_file
		character (len=100) :: output_file
		character (len=100) :: latvar,lonvar
		integer :: ntimesteps,timestep,outputinterval
		real :: dx,dt
		type(physics_type)::physics
	end type options_type
end module data_structures	