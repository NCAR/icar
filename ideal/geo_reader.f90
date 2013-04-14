module geo_reader
	implicit none
	type geolut
		integer,allocatable,dimension(:,:,:)::x,y
		real,allocatable,dimension(:,:,:)::w
	end type geolut
	type geo_info
		real,allocatable,dimension(:,:)::lat,lon
	end type geo_info
	type reader_type
		real, allocatable, dimension(:,:,:) :: datavar
		real, allocatable, dimension(:,:,:) :: qv,cloud,ice,nice,qrain,nrain,qsnow,qgrau
		real, allocatable, dimension(:,:) :: terrain,rain,snow,graupel,dzdx,dzdy
		type(geolut)::geolut
		real::dx,dt
	end type reader_type
	
contains
	subroutine geo_LUT(geohi, geolo)
		implicit none
		type(domain_type),intent(in)::domain
		type(bc_type),intent(in)::boundary
		integer :: nx,ny
		
		ny=size(geolo%lat,1)
		nx=size(geolo%lat,2)
		
		allocate(geo_LUT%x(4,ny,nx))
		allocate(geo_LUT%y(4,ny,nx))
		allocate(geo_LUT%w(4,ny,nx))
		
		geo_LUT%x=0
		geo_LUT%y=0
		geo_LUT%w=0.0
	end function geo_LUT
		
end module geo_reader