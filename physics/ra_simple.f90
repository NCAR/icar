!----------------------------------------------------------
!
! Very simple radiation code modeled after the description
! in Finch and Best(2004?) of Reiff et al. 1984 shortwave and Idso and Jackson (1969) longwave. 
! 
! Clearsky Shortwave radiation is calculated as a function of day of year and time of day. 
! Cloudy Shortwave is calculated as clearsky SW * f(cloud cover) [0.25-1]
! 
! Cloud cover is calculated as in Xu and Randal (1996) as f(surface_RH, qc+qs+qr)
! 
! Clearsky Longwave radiation is calculated as f(Tair)
! Cloudy longwave is scaled up by a f(cloud cover) [1-1.2]
!
! The entry point to the code is ra_simple. 
!
! Call tree graph :
! ra_simple->
! 	[cloudfrac->],
! 	[shortwave->],
! 	[longwave->]
! 
! High level routine descriptions / purpose
!   ra_simple    		- loops over X,Y grid cells, calls cloudfrac, shortwave,longwave on columns
!   cloudfrac           - calculates the cloud fraction following Xu and Randall (1996)
!   shortwave           - calculates shortwave at the surface following Reiff et al (1984)
!   longwave           	- calculates longwave at the surface following Idso and Jackson (1969)
! 
! Driver inputs: p,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz
!   p   = pressure                      - 3D - input  - Pa		- (nx,nz,ny)
!   th  = potential temperature         - 3D - in/out - K		- (nx,nz,ny)
!   pii = inverse exner function        - 3D - input  - []		- (nx,nz,ny)
!   rho = air density                   - 3D - input  - kg/m^3	- (nx,nz,ny)
!   qv  = specific humidity             - 3D - input  - kg/kg	- (nx,nz,ny)
!   qc  = cloud water content           - 3D - input  - kg/kg	- (nx,nz,ny)
!   qr  = rain water content            - 3D - input  - kg/kg	- (nx,nz,ny)
!   qs  = snow water content            - 3D - input  - kg/kg	- (nx,nz,ny)
!   swdown = shortwave down at surface	- 2D - output - W/m^2	- (nx,ny)
!   lwdown = longwave down at surface	- 2D - output - W/m^2	- (nx,ny)
!   dt = time step                      - 0D - input  - seconds	- scalar
!
!----------------------------------------------------------
module module_ra_simple
	use data_structures
	implicit none
	real, allocatable, dimension(:,:) :: cos_lat_m,sin_lat_m
	real, parameter :: So=1367.0 ! Solar "constant" W/m^2
contains
	subroutine ra_simple_init(domain,options)
		type(domain_type), intent(in) :: domain
		type(options_type),intent(in)    :: options
		integer :: nx,ny
		
		write(*,*) "Initializing simple radiation"
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
		allocate(cos_lat_m(nx,ny))
		allocate(sin_lat_m(nx,ny))
		
		cos_lat_m=cos(domain%lat)
		sin_lat_m=sin(domain%lat)
		
	end subroutine ra_simple_init

! 	day_of_year=calc_day_of_year(date)
! 	time_of_day=calc_time_of_day(date)

	
	function relative_humidity(t,qv,p,j,nx)
		real,intent(out),dimension(nx) :: relative_humidity
		real,intent(in),dimension(nx)::t
		real,intent(in),dimension(:,:,:)::qv
		integer,intent(in)::j,nx
		real,dimension(nx) :: mr, e, es
		
		! convert sensible humidity to mixing ratio
		mr=qv(:,1,j)/(1-qv(:,1,j))
		! convert mixing ratio to vapor pressure
		e=mr*p(:,1,j)/(0.62197+mr)
		! convert temperature to saturated vapor pressure
		es=6.112*exp(17.67*(t-273.15)/(t-29.65))
		! finally return relative humidity
		relative_humidity= e/es
	end function relative_humidity
	
	function shortwave(rh, cloud_cover, solar_elevation)
		real, intent(inout) :: rh,cloud_cover,solar_elevation
		
	end subroutine shortwave
	
	subroutine ra_simple(theta,pii,qv,qc,qs,qr,p,swdown,lwdown,cloudfrac,lat,lon,date,options,dt)
		real,dimension(:,:,:), intent(inout) :: theta
		real,dimension(:,:,:), intent(in) :: pii,qv,qc,qs,qr,p
		real,dimension(:,:), intent(out) :: swdown,lwdown,cloudfrac
		real,dimension(:,:), intent(in) :: lat,lon
		double precision, intent(in) :: date
		type(options_type),intent(in)    :: options
		real, intent(in) :: dt
		integer :: nx,ny,j
		real, allocatable, dimension(:) :: rh,T_air,solar_elevation
		real :: solar_elevation
		
		
		nx=size(lat,1)
		ny=size(lat,2)
		
		allocate(rh(nx))
		allocate(T_air(nx))
		allocate(solar_elevation(nx))
		do j=2,ny-1
			solar_elevation=calc_solar_elevation(date,lon,j,nx)
			T_air=theta(:,1,j)*pii(:,1,j)
			rh=calc_rh(T_air,qv,p,j,nx)
			
			swdown(:,j) = shortwave(rh,cloud_cover,solar_elevation)
			call longwave(T_air)
		end do
		
	end subroutine ra_simple
end module module_ra_simple
