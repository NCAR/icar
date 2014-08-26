!----------------------------------------------------------
!
! Very simple land surface model code
! 
! Rain is partitioned into infiltration and runoff
! Snow is accumulated on the surface, then melts, runsoff, or sublimates
! Soil moisture is permitted to be lost to ET or subsurface flow
! 
! ET, Sensible Heat Flux, and Longwave are partitioned using Penman Monteith.
!
! The entry point to the code is lsm_simple. 
!
! Call tree graph :
! lsm_simple->
! 	[->],
! 	[->],
! 	[->]
! 
! High level routine descriptions / purpose
!   ls,_simple    		- loops over X,Y grid cells, calls a, b, c
!   a           - calculates 
!   b           - calculates 
!   x           	- calculates 
! 
! Driver inputs: p,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz
!   p   = pressure                      - 3D - input  - Pa		- (nx,nz,ny)
!   th  = potential temperature         - 3D - in/out - K		- (nx,nz,ny)
!   pii = inverse exner function        - 3D - input  - []		- (nx,nz,ny)
!   rho = air density                   - 3D - input  - kg/m^3	- (nx,nz,ny)
!   qv  = specific humidity             - 3D - in/out - kg/kg	- (nx,nz,ny)
!   rain= rainfall                      - 2D - input  - mm  	- (nx,ny)
!   snow= snowfall                      - 2D - input  - mm	    - (nx,ny)
!   wind= wind speed                    - 2D - input  - m/s  	- (nx,ny)
!   swdown = shortwave down at surface	- 2D - input  - W/m^2	- (nx,ny)
!   lwdown = longwave down at surface	- 2D - input  - W/m^2	- (nx,ny)
!   dt = time step                      - 0D - input  - seconds	- scalar
!
!----------------------------------------------------------
module module_lsm_simple
	use data_structures
	implicit none
contains
	subroutine lsm_simple_init(domain,options)
		implicit none
		type(domain_type), intent(in) :: domain
		type(options_type),intent(in)    :: options
		integer :: nx,ny
		
		write(*,*) "Initializing simple land surface model"
		
	end subroutine lsm_simple_init
	
	
	subroutine lsm_simple(theta,pii,qv,rain,snow,p,swdown,lwdown, wind, &
						sensible_heat, latent_heat, tskin, tsoil, vwc, swe, options,dt)
		implicit none
		real,dimension(:,:,:), intent(inout) :: theta, qv
		real,dimension(:,:,:), intent(in) :: pii,p
		real,dimension(:,:), intent(in) :: rain,snow,swdown,lwdown,wind
		real,dimension(:,:), intent(inout) :: sensible_heat, latent_heat, tskin, tsoil, vwc, swe
		type(options_type),intent(in)    :: options
		real, intent(in) :: dt
		integer :: nx,ny,j,k,nz
		
		
		
	end subroutine lsm_simple
end module module_lsm_simple
