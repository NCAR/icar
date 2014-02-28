module planetary_boundary_layer
	use data_structures
	use pbl_simple,    only : simple_pbl, finalize_simple_pbl, init_simple_pbl
	implicit none
	
	private
	public :: pbl_init, pbl, pbl_finalize
! 	these are now defined in data_structures.f90
! 	real, parameter :: LH_vaporization=2260000.0 ! J/kg
! 	real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
! 	real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
! 	real, parameter :: g=9.81 ! gravity m/s^2

contains
	subroutine pbl_init(domain,options)
		implicit none
		type(domain_type),intent(in)::domain
		type(options_type),intent(in)::options
		if (options%physics%boundarylayer==2) then
			call init_simple_pbl(domain,options)
		endif
	end subroutine pbl_init
	
	subroutine pbl(domain,options,dt_in)
		implicit none
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		real,intent(in)::dt_in
		
		if (options%physics%boundarylayer==2) then
			call simple_pbl(domain,dt_in)
		endif
						
	end subroutine pbl
	
	subroutine pbl_finalize(options)
		implicit none
		type(options_type),intent(in)::options
		if (options%physics%boundarylayer==2) then
			call finalize_simple_pbl()
		endif
	end subroutine pbl_finalize
end module planetary_boundary_layer