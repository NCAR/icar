!> ----------------------------------------------------------------------------
!!
!!  Driver to call different advection schemes
!!
!!	Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module advection
	use data_structures
	use adv_upwind, only : upwind
	use adv_mpdata, only : mpdata

	implicit none
	private
	public::advect
contains
	subroutine advect(domain,options,dt)
		type(domain_type), intent(inout) :: domain
		type(options_type),intent(in) :: options
		real,intent(in) :: dt
		
		if (options%physics%advection==1) then
			call upwind(domain,options,dt)
		elseif(options%physics%advection==2) then
			call mpdata(domain,options,dt)
		endif
! 		NOTE if advection==0 (e.g.) do nothing!
		
	end subroutine advect

end module advection
