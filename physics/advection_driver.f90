!> ----------------------------------------------------------------------------
!!
!!  Driver to call different advection schemes
!!
!!  Author: Ethan Gutmann (gutmann@ucar.edu)
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
        
        if (options%physics%advection==ADV_UPWIND) then
            call upwind(domain,options,dt)
        elseif(options%physics%advection==ADV_MPDATA) then
            call mpdata(domain,options,dt)
        endif
        
    end subroutine advect

end module advection
