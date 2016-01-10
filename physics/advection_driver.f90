!> ----------------------------------------------------------------------------
!!  Driver to call different advection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module advection
    use data_structures
    use adv_upwind, only : upwind
    use adv_mpdata, only : mpdata, mpdata_init


    implicit none
    private
    public::advect,adv_init
contains
    
    subroutine adv_init(domain,options)
        implicit none
        type(domain_type), intent(inout) :: domain
        type(options_type),intent(in) :: options

!         if (options%physics%advection==kADV_UPWIND) then
!             call upwind(domain,options,dt)
        if(options%physics%advection==kADV_MPDATA) then
            call mpdata_init(domain,options)
        endif
        
    end subroutine adv_init
    
    subroutine advect(domain,options,dt)
        implicit none
        type(domain_type), intent(inout) :: domain
        type(options_type),intent(in) :: options
        real,intent(in) :: dt
        integer :: nx, nz, ny
        
        nx=size(domain%p,1)
        nz=size(domain%p,2)
        ny=size(domain%p,3)
        
        if (.not.allocated(domain%tend%qv_adv)) then
            allocate(domain%tend%qv_adv(nx,nz,ny))
            domain%tend%qv_adv=0
        endif
            
        
        if (options%physics%advection==kADV_UPWIND) then
            call upwind(domain,options,dt)
        elseif(options%physics%advection==kADV_MPDATA) then
            call mpdata(domain,options,dt)
        endif
        
    end subroutine advect

end module advection
