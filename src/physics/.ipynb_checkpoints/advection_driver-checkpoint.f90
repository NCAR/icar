!> ----------------------------------------------------------------------------
!!  Driver to call different advection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module advection
    use data_structures
    use icar_constants
    use adv_upwind,                 only : upwind, upwind_init, upwind_var_request
    use adv_mpdata,                 only : mpdata, mpdata_init
    ! use debug_module,               only : domain_fix
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t

    implicit none
    private
    public :: advect, adv_init, adv_var_request
contains

    subroutine adv_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options

        if (this_image()==1) write(*,*) ""
        if (this_image()==1) write(*,*) "Initializing Advection"
        if (options%physics%advection==kADV_UPWIND) then
            if (this_image()==1) write(*,*) "    Upwind"
            call upwind_init(domain,options)
        else if(options%physics%advection==kADV_MPDATA) then
            if (this_image()==1) write(*,*) "    MP-DATA"
            call mpdata_init(domain,options)
        endif

    end subroutine adv_init

    subroutine adv_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%advection==kADV_UPWIND) then
            call upwind_var_request(options)
        else if (options%physics%advection==kADV_MPDATA) then
            call upwind_var_request(options)
        endif

    end subroutine

    subroutine advect(domain,options,dt)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,intent(in) :: dt
        ! integer :: nx, nz, ny
        !
        ! nx=size(domain%p,1)
        ! nz=size(domain%p,2)
        ! ny=size(domain%p,3)
        !
        ! if (.not.allocated(domain%tend%qv_adv)) then
        !     allocate(domain%tend%qv_adv(nx,nz,ny))
        !     domain%tend%qv_adv=0
        ! endif


        if (options%physics%advection==kADV_UPWIND) then
            call upwind(domain,options,dt)
        else if(options%physics%advection==kADV_MPDATA) then
            call mpdata(domain,options,dt)
        endif

        ! if (options%advect_density) then
        !     call domain_fix(domain)
        ! endif
    end subroutine advect

end module advection
