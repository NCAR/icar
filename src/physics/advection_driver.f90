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
    use adv_std,                    only : adv_std_init, adv_std_var_request, adv_std_advect3d, adv_std_compute_wind
    use adv_mpdata,                 only : mpdata_init, mpdata_advect3d, mpdata_compute_wind
    use adv_fluxcorr,               only : init_fluxcorr
    ! use debug_module,               only : domain_fix
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t

    implicit none
    private

    public :: advect, adv_init, adv_var_request
contains

    subroutine adv_init(domain,options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in) :: options

        if (this_image()==1) write(*,*) ""
        if (this_image()==1) write(*,*) "Initializing Advection"
        if (options%physics%advection==kADV_STD) then
            if (this_image()==1) write(*,*) "    Standard"
            call adv_std_init(domain,options)
        else if(options%physics%advection==kADV_MPDATA) then
            if (this_image()==1) write(*,*) "    MP-DATA"
            call mpdata_init(domain,options)
        endif
        
        if (options%adv_options%flux_corr > 0) call init_fluxcorr(domain)

        
    end subroutine adv_init

    subroutine adv_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        !if (options%physics%advection==kADV_UPWIND) then
        !    call upwind_var_request(options)
        if (options%physics%advection==kADV_MPDATA) then
            call adv_std_var_request(options)
        else
            call adv_std_var_request(options)
        endif
    end subroutine
    
    subroutine advect(domain, options, dt)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,intent(in) :: dt

        type(variable_t) :: var_to_advect
        real, allocatable :: temp(:,:,:)

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
        
        !Allocate storage variable for temp-quantities
        if (options%time_options%RK3) allocate(temp(domain%ims:domain%ime,domain%kms:domain%kme,domain%jms:domain%jme))

        if (options%physics%advection==kADV_STD) then
            call adv_std_compute_wind(domain,options,dt)
        else if(options%physics%advection==kADV_MPDATA) then
            call mpdata_compute_wind(domain,options,dt)
        endif
        

        !Loop through all vars to advect
                
        ! make sure the dictionary is reset to point to the first variable
        call domain%adv_vars%reset_iterator()
        ! Now iterate through the dictionary as long as there are more elements present
        do while (domain%adv_vars%has_more_elements())
            ! get the next variable
            var_to_advect = domain%adv_vars%next()

            if (options%time_options%RK3) then
            
                if (options%physics%advection==kADV_STD) then
                
                    !Initial advection-tendency calculations
                    temp = var_to_advect%data_3d
                    call adv_std_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, domain%jacobian,t_factor_in=0.333)
                    call adv_std_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, domain%jacobian,t_factor_in=0.5)
                                            
                    !final advection call with tendency-fluxes
                    call adv_std_advect3d(temp,var_to_advect%data_3d, domain%advection_dz, domain%jacobian,flux_corr=options%adv_options%flux_corr)
                                            
                    var_to_advect%data_3d = temp                             
                                                                 
                else if(options%physics%advection==kADV_MPDATA) then
                
                    ! Not yet implemented (is it compatable w/ RK3?)
                endif
            else
                if (options%physics%advection==kADV_STD) then
                    call adv_std_advect3d(var_to_advect%data_3d,var_to_advect%data_3d, domain%advection_dz, domain%jacobian,flux_corr=options%adv_options%flux_corr)
                else if(options%physics%advection==kADV_MPDATA) then                                    
                    call mpdata_advect3d(var_to_advect%data_3d, domain%density%data_3d, domain%jacobian, domain%advection_dz, options)
                endif
            endif

        enddo

    end subroutine advect

end module advection
