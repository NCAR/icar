!>----------------------------------------------------------
!! This module provides a wrapper to call various microphysics models
!! It sets up variables specific to the physics package to be used including 
!! history variables not currently stored in the domain level data structure
!!
!! The main entry point to the code is mp(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  mp_init->[ external initialization routines]
!!  mp->[   external microphysics routines]
!!  mp_finish
!! 
!! High level routine descriptions / purpose
!!   mp_init            - allocates module data and initializes physics package
!!   mp                 - sets up and calls main physics package
!!   mp_finish          - deallocates module memory, place to do the same for physics
!! 
!! Inputs: domain, options, dt
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module microphysics
    use data_structures
    use module_mp_thompson, only: mp_gt_driver,thompson_init
    use module_mp_simple, only:mp_simple_driver
    implicit none
!   these are now defined in data_structures.f90
!   real, parameter :: LH_vaporization=2260000.0 ! J/kg
!   real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
!   real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
!   real, parameter :: g=9.81 ! gravity m/s^2

    real,allocatable,dimension(:,:) :: SR, last_rain, last_snow, this_precip
    integer, parameter :: npoints=8
    real,    dimension(npoints) :: dist_fraction = [ 0.1,0.15,0.1, 0.15,0.15, 0.1,0.15,0.1]
    integer, dimension(npoints) :: x_list = [ -1,0,1, -1,1, -1,0,1]
    integer, dimension(npoints) :: y_list = [ 1,1,1, 0,0, -1,-1,-1]
contains
    subroutine mp_init(options)
        implicit none
        type(options_type), intent(in)::options
        
        write(*,*) "Initializing Microphysics"
        if (options%physics%microphysics==kMP_THOMPSON) then
            write(*,*) "    Thompson Microphysics"
            call thompson_init(options%mp_options)
        endif
        
    end subroutine mp_init
    
    subroutine distribute_precip(current_precip, last_precip, local_fraction)
        real, dimension(:,:), intent(inout) :: current_precip, last_precip
        real, intent(in) :: local_fraction
        ! relies on module variable this_precip as a temporary array
        
        integer :: i,j, nx,ny
        integer :: x,y, point
        
        nx=size(current_precip,1)
        ny=size(current_precip,2)
        
        do j=2,ny-1
            do i=2,nx-1
                this_precip(i,j) = current_precip(i,j)-last_precip(i,j)
                current_precip(i,j) = last_precip(i,j)
            end do
        end do
        do j=2,ny-1
            do i=2,nx-1
                current_precip(i,j) = current_precip(i,j)+this_precip(i,j)*local_fraction
                do point=1,npoints
                    x = i + x_list(point)
                    y = j + y_list(point)
                    current_precip(i,j) = current_precip(i,j) + this_precip(x,y) * (1-local_fraction) * dist_fraction(point)
                end do
            end do
        end do
                
        
        
    end subroutine distribute_precip
    
    subroutine mp(domain,options,dt_in)
        implicit none
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        real,intent(in)::dt_in
        integer ::ids,ide,jds,jde,kds,kde,itimestep=1
        integer ::its,ite,jts,jte,kts,kte, nx,ny
        
        ids=1
        ide=size(domain%qv,1)
        nx=ide
        kds=1
        kde=size(domain%qv,2)
        jds=1
        jde=size(domain%qv,3)
        ny=jde
        
        ! snow rain ratio
        if (.not.allocated(SR)) then
            allocate(SR(nx,ny))
            SR=0
        endif
        ! last snow amount
        if (.not.allocated(last_snow)) then
            allocate(last_snow(nx,ny))
            last_snow=0
        endif
        ! last rain amount
        if (.not.allocated(last_rain)) then
            allocate(last_rain(nx,ny))
            last_rain=0
        endif
        ! temporary precip amount
        if (.not.allocated(this_precip)) then
            allocate(this_precip(nx,ny))
            this_precip=0
        endif
        
        if (options%physics%microphysics==kMP_THOMPSON) then
            kts=kds
            kte=kde
            if (options%mp_options%top_mp_level>0) then
                kte=min(kte, options%mp_options%top_mp_level)
            endif
            if (options%ideal) then
                ! for ideal runs process the boundaries as well to be consistent with WRF
                its=ids;ite=ide
                jts=jds;jte=jde
            else
                its=ids+1;ite=ide-1
                jts=jds+1;jte=jde-1
            endif
            last_rain=domain%rain
            last_snow=domain%snow
            call mp_gt_driver(domain%qv, domain%cloud, domain%qrain, domain%ice, &
                            domain%qsnow, domain%qgrau, domain%nice, domain%nrain, &
                            domain%th, domain%pii, domain%p, domain%dz_inter, dt_in, itimestep, &
                            domain%rain, domain%rain, &
                            domain%snow, domain%snow, &
                            domain%graupel, domain%graupel, &
                            SR, &
                            ids,ide, jds,jde, kds,kde, &    ! domain dims
                            ids,ide, jds,jde, kds,kde, &    ! memory dims
                            its,ite, jts,jte, kts,kte)      ! tile dims
        elseif (options%physics%microphysics==kMP_SB04) then
            call mp_simple_driver(domain%p,domain%th,domain%pii,domain%rho,domain%qv,domain%cloud, &
                            domain%qrain,domain%qsnow,domain%rain,domain%snow,&
                            dt_in,domain%dz,ide,jde,kde)
        endif
        
        if (options%mp_options%local_precip_fraction<1) then
            call distribute_precip(domain%rain, last_rain, options%mp_options%local_precip_fraction)
            call distribute_precip(domain%snow, last_snow, options%mp_options%local_precip_fraction)
        endif
                        
    end subroutine mp
    
    subroutine mp_finish()
        implicit none
        if (allocated(SR)) then
            deallocate(SR)
        endif
    end subroutine mp_finish
end module microphysics
