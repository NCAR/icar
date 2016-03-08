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
    use module_mp_thompson, only:   mp_gt_driver, thompson_init
    use module_mp_simple,   only:   mp_simple_driver
    implicit none

    ! permit the microphysics to update on a longer time step than the advection
    integer :: update_interval
    real*8 :: last_model_time
    ! temporary variables
    real,allocatable,dimension(:,:) :: SR, last_rain, last_snow, this_precip
    integer, parameter :: npoints=8
    real,    dimension(npoints) :: dist_fraction = [ 0.1,0.15,0.1, 0.15,0.15, 0.1,0.15,0.1]
    integer, dimension(npoints) :: x_list = [ -1,0,1, -1,1, -1,0,1]
    integer, dimension(npoints) :: y_list = [ 1,1,1, 0,0, -1,-1,-1]
contains
    
    
    !>----------------------------------------------------------
    !! Initialize microphysical routines
    !!
    !! This routine will call the initialization routines for the specified microphysics packages. 
    !! It also initializes any module level variables, e.g. update_interval
    !! 
    !! @param   options     ICAR model options to specify required initializations
    !!
    !!----------------------------------------------------------
    subroutine mp_init(options)
        implicit none
        type(options_type), intent(in)::options
        
        write(*,*) "Initializing Microphysics"
        if (options%physics%microphysics==kMP_THOMPSON) then
            write(*,*) "    Thompson Microphysics"
            call thompson_init(options%mp_options)
        endif

        update_interval = options%mp_options%update_interval
        last_model_time = -999
    end subroutine mp_init
    
    !>----------------------------------------------------------
    !! Distribute the microphysics precipitation to neighboring grid cells
    !!
    !! Because ICAR can be too aggressive at putting precip on mountain tops, this
    !! routine smooths out the precip by keeping only a fraction of it locally, and 
    !! distributing the rest to the neighboring grid cells, weighted by distance. 
    !! 
    !! @param   [inout]current_precip   accumulated model precip at this time step
    !! @param   [in]last_precip         accumulated model precip prior to microphysics call
    !! @param   [in]local_fraction      fraction of precip to maintain in the local gridcell
    !!
    !!----------------------------------------------------------
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
    
    !>----------------------------------------------------------
    !! Microphysical driver
    !!
    !! This routine handles calling the individual microphysics routine specified, that 
    !! includes creating and passing any temporary variables, and checking when to update
    !! the microphysics based on the specified update_interval. 
    !! 
    !! @param   domain      ICAR model domain structure
    !! @param   options     ICAR model options structure
    !! @param   dt_in       Current driving time step (this is the advection step)
    !! @param   model_time  Current model time (to check if it will exceed the update_interval)
    !!
    !!----------------------------------------------------------
    subroutine mp(domain,options,dt_in, model_time)
        implicit none
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        real,intent(in)::dt_in
        double precision, intent(in) :: model_time
        real :: mp_dt
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

        ! if this is the first time mp is called, set last time such that mp will update
        if (last_model_time==-999) then
            last_model_time = (model_time-max(real(update_interval),dt_in))
        endif
        
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

        ! only run the microphysics if the next time step would put it over the update_interval time
        if (((model_time+dt_in)-last_model_time)>=update_interval) then
            ! calculate the actual time step for the microphysics
            mp_dt = model_time-last_model_time
            ! reset the counter so we know that *this* is the last time we've run the microphysics
            last_model_time = model_time
            
            ! If we are going to distribute the current precip over a few grid cells, we need to keep track of
            ! the last_precip so we know how much fell
            if (options%mp_options%local_precip_fraction<1) then
                last_rain=domain%rain
                last_snow=domain%snow
            endif
            
            ! run the thompson microphysics
            if (options%physics%microphysics==kMP_THOMPSON) then
                kts=kds
                kte=kde
                ! set the current tile to the top layer to process microphysics for
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
                ! call the thompson microphysics
                write(*,*) "Befor call thomspon mp"
                call mp_gt_driver(domain%qv, domain%cloud, domain%qrain, domain%ice, &
                                domain%qsnow, domain%qgrau, domain%nice, domain%nrain, &
                                domain%th, domain%pii, domain%p, domain%dz_inter, mp_dt, itimestep, &
                                domain%rain, last_rain, &       ! last_rain is not used in thompson
                                domain%snow, last_snow, &       ! last_snow is not used in thompson
                                domain%graupel, this_precip, &  ! this_precip is not used in thompson
                                SR, &
                                ids,ide, jds,jde, kds,kde, &    ! domain dims
                                ids,ide, jds,jde, kds,kde, &    ! memory dims
                                its,ite, jts,jte, kts,kte)      ! tile dims
                write(*,*) "After call thomspon mp"
                                
            elseif (options%physics%microphysics==kMP_SB04) then
                ! call the simple microphysics routine of SB04
                call mp_simple_driver(domain%p,domain%th,domain%pii,domain%rho,domain%qv,domain%cloud, &
                                domain%qrain,domain%qsnow,domain%rain,domain%snow,&
                                mp_dt,domain%dz,ide,jde,kde)
            endif
            
            if (options%mp_options%local_precip_fraction<1) then
                call distribute_precip(domain%rain, last_rain, options%mp_options%local_precip_fraction)
                call distribute_precip(domain%snow, last_snow, options%mp_options%local_precip_fraction)
            endif
        endif
        
    end subroutine mp
    
    !>----------------------------------------------------------
    !! Finalize microphysical routines
    !!
    !! This routine will call the finalization routines (if any) for the specified microphysics packages. 
    !! It also deallocates any module level variables, e.g. SR
    !! 
    !! @param   options     ICAR model options to specify required initializations
    !!
    !!----------------------------------------------------------
    subroutine mp_finish(options)
        implicit none
        type(options_type),intent(in)::options
        
        if (allocated(SR)) then
            deallocate(SR)
        endif
        if (allocated(last_snow)) then
            deallocate(last_snow)
        endif
        if (allocated(last_rain)) then
            deallocate(last_rain)
        endif
        if (allocated(this_precip)) then
            deallocate(this_precip)
        endif
    end subroutine mp_finish
end module microphysics
