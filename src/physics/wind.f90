!>------------------------------------------------------------
!! Module to manage the ICAR wind field, including calls to linear winds
!! importantly it also rotates the wind field into the ICAR grid and
!! balances the U, V, and W fields for "mass" conservation
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module wind
    use icar_constants
    use data_structures

    use linear_theory_winds, only : linear_perturb
    ! use mod_blocking,        only : add_blocked_flow

    use domain_interface,  only : domain_t
    use options_interface, only : options_t

    use io_routines, only : io_read

    implicit none
    private
    public::update_winds, init_winds, wind_var_request
    real, parameter::deg2rad=0.017453293 !2*pi/360
contains



        subroutine wind_linear_var_request(options)
            implicit none
            type(options_t), intent(inout) :: options

            ! List the variables that are required to be allocated for the linear wind solution
            call options%alloc_vars( &
                            [kVARS%nsquared,    kVARS%potential_temperature,   kVARS%exner,            &
                             kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,      &
                             kVARS%u,           kVARS%v,                       kVARS%w,                &
                             kVARS%dz ])

            ! List the variables that are required to be advected
            ! call options%advect_vars( &
            !               [, &
            !                , &
            !                ] )

            ! List the variables that are required for restarts with the linear wind solution
            call options%restart_vars( &
                            [kVARS%nsquared,    kVARS%potential_temperature,                           &
                             kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,      &
                             kVARS%u,           kVARS%v,                       kVARS%w,                &
                             kVARS%dz ])

        end subroutine

        subroutine wind_var_request(options)
            implicit none
            type(options_t), intent(inout) :: options

            if (options%physics%windtype == kWIND_LINEAR) then
                call wind_linear_var_request(options)
            endif

        end subroutine wind_var_request




    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx + dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------
    subroutine balance_uvw(u,v,w, dz, dx, options)
        implicit none
        real,           intent(inout) :: u(:,:,:), v(:,:,:), w(:,:,:)
        real,           intent(in)    :: dz(:,:,:)
        real,           intent(in)    :: dx
        type(options_t),intent(in)    :: options

        real, allocatable, dimension(:,:) :: du, dv, divergence, rhou, rhov,rhow
        real, allocatable, dimension(:,:,:) :: dzu, dzv

        integer :: k, ims, ime, jms, jme, kms, kme, i,j

        ! associate(u => domain%u%data_3d,  &
        !           v => domain%v%data_3d,  &
        !           w => domain%w%data_3d  )

        ims = lbound(w,1)
        ime = ubound(w,1)
        kms = lbound(w,2)
        kme = ubound(w,2)
        jms = lbound(w,3)
        jme = ubound(w,3)

        w = 0

        !------------------------------------------------------------
        ! These could be module level to prevent lots of allocation/deallocation/reallocations
        ! but these are relatively small, and should be allocated efficiently on the stack
        !------------------------------------------------------------
        ! if (options%advect_density) then
        !     allocate(rhou(nx-1,ny-2))
        !     allocate(rhov(nx-2,ny-1))
        !     allocate(rhow(nx-2,ny-2))
        ! endif

        allocate(du(ims+1:ime-1,jms+1:jme-1))
        allocate(dv(ims+1:ime-1,jms+1:jme-1))
        allocate(divergence(ims+1:ime-1,jms+1:jme-1))
        allocate(dzv(ims:ime,kms:kme,jms:jme))
        allocate(dzu(ims:ime,kms:kme,jms:jme))

        dzv = 0
        dzu = 0

        dzv(:,:,jms+1:jme) = (dz(:,:,jms+1:jme) + dz(:,:,jms:jme-1)) / 2
        dzu(ims+1:ime,:,:) = (dz(ims+1:ime,:,:) + dz(ims:ime-1,:,:)) / 2

        ! If this becomes a bottle neck in the code it could be parallelized over y
        ! loop over domain levels

        do k = kms, kme
            !------------------------------------------------------------
            ! If we are incorporating density into the advection equation
            ! then it needs to be incorporated when balancing the wind field
            !
            ! Note that the "else" case below does the same thing without density
            ! and is much easier to understand
            !------------------------------------------------------------

            ! this is the else, advect density is not supported at the moment
            !------------------------------------------------------------
            ! If we are not incorporating density this is simpler
            !------------------------------------------------------------
            ! calculate horizontal divergence
            !   in the North-South direction
            dv =  v(ims+1:ime-1, k, jms+2:jme)   * dzv(ims+1:ime-1, k, jms+2:jme)   &
                - v(ims+1:ime-1, k, jms+1:jme-1) * dzv(ims+1:ime-1, k, jms+1:jme-1)
            !   in the East-West direction
            du =  u(ims+2:ime, k, jms+1:jme-1)   * dzu(ims+2:ime, k, jms+1:jme-1)   &
                - u(ims+1:ime-1, k, jms+1:jme-1) * dzu(ims+1:ime-1, k, jms+1:jme-1)
            !   in net
            divergence = du + dv

            ! Then calculate w to balance
            if (k==kms) then
                ! if this is the first model level start from 0 at the ground
                ! note the are out for w is dx^2, but there is a dx in the divergence term that is dropped to balance
                w(ims+1:ime-1,k,jms+1:jme-1) = 0 - divergence / dx
            else
                ! else calculate w as a change from w at the level below
                w(ims+1:ime-1,k,jms+1:jme-1) = w(ims+1:ime-1,k-1,jms+1:jme-1) - divergence / dx
            endif

            !------------------------------------------------------------
            ! Now do the same for the convective wind field if needed
            !------------------------------------------------------------
            ! if (options%physics%convection > 0) then
            !     ! calculate horizontal divergence
            !     dv = domain%v_cu(2:nx-1,i,3:ny) - domain%v_cu(2:nx-1,i,2:ny-1)
            !     du = domain%u_cu(3:nx,i,2:ny-1) - domain%u_cu(2:nx-1,i,2:ny-1)
            !     divergence = du + dv
            !     ! Then calculate w to balance
            !     if (i==1) then
            !         ! if this is the first model level start from 0 at the ground
            !         domain%w_cu(2:nx-1,i,2:ny-1) = 0 - divergence
            !     else
            !         ! else calculate w as a change from w at the level below
            !         domain%w_cu(2:nx-1,i,2:ny-1) = domain%w_cu(2:nx-1,i-1,2:ny-1) - divergence
            !     endif
            ! endif

        enddo
        ! end associate

    end subroutine balance_uvw

    !>------------------------------------------------------------
    !! Correct for a grid that is locally rotated with respect to EW,NS
    !!
    !! Assumes forcing winds are EW, NS relative, not grid relative.
    !!
    !!------------------------------------------------------------
    subroutine make_winds_grid_relative(u, v, w, sintheta, costheta)
        real, intent(inout) :: u(:,:,:), v(:,:,:), w(:,:,:)
        double precision, intent(in)    :: sintheta(:,:), costheta(:,:)

        real, dimension(:), allocatable :: u_local,v_local

        integer :: k, j, ims, ime, jms, jme, kms, kme

        ims = lbound(w,1)
        ime = ubound(w,1)
        kms = lbound(w,2)
        kme = ubound(w,2)
        jms = lbound(w,3)
        jme = ubound(w,3)

        allocate(u_local(ims:ime))
        allocate(v_local(ims:ime))
        !assumes u and v come in on a staggered Arakawa C-grid with one additional grid point in x/y for u/v respectively
        ! destagger to a centered grid (the mass grid)
        u(:ime,:,:) = (u(:ime,:,:) + u(ims+1:,:,:))/2
        v(:,:,:jme) = (v(:,:,:jme) + v(:,:,jms+1:))/2

        do j = jms, jme
            do k = kms, kme
                ! rotate wind field to the real grid
                u_local = u(ims:ime,k,j) * costheta(:,j) + v(ims:ime,k,j) * sintheta(ims:ime,j)
                v_local = v(ims:ime,k,j) * costheta(:,j) + u(ims:ime,k,j) * sintheta(ims:ime,j)
                u(:ime,k,j) = u_local
                v(:ime,k,j) = v_local
            enddo
        enddo
        deallocate(u_local, v_local)

        ! put the fields back onto a staggered grid, having effectively lost two grid cells in the staggered directions
        ! estimate the "lost" grid cells by extrapolating beyond the remaining
        u(ims+1:ime,:,:) =  (u(ims:ime-1,:,:) + u(ims+1:ime,:,:))/2
        u(ims,:,:)       = 2*u(ims,:,:)       - u(ims+1,:,:)
        u(ime+1,:,:)     = 2*u(ime,:,:)       - u(ime-1,:,:)

        v(:,:,jms+1:jme) =  (v(:,:,jms:jme-1) + v(:,:,jms+1:jme))/2
        v(:,:,jms)       = 2*v(:,:,jms)       - v(:,:,jms+1)
        v(:,:,jme+1)     = 2*v(:,:,jme)       - v(:,:,jme-1)

    end subroutine


    !>------------------------------------------------------------
    !! Apply wind field physics and adjustments
    !!
    !! This will call the linear wind module if necessary, otherwise it just updates for
    !! This should ONLY be called once for each forcing step, otherwise effects will be additive.
    !!
    !!------------------------------------------------------------
    subroutine update_winds(domain, options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options

        real, allocatable, dimension(:,:,:) :: temparray
        integer :: nx, ny, nz, i, j

        if (.not.allocated(domain%sintheta)) then
            call init_winds(domain, options)

            ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
            call make_winds_grid_relative(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%sintheta, domain%costheta)

            ! flow blocking parameterization
            ! if (options%block_options%block_flow) then
            !     call add_blocked_flow(domain, options)
            ! endif

            ! linear winds
            if (options%physics%windtype==kWIND_LINEAR) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density, update=.False.)
            ! simple acceleration over topography
            elseif (options%physics%windtype==kCONSERVE_MASS) then
                call mass_conservative_acceleration(domain%u%data_3d, domain%v%data_3d, domain%zr_u, domain%zr_v)
            endif
            ! else assumes even flow over the mountains

            ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)
            call balance_uvw(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%advection_dz, domain%dx, options)

        else

            ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
            call make_winds_grid_relative(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%w%meta_data%dqdt_3d, domain%sintheta, domain%costheta)

            ! linear winds
            if (options%physics%windtype==kWIND_LINEAR) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density, update=.True.)
            ! simple acceleration over topography
            elseif (options%physics%windtype==kCONSERVE_MASS) then
                call mass_conservative_acceleration(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%zr_u, domain%zr_v)
            endif
            ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)
            call balance_uvw(domain% u %meta_data%dqdt_3d,      &
                             domain% v %meta_data%dqdt_3d,      &
                             domain% w %meta_data%dqdt_3d,      &
                             domain% advection_dz,              &
                             domain% dx,                        &
                             options)

        endif

    end subroutine update_winds

    subroutine mass_conservative_acceleration(u, v, u_accel, v_accel)
        implicit none
        real, intent(inout) :: u(:,:,:)
        real, intent(inout) :: v(:,:,:)
        real, intent(in)    :: u_accel(:,:,:)
        real, intent(in)    :: v_accel(:,:,:)

        u = u / u_accel
        v = v / v_accel

    end subroutine mass_conservative_acceleration

    !>------------------------------------------------------------
    !! Setup initial fields (i.e. grid relative rotation fields)
    !!
    !!------------------------------------------------------------
    subroutine init_winds(domain,options)
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options

        integer :: i, j, ims, ime, jms, jme, kms, kme
        integer :: starti, endi
        double precision :: dist, dlat, dlon

        real, allocatable :: temporary_2d(:,:)

        call allocate_winds(domain)

        if (options%parameters%fixed_dz_advection) then
            do i=domain%grid%kms, domain%grid%kme
                domain%advection_dz(:,i,:) = options%parameters%dz_levels(i)
            enddo
        else
            domain%advection_dz = domain%dz_interface%data_3d
        endif



        if (options%parameters%sinalpha_var /= "") then
            ims = lbound(domain%latitude%data_2d, 1)
            ime = ubound(domain%latitude%data_2d, 1)
            jms = lbound(domain%latitude%data_2d, 2)
            jme = ubound(domain%latitude%data_2d, 2)

            if (this_image()==1) print*, "Reading Sinalph/cosalpha"

            call io_read(options%parameters%init_conditions_file, options%parameters%sinalpha_var, temporary_2d)
            domain%sintheta = temporary_2d(ims:ime, jms:jme)

            call io_read(options%parameters%init_conditions_file, options%parameters%cosalpha_var, temporary_2d)
            domain%costheta = temporary_2d(ims:ime, jms:jme)

            deallocate(temporary_2d)
        else

            associate(lat => domain%latitude%data_2d,  &
                      lon => domain%longitude%data_2d  )

            ims = lbound(lat,1)
            ime = ubound(lat,1)
            jms = lbound(lat,2)
            jme = ubound(lat,2)
            do j = jms, jme
                do i = ims, ime
                    ! in case we are in the first or last grid, reset boundaries
                    starti = max(ims, i-2)
                    endi   = min(ime, i+2)

                    ! change in latitude
                    dlat = DBLE(lat(endi,j)) - lat(starti,j)
                    ! change in longitude
                    dlon = DBLE(lon(endi,j) - lon(starti,j)) * cos(deg2rad*DBLE(lat(i,j)))
                    ! distance between two points
                    dist = sqrt(DBLE(dlat)**2 + DBLE(dlon)**2)

                    ! sin/cos of angles for use in rotating fields later
                    domain%costheta(i, j) = abs(dlon / dist)
                    domain%sintheta(i, j) = (-1) * dlat / dist

                enddo
            enddo

            end associate
        endif
        if (options%parameters%debug .and.(this_image()==1)) then
            print*, ""
            print*, "Domain Geometry"
            print*, "MAX / MIN SIN(theta) (ideally 0)"
            print*, "   ", maxval(domain%sintheta), minval(domain%sintheta)
            print*, "MAX / MIN COS(theta) (ideally 1)"
            print*, "   ", maxval(domain%costheta), minval(domain%costheta)
            print*, ""
        endif


    end subroutine init_winds

    !>------------------------------------------------------------
    !! Allocate memory used in various wind related routines
    !!
    !!------------------------------------------------------------
    subroutine allocate_winds(domain)
        type(domain_t), intent(inout) :: domain
        integer :: ims, ime, jms, jme, kms, kme

        ims = lbound(domain%latitude%data_2d, 1)
        ime = ubound(domain%latitude%data_2d, 1)
        jms = lbound(domain%latitude%data_2d, 2)
        jme = ubound(domain%latitude%data_2d, 2)
        kms = lbound(domain%w%data_3d, 2)
        kme = ubound(domain%w%data_3d, 2)

        if (.not.allocated(domain%sintheta)) then
            allocate(domain%sintheta(ims:ime, jms:jme))
            domain%sintheta = 0
        endif
        if (.not.allocated(domain%costheta)) then
            allocate(domain%costheta(ims:ime, jms:jme))
            domain%costheta = 0
        endif

        if (.not.allocated(domain%advection_dz)) then
            allocate(domain%advection_dz(ims:ime,kms:kme,jms:jme))
        endif

        ! note w is special cased because it does not have a forcing variable, so it is not necessarily allocated automatically
        if (.not.associated(domain%w%meta_data%dqdt_3d)) then
            allocate(domain%w%meta_data%dqdt_3d(ims:ime,kms:kme,jms:jme))
            domain%w%meta_data%dqdt_3d = 0
        endif

        ! if (.not.allocated(domain%dzdx)) then
        !     allocate(domain%dzdx(nx-1,ny))
        ! endif
        ! if (.not.allocated(domain%dzdy)) then
        !     allocate(domain%dzdy(nx,ny-1))
        ! endif

    end subroutine allocate_winds

    !>------------------------------------------------------------
    !! Provides a routine to deallocate memory allocated in allocate_winds
    !!
    !!------------------------------------------------------------
    ! subroutine finalize_winds(domain)
    !     type(domain_t), intent(inout) :: domain
    !
    !     if (allocated(domain%sintheta)) then
    !         deallocate(domain%sintheta)
    !     endif
    !     if (allocated(domain%costheta)) then
    !         deallocate(domain%costheta)
    !     endif
    !     if (allocated(domain%dzdx)) then
    !         deallocate(domain%dzdx)
    !     endif
    !     if (allocated(domain%dzdy)) then
    !         deallocate(domain%dzdy)
    !     endif
    !
    ! end subroutine finalize_winds
end module wind
