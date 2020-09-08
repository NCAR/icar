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
    use linear_theory_winds, only : linear_perturb
    ! use mod_blocking,        only : add_blocked_flow
    use data_structures
    use exchangeable_interface,   only : exchangeable_t
    use domain_interface,  only : domain_t
    use options_interface, only : options_t

    use io_routines, only : io_read

    implicit none
    private
    public::update_winds, init_winds
    real, parameter::deg2rad=0.017453293 !2*pi/360
contains

    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx + dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------
    subroutine balance_uvw(u,v,w, dzdx,dzdy,dz,dx,jaco,smooth_height, options)
        implicit none
        real,           intent(inout) :: u(:,:,:), v(:,:,:), w(:,:,:)
        real,           intent(in)    :: dzdx(:,:,:), dzdy(:,:,:), dz(:,:,:), jaco(:,:,:)
        real,           intent(in)    :: dx, smooth_height
        type(options_t),intent(in)    :: options

        real, allocatable, dimension(:,:) :: rhou, rhov, rhow
        real, allocatable, dimension(:,:,:) :: divergence

        integer :: ims, ime, jms, jme, kms, kme, k

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

        allocate(divergence(ims:ime,kms:kme,jms:jme))

        call calc_divergence(divergence,u,v,w,dzdx,dzdy,dz,dx,jaco,smooth_height,horz_only=.True.)

        write(*,*)"abs of div: ",maxval(ABS(divergence(ims+1:ime-1,:,jms+1:jme-1)))

        ! If this becomes a bottle neck in the code it could be parallelized over y
        ! loop over domain levels
        do k = kms,kme
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
            !   in the East-West direction
            !   in net

            ! Then calculate w to balance
                ! if this is the first model level start from 0 at the ground
                ! note the are out for w is dx^2, but there is a dx in the divergence term that is dropped to balance
                if (k==kms) then
                    w(:,k,:) = 0 - divergence(:,k,:) * dz(:,k,:)
                else
                    w(:,k,:) = w(:,k-1,:) - divergence(:,k,:) * dz(:,k,:)
                endif
            end do
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

        ! end associate

    end subroutine balance_uvw


    subroutine calc_divergence(div, u, v, w, dzdx, dzdy, dz, dx, jaco, smooth_height,horz_only)
        implicit none
        real,           intent(inout) :: div(:,:,:)
        real,           intent(in)    :: u(:,:,:), v(:,:,:), w(:,:,:), dz(:,:,:), dzdx(:,:,:), dzdy(:,:,:), jaco(:,:,:)
        real,           intent(in)    :: dx, smooth_height
        logical, optional, intent(in)  :: horz_only
        ! type(options_t),intent(in)    :: options

        real, allocatable, dimension(:,:,:) :: diff_U, diff_V, slope_U, slope_V
        integer :: ims, ime, jms, jme, kms, kme, k
        logical :: horz

        horz = .False.
        if (present(horz_only)) horz=horz_only

        ims = lbound(w,1)
        ime = ubound(w,1)
        kms = lbound(w,2)
        kme = ubound(w,2)
        jms = lbound(w,3)
        jme = ubound(w,3)

        allocate(diff_U(ims:ime,kms:kme,jms:jme))
        allocate(diff_V(ims:ime,kms:kme,jms:jme))
        allocate(slope_U(ims:ime,kms:kme,jms:jme))
        allocate(slope_V(ims:ime,kms:kme,jms:jme))

        diff_U = u(ims+1:ime+1, :, jms:jme) - u(ims:ime, :, jms:jme)
        diff_V = v(ims:ime, :, jms+1:jme+1) - v(ims:ime, :, jms:jme)

        do k = kms,kme
            slope_U(:,k,:) = (u(ims+1:ime+1, k, jms:jme)*dzdx(ims+1:ime+1, kms, jms:jme) + &
                              u(ims:ime, k, jms:jme)*dzdx(ims:ime, kms, jms:jme))/2
            slope_V(:,k,:) = (v(ims:ime, k, jms+1:jme+1)*dzdy(ims:ime, kms, jms+1:jme+1) + &
                              v(ims:ime, k, jms:jme)*dzdy(ims:ime, kms, jms:jme))/2
                              
            !slope_U(:,k,:) = (u(ims+1:ime+1, k, jms:jme)*0.2 + &
            !                  u(ims:ime, k, jms:jme)*0.2)/2
            !slope_V(:,k,:) = (v(ims:ime, k, jms+1:jme+1)*0.2 + &
            !                  v(ims:ime, k, jms:jme)*0.2)/2
        enddo
        

        div(ims:ime,kms:kme,jms:jme) = (diff_U+diff_V)/(dx) - (1/(smooth_height*jaco(ims:ime,kms:kme,jms:jme)))*(slope_U+slope_V)

        if (.NOT.(horz)) then
            do k = kms,kme
                if (k == kms) then
                    div(ims:ime, k, jms:jme) = div(ims:ime, k, jms:jme) + w(ims:ime, k, jms:jme)/dz(ims:ime, k, jms:jme)
                else
                    div(ims:ime, k, jms:jme) = div(ims:ime, k, jms:jme) + &
                                   (w(ims:ime,k,jms:jme)-w(ims:ime,k-1,jms:jme))/dz(ims:ime,k,jms:jme)
                endif
            enddo
        endif
        

    end subroutine calc_divergence

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

            write(*,*) 'Shape of U: ',SHAPE(domain%u%data_3d)
            write(*,*) 'Shape of V: ',SHAPE(domain%v%data_3d)
            write(*,*) 'Shape of W: ',SHAPE(domain%w%data_3d)


            ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
            call make_winds_grid_relative(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%sintheta, domain%costheta)

            ! flow blocking parameterization
            ! if (options%block_options%block_flow) then
            !     call add_blocked_flow(domain, options)
            ! endif

            ! linear winds
            if (options%physics%windtype==kWIND_LINEAR) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density)
            ! simple acceleration over topography
            elseif (options%physics%windtype==kCONSERVE_MASS) then
                call mass_conservative_acceleration(domain%u%data_3d, domain%v%data_3d, domain%zr_u, domain%zr_v)
            elseif (options%physics%windtype==kITERATIVE_WINDS) then
                call iterative_winds(domain, options)

            endif
            ! else assumes even flow over the mountains

            ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)
            call balance_uvw(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%dzdx, domain%dzdy, domain%advection_dz, domain%dx, domain%jacobian, domain%smooth_height, options)

        else

            write(*,*) 'Shape of U: ',SHAPE(domain%u%meta_data%dqdt_3d)
            write(*,*) 'Shape of V: ',SHAPE(domain%v%meta_data%dqdt_3d)
            write(*,*) 'Shape of W: ',SHAPE(domain%w%meta_data%dqdt_3d)


            ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
            call make_winds_grid_relative(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%w%meta_data%dqdt_3d, domain%sintheta, domain%costheta)

            ! linear winds
            if (options%physics%windtype==kWIND_LINEAR) then
                call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%parameters%advect_density, update=.True.)
            ! simple acceleration over topography
            elseif (options%physics%windtype==kCONSERVE_MASS) then
                call mass_conservative_acceleration(domain%u%meta_data%dqdt_3d, domain%v%meta_data%dqdt_3d, domain%zr_u, domain%zr_v)
            elseif (options%physics%windtype==kITERATIVE_WINDS) then
                call iterative_winds(domain, options, update_in=.True.)
            endif
            ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)

            write(*,*) 'max of U: ',maxval(domain%u%meta_data%dqdt_3d)
            write(*,*) 'max of V: ',maxval(domain%v%meta_data%dqdt_3d)
            write(*,*) 'max of W: ',maxval(domain%w%meta_data%dqdt_3d)

            call balance_uvw(domain% u %meta_data%dqdt_3d,      &
                             domain% v %meta_data%dqdt_3d,      &
                             domain% w %meta_data%dqdt_3d,      &
                             domain%dzdx, domain%dzdy,          &
                             domain%advection_dz, domain%dx,    &
                             domain%jacobian, domain%smooth_height, options)

        endif


    end subroutine update_winds

    subroutine iterative_winds(domain, options, update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: update_in
        
        ! interal parameters
        real, allocatable, dimension(:,:,:) :: div, ADJ,ADJ_coef, U_cor, V_cor, current_u, current_v, current_w
        real    :: corr_factor
        integer :: it, k, j, i, ims, ime, jms, jme, kms, kme, wind_k
        logical :: update
        
        update=.False.
        if (present(update_in)) update=update_in

        ims = lbound(domain%w%data_3d,1)
        ime = ubound(domain%w%data_3d,1)
        kms = lbound(domain%w%data_3d,2)
        kme = ubound(domain%w%data_3d,2)
        jms = lbound(domain%w%data_3d,3)
        jme = ubound(domain%w%data_3d,3)

        !If we are doing an update, we need to swap meta data into data_3d fields so it can be exchanged while balancing
        !First, we save a copy of the current data_3d so that we can substitute it back in later
        if (update) then
             current_u = domain%u%data_3d
             current_v = domain%v%data_3d
             current_w = domain%w%data_3d
             
             domain%u%data_3d = domain%u%meta_data%dqdt_3d
             domain%v%data_3d = domain%v%meta_data%dqdt_3d
             domain%w%data_3d = domain%w%meta_data%dqdt_3d
        endif
        
        !Do an initial exchange to make sure the U and V grids are similar for calculating w
        call domain%u%exchange_u()
        call domain%v%exchange_v()

        !First call bal_uvw to generate an initial-guess for vertical winds
        call balance_uvw(domain%u%data_3d,domain%v%data_3d,domain%w%data_3d,domain%dzdx,domain%dzdy,domain%advection_dz,&
                         domain%dx,domain%jacobian,domain%smooth_height,options)
  
        allocate(div(ims:ime,kms:kme,jms:jme))
        allocate(ADJ_coef(ims:ime,kms:kme,jms:jme))
        allocate(ADJ(ims:ime,kms:kme,jms:jme))
        allocate(U_cor(ims:ime,kms:kme,jms:jme))
        allocate(V_cor(ims:ime,kms:kme,jms:jme))

        ! First calculate and apply correction to w-winds
        ! !$omp parallel shared(w,dz,dzdx,dzdy,jaco,dx,smooth_height) firstprivate(nx,ny,nz) private(i,j,k,div,ADJ,ADJ_coef,corr_factor)
        ! !$omp do schedule(static)
        
        wind_k = kme
        do k = kms,kme
            if (sum(domain%advection_dz(ims,1:k,jms)) > domain%smooth_height) then
                wind_k = k
                exit
            endif
        enddo
        
        do k = kms,kme
            corr_factor = ((sum(domain%advection_dz(ims,1:k,jms)))/domain%smooth_height)
            !corr_factor = (k*1.0)/wind_k
            corr_factor = min(corr_factor,1.0)
            do i = ims,ime
                do j = jms,jme
                    domain%w%data_3d(i,k,j) = domain%w%data_3d(i,k,j) - corr_factor * (domain%w%data_3d(i,wind_k,j))
                enddo
            enddo
            ! Compute this now, since it wont change in the loop
            ADJ_coef(:,k,:) = -2/domain%dx 
        enddo
        ! !$omp end do
        ! !$omp barrier
        ! !$omp end parallel
        
        !Compute relative correction factors for U and V based on input speeds
        U_cor = ABS(domain%u%data_3d(ims:ime,:,jms:jme))/ &
                (ABS(domain%u%data_3d(ims:ime,:,jms:jme))+ABS(domain%v%data_3d(ims:ime,:,jms:jme)))            
        V_cor = 1 - U_cor !ABS(domain%v%data_3d(ims:ime,:,jms:jme))/ &
        !        (ABS(domain%u%data_3d(ims:ime,:,jms:jme))+ABS(domain%v%data_3d(ims:ime,:,jms:jme)))      
        
        U_cor = 0.5
        V_cor = 0.5
        ! Now, fixing w-winds, iterate over U/V to reduce divergence with new w-winds
        ! !$omp do schedule(static)

        do it = 0,options%parameters%wind_iterations
        
            !Compute divergence in new wind field
            call calc_divergence(div,domain%u%data_3d,domain%v%data_3d,domain%w%data_3d,domain%dzdx,domain%dzdy, &
                                domain%advection_dz,domain%dx,domain%jacobian,domain%smooth_height)

            !Compute adjustment based on divergence
            ADJ = div/ADJ_coef
        
            !Distribute divergence among the U and V fields                                                                    
            domain%u%data_3d(ims+2:ime,:,jms+1:jme-1) = domain%u%data_3d(ims+2:ime,:,jms+1:jme-1) + &
                                                        (ADJ(ims+1:ime-1,:,jms+1:jme-1) * U_cor(ims+2:ime,:,jms+1:jme-1))
                                                        
            domain%u%data_3d(ims+2:ime,:,jms+1:jme-1) = domain%u%data_3d(ims+2:ime,:,jms+1:jme-1) - &
                                                        (ADJ(ims+2:ime,:,jms+1:jme-1) * U_cor(ims+2:ime,:,jms+1:jme-1))
                                                        
            domain%v%data_3d(ims+1:ime-1,:,jms+2:jme) = domain%v%data_3d(ims+1:ime-1,:,jms+2:jme) + &
                                                        (ADJ(ims+1:ime-1,:,jms+1:jme-1) * V_cor(ims+1:ime-1,:,jms+2:jme))
                                                        
            domain%v%data_3d(ims+1:ime-1,:,jms+2:jme) = domain%v%data_3d(ims+1:ime-1,:,jms+2:jme) - &
                                                        (ADJ(ims+1:ime-1,:,jms+2:jme) * V_cor(ims+1:ime-1,:,jms+2:jme))
            call domain%u%exchange_u()
            call domain%v%exchange_v()
            
        enddo

        !If an update loop, swap meta_data and data_3d fields back
        if (update) then
            domain%u%meta_data%dqdt_3d = domain%u%data_3d
            domain%v%meta_data%dqdt_3d = domain%v%data_3d
            domain%w%meta_data%dqdt_3d = domain%w%data_3d
            
            domain%u%data_3d = current_u
            domain%v%data_3d = current_v
            domain%w%data_3d = current_w
        endif
        
        ! !$omp end do
        ! !$omp end parallel

    end subroutine iterative_winds

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
