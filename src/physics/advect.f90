!> ----------------------------------------------------------------------------
!!  A simple upwind advection scheme
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module adv_upwind
    use data_structures
    use options_interface, only: options_t
    use domain_interface,  only: domain_t

    implicit none
    private
    real,dimension(:,:,:),allocatable :: U_m, V_m, W_m, lastqv_m

    ! For use advecting a (convective?) wind field
    ! real,dimension(:,:,:),allocatable :: U_4cu_u, V_4cu_u, W_4cu_u
    ! real,dimension(:,:,:),allocatable :: U_4cu_v, V_4cu_v, W_4cu_v

    public :: upwind, upwind_init, upwind_var_request

contains

    subroutine upwind_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options

        integer :: nx, nz, ny

        nx = domain%grid%nx
        nz = domain%grid%nz
        ny = domain%grid%ny

        ! if module level arrays are already allocated for some reason, deallocate them first
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        if (allocated(lastqv_m)) deallocate(lastqv_m)

        ! allocate the module level arrays
        allocate(U_m     (nx-1,nz,ny  ))
        allocate(V_m     (nx,  nz,ny-1))
        allocate(W_m     (nx,  nz,ny  ))
        allocate(lastqv_m(nx,  nz,ny  ))

        !     if (.not.allocated(U_4cu_u)) then
        !         allocate(U_4cu_u(nx,  nz, ny))
        !         U_4cu_u = 0
        !         allocate(V_4cu_u(nx+1,nz, ny-1))
        !         V_4cu_u = 0
        !         allocate(W_4cu_u(nx+1,nz, ny))
        !         W_4cu_u = 0
        !
        !         allocate(U_4cu_v(nx-1,nz, ny+1))
        !         U_4cu_v = 0
        !         allocate(V_4cu_v(nx,  nz, ny))
        !         V_4cu_v = 0
        !         allocate(W_4cu_v(nx,  nz, ny+1))
        !         W_4cu_v = 0
        !     endif
    end subroutine

    subroutine upwind_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for upwind advection
        call options%alloc_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface])

        ! List the variables that are required for restarts with upwind advection
        call options%restart_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface])

    end subroutine

!     Note this routine has been manually inlined because the compiler didn't seem to optimize it well and made the array copies
!     the routine is left in for documentation.
!     subroutine flux2(l,r,U,nx,nz,ny,f)
!     !     Calculate the donor cell flux function
!     !     l = left gridcell scalar
!     !     r = right gridcell scalar
!     !     U = Courant number (u*dt/dx)
!     !
!     !     If U is positive, return l*U if U is negative return r*U
!     !     By using the mathematical form instead of the logical form,
!     !     we can run on the entire grid simultaneously, and avoid branches
!
!     !   arguments
!         implicit none
!         real, dimension(1:nx,1:nz,1:ny), intent(in) :: l,r,U
!         real, dimension(1:nx,1:nz,1:ny), intent(inout) :: f
!         integer,intent(in) :: ny,nz,nx
!         !   internal parameter
!         integer ::  err,i!,j,Ny,Nz,Nx
!         !   main code
!         f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2
!
!     end subroutine flux2

    subroutine advect3d(q,u,v,w,rho,dz,dx,nx,nz,ny,jaco,options)
        implicit none
        real, dimension(1:nx,  1:nz,1:ny),  intent(inout)   :: q
        real, dimension(1:nx,  1:nz,1:ny),  intent(in)      :: w
        real, dimension(1:nx-1,1:nz,1:ny),  intent(in)      :: u
        real, dimension(1:nx,  1:nz,1:ny-1),intent(in)      :: v
        real, dimension(:, :, :),           intent(in)      :: rho
        real, dimension(:, :, :),           intent(in)      :: dz
        real,                               intent(in)      :: dx
        integer,                            intent(in)      :: ny, nz, nx
        real, dimension(:, :, :),           intent(in)      :: jaco
        type(options_t), intent(in)::options

        ! interal parameters
        integer                         :: err, i, j, k
        real, dimension(1:nx,1:nz,1:ny) :: qin
        real, dimension(1:nx-1,1:nz)    :: f1   ! historical note, there used to be an f2 to store f[x+1]
        real, dimension(1:nx-2,1:nz)    :: f3,f4
        real, dimension(1:nx-2,1:nz-1)  :: f5

        ! Multiply dz by jacobian so that all advection terms will be divided by mass-centered jacobian

        !$omp parallel shared(qin,q,u,v,w,rho,dz,jaco) firstprivate(nx,ny,nz) private(i,f1,f3,f4,f5)
        !$omp do schedule(static)
        do i=1,ny
            qin(:,:,i)=q(:,:,i)
        enddo
        !$omp end do
        !$omp barrier
        !$omp do schedule(static)
        do i=2,ny-1
            ! by manually inlining the flux2 call we should remove extra array copies that the compiler doesn't remove.
            ! equivalent flux2 calls are left in for reference (commented) to restore recall that f1,f3,f4... arrays should be 3D : n x m x 1
            ! calculate fluxes between grid cells
            ! call flux2(qin(1:nx-1,:,i),     qin(2:nx,:,i),     u(1:nx-1,:,i),     nx-1,nz,  1,f1)  ! f1 = Ux0 and Ux1
            ! call flux2(qin(2:nx-1,:,i),     qin(2:nx-1,:,i+1), v(2:nx-1,:,i),     nx-2,nz,  1,f3)  ! f3 = Vy1
            ! call flux2(qin(2:nx-1,:,i-1),   qin(2:nx-1,:,i),   v(2:nx-1,:,i-1),   nx-2,nz,  1,f4)  ! f4 = Vy0
            ! call flux2(qin(2:nx-1,1:nz-1,i),qin(2:nx-1,2:nz,i),w(2:nx-1,1:nz-1,i),nx-2,nz-1,1,f5)  ! f5 = Wz0 and Wz1
            f1= ((u(1:nx-1,:,i)      + ABS(u(1:nx-1,:,i)))      * qin(1:nx-1,:,i)    + &
                 (u(1:nx-1,:,i)      - ABS(u(1:nx-1,:,i)))      * qin(2:nx,:,i))     / 2

            f3= ((v(2:nx-1,:,i)      + ABS(v(2:nx-1,:,i)))      * qin(2:nx-1,:,i)    + &
                 (v(2:nx-1,:,i)      - ABS(v(2:nx-1,:,i)))      * qin(2:nx-1,:,i+1)) / 2

            f4= ((v(2:nx-1,:,i-1)    + ABS(v(2:nx-1,:,i-1)))    * qin(2:nx-1,:,i-1)  + &
                 (v(2:nx-1,:,i-1)    - ABS(v(2:nx-1,:,i-1)))    * qin(2:nx-1,:,i))   / 2

            f5= ((w(2:nx-1,1:nz-1,i) + ABS(w(2:nx-1,1:nz-1,i))) * qin(2:nx-1,1:nz-1,i) + &
                 (w(2:nx-1,1:nz-1,i) - ABS(w(2:nx-1,1:nz-1,i))) * qin(2:nx-1,2:nz,i))  / 2

           ! if (options%parameters%advect_density) then
           !     ! perform horizontal advection
           !     q(2:nx-1,:,i)      = q(2:nx-1,:,i)      - ((f1(2:nx-1,:) - f1(1:nx-2,:)) + (f3 - f4)) &
           !                          / rho(2:nx-1,:,i) / dz(2:nx-1,:,i)
           !     ! then vertical
           !     ! (order doesn't matter because fluxes f1-6 are calculated before applying them)
           !     ! add fluxes to middle layers
           !     q(2:nx-1,2:nz-1,i) = q(2:nx-1,2:nz-1,i) - (f5(:,2:nz-1) - f5(:,1:nz-2))                       &
           !                          / rho(2:nx-1,2:nz-1,i) / dz(2:nx-1,2:nz-1,i)
           !     ! add fluxes to bottom layer
           !     q(2:nx-1,1,i)      = q(2:nx-1,1,i)      - f5(:,1)                                             &
           !                          / rho(2:nx-1,1,i) / dz(2:nx-1,1,i)
           !     ! add fluxes to top layer
           !     q(2:nx-1,nz,i)     = q(2:nx-1,nz,i)     - (qin(2:nx-1,nz,i) * w(2:nx-1,nz,i)-f5(:,nz-1))      &
           !                          / rho(2:nx-1,nz,i) / dz(2:nx-1,nz,i)
           ! else
               ! perform horizontal advection, from difference terms
               q(2:nx-1,:,i)      = q(2:nx-1,:,i)       - ((f1(2:nx-1,:) - f1(1:nx-2,:)) + (f3 - f4)) /(dx*dz(2:nx-1,:,i)*jaco(2:nx-1,:,i))
               ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
               q(2:nx-1,2:nz-1,i) = q(2:nx-1,2:nz-1,i)  - (f5(:,2:nz-1) - f5(:,1:nz-2)) / (dz(2:nx-1,2:nz-1,i)*jaco(2:nx-1,2:nz-1,i))
               ! add fluxes to bottom layer
               q(2:nx-1,1,i)      = q(2:nx-1,1,i)       - f5(:,1) / (dz(2:nx-1,1,i)*jaco(2:nx-1,1,i))
               ! add fluxes to top layer
               q(2:nx-1,nz,i)     = q(2:nx-1,nz,i)      - (qin(2:nx-1,nz,i) * w(2:nx-1,nz,i) - f5(:,nz-1)) / (dz(2:nx-1,nz,i)*jaco(2:nx-1,nz,i))
           ! endif
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine advect3d

    ! subroutine setup_cu_winds(domain, options, dt)
    !     implicit none
    !     type(domain_type),  intent(in) :: domain
    !     type(options_type), intent(in) :: options
    !     real,               intent(in) :: dt
    !
    !     real    :: dx
    !     integer :: nx,nz,ny
    !
    !     dx = domain%dx
    !     nx = size(domain%dz,1)
    !     nz = size(domain%dz,2)
    !     ny = size(domain%dz,3)
    !
    !
    !     U_4cu_u           =  (dt/dx) * (domain%u(1:nx,:,:)      + domain%u(2:nx+1,:,:)) / 2
    !     V_4cu_u(2:nx,:,:) =  (dt/dx) * (domain%v(1:nx-1,:,2:ny) + domain%v(2:nx,:,2:ny)) / 2
    !     W_4cu_u(2:nx,:,:) =  (dt/dx) * (domain%w(1:nx-1,:,:)    + domain%w(2:nx,:,:)) / 2
    !     call rebalance_cu_winds(U_4cu_u, V_4cu_u, W_4cu_u)
    !
    !     U_4cu_v(:,:,2:ny) =  (dt/dx) * (domain%u(2:nx,:,1:ny-1) + domain%u(2:nx,:,2:ny)) / 2
    !     V_4cu_v           =  (dt/dx) * (domain%v(:,:,1:ny)      + domain%v(:,:,2:ny+1)) / 2
    !     W_4cu_v(:,:,2:ny) =  (dt/dx) * (domain%w(:,:,1:ny-1)    + domain%w(:,:,2:ny)) / 2
    !     call rebalance_cu_winds(U_4cu_v, V_4cu_v, W_4cu_v)
    !
    ! end subroutine setup_cu_winds

    ! subroutine rebalance_cu_winds(u,v,w)
    !     implicit none
    !     ! u, v, w 3D east-west, south-north, and up-down winds repsectively
    !     ! note for this code, u is [nx-1,nz,ny] and v is [nx,nz,ny-1]
    !     real, dimension(:,:,:), intent(inout) :: u, v, w
    !
    !     real, allocatable, dimension(:,:) :: divergence, du, dv
    !     integer :: i,nx,ny,nz
    !
    !     nx = size(w,1)
    !     nz = size(w,2)
    !     ny = size(w,3)
    !
    !     allocate(divergence(nx-2,ny-2))
    !     allocate(du(nx-2,ny-2))
    !     allocate(dv(nx-2,ny-2))
    !
    !     do i=1,nz
    !         ! calculate horizontal divergence
    !         dv = v(2:nx-1,i,2:ny-1) - v(2:nx-1,i,1:ny-2)
    !         du = u(2:nx-1,i,2:ny-1) - u(1:nx-2,i,2:ny-1)
    !         divergence = du + dv
    !         ! Then calculate w to balance
    !         if (i==1) then
    !             ! if this is the first model level start from 0 at the ground
    !             w(2:nx-1,i,2:ny-1) = 0 - divergence
    !         else
    !             ! else calculate w as a change from w at the level below
    !             w(2:nx-1,i,2:ny-1) = w(2:nx-1,i-1,2:ny-1) - divergence
    !         endif
    !     enddo
    !
    ! end subroutine rebalance_cu_winds
    !
    ! subroutine advect_cu_winds(domain, options, dt)
    !     implicit none
    !     type(domain_type),  intent(inout) :: domain
    !     type(options_type), intent(in)    :: options
    !     real,               intent(in)    :: dt
    !
    !     integer :: nx,nz,ny
    !
    !     nx = size(domain%dz,1)
    !     nz = size(domain%dz,2)
    !     ny = size(domain%dz,3)
    !
    !     ! first put the background u,v,w winds on a staggered grid with respect to the u grid
    !     ! then advect the u winds
    !     if (options%advect_density) then
    !         print*, "ERROR: Density advection not enabled when using convective winds"
    !         print*, "   Requires update to wind.f90 balance_uvw and advect.f90 (at least)"
    !         stop
    !     endif
    !
    !     call setup_cu_winds(domain, options, dt)
    !
    !     ! set the top boundary condition for CU winds to 0 to prevent artifacts coming in from the "top"
    !     domain%u_cu(:,nz,:) = 0
    !     domain%v_cu(:,nz,:) = 0
    !
    !     call advect3d(domain%u_cu, U_4cu_u,V_4cu_u,W_4cu_u, domain%rho, domain%dz_inter, nx+1,nz,ny, options)
    !     call advect3d(domain%v_cu, U_4cu_v,V_4cu_v,W_4cu_v, domain%rho, domain%dz_inter, nx,nz,ny+1, options)
    !
    ! end subroutine advect_cu_winds


    subroutine test_divergence(u, v, w, dz, dx, jaco_u, jaco_v, jaco_w, ims, ime, jms, jme, kms, kme)
        implicit none
        real, intent(in) :: u(ims:ime+1,kms:kme,jms:jme), jaco_u(ims:ime+1,kms:kme,jms:jme)
        real, intent(in) :: v(ims:ime,kms:kme,jms:jme+1), jaco_v(ims:ime,kms:kme,jms:jme+1)
        real, intent(in) :: w(ims:ime,kms:kme,jms:jme), jaco_w(ims:ime,kms:kme,jms:jme)
        real, intent(in) :: dz(ims:ime,kms:kme,jms:jme)
        real, intent(in) :: dx
        integer, intent(in) :: ims, ime, jms, jme, kms, kme

        real, allocatable :: du(:,:), dv(:,:), dzu(:,:,:), dzv(:,:,:)
        integer :: i,j,k

        allocate(du(ims+1:ime-1,jms+1:jme-1))
        allocate(dv(ims+1:ime-1,jms+1:jme-1))
        allocate(dzv(ims:ime,kms:kme,jms:jme))
        allocate(dzu(ims:ime,kms:kme,jms:jme))

        dzv = 0
        dzu = 0

        dzv(:,:,jms+1:jme) = (dz(:,:,jms+1:jme) + dz(:,:,jms:jme-1)) / 2
        dzu(ims+1:ime,:,:) = (dz(ims+1:ime,:,:) + dz(ims:ime-1,:,:)) / 2

        do i=ims+1,ime-1
            do j=jms+1,jme-1
                do k=kms+1,kme
                    dv(i,j) = (v(i,k,j+1) * jaco_v(i,k,j+1) - v(i,k,j) * jaco_v(i,k,j))/dx
                    du(i,j) = (u(i+1,k,j) * jaco_u(i+1,k,j) - u(i,k,j) * jaco_u(i,k,j))/dx

                    if (abs(du(i,j) + dv(i,j) + (w(i,k,j)*jaco_w(i,k,j)-w(i,k-1,j)*jaco_w(i,k-1,j))/(dz(i,k,j))) > 1e-3) then
                        print*, this_image(), i,j,k , abs(du(i,j) + dv(i,j) + (w(i,k,j)*jaco_w(i,k,j)-w(i,k-1,j)*jaco_w(i,k-1,j))/(dz(i,k,j)))
                        print*, "Winds are not balanced on entry to advect"
                        !error stop
                    endif
                enddo
            enddo
        enddo

    end subroutine test_divergence

    subroutine setup_module_winds(u,v,w, dz, dx, options, dt,jaco,jaco_u,jaco_v, jaco_w)
        implicit none
        real,               intent(in)  :: u(:,:,:)
        real,               intent(in)  :: v(:,:,:)
        real,               intent(in)  :: w(:,:,:)
        real,               intent(in)  :: dz(:,:,:)
        real,               intent(in)  :: dx
        type(options_t),    intent(in)  :: options
        real,               intent(in)  :: dt
        real,               intent(in)  :: jaco(:, :, :)
        real,               intent(in)  :: jaco_u(:, :, :)
        real,               intent(in)  :: jaco_v(:, :, :)
        real,               intent(in)  :: jaco_w(:, :, :)

        integer :: nx, nz, ny, i

        nx = size(w,1)
        nz = size(w,2)
        ny = size(w,3)

        ! if this if the first time we are called, we need to allocate the module level arrays
        ! Could/should be put in an init procedure
        if (.not.allocated(U_m)) then
            allocate(U_m     (nx-1,nz,ny  ))
            allocate(V_m     (nx,  nz,ny-1))
            allocate(W_m     (nx,  nz,ny  ))
            allocate(lastqv_m(nx,  nz,ny  ))
        endif

        ! calculate U,V,W normalized for dt/dx (dx**2 for density advection so we can skip a /dx in the actual advection code)
        ! if (options%parameters%advect_density) then
        !     U_m = domain%ur(2:nx,:,:) * (dt/dx**2)
        !     V_m = domain%vr(:,:,2:ny) * (dt/dx**2)
        !     W_m = domain%wr           * (dt/dx**2)
        ! else
            ! if (options%physics%convection > 0) then
                ! print*, "Advection of convective winds not enabled in ICAR >=1.5 yet"
                ! stop
                ! U_m = (domain%u_cu(2:nx,:,:) + domain%u(2:nx,:,:)) * (dt/dx)
                ! V_m = (domain%v_cu(:,:,2:ny) + domain%v(:,:,2:ny)) * (dt/dx)
                ! W_m = (domain%w_cu + domain%w)                     * (dt/dx)
                ! call rebalance_cu_winds(U_m,V_m,W_m)
            ! else
                U_m = u(2:nx,:,:) * ((dz(1:nx-1,:,:)+dz(2:nx,:,:))/2 * dt) * jaco_u(2:nx,:,:)
                V_m = v(:,:,2:ny) * ((dz(:,:,1:ny-1)+dz(:,:,2:ny))/2 * dt) * jaco_v(:,:,2:ny)
                W_m = w * dt * jaco_w
            ! endif
        ! endif

    end subroutine setup_module_winds

    subroutine setup_advection_dz(domain, options, nx, ny, nz)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        integer,         intent(in)    :: nx, ny, nz
        integer :: i, ims, ime, jms, jme, kms, kme

        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme

        if (.not.allocated(domain%advection_dz)) then
            allocate(domain%advection_dz(ims:ime,kms:kme,jms:jme))
        else
            return
        endif

        if (options%parameters%fixed_dz_advection) then
            do i=kms,kme
                domain%advection_dz(:,i,:) = options%parameters%dz_levels(i)
            enddo
        else
            domain%advection_dz = domain%dz_interface%data_3d
        endif

    end subroutine setup_advection_dz


    ! primary entry point, advect all scalars in domain
    subroutine upwind(domain,options,dt)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        real,            intent(in)    :: dt

        real    :: dx
        integer :: nx, nz, ny, i

        nx = domain%grid%nx
        nz = domain%grid%nz
        ny = domain%grid%ny

        call setup_advection_dz(domain, options, nx,ny,nz)

        ! calculate U,V,W normalized for dt/dx (dx**2 for density advection so we can skip a /dx in the actual advection code)
        call setup_module_winds(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%advection_dz, domain%dx, options, dt,domain%jacobian,domain%jacobian_u,domain%jacobian_v,domain%jacobian_w)

        ! lastqv_m=domain%qv

        if (options%parameters%debug) then
            call test_divergence(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%advection_dz, domain%dx, domain%jacobian_u, domain%jacobian_v, domain%jacobian_w, domain%ims, domain%ime, domain%jms, domain%jme, domain%kms, domain%kme)
        endif

        if (options%vars_to_advect(kVARS%water_vapor)>0)                  call advect3d(domain%water_vapor%data_3d,             U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%cloud_water)>0)                  call advect3d(domain%cloud_water_mass%data_3d,        U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%rain_in_air)>0)                  call advect3d(domain%rain_mass%data_3d,               U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%snow_in_air)>0)                  call advect3d(domain%snow_mass%data_3d,               U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%potential_temperature)>0)        call advect3d(domain%potential_temperature%data_3d,   U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%cloud_ice)>0)                    call advect3d(domain%cloud_ice_mass%data_3d,          U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%graupel_in_air)>0)               call advect3d(domain%graupel_mass%data_3d,            U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%ice_number_concentration)>0)     call advect3d(domain%cloud_ice_number%data_3d,        U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%rain_number_concentration)>0)    call advect3d(domain%rain_number%data_3d,             U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%snow_number_concentration)>0)    call advect3d(domain%snow_number%data_3d,             U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%graupel_number_concentration)>0) call advect3d(domain%graupel_number%data_3d,          U_m,V_m,W_m, domain%density%data_3d, domain%advection_dz, domain%dx, nx,nz,ny, domain%jacobian, options)

        ! if (options%physics%convection > 0) then
        !     call advect_cu_winds(domain, options, dt)
        ! endif

        ! used in some physics routines
        ! domain%tend%qv_adv = (domain%qv - lastqv_m) / dt
    end subroutine upwind

end module adv_upwind
