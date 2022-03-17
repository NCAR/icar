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

        integer :: ims, ime, jms, jme, kms, kme

        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme

        ! if module level arrays are already allocated for some reason, deallocate them first
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        if (allocated(lastqv_m)) deallocate(lastqv_m)

        ! allocate the module level arrays
        allocate(U_m     (ims+1:ime,kms:kme,jms:jme  ))
        allocate(V_m     (ims:ime,  kms:kme,jms+1:jme))
        allocate(W_m     (ims:ime,  kms:kme,jms:jme  ))
        allocate(lastqv_m(ims:ime,  kms:kme,jms:jme  ))

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

    subroutine advect3d(q,rho_in,dz,ims,ime,kms,kme,jms,jme,jaco,options)

        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(inout)   :: q
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: rho_in
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: dz
        integer,                            intent(in)      :: ims, ime, kms, kme, jms, jme
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: jaco
        type(options_t), intent(in)::options

        ! interal parameters
        integer                         :: err, i, j, k
        real, dimension(ims:ime-1,kms:kme)    :: f1   ! historical note, there used to be an f2 to store f[x+1]
        real, dimension(ims:ime-2,kms:kme)    :: f3,f4
        real, dimension(ims:ime-2,kms:kme-1)  :: f5
        real, dimension(ims:ime,kms:kme,jms:jme) :: qin, rho
        
        rho = 1
        
        do i=jms,jme
            qin(:,:,i)=q(:,:,i)
        enddo
        
        if (options%parameters%advect_density) rho = rho_in
        
        
        ! !$omp parallel shared(qin,q,u,v,w,rho,dz,jaco) firstprivate(ims,ime,kms,kme,jms,jme) private(i,f1,f3,f4,f5)
        ! !$omp do schedule(static)


        ! !$omp end do
        ! !$omp barrier
        ! !$omp do schedule(static)
        do i=jms+1,jme-1
            ! by manually inlining the flux2 call we should remove extra array copies that the compiler doesn't remove.
            ! equivalent flux2 calls are left in for reference (commented) to restore recall that f1,f3,f4... arrays should be 3D : n x m x 1
            ! calculate fluxes between grid cells
            ! call flux2(qin(1:nx-1,:,i),     qin(2:nx,:,i),     u(1:nx-1,:,i),     nx-1,nz,  1,f1)  ! f1 = Ux0 and Ux1
            ! call flux2(qin(2:nx-1,:,i),     qin(2:nx-1,:,i+1), v(2:nx-1,:,i),     nx-2,nz,  1,f3)  ! f3 = Vy1
            ! call flux2(qin(2:nx-1,:,i-1),   qin(2:nx-1,:,i),   v(2:nx-1,:,i-1),   nx-2,nz,  1,f4)  ! f4 = Vy0
            ! call flux2(qin(2:nx-1,1:nz-1,i),qin(2:nx-1,2:nz,i),w(2:nx-1,1:nz-1,i),nx-2,nz-1,1,f5)  ! f5 = Wz0 and Wz1
            f1= ((U_m(ims+1:ime,:,i)      + ABS(U_m(ims+1:ime,:,i)))      * qin(ims:ime-1,:,i)    + &
                 (U_m(ims+1:ime,:,i)      - ABS(U_m(ims+1:ime,:,i)))      * qin(ims+1:ime,:,i))     / 2

            f3= ((V_m(ims+1:ime-1,:,i+1)      + ABS(V_m(ims+1:ime-1,:,i+1)))      * qin(ims+1:ime-1,:,i)  + &
                 (V_m(ims+1:ime-1,:,i+1)      - ABS(V_m(ims+1:ime-1,:,i+1)))      * qin(ims+1:ime-1,:,i+1)) / 2

            f4= ((V_m(ims+1:ime-1,:,i)    + ABS(V_m(ims+1:ime-1,:,i)))    * qin(ims+1:ime-1,:,i-1) + &
                (V_m(ims+1:ime-1,:,i)    - ABS(V_m(ims+1:ime-1,:,i)))    * qin(ims+1:ime-1,:,i))  / 2

            f5= ((W_m(ims+1:ime-1,kms:kme-1,i) + ABS(W_m(ims+1:ime-1,kms:kme-1,i))) * qin(ims+1:ime-1,kms:kme-1,i) + &
                 (W_m(ims+1:ime-1,kms:kme-1,i) - ABS(W_m(ims+1:ime-1,kms:kme-1,i))) * qin(ims+1:ime-1,kms+1:kme,i))  / 2

                
               ! perform horizontal advection, from difference terms
               q(ims+1:ime-1,:,i)      = q(ims+1:ime-1,:,i)       - ((f1(ims+1:ime-1,:) - f1(ims:ime-2,:)) + (f3 - f4)) &
                                   / (jaco(ims+1:ime-1,:,i)*rho(ims+1:ime-1,:,i))                      
               ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
               q(ims+1:ime-1,kms+1:kme-1,i) = q(ims+1:ime-1,kms+1:kme-1,i)  - (f5(:,kms+1:kme-1) - f5(:,kms:kme-2)) &
                                   / (dz(ims+1:ime-1,kms+1:kme-1,i)*jaco(ims+1:ime-1,kms+1:kme-1,i)*rho(ims+1:ime-1,kms+1:kme-1,i))
               ! add fluxes to bottom layer
               q(ims+1:ime-1,kms,i)      = q(ims+1:ime-1,kms,i)       - f5(:,kms) &
                                   / (dz(ims+1:ime-1,kms,i)*jaco(ims+1:ime-1,kms,i) * rho(ims+1:ime-1,kms,i) )
               ! add fluxes to top layer
               q(ims+1:ime-1,kme,i)     = q(ims+1:ime-1,kme,i)     - (qin(ims+1:ime-1,kme,i) * W_m(ims+1:ime-1,kme,i) - f5(:,kme-1)) &
                                   / (dz(ims+1:ime-1,kme,i)*jaco(ims+1:ime-1,kme,i)*rho(ims+1:ime-1,kme,i))


        enddo
        ! !$omp end do
        ! !$omp end parallel
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


    subroutine test_divergence(dz, ims, ime, kms, kme, jms, jme)
        implicit none
        real, intent(in) :: dz(ims:ime,kms:kme,jms:jme)
        integer, intent(in) :: ims, ime, jms, jme, kms, kme

        real, allocatable :: du(:,:), dv(:,:), dw(:,:)
        integer :: i,j,k

        allocate(du(ims+1:ime-1,jms+1:jme-1))
        allocate(dv(ims+1:ime-1,jms+1:jme-1))
        allocate(dw(ims+1:ime-1,jms+1:jme-1))

        do i=ims+1,ime-1
            do j=jms+1,jme-1
                do k=kms,kme
                    du(i,j) = (U_m(i+1,k,j)-U_m(i,k,j))
                    dv(i,j) = (V_m(i,k,j+1)-V_m(i,k,j))
                    if (k==kms) then
                        dw(i,j) = (W_m(i,k,j))/dz(i,k,j)
                    else
                        dw(i,j) = (W_m(i,k,j)-W_m(i,k-1,j))/dz(i,k,j)
                    endif
                    if (abs(du(i,j) + dv(i,j) + dw(i,j)) > 1e-3) then
                        print*, this_image(), i,k,j , abs(du(i,j) + dv(i,j) + dw(i,j))
                        print*, "Winds are not balanced on entry to advect"
                        !error stop
                    endif
                enddo
            enddo
        enddo

    end subroutine test_divergence

    subroutine setup_module_winds(u,v,w, dx, options, dt,jaco_u,jaco_v, jaco_w,rho_in,ims, ime, kms, kme, jms, jme)
        implicit none
        real,               intent(in)  :: u(ims:ime+1,kms:kme,jms:jme), v(ims:ime,kms:kme,jms:jme+1)
        real,               intent(in)  :: w(ims:ime,kms:kme,jms:jme)
        real,               intent(in)  :: dx
        type(options_t),    intent(in)  :: options
        real,               intent(in)  :: dt
        real,               intent(in)  :: jaco_u(ims:ime+1,kms:kme,jms:jme)
        real,               intent(in)  :: jaco_v(ims:ime,kms:kme,jms:jme+1), jaco_w(ims:ime,kms:kme,jms:jme)
        real,               intent(in)  :: rho_in(ims:ime,kms:kme,jms:jme)
        integer, intent(in) :: ims, ime, jms, jme, kms, kme

        real, dimension(ims:ime,kms:kme,jms:jme) :: rho
        integer :: i
        

        ! if this if the first time we are called, we need to allocate the module level arrays
        ! Could/should be put in an init procedure
        if (.not.allocated(U_m)) then
            allocate(U_m     (ims+1:ime,kms:kme,jms:jme  ))
            allocate(V_m     (ims:ime,  kms:kme,jms+1:jme))
            allocate(W_m     (ims:ime,  kms:kme,jms:jme  ))
            allocate(lastqv_m(ims:ime,  kms:kme,jms:jme  ))
        endif
        
        rho = 1
        if (options%parameters%advect_density) rho = rho_in

        ! if (options%physics%convection > 0) then
            ! print*, "Advection of convective winds not enabled in ICAR >=1.5 yet"
            ! stop
            ! U_m = (domain%u_cu(2:nx,:,:) + domain%u(2:nx,:,:)) * (dt/dx)
            ! V_m = (domain%v_cu(:,:,2:ny) + domain%v(:,:,2:ny)) * (dt/dx)
            ! W_m = (domain%w_cu + domain%w)                     * (dt/dx)
            ! call rebalance_cu_winds(U_m,V_m,W_m)
        ! else
             ! Divide only U and V by dx. This minimizes the number of operations per advection step. W cannot be divided by dz,
             ! since non-uniform dz spacing does not allow for the same spacing to be assumed on either side of a k+1/2 interface,
             ! as is required for the upwind scheme.
             U_m = u(ims+1:ime,:,:) * dt * jaco_u(ims+1:ime,:,:) * (rho(ims+1:ime,:,:)+rho(ims:ime-1,:,:))*0.5 / dx
             V_m = v(:,:,jms+1:jme) * dt * jaco_v(:,:,jms+1:jme) * (rho(:,:,jms+1:jme)+rho(:,:,jms:jme-1))*0.5 / dx
             W_m(:,kms:kme-1,:) = w(:,kms:kme-1,:) * dt * jaco_w(:,kms:kme-1,:) * (rho(:,kms+1:kme,:)+rho(:,kms:kme-1,:))*0.5
             W_m(:,kme,:) = w(:,kme,:) * dt * jaco_w(:,kme,:) * rho(:,kme,:)
         ! endif

    end subroutine setup_module_winds

    subroutine setup_advection_dz(domain, options)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
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

        do i=kms,kme
            domain%advection_dz(:,i,:) = options%parameters%dz_levels(i)
        enddo

    end subroutine setup_advection_dz


    ! primary entry point, advect all scalars in domain
    subroutine upwind(domain,options,dt)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        real,            intent(in)    :: dt

        real    :: dx
        integer :: i

        call setup_advection_dz(domain, options)

        ! calculate U,V,W normalized for dt/dx (dx**2 for density advection so we can skip a /dx in the actual advection code)
        call setup_module_winds(domain%u%data_3d, domain%v%data_3d, domain%w%data_3d, domain%dx, options, dt,domain%jacobian_u,domain%jacobian_v,domain%jacobian_w,domain%density%data_3d, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme)

        ! lastqv_m=domain%qv

        if (options%parameters%debug) then
            call test_divergence(domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme)
        endif

        if (options%vars_to_advect(kVARS%water_vapor)>0)                  call advect3d(domain%water_vapor%data_3d,             domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%cloud_water)>0)                  call advect3d(domain%cloud_water_mass%data_3d,        domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%rain_in_air)>0)                  call advect3d(domain%rain_mass%data_3d,               domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%snow_in_air)>0)                  call advect3d(domain%snow_mass%data_3d,               domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%potential_temperature)>0)        call advect3d(domain%potential_temperature%data_3d,   domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%cloud_ice)>0)                    call advect3d(domain%cloud_ice_mass%data_3d,          domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%graupel_in_air)>0)               call advect3d(domain%graupel_mass%data_3d,            domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%ice_number_concentration)>0)     call advect3d(domain%cloud_ice_number%data_3d,        domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%rain_number_concentration)>0)    call advect3d(domain%rain_number%data_3d,             domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%snow_number_concentration)>0)    call advect3d(domain%snow_number%data_3d,             domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)
        if (options%vars_to_advect(kVARS%graupel_number_concentration)>0) call advect3d(domain%graupel_number%data_3d,          domain%density%data_3d, domain%advection_dz, domain%ims, domain%ime, domain%kms, domain%kme, domain%jms, domain%jme, domain%jacobian, options)

        ! if (options%physics%convection > 0) then
        !     call advect_cu_winds(domain, options, dt)
        ! endif

        ! used in some physics routines
        ! domain%tend%qv_adv = (domain%qv - lastqv_m) / dt
    end subroutine upwind

end module adv_upwind
