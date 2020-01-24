!> ----------------------------------------------------------------------------
!!  A simple upwind advection scheme
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module adv_upwind
    use data_structures

    implicit none
    private
    real,dimension(:,:,:),allocatable :: U_m, V_m, W_m, lastqv_m
    real,dimension(:,:,:),allocatable :: U_4cu_u, V_4cu_u, W_4cu_u
    real,dimension(:,:,:),allocatable :: U_4cu_v, V_4cu_v, W_4cu_v
    real,dimension(:,:,:),allocatable :: Uc_m, Vc_m
    public::upwind

contains
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
!
!    modifications:
!
!    20190402 jh ... added integer that indicates which boundary condition to use for the quantity.
!                    0 ... zero gradient (default)
!                    1 ... constant gradient
!                    2 ... zero value
!                    3 ... zero value but for the flux gradient over the topmost vertical layer
!
    subroutine advect3d(q,u,v,w,rho,dz,nx,nz,ny,options,bctop)
        implicit none
        real, dimension(1:nx,  1:nz,1:ny),  intent(inout)   :: q
        real, dimension(1:nx,  1:nz,1:ny),  intent(in)      :: w
        real, dimension(1:nx-1,1:nz,1:ny),  intent(in)      :: u
        real, dimension(1:nx,  1:nz,1:ny-1),intent(in)      :: v
        real, dimension(:, :, :),           intent(in)      :: rho
        real, dimension(:, :, :),           intent(in)      :: dz
        integer,                            intent(in)      :: ny, nz, nx
        integer,                            intent(in)      :: bctop
        type(options_type), intent(in)::options


        ! interal parameters
        integer                         :: err, i
        real, dimension(1:nx,1:nz,1:ny) :: qin
        real, dimension(1:nx,1:nz,1:ny) :: qinbc ! equal to qin if bc_top is not set or set to zero
        real, dimension(1:nx-1,1:nz)    :: f1    ! historical note, there used to be an f2 to store f[x+1]
        real, dimension(1:nx-2,1:nz)    :: f3,f4
        real, dimension(1:nx-2,1:nz-1)  :: f5

        !$omp parallel shared(qin,qinbc,q,u,v,w,rho,dz) firstprivate(nx,ny,nz) private(i,f1,f3,f4,f5)
        !$omp do schedule(static)
        do i=1,ny
            qin(:,:,i)=q(:,:,i)
            qinbc(:,:,i)=q(:,:,i)
        enddo
        !$omp end do

        !$omp barrier

        ! ---------------
        ! jhorak:
        ! if alternative boundary conditions for the model top are set, this if clause triggers.
        ! the standard bc (bctop = 0) currently is a zero gradient boundary condition for all quantities
        ! in this case the if clause is ignored.
        ! ---------------
        if (bctop /= 0) then
            if (bctop == 1) then
                ! constant gradient BC. If q_{n+1}<= 0 it is set to zero.
                ! we basically extrapolate q_{n+1} linearly from q_{n} and q_{n-1}.
                ! here n is the number of vertical levels
                qinbc(:,nz,:) = qin(:,nz,:)+( qin(:,nz,:)-qin(:,nz-1,:) )
                ! if one of the quantities is now negative we correct that to zero.
                where (qinbc(:,nz,:) < 0) qinbc(:,nz,:) = 0
                
            elseif (bctop == 2) then
                ! zero value BC. We set q_{n+1} = 0 where n is the number of vertical levels
                ! only inflow (downdraft) needs adjusting since for updrafts no quantity outside
                ! the domain is necessary for the calculation
                where (w(:,nz,:) < 0) qinbc(:,nz,:) = 0
                
            elseif (bctop == 3) then
                ! zero value BC for flux gradient of the topmost level for downdrafts
                ! this ensures that a downdraft at the model top wont change the content
                ! of the topmost grid cell. Again only needed for downdrafts since
                ! updraft advection doesn't required an assumption about the quantity
                ! outside of the model domain.
                where (w(:,nz,:) < 0) qinbc(:,nz,:) = (w(:,nz-1,:)/w(:,nz,:)) * qin(:,nz,:)
            endif
        endif

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
                 (w(2:nx-1,1:nz-1,i) - ABS(w(2:nx-1,1:nz-1,i))) * qin(2:nx-1,2:nz,i)) / 2


           if (options%advect_density) then
               ! perform horizontal advection
               q(2:nx-1,:,i)      = q(2:nx-1,:,i)      - ((f1(2:nx-1,:) - f1(1:nx-2,:)) + (f3 - f4)) &
                                    / rho(2:nx-1,:,i) / dz(2:nx-1,:,i)
               ! then vertical
               ! (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
               q(2:nx-1,2:nz-1,i) = q(2:nx-1,2:nz-1,i) - (f5(:,2:nz-1) - f5(:,1:nz-2))                       &
                                    / rho(2:nx-1,2:nz-1,i) / dz(2:nx-1,2:nz-1,i)
               ! add fluxes to bottom layer
               q(2:nx-1,1,i)      = q(2:nx-1,1,i)      - f5(:,1)                                             &
                                    / rho(2:nx-1,1,i) / dz(2:nx-1,1,i)
               ! add fluxes to top layer
               q(2:nx-1,nz,i)     = q(2:nx-1,nz,i) - (qinbc(2:nx-1,nz,i) * w(2:nx-1,nz,i) - f5(:,nz-1))      &
                                    / rho(2:nx-1,nz,i) / dz(2:nx-1,nz,i)
           else
               ! perform horizontal advection
               q(2:nx-1,:,i)      = q(2:nx-1,:,i)       - ((f1(2:nx-1,:) - f1(1:nx-2,:)) + (f3 - f4))
               ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
               q(2:nx-1,2:nz-1,i) = q(2:nx-1,2:nz-1,i)  - (f5(:,2:nz-1) - f5(:,1:nz-2))
               ! add fluxes to bottom layer
               q(2:nx-1,1,i)      = q(2:nx-1,1,i)       - f5(:,1)
               ! add fluxes to top layer
               q(2:nx-1,nz,i)     = q(2:nx-1,nz,i) - (qinbc(2:nx-1,nz,i) * w(2:nx-1,nz,i) - f5(:,nz-1))               
           endif
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine advect3d

    subroutine setup_cu_winds(domain, options, dt)
        implicit none
        type(domain_type),  intent(in) :: domain
        type(options_type), intent(in) :: options
        real,               intent(in) :: dt

        real    :: dx
        integer :: nx,nz,ny

        dx = domain%dx
        nx = size(domain%dz,1)
        nz = size(domain%dz,2)
        ny = size(domain%dz,3)

        if (.not.allocated(U_4cu_u)) then
            allocate(U_4cu_u(nx,  nz, ny))
            U_4cu_u = 0
            allocate(V_4cu_u(nx+1,nz, ny-1))
            V_4cu_u = 0
            allocate(W_4cu_u(nx+1,nz, ny))
            W_4cu_u = 0

            allocate(U_4cu_v(nx-1,nz, ny+1))
            U_4cu_v = 0
            allocate(V_4cu_v(nx,  nz, ny))
            V_4cu_v = 0
            allocate(W_4cu_v(nx,  nz, ny+1))
            W_4cu_v = 0
        endif

        U_4cu_u           =  (dt/dx) * (domain%u(1:nx,:,:)      + domain%u(2:nx+1,:,:)) / 2
        V_4cu_u(2:nx,:,:) =  (dt/dx) * (domain%v(1:nx-1,:,2:ny) + domain%v(2:nx,:,2:ny)) / 2
        W_4cu_u(2:nx,:,:) =  (dt/dx) * (domain%w(1:nx-1,:,:)    + domain%w(2:nx,:,:)) / 2
        call rebalance_cu_winds(U_4cu_u, V_4cu_u, W_4cu_u)

        U_4cu_v(:,:,2:ny) =  (dt/dx) * (domain%u(2:nx,:,1:ny-1) + domain%u(2:nx,:,2:ny)) / 2
        V_4cu_v           =  (dt/dx) * (domain%v(:,:,1:ny)      + domain%v(:,:,2:ny+1)) / 2
        W_4cu_v(:,:,2:ny) =  (dt/dx) * (domain%w(:,:,1:ny-1)    + domain%w(:,:,2:ny)) / 2
        call rebalance_cu_winds(U_4cu_v, V_4cu_v, W_4cu_v)

    end subroutine setup_cu_winds

    subroutine rebalance_cu_winds(u,v,w)
        implicit none
        ! u, v, w 3D east-west, south-north, and up-down winds repsectively
        ! note for this code, u is [nx-1,nz,ny] and v is [nx,nz,ny-1]
        real, dimension(:,:,:), intent(inout) :: u, v, w

        real, allocatable, dimension(:,:) :: divergence, du, dv
        integer :: i,nx,ny,nz

        nx = size(w,1)
        nz = size(w,2)
        ny = size(w,3)

        allocate(divergence(nx-2,ny-2))
        allocate(du(nx-2,ny-2))
        allocate(dv(nx-2,ny-2))

        do i=1,nz
            ! calculate horizontal divergence
            dv = v(2:nx-1,i,2:ny-1) - v(2:nx-1,i,1:ny-2)
            du = u(2:nx-1,i,2:ny-1) - u(1:nx-2,i,2:ny-1)
            divergence = du + dv
            ! Then calculate w to balance
            if (i==1) then
                ! if this is the first model level start from 0 at the ground
                w(2:nx-1,i,2:ny-1) = 0 - divergence
            else
                ! else calculate w as a change from w at the level below
                w(2:nx-1,i,2:ny-1) = w(2:nx-1,i-1,2:ny-1) - divergence
            endif
        enddo

    end subroutine rebalance_cu_winds

    subroutine advect_cu_winds(domain, options, dt)
        implicit none
        type(domain_type),  intent(inout) :: domain
        type(options_type), intent(in)    :: options
        real,               intent(in)    :: dt

        integer :: nx,nz,ny

        nx = size(domain%dz,1)
        nz = size(domain%dz,2)
        ny = size(domain%dz,3)

        ! first put the background u,v,w winds on a staggered grid with respect to the u grid
        ! then advect the u winds
        if (options%advect_density) then
            print*, "ERROR: Density advection not enabled when using convective winds"
            print*, "   Requires update to wind.f90 balance_uvw and advect.f90 (at least)"
            stop
        endif

        call setup_cu_winds(domain, options, dt)

        ! set the top boundary condition for CU winds to 0 to prevent artifacts coming in from the "top"
        domain%u_cu(:,nz,:) = 0
        domain%v_cu(:,nz,:) = 0

        call advect3d(domain%u_cu, U_4cu_u,V_4cu_u,W_4cu_u, domain%rho, domain%dz_inter, nx+1,nz,ny, options, 0) ! the last parameter (0) indicates the boundary condition to be used for advection
        call advect3d(domain%v_cu, U_4cu_v,V_4cu_v,W_4cu_v, domain%rho, domain%dz_inter, nx,nz,ny+1, options, 0) ! at the model top. 0 corresponds to the default zero gradient BC

    end subroutine advect_cu_winds

    subroutine setup_module_winds(domain, options, dt)
        implicit none
        type(domain_type),  intent(in)  :: domain
        type(options_type), intent(in)  :: options
        real,               intent(in)  :: dt

        real    :: dx
        integer :: nx, nz, ny, i

        dx = domain%dx
        nx = size(domain%dz,1)
        nz = size(domain%dz,2)
        ny = size(domain%dz,3)

        ! if this if the first time we are called, we need to allocate the module level arrays
        ! Could/should be put in an init procedure
        if (.not.allocated(U_m)) then
            allocate(U_m     (nx-1,nz,ny  ))
            allocate(V_m     (nx,  nz,ny-1))
            allocate(W_m     (nx,  nz,ny  ))
            allocate(lastqv_m(nx,  nz,ny  ))
        endif

        ! calculate U,V,W normalized for dt/dx (dx**2 for density advection so we can skip a /dx in the actual advection code)
        if (options%advect_density) then
            U_m = domain%ur(2:nx,:,:) * (dt/dx**2)
            V_m = domain%vr(:,:,2:ny) * (dt/dx**2)
            W_m = domain%wr           * (dt/dx**2)
        else
            if (options%physics%convection > 0) then
                U_m = (domain%u_cu(2:nx,:,:) + domain%u(2:nx,:,:)) * (dt/dx)
                V_m = (domain%v_cu(:,:,2:ny) + domain%v(:,:,2:ny)) * (dt/dx)
                W_m = (domain%w_cu + domain%w)                     * (dt/dx)
                call rebalance_cu_winds(U_m,V_m,W_m)
            else
                U_m = domain%u(2:nx,:,:) * (dt/dx)
                V_m = domain%v(:,:,2:ny) * (dt/dx)
                ! note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
                W_m = domain%w           * (dt/dx)
            endif
        endif

    end subroutine setup_module_winds

    ! primary entry point, advect all scalars in domain
    subroutine upwind(domain,options,dt)
        implicit none
        type(domain_type),  intent(inout) :: domain
        type(options_type), intent(in)    :: options
        real,               intent(in)    :: dt

        real    :: dx
        integer :: nx, nz, ny, i

        nx = size(domain%dz,1)
        nz = size(domain%dz,2)
        ny = size(domain%dz,3)

        ! calculate U,V,W normalized for dt/dx (dx**2 for density advection so we can skip a /dx in the actual advection code)
        call setup_module_winds(domain, options, dt)

        lastqv_m=domain%qv

        ! 20190402 jh ... added the boundary condition to be used at the model top to the function call of
        !                 advect3d. depending on the option set the behaviour is modified accordingly.
        !                 the boundary condition settings currently implemented are
        !
        !                 0 ... zero gradient on the quantitiy (default)
        !                 1 ... constant gradient on the quantitiy
        !                 2 ... zero value on the quantity
        !                 3 ... zero value on the flux gradient over the topmost vertical layer
        !                   
        !                 currently implemented for upwind advection
        
        call advect3d(domain%qv,          U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
        call advect3d(domain%cloud,       U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
        call advect3d(domain%qrain,       U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
        call advect3d(domain%qsnow,       U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
        call advect3d(domain%th,          U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_th_top)
        if (options%physics%microphysics == kMP_THOMPSON) then
            call advect3d(domain%ice,     U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%qgrau,   U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%nice,    U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%nrain,   U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
        elseif (options%physics%microphysics == kMP_MORRISON) then
            call advect3d(domain%ice,     U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%qgrau,   U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%nice,    U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%nrain,   U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%nsnow,   U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%ngraupel,U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
        elseif (options%physics%microphysics == kMP_WSM6) then
            call advect3d(domain%ice,     U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
            call advect3d(domain%qgrau,   U_m,V_m,W_m, domain%rho, domain%dz_inter, nx,nz,ny, options,options%adv_options%bc_top)
        endif

        ! if (options%physics%convection > 0) then
        !     call advect_cu_winds(domain, options, dt)
        ! endif

        ! used in some physics routines
        domain%tend%qv_adv = (domain%qv - lastqv_m) / dt
    end subroutine upwind

end module adv_upwind
