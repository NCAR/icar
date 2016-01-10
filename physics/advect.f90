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
    real,dimension(:,:,:),allocatable::U_m,V_m,W_m,lastqv_m
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

    subroutine advect3d(q,u,v,w,rho,dz,nx,nz,ny,debug,options)
        implicit none
        real,dimension(1:nx,1:nz,1:ny), intent(inout) :: q
        real,dimension(1:nx,1:nz,1:ny), intent(in) :: w
        real,dimension(1:nx-1,1:nz,1:ny),intent(in) :: u
        real,dimension(1:nx,1:nz,1:ny-1),intent(in) :: v
        real,dimension(1:nx,1:nz,1:ny), intent(in) :: rho
        real,dimension(1:nx,1:nz,1:ny), intent(in) :: dz
        integer, intent(in) :: ny,nz,nx,debug
        type(options_type), intent(in)::options
        ! interal parameters
        integer :: err,i
        real,dimension(1:nx,1:nz,1:ny) :: qin
        real, dimension(1:nx-1,1:nz) :: f1 ! there used to be an f2 to store f[x+1]
        real, dimension(1:nx-2,1:nz) :: f3,f4
        real, dimension(1:nx-2,1:nz-1) ::f5
        !$omp parallel shared(qin,q,u,v,w,rho,dz) firstprivate(nx,ny,nz) private(i,f1,f3,f4,f5)
        !$omp do schedule(static)
        do i=1,ny
            qin(:,:,i)=q(:,:,i)
        enddo
        !$omp end do
        !$omp barrier
        !$omp do schedule(static)
        do i=2,ny-1
!           by manually inlining the flux2 call we should remove extra array copies that the compiler doesn't remove. 
!           equivalent flux2 calls are left in for reference (commented) to restore recall that f1,f3,f4... arrays should be 3D : n x m x 1
!           calculate fluxes between grid cells
!            call flux2(qin(1:nx-1,:,i),     qin(2:nx,:,i),     u(1:nx-1,:,i),     nx-1,nz,  1,f1)  ! f1 = Ux0 and Ux1
!            call flux2(qin(2:nx-1,:,i),     qin(2:nx-1,:,i+1), v(2:nx-1,:,i),     nx-2,nz,  1,f3)  ! f3 = Vy1
!            call flux2(qin(2:nx-1,:,i-1),   qin(2:nx-1,:,i),   v(2:nx-1,:,i-1),   nx-2,nz,  1,f4)  ! f4 = Vy0
!            call flux2(qin(2:nx-1,1:nz-1,i),qin(2:nx-1,2:nz,i),w(2:nx-1,1:nz-1,i),nx-2,nz-1,1,f5)  ! f5 = Wz0 and Wz1
           f1= ((u(1:nx-1,:,i)      + ABS(u(1:nx-1,:,i)))      * qin(1:nx-1,:,i) + &
                (u(1:nx-1,:,i)      - ABS(u(1:nx-1,:,i)))      * qin(2:nx,:,i))/2
           f3= ((v(2:nx-1,:,i)      + ABS(v(2:nx-1,:,i)))      * qin(2:nx-1,:,i) + &
                (v(2:nx-1,:,i)      - ABS(v(2:nx-1,:,i)))      * qin(2:nx-1,:,i+1))/2
           f4= ((v(2:nx-1,:,i-1)    + ABS(v(2:nx-1,:,i-1)))    * qin(2:nx-1,:,i-1) + &
                (v(2:nx-1,:,i-1)    - ABS(v(2:nx-1,:,i-1)))    * qin(2:nx-1,:,i))/2
           f5= ((w(2:nx-1,1:nz-1,i) + ABS(w(2:nx-1,1:nz-1,i))) * qin(2:nx-1,1:nz-1,i) + &
                (w(2:nx-1,1:nz-1,i) - ABS(w(2:nx-1,1:nz-1,i))) * qin(2:nx-1,2:nz,i))/2
           
           if (options%advect_density) then
               ! perform horizontal advection
               q(2:nx-1,:,i)=q(2:nx-1,:,i) - ((f1(2:nx-1,:)-f1(1:nx-2,:)) + (f3(:,:)-f4(:,:)))/rho(2:nx-1,:,i)/dz(2:nx-1,:,i)
               ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
               q(2:nx-1,2:nz-1,i)=q(2:nx-1,2:nz-1,i)-(f5(:,2:nz-1)-f5(:,1:nz-2))/rho(2:nx-1,2:nz-1,i)/dz(2:nx-1,2:nz-1,i)
               ! add fluxes to bottom layer
               q(2:nx-1,1,i)=q(2:nx-1,1,i)-f5(:,1)/rho(2:nx-1,1,i)/dz(2:nx-1,1,i)
               ! add fluxes to top layer
               q(2:nx-1,nz,i)=q(2:nx-1,nz,i)-(qin(2:nx-1,nz,i)*w(2:nx-1,nz,i)-f5(:,nz-1))/dz(2:nx-1,nz,i)/rho(2:nx-1,nz,i)
           else
               ! perform horizontal advection
               q(2:nx-1,:,i)=q(2:nx-1,:,i) - ((f1(2:nx-1,:)-f1(1:nx-2,:)) + (f3(:,:)-f4(:,:)))
               ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
               q(2:nx-1,2:nz-1,i)=q(2:nx-1,2:nz-1,i)-(f5(:,2:nz-1)-f5(:,1:nz-2))
               ! add fluxes to bottom layer
               q(2:nx-1,1,i)=q(2:nx-1,1,i)-f5(:,1)
               ! add fluxes to top layer
               q(2:nx-1,nz,i)=q(2:nx-1,nz,i)-(qin(2:nx-1,nz,i)*w(2:nx-1,nz,i)-f5(:,nz-1))
           endif
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine advect3d

!   primary entry point, advect all scalars in domain
    subroutine upwind(domain,options,dt)
        implicit none
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        real,intent(in)::dt
        
        real::dx
        integer::nx,nz,ny,i
            
        dx=domain%dx
        nx=size(domain%dz,1)
        nz=size(domain%dz,2)
        ny=size(domain%dz,3)
        
!       if this if the first time we are called, we need to allocate the module level arrays
        if (.not.allocated(U_m)) then
            allocate(U_m(nx-1,nz,ny))
            allocate(V_m(nx,nz,ny-1))
            allocate(W_m(nx,nz,ny))
            allocate(lastqv_m(nx,nz,ny))
        endif
        lastqv_m=domain%qv
!       calculate U,V,W normalized for dt/dx (dx**2 for density advection so we can skip a /dx in the actual advection code)
        if (options%advect_density) then
            U_m=domain%ur(2:nx,:,:)*(dt/dx**2)
            V_m=domain%vr(:,:,2:ny)*(dt/dx**2)
            W_m=domain%wr*(dt/dx**2)
        else
            U_m=domain%u(2:nx,:,:)*(dt/dx)
            V_m=domain%v(:,:,2:ny)*(dt/dx)
    !       note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
            W_m=domain%w*(dt/dx)
        endif

        call advect3d(domain%qv,   U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
        call advect3d(domain%cloud,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
        call advect3d(domain%qrain,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
        call advect3d(domain%qsnow,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
        call advect3d(domain%th,   U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
        if (options%physics%microphysics==kMP_THOMPSON) then
            call advect3d(domain%ice,  U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
            call advect3d(domain%qgrau,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
            call advect3d(domain%nice, U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
            call advect3d(domain%nrain,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
        endif
        
        domain%tend%qv_adv=(domain%qv-lastqv_m)/dt
    end subroutine upwind
    
end module adv_upwind
