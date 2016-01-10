!> ----------------------------------------------------------------------------
!!  The MPDATA advection scheme
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module adv_mpdata
    use data_structures
    
    implicit none
    private
    real,dimension(:,:,:),allocatable::U_m,V_m,W_m
    integer :: order
    
    public:: mpdata, mpdata_init
    public:: advect3d ! for test_mpdata testing only!
    

contains

    subroutine flux1(l,r,U,f)
    !     Calculate the donor cell flux function
    !     l = left gridcell scalar 
    !     r = right gridcell scalar
    !     U = Courant number (u*dt/dx)
    !     
    !     If U is positive, return l*U if U is negative return r*U
    !     By using the mathematical form instead of the logical form, 
    !     we can run on the entire grid simultaneously, and avoid branches

    !   arguments
        implicit none
        real, dimension(:), intent(in) :: l,r,U
        real, dimension(:), intent(inout) :: f
        
        !   main code
        f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2

    end subroutine flux1
    
    subroutine upwind_advection(q, u, v, w, q2, nx,nz,ny)
        implicit none
        real,dimension(1:nx,1:nz,1:ny),  intent(in) :: q
        real,dimension(1:nx-1,1:nz,1:ny),intent(in) :: u
        real,dimension(1:nx,1:nz,1:ny-1),intent(in) :: v
        real,dimension(1:nx,1:nz,1:ny),  intent(in) :: w
        real,dimension(1:nx,1:nz,1:ny),  intent(inout) :: q2
        integer, intent(in) :: ny,nz,nx
        
        ! interal parameters
        integer :: i
        real, dimension(1:nx-1,1:nz) :: f1 ! there used to be an f2 to store f[x+1]
        real, dimension(1:nx-2,1:nz) :: f3,f4
        real, dimension(1:nx-2,1:nz-1) ::f5
        !$omp parallel shared(q2,q,u,v,w) firstprivate(nx,ny,nz) private(i,f1,f3,f4,f5)
        !$omp do schedule(static)
        do i=1,ny
            if ((i==1).or.(i==ny)) then
                q2(:,:,i)=q(:,:,i)
            else
                q2(1,:,i)=q(1, :,i)
                q2(nx,:,i)=q(nx, :,i)
                
!           by manually inlining the flux2 call we should remove extra array copies that the compiler doesn't remove. 
!           equivalent flux2 calls are left in for reference (commented) to restore recall that f1,f3,f4... arrays should be 3D : n x m x 1
!           calculate fluxes between grid cells
!            call flux2(q(1:nx-1,:,i),     q(2:nx,:,i),     u(1:nx-1,:,i),     nx-1,nz,  1,f1)  ! f1 = Ux0 and Ux1
!            call flux2(q(2:nx-1,:,i),     q(2:nx-1,:,i+1), v(2:nx-1,:,i),     nx-2,nz,  1,f3)  ! f3 = Vy1
!            call flux2(q(2:nx-1,:,i-1),   q(2:nx-1,:,i),   v(2:nx-1,:,i-1),   nx-2,nz,  1,f4)  ! f4 = Vy0
!            call flux2(q(2:nx-1,1:nz-1,i),q(2:nx-1,2:nz,i),w(2:nx-1,1:nz-1,i),nx-2,nz-1,1,f5)  ! f5 = Wz0 and Wz1
               f1= ((u(1:nx-1,:,i)      + ABS(u(1:nx-1,:,i)))      * q(1:nx-1,:,i) + &
                    (u(1:nx-1,:,i)      - ABS(u(1:nx-1,:,i)))      * q(2:nx,:,i))/2
               f3= ((v(2:nx-1,:,i)      + ABS(v(2:nx-1,:,i)))      * q(2:nx-1,:,i) + &
                    (v(2:nx-1,:,i)      - ABS(v(2:nx-1,:,i)))      * q(2:nx-1,:,i+1))/2
               f4= ((v(2:nx-1,:,i-1)    + ABS(v(2:nx-1,:,i-1)))    * q(2:nx-1,:,i-1) + &
                    (v(2:nx-1,:,i-1)    - ABS(v(2:nx-1,:,i-1)))    * q(2:nx-1,:,i))/2
               f5= ((w(2:nx-1,1:nz-1,i) + ABS(w(2:nx-1,1:nz-1,i))) * q(2:nx-1,1:nz-1,i) + &
                    (w(2:nx-1,1:nz-1,i) - ABS(w(2:nx-1,1:nz-1,i))) * q(2:nx-1,2:nz,i))/2
           
               ! perform horizontal advection
               q2(2:nx-1,:,i)=q(2:nx-1,:,i) - ((f1(2:nx-1,:)-f1(1:nx-2,:)) + (f3(:,:)-f4(:,:)))
               ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
               q2(2:nx-1,2:nz-1,i)=q2(2:nx-1,2:nz-1,i)-(f5(:,2:nz-1)-f5(:,1:nz-2))
               ! add fluxes to bottom layer
               q2(2:nx-1,1,i)=q2(2:nx-1,1,i)-f5(:,1)
               ! add fluxes to top layer
               q2(2:nx-1,nz,i)=q2(2:nx-1,nz,i)-(q(2:nx-1,nz,i)*w(2:nx-1,nz,i)-f5(:,nz-1))
           endif
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine upwind_advection

    subroutine mpdata_fluxes(q,u,v,w,u2,v2,w2, nx,nz,ny)
        implicit none
        real, dimension(nx,nz,ny),   intent(in) :: q,w
        real, dimension(nx-1,nz,ny), intent(in) :: u
        real, dimension(nx,nz,ny-1), intent(in) :: v
        real, dimension(nx-1,nz,ny), intent(out) :: u2
        real, dimension(nx,nz,ny-1), intent(out) :: v2
        real, dimension(nx,nz,ny),   intent(out) :: w2
        integer, intent(in) :: nx,ny,nz
        
        real, dimension(nx-1) :: rx, lx, denomx
        real, dimension(nx) :: r, l, denom
        integer :: i, j
        
        ! This might run faster if tiled over x and y to be more cache friendly. 
        !$omp parallel shared(q,u,v,w,u2,v2,w2) firstprivate(nx,ny,nz) &
        !$omp private(i,j, rx,lx,r,l, denomx,denom)
        !$omp do schedule(static)
        do i=1,ny
            do j=1,nz
                ! -----------------------
                ! First compute the U component
                ! -----------------------
                if ((i>1).and.(i<ny)) then
                    rx=q(2:nx,j,i)
                    lx=q(1:nx-1,j,i)
                    ! In MPDATA papers (r-l)/(r+l) is usually refered to as "A"
                    ! compute the denomenator first so we can check that it is not zero
                    denomx=(rx + lx)
                    where(denomx==0) denomx=1e-10
                    ! U2 is the diffusive pseudo-velocity
                    u2(:,j,i) = abs(u(:,j,i)) - u(:,j,i)**2
                    u2(:,j,i) = u2(:,j,i) * (rx-lx) / denomx
                else
                    u2(:,j,i)=0
                endif
                

                ! next compute the V and W components
                if (i==1) then
                    w2(:,j,i)=0
                else
                    ! -----------------------
                    ! compute the V component
                    ! -----------------------
                    r=q(:,j,i)
                    l=q(:,j,i-1)
                    ! In MPDATA papers A = (r-l)/(r+l)
                    ! compute the denomenator first so we can check that it is not zero
                    denom=(r + l)
                    where(denom==0) denom=1e-10
                    ! U2 is the diffusive pseudo-velocity
                    v2(:,j,i-1) = abs(v(:,j,i-1)) - v(:,j,i-1)**2
                    v2(:,j,i-1) = v2(:,j,i-1) * (r-l) / denom
                    
                    
                    ! -----------------------
                    ! compute the w component
                    ! -----------------------
                    if (i==ny) then
                        w2(:,j,i)=0
                    else
                        if (j<nz) then
                            r=q(:,j+1,i)
                            l=q(:,j,i)
                            ! In MPDATA papers A = (r-l)/(r+l)
                            ! compute the denomenator first so we can check that it is not zero
                            denom=(r + l)
                            where(denom==0) denom=1e-10
                            ! U2 is the diffusive pseudo-velocity
                            w2(:,j,i) = abs(w(:,j,i)) - w(:,j,i)**2
                            w2(:,j,i) = w2(:,j,i) * (r-l) / denom
                        else
                            w2(:,j,i) = 0
                        endif
                    endif
                    
                endif
            end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine mpdata_fluxes

    subroutine flux_limiter(q, q2, u,v,w, nx,nz,ny)
        implicit none
        real,dimension(1:nx,1:nz,1:ny),  intent(in)    :: q, q2
        real,dimension(1:nx-1,1:nz,1:ny),intent(inout) :: u
        real,dimension(1:nx,1:nz,1:ny-1),intent(inout) :: v
        real,dimension(1:nx,1:nz,1:ny),  intent(inout) :: w
        integer, intent(in) :: nx,nz,ny
        
        integer :: i,j,k,n
        real, dimension(:), pointer :: q1, U2, l, f
        ! q1 = q after applying previous iteration advection
        ! l  = q before applying previous iteration
        ! U2 is the anti-diffusion pseudo-velocity
        ! f is the first pass calculation of MPDATA fluxes
        real, dimension(nx),   target :: q1x,lx
        real, dimension(nx-1), target :: fx, U2x
        real, dimension(ny),   target :: q1y,ly
        real, dimension(ny-1), target :: fy, U2y
        real, dimension(nz),   target :: q1z,lz
        real, dimension(nz-1), target :: fz, U2z
        logical :: flux_is_w
        
        real :: qmax_i,qmin_i,qmax_i2,qmin_i2
        real :: beta_in_i, beta_out_i, beta_in_i2, beta_out_i2
        real :: fin_i, fout_i, fin_i2, fout_i2
        
        ! NOTE: before inclusion of FCT_core the following variables must be setup: 
        ! q1 and l (l=q0)
        !$omp parallel shared(q2,q,u,v,w) firstprivate(nx,ny,nz) default(private)
        !$omp do schedule(static)
        do j=2,ny-1
            flux_is_w=.False.
            n=nx
            q1=>q1x
            l =>lx
            U2=>U2x
            f=>fx
            do k=1,nz
                ! setup u
                q1=q2(:,k,j)
                U2=u(:,k,j)
                l =q(:,k,j)
                call flux1(q1(1:n-1),q1(2:n),U2,f)
                
                include "adv_mpdata_FCT_core.f90"
                u(:,k,j)=U2
            end do
            
            n=nz
            q1=>q1z
            l =>lz
            U2=>U2z
            f=>fz
            flux_is_w=.True.
            do k=2,nx-1
                ! setup w
                q1=q2(k,:,j)
                U2=w(k,1:n-1,j)
                l =q(k,:,j)
                call flux1(q1(1:n-1),q1(2:n),U2,f)
                ! NOTE: need to check this a little more
                include "adv_mpdata_FCT_core.f90"
                w(k,1:n-1,j)=U2
                w(k,n,j)=0
            end do
            
        end do
        !$omp end do
        
        flux_is_w=.False.
        n=ny
        q1=>q1y
        l =>ly
        U2=>U2y
        f=>fy
        ! NOTE: This it typically not the correct order for the loop variables
        ! but in this case it permits parallelization over a larger number (nx instead of nz)
        ! and because all data are copied from an oddly spaced grid regardless, it *probably* doesn't slow it down
        ! I'd like to re-write the v-flux delimiter to operate on all x simulataneously at some point...
        !$omp do
        do j=1,nx
            do k=1,nz
                q1=q2(j,k,:)
                U2=v(j,k,:)
                l =q(j,k,:)
                call flux1(q1(1:n-1),q1(2:n),U2,f)
                
                include "adv_mpdata_FCT_core.f90"
                v(j,k,:)=U2
            end do
        end do
        !$omp end do
        !$omp end parallel
        
    end subroutine flux_limiter

    subroutine advect3d(q,u,v,w,rho,dz,nx,nz,ny,options, err)
        implicit none
        real,dimension(1:nx,1:nz,1:ny), intent(inout) :: q
        real,dimension(1:nx-1,1:nz,1:ny),intent(in) :: u
        real,dimension(1:nx,1:nz,1:ny-1),intent(in) :: v
        real,dimension(1:nx,1:nz,1:ny), intent(in) :: w
        real,dimension(1:nx,1:nz,1:ny), intent(in) :: rho
        real,dimension(1:nx,1:nz,1:ny), intent(in) :: dz
        integer, intent(in) :: ny,nz,nx
        type(options_type), intent(in)::options
        integer, intent(inout) :: err

        ! used for intermediate values in the mpdata calculation
        real,dimension(1:nx,1:nz,1:ny)   :: q2
        real,dimension(1:nx-1,1:nz,1:ny) :: u2
        real,dimension(1:nx,1:nz,1:ny-1) :: v2
        real,dimension(1:nx,1:nz,1:ny)   :: w2
        
        integer :: iord, i
        
        do iord=1,options%adv_options%mpdata_order
            if (iord==1) then
                call upwind_advection(q, u, v, w, q2, nx,nz,ny)
            else
                call mpdata_fluxes(q2, u, v, w, u2,v2,w2, nx,nz,ny)
        
                if (options%adv_options%flux_corrected_transport) then
                    call flux_limiter(q, q2, u2,v2,w2, nx,nz,ny)
                endif
                
                call upwind_advection(q2, u2, v2, w2, q, nx,nz,ny)
            endif
            
            ! 
            if (iord/=options%adv_options%mpdata_order) then
                if (iord>1) then
                    !$omp parallel shared(q,q2) firstprivate(ny) private(i)
                    !$omp do schedule(static)
                    do i=1,ny
                        q2(:,:,i)=q(:,:,i)
                    enddo
                    !$omp end do
                    !$omp end parallel
                endif
            else 
                if (iord==1) then
                    !$omp parallel shared(q,q2) firstprivate(ny) private(i)
                    !$omp do schedule(static)
                    do i=1,ny
                        q(:,:,i)=q2(:,:,i)
                    enddo
                    !$omp end do
                    !$omp end parallel
                endif
                
            endif
        end do
        

    end subroutine advect3d
    
    subroutine mpdata_init(domain,options)
        type(domain_type), intent(inout) :: domain
        type(options_type), intent(in) :: options
        
        order    = 0
    end subroutine mpdata_init
    
!   primary entry point, advect all scalars in domain
    subroutine mpdata(domain,options,dt)
        implicit none
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        real,intent(in)::dt
        
        real::dx
        integer::nx,nz,ny,i, error
        
        dx=domain%dx
        nx=size(domain%dz,1)
        nz=size(domain%dz,2)
        ny=size(domain%dz,3)
        
!       if this if the first time we are called, we need to allocate the module level arrays
        if (.not.allocated(U_m)) then
            allocate(U_m(nx-1,nz,ny))
            allocate(V_m(nx,nz,ny-1))
            allocate(W_m(nx,nz,ny))
        endif
        
!       calculate U,V,W normalized for dt/dx
        if (options%advect_density) then
            U_m=domain%ur(2:nx,:,:)*(dt/dx**2)
            V_m=domain%vr(:,:,2:ny)*(dt/dx**2)
!           note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
            W_m=domain%wr*(dt/dx**2)
        else
            U_m=domain%u(2:nx,:,:)*(dt/dx)
            V_m=domain%v(:,:,2:ny)*(dt/dx)
!           note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
            W_m=domain%w*(dt/dx)
        endif
        
        error=0
        
!         print*, "qv"
        call advect3d(domain%qv,   U_m,V_m,W_m, domain%rho,domain%dz,nx,nz,ny,options, error)
        if (error/=0) then
            print*, "qv", error
            if (error>0) then
                stop "MPDATA Error"
            endif
        endif
        
        error=0
!         print*, "qc"
        call advect3d(domain%cloud,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,options, error)
        if (error/=0) then
            print*, "qc", error
            if (error>0) then
                stop "MPDATA Error"
            endif
        endif

        error=0
!         print*, "qr"
        call advect3d(domain%qrain,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,options, error)
        if (error/=0) then
            print*, "qr", error
            if (error>0) then
                stop "MPDATA Error"
            endif
        endif

        error=0
!         print*, "qs"
        call advect3d(domain%qsnow,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,options, error)
        if (error/=0) then
            print*, "qs", error
            if (error>0) then
                stop "MPDATA Error"
            endif
        endif

        error=0
!         print*, "th"
        call advect3d(domain%th,   U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,options, error)
        if (error/=0) then
            print*, "th", error
            if (error>0) then
                stop "MPDATA Error"
            endif
        endif
        if (options%physics%microphysics==kMP_THOMPSON) then
            error=0
!             print*, "qi"
            call advect3d(domain%ice,  U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,options, error)
            if (error/=0) then
                print*, "qi", error
                if (error>0) then
                    stop "MPDATA Error"
                endif
            endif

            error=0
!             print*, "qg"
            call advect3d(domain%qgrau,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,options, error)
            if (error/=0) then
                print*, "qg", error
                if (error>0) then
                    stop "MPDATA Error"
                endif
            endif

            error=0
!             print*, "ni"
            call advect3d(domain%nice, U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,options, error)
            if ((error/=0).and.(error/=2)) then
                print*, "ni", error
!                 if (error>0) then
!                     stop "MPDATA Error"
!                 endif
            endif

            error=0
!             print*, "nr"
            call advect3d(domain%nrain,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,options, error)
            if ((error/=0).and.(error/=2)) then
                print*, "nr", error
!                 if (error>0) then
!                     stop "MPDATA Error"
!                 endif
            endif
        endif
        order=mod(order+1,3)
        
    end subroutine mpdata
end module adv_mpdata
