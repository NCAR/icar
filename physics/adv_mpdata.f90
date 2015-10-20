!> ----------------------------------------------------------------------------
!!
!!  The MPDATA advection scheme
!!  Not fully implemented / translated from python yet
!!
!!  Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module adv_mpdata
    use data_structures
    
    implicit none
    private
    real,dimension(:,:,:),allocatable::U_m,V_m,W_m
    integer :: order
    public:: mpdata, mpdata_init
    public:: advect_u, advect_w, advect_v ! for test_mpdata testing only!

contains
    subroutine flux2(l,r,U,nx,nz,ny,f)
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
        real, dimension(1:nx,1:nz,1:ny), intent(in) :: l,r,U
        real, dimension(1:nx,1:nz,1:ny), intent(inout) :: f
        integer,intent(in) :: ny,nz,nx
        !   internal parameter
        integer ::  err,i!,j,Ny,Nz,Nx
        !   main code
        f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2

    end subroutine flux2

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

    subroutine advect2d(q,u,v,rho,dz,nx,nz,ny,options)
!     horizontally advect a scalar q by wind field (u,v)
!     q = input scalar field (e.g. cloud water [kg])
!     u,v = horizontal and vertical wind speeds [m/s] on a staggered grid. normalized by dt/dx
!     
!     Algorithm from MPDATA: 
!         Smolarkiewicz, PK and Margolin, LG (1998)
!             MPDATA: A Finite-Differnce Solver for Geophysical Flows
!             Journal of Computational Physics 140, p459-480. CP985901
        real,dimension(:,:,:), intent(inout) :: q,u,v,rho,dz
        integer,intent(in)::nx,nz,ny
        type(options_type),intent(in) :: options
        
!     # The above calculation is just the standard upwind formulations
!     # below, calculate the diffusivity correction term of MPDATA
!     # f=np.double(0.5)
!     # 
!     # # define A,B [l,r,u,b] (see MPDATA review reference)
!     # # l,r,u,b = left, right, upper, bottom edges of the grid cells
!     # Al=(q1(1:-1,:-2)-q1(1:-1,1:-1))/(q1(1:-1,:-2)+q1(1:-1,1:-1))
!     # Ar=(q1(1:-1,2:)-q1(1:-1,1:-1))/(q1(1:-1,2:)+q1(1:-1,1:-1))
!     # Au=(q1(:-2,1:-1)-q1(1:-1,1:-1))/(q1(:-2,1:-1)+q1(1:-1,1:-1))
!     # Ab=(q1(2:,1:-1)-q1(1:-1,1:-1))/(q1(2:,1:-1)+q1(1:-1,1:-1))
!     # 
!     # q11=q1(2:,2:)+q1(2:,1:-1)
!     # Br=0.5*((q11-q1(:-2,2:)-q1(:-2,1:-1))
!     #     /(q11+q1(:-2,2:)+q1(:-2,1:-1)))
!     # q11=q1(2:,:-2)+q1(2:,1:-1)
!     # Bl=0.5*((q11-q1(:-2,:-2)-q1(:-2,1:-1))
!     #     /(q11+q1(:-2,:-2)+q1(:-2,1:-1)))
!     # q11=q1(:-2,2:)+q1(1:-1,2:)
!     # Bu=0.5*((q11-q1(:-2,:-2)-q1(1:-1,:-2))
!     #     /(q11+q1(:-2,:-2)+q1(1:-1,:-2)))
!     # q11=q1(2:,2:)+q1(1:-1,2:)
!     # Bb=0.5*((q11-q1(2:,:-2)-q1(1:-1,:-2))
!     #     /(q11+q1(2:,:-2)+q1(1:-1,:-2)))
!     # 
!     # # compute diffusion correction U/V terms (see MPDATA review reference)
!     # Uabs=np.abs(U)
!     # # first find U/V terms on the grid cell borders
!     # curUabs=(Uabs[1:-1,:-2)+Uabs[1:-1,1:-1))/2
!     # curU=Ux0
!     # curV=(V[1:-1,:-2)+V[1:-1,1:-1))/2
!     # # then compute Ul
!     # Ul=curUabs*(1-curUabs)*Al - 2*f*curU*curV*Bl
!     # # compute Ur using the same two steps
!     # curUabs=(Uabs[1:-1,2:)+Uabs[1:-1,1:-1))/2
!     # curU=Ux1
!     # curV=(V[1:-1,2:)+V[1:-1,1:-1))/2
!     # Ur=curUabs*(1-curUabs)*Ar - 2*f*curU*curV*Br
!     # # compute Vu
!     # Vabs=np.abs(V)
!     # curVabs=(Vabs[:-2,1:-1)+Vabs[1:-1,1:-1))/2
!     # curV=Vy0
!     # curU=(U[:-2,1:-1)+U[1:-1,1:-1))/2
!     # Vu=curVabs*(1-curVabs)*Bu - 2*f*curU*curV*Bu
!     # # compute Vb
!     # curVabs=(Vabs[2:,1:-1)+Vabs[1:-1,1:-1))/2
!     # curV=Vy1
!     # curU=(U[2:,1:-1)+U[1:-1,1:-1))/2
!     # Vb=curVabs*(1-curVabs)*Bb - 2*f*curU*curV*Bb
!     # 
!     # q[1:-1,1:-1)=(q1(1:-1,1:-1)
!     #     -(F(q1(1:-1,1:-1),q1(1:-1,2:),Ul)
!     #     -F(q1(1:-1,:-2),q1(1:-1,1:-1),Ur))
!     #     -(F(q1(1:-1,1:-1),q1(2:,1:-1),Vu)
!     #     -F(q1(:-2,1:-1),q1(1:-1,1:-1),Vb)))
    end subroutine advect2d

    subroutine advect_w(q,u,nx,n,ny)
        real,dimension(:,:,:), intent(inout) :: q
        real,dimension(:,:,:), intent(in) :: u
        integer,intent(in)::n,nx,ny
        real,dimension(n)::q1
        real,dimension(n-1)::f,l,r,U2, denom
!         real,dimension(n-2)::fl,fr, l2,c2,r2, U3, denom, Vl, Vr
        integer::y, x
        
        do y=2,ny-1
            do x=2,nx-1
                l  = q(x,1:n-1,y)
                r  = q(x,2:n,y)
                U2 = u(x,1:n-1,y)
                
                ! This is the core 1D implementation of mpdata
                ! advect_w requires its own implementation to handle the top and bottom boundaries
                call flux1(l,r,U2,f)
                
                ! for advect_w we need to add the fluxes to the top and bottom boundaries too. 
                q1(1)=l(1)-f(1)
                ! note inlined flux1 call for the top boundary condition
                q1(n)=r(n-1)+f(n-1) - &
                        ((u(x,n,y)+ABS(u(x,n,y))) * r(n-1) + (u(x,n,y)-ABS(u(x,n,y))) * r(n-1))/2
                
                q1(2:n-1) = l(2:n-1) + (f(1:n-2) - f(2:n-1))

                ! This is the core 1D implementation of MPDATA correction
                
                ! This is the MPDATA diffusion correction term for 1D flow
                ! U is defined on the mass grid for the pseudo-velocities?
                ! left and right pseudo (diffusive) velocities
    
                ! we will copy the q1 data into r to potentially minimize aliasing problems 
                ! for the compiler, and improve memory alignment for vectorization
                r  = q1(2:n) 
                ! l  = q1(1:n-1) ! no need to copy these data over again
    
                ! In MPDATA papers (r-l)/(r+l) is usually refered to as "A"
                ! compute the denomenator first so we can check that it is not zero
                denom=(r + q1(1:n-1))
                where(denom==0) denom=1e-10
                ! U2 is the diffusive pseudo-velocity
                U2 = abs(U2) - U2**2
                U2 = U2 * (r-q1(1:n-1)) / denom
    
                ! now calculate the MPDATA flux term
                call flux1(q1(1:n-1),r,U2,f)
                
                q(x,2:n-1,y) = q1(2:n-1) + (f(1:n-2) - f(2:n-1))
                q(x,1,y) = q1(1) - f(1)
                q(x,n,y) = q1(n) + f(n-1) 
                ! note we don't subtract diffusive fluxes out the top as the layer above it is presumed to be the same
                ! as a result the diffusion pseudo-velocity term would be 0
            enddo
        enddo
    end subroutine advect_w

    subroutine advect_v(q,u,nx,nz,n)
        real,dimension(:,:,:), intent(inout) :: q
        real,dimension(:,:,:), intent(in) :: u
        integer,intent(in)::n,nz,nx
        real,dimension(n)::q1
        real,dimension(n-1)::f,l,r,U2, denom
        integer::x, z
        
        ! might be more cache friendly if we tile over x? 
        ! though that won't work with including mpdata_core
        do z=1,nz
            do x=2,nx-1
                l  = q(x,z,1:n-1)
                r  = q(x,z,2:n)
                U2 = u(x,z,1:n-1)
                
                include 'adv_mpdata_core.f90'
                
                q(x,z,2:n-1) = q1(2:n-1) + (f(:n-2) - f(2:n-1))
            enddo
        enddo
    end subroutine advect_v

    subroutine advect_u(q,u,n,nz,ny)
        real,dimension(:,:,:), intent(inout) :: q
        real,dimension(:,:,:), intent(in) :: u
        integer,intent(in)::n,nz,ny
        real,dimension(n)::q1
        real,dimension(n-1)::f,l,r,U2, denom
        integer::y, z
        
        ! loop over internal y columns
        do y=2,ny-1
            ! loop over all z layers
            do z=1,nz
                
                ! copy the data into local (cached) variables.  Makes the include mpdata_core possible
                l  = q(1:n-1,z,y)
                r  = q(2:n,z,y)
                U2 = u(1:n-1,z,y)

                include 'adv_mpdata_core.f90'
                
                q(2:n-1,z,y) = q1(2:n-1) + (f(:n-2) - f(2:n-1))
            enddo
        enddo
    end subroutine advect_u


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
        
        ! perform an Alternating Direction Explicit MP-DATA time step
        ! but swap the order of the alternating directions with every call
        if (order==0) then
!             if (any(isnan(q))) then
!                 print*, "pre"
!                 write(*,*) maxval(q)
!                 stop "MPDATA:NaN error"
!             endif
            
            call advect_u(q,u,nx,nz,ny)
!             if (any(isnan(q))) then
!                 print*, "post-u pre-v"
!                 write(*,*) maxval(q)
!                 stop "MPDATA:NaN error"
!             endif
            call advect_v(q,v,nx,nz,ny)
!             if (any(isnan(q))) then
!                 print*, "post-v pre-w"
!                 write(*,*) maxval(q)
!                 stop "MPDATA:NaN error"
!             endif
            call advect_w(q,w,nx,nz,ny)
!             if (any(isnan(q))) then
!                 print*, "post-w"
!                 write(*,*) maxval(q)
!                 stop "MPDATA:NaN error"
!             endif
        elseif (order==1) then
            call advect_v(q,v,nx,nz,ny)
            call advect_w(q,w,nx,nz,ny)
            call advect_u(q,u,nx,nz,ny)
        elseif (order==2) then
            call advect_w(q,w,nx,nz,ny)
            call advect_u(q,u,nx,nz,ny)
            call advect_v(q,v,nx,nz,ny)
        endif
        
        if (options%debug) then
            if (minval(q)<0) then
                write(*,*) minval(q)
                stop "MPDATA:minval error"
            endif
            if (maxval(q)>600) then
                write(*,*) maxval(q)
                stop "MPDATA:maxval error"
            endif
            if (any(isnan(q))) then
                write(*,*) maxval(q)
                stop "MPDATA:NaN error"
            endif
        endif
    end subroutine advect3d
    
    subroutine mpdata_init(domain,options)
        type(domain_type), intent(inout) :: domain
        type(options_type), intent(in) :: options
        
        order=0
    end subroutine mpdata_init
    
!   primary entry point, advect all scalars in domain
    subroutine mpdata(domain,options,dt)
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
        endif
        
!       calculate U,V,W normalized for dt/dx
        if (options%advect_density) then
            U_m=domain%ur(1:nx-1,:,:)*(dt/dx**2)
            V_m=domain%vr(:,:,1:ny-1)*(dt/dx**2)
!           note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
            W_m=domain%wr*(dt/dx**2)
        else
            U_m=domain%u(1:nx-1,:,:)*(dt/dx)
            V_m=domain%v(:,:,1:ny-1)*(dt/dx)
!           note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
            W_m=domain%w*(dt/dx)
        endif

!         print*, "qv"
        call advect3d(domain%qv,   U_m,V_m,W_m, domain%rho,domain%dz,nx,nz,ny,0,options)
!         print*, "qc"
        call advect3d(domain%cloud,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
!         print*, "qr"
        call advect3d(domain%qrain,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
!         print*, "qs"
        call advect3d(domain%qsnow,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
!         print*, "th"
        call advect3d(domain%th,   U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
        if (options%physics%microphysics==kMP_THOMPSON) then
!             print*, "qi"
            call advect3d(domain%ice,  U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
!             print*, "qg"
            call advect3d(domain%qgrau,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
!             print*, "ni"
            call advect3d(domain%nice, U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
!             print*, "nr"
            call advect3d(domain%nrain,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
        endif
!         print*, "CLOUD", minval(domain%cloud), maxval(domain%cloud)
        order=mod(order+1,3)
    end subroutine mpdata
end module adv_mpdata
