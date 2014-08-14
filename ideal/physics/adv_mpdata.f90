module adv_mpdata
	use data_structures
    
    implicit none
	private
	real,dimension(:,:,:),allocatable::U_m,V_m,W_m
	public::mpdata
	
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
!     # Al=(q1[1:-1,:-2]-q1[1:-1,1:-1])/(q1[1:-1,:-2]+q1[1:-1,1:-1])
!     # Ar=(q1[1:-1,2:]-q1[1:-1,1:-1])/(q1[1:-1,2:]+q1[1:-1,1:-1])
!     # Au=(q1[:-2,1:-1]-q1[1:-1,1:-1])/(q1[:-2,1:-1]+q1[1:-1,1:-1])
!     # Ab=(q1[2:,1:-1]-q1[1:-1,1:-1])/(q1[2:,1:-1]+q1[1:-1,1:-1])
!     # 
!     # q11=q1[2:,2:]+q1[2:,1:-1]
!     # Br=0.5*((q11-q1[:-2,2:]-q1[:-2,1:-1])
!     #     /(q11+q1[:-2,2:]+q1[:-2,1:-1]))
!     # q11=q1[2:,:-2]+q1[2:,1:-1]
!     # Bl=0.5*((q11-q1[:-2,:-2]-q1[:-2,1:-1])
!     #     /(q11+q1[:-2,:-2]+q1[:-2,1:-1]))
!     # q11=q1[:-2,2:]+q1[1:-1,2:]
!     # Bu=0.5*((q11-q1[:-2,:-2]-q1[1:-1,:-2])
!     #     /(q11+q1[:-2,:-2]+q1[1:-1,:-2]))
!     # q11=q1[2:,2:]+q1[1:-1,2:]
!     # Bb=0.5*((q11-q1[2:,:-2]-q1[1:-1,:-2])
!     #     /(q11+q1[2:,:-2]+q1[1:-1,:-2]))
!     # 
!     # # compute diffusion correction U/V terms (see MPDATA review reference)
!     # Uabs=np.abs(U)
!     # # first find U/V terms on the grid cell borders
!     # curUabs=(Uabs[1:-1,:-2]+Uabs[1:-1,1:-1])/2
!     # curU=Ux0
!     # curV=(V[1:-1,:-2]+V[1:-1,1:-1])/2
!     # # then compute Ul
!     # Ul=curUabs*(1-curUabs)*Al - 2*f*curU*curV*Bl
!     # # compute Ur using the same two steps
!     # curUabs=(Uabs[1:-1,2:]+Uabs[1:-1,1:-1])/2
!     # curU=Ux1
!     # curV=(V[1:-1,2:]+V[1:-1,1:-1])/2
!     # Ur=curUabs*(1-curUabs)*Ar - 2*f*curU*curV*Br
!     # # compute Vu
!     # Vabs=np.abs(V)
!     # curVabs=(Vabs[:-2,1:-1]+Vabs[1:-1,1:-1])/2
!     # curV=Vy0
!     # curU=(U[:-2,1:-1]+U[1:-1,1:-1])/2
!     # Vu=curVabs*(1-curVabs)*Bu - 2*f*curU*curV*Bu
!     # # compute Vb
!     # curVabs=(Vabs[2:,1:-1]+Vabs[1:-1,1:-1])/2
!     # curV=Vy1
!     # curU=(U[2:,1:-1]+U[1:-1,1:-1])/2
!     # Vb=curVabs*(1-curVabs)*Bb - 2*f*curU*curV*Bb
!     # 
!     # q[1:-1,1:-1]=(q1[1:-1,1:-1]
!     #     -(F(q1[1:-1,1:-1],q1[1:-1,2:],Ul)
!     #     -F(q1[1:-1,:-2],q1[1:-1,1:-1],Ur))
!     #     -(F(q1[1:-1,1:-1],q1[2:,1:-1],Vu)
!     #     -F(q1[:-2,1:-1],q1[1:-1,1:-1],Vb)))
	end subroutine advect2d


	subroutine advect1d(q,u,rho,dz,nx,nz,ny,options)
		real,dimension(:,:,:), intent(inout) :: q,u,rho,dz
		integer,intent(in)::nx,nz,ny
		type(options_type),intent(in) :: options
		real,dimension(nx,nz,ny)::q1
		real,dimension(nx,nz-1,ny)::U2,fluxes
		integer::i
		
		q1=q

!       This is the MPDATA diffusion correction term for 1D flow
! 		U2 is the diffusion correction pseudo-velocity
		U2=abs(u(:,1:nz-1,:))-u(:,1:nz-1,:)**2
! 		note U2(nz)=0 because we assume q(nz+1)=q(nz), thus diffusion = 0
		U2= U2 * (q1(:,2:nz,:)-q1(:,1:nz-1,:))/(q1(:,1:nz-1,:)+q1(:,2:nz,:))
! 	    # Vl=(U2[:-2]+U2[1:-1])/2*(q1[:-2]-q1[1:-1])/(q1[:-2]+q1[1:-1])
! 		fluxout=vr,fluxin=vl
!         - (F(q1[1:-1],q1[2:],Vr) - F(q1[1:-1],q1[:-2],Vl)))

! 	    # Vr=(U2[2:]+U2[1:-1])/2*(q1[2:]-q1[1:-1])/(q1[2:]+q1[1:-1])
		
! 		fluxes= ((U2+ABS(U2)) * q(:,1:nz-1,:)  +  (Ur-ABS(Ur)) * q(:,2:nz,:))/2

! 		q(2:nx-1,2:nz-1,2:ny-1)=q1(2:nx-1,2:nz-1,2:ny-1) &
! 				- (fluxes(2:nx-1,2:nz-1,2:ny-1) - fluxes(2:nx-1,1:nz-2,2:ny-1))
! 		
! 		q(2:nx-1,1,2:ny-1)=q(2:nx-1,1,2:ny-1) - fluxes(2:nx-1,1,2:ny-1)
	end subroutine advect1d


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
! 			by manually inlining the flux2 call we should remove extra array copies that the compiler didn't seem to be able to remove. 
! 			equivalent flux2 calls are left in for reference (commented) to restore recall that fx arrays should be 3D : n x m x 1
! 			calculate fluxes between grid cells
           f1= ((u(1:nx-1,:,i)+ABS(u(1:nx-1,:,i))) * qin(1:nx-1,:,i) + &
		    (u(1:nx-1,:,i)-ABS(u(1:nx-1,:,i))) * qin(2:nx,:,i))/2
!            call flux2(qin(1:nx-1,:,i),qin(2:nx,:,i),u(1:nx-1,:,i),nx-1,nz,1,f1)  !Ux1
           f3= ((v(2:nx-1,:,i)+ABS(v(2:nx-1,:,i))) * qin(2:nx-1,:,i) + &
		    (v(2:nx-1,:,i)-ABS(v(2:nx-1,:,i))) * qin(2:nx-1,:,i+1))/2
! 		   call flux2(qin(2:nx-1,:,i),qin(2:nx-1,:,i+1),v(2:nx-1,:,i),nx-2,nz,1,f3)  !Vy1
		   f4= ((v(2:nx-1,:,i-1)+ABS(v(2:nx-1,:,i-1))) * qin(2:nx-1,:,i-1) + &
		    (v(2:nx-1,:,i-1)-ABS(v(2:nx-1,:,i-1))) * qin(2:nx-1,:,i))/2
!            call flux2(qin(2:nx-1,:,i-1),qin(2:nx-1,:,i),v(2:nx-1,:,i-1),nx-2,nz,1,f4)  !Vy0
		   f5= ((w(2:nx-1,1:nz-1,i)+ABS(w(2:nx-1,1:nz-1,i))) * qin(2:nx-1,1:nz-1,i) + &
		    (w(2:nx-1,1:nz-1,i)-ABS(w(2:nx-1,1:nz-1,i))) * qin(2:nx-1,2:nz,i))/2
!            call flux2(qin(2:nx-1,1:nz-1,i),qin(2:nx-1,2:nz,i),w(2:nx-1,1:nz-1,i),nx-2,nz-1,1,f5)
		   
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

! 	primary entry point, advect all scalars in domain
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
		
! 		if this if the first time we are called, we need to allocate the module level arrays
		if (.not.allocated(U_m)) then
			allocate(U_m(nx-1,nz,ny))
			allocate(V_m(nx,nz,ny-1))
			allocate(W_m(nx,nz,ny))
		endif
		
! 		calculate U,V,W normalized for dt/dx
		if (options%advect_density) then
			U_m=domain%ur(1:nx-1,:,:)*(dt/dx**2)
			V_m=domain%vr(:,:,1:ny-1)*(dt/dx**2)
! 			note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
			W_m=domain%wr*(dt/dx**2)
		else
			U_m=domain%u(1:nx-1,:,:)*(dt/dx)
			V_m=domain%v(:,:,1:ny-1)*(dt/dx)
! 			note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
			W_m=domain%w*(dt/dx)
		endif

		call advect3d(domain%qv,   U_m,V_m,W_m,	domain%rho,domain%dz,nx,nz,ny,0,options)
		!MPDATA advection correction
! 		call advect2d(domain%qv,   U_m,V_m,		domain%rho,domain%dz,nx,nz,ny,options)
! 		call advect1d(domain%qv,   W_m,			domain%rho,domain%dz,nx,nz,ny,options)
		
		call advect3d(domain%cloud,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
		call advect3d(domain%qrain,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
		call advect3d(domain%qsnow,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
		call advect3d(domain%th,   U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
		if (options%physics%microphysics==1) then
			call advect3d(domain%ice,  U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
			call advect3d(domain%qgrau,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
			call advect3d(domain%nice, U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
			call advect3d(domain%nrain,U_m,V_m,W_m,domain%rho,domain%dz,nx,nz,ny,0,options)
		endif
	end subroutine mpdata
end module adv_mpdata
