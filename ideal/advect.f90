module advection
	use data_structures
    
    implicit none
	private
	public::advect
	
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
        real, dimension(1:nx,1:nz,1:ny), intent(in) :: l,r,U
        real, dimension(1:nx,1:nz,1:ny), intent(inout) :: f
        integer,intent(in) :: ny,nz,nx
        !   internal parameter
        integer ::  err,i!,j,Ny,Nz,Nx
        !   main code
        f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2

    end subroutine flux2

    subroutine advect3d(q,u,v,w,nx,nz,ny,debug)
	    real,dimension(1:nx,1:nz,1:ny), intent(inout) :: q
	    real,dimension(1:nx,1:nz,1:ny), intent(in) :: w
        real,dimension(1:nx-1,1:nz,1:ny),intent(in) :: u
        real,dimension(1:nx,1:nz,1:ny-1),intent(in) :: v
		integer, intent(in) :: ny,nz,nx,debug
        ! interal parameters
        integer :: err,i
	    real,dimension(1:nx,1:nz,1:ny) :: qin
        real, dimension(1:nx-2,1:nz,1) :: f1,f2,f3,f4
        real, dimension(1:nx-2,1:nz-2,1) ::f5,f6
        !$omp parallel shared(q,u,v,w) firstprivate(nx,ny,nz) private(i,f1,f2,f3,f4,f5,f6)
        !$omp do
        do i=1,ny
			qin(:,:,i)=q(:,:,i)
		enddo
        !$omp end do
        !$omp do
        do i=2,ny-1
! 			note f1,f2 are largely duplicates (ditto f5,f6 and f3,f4 between loops (e.g. f3=previous f4))
! 			this is just a convenient notation for now. 
! 			ideally this redundancy should be used to decrease the number of flux2 calls (cut in 1/2)
           call flux2(qin(2:nx-1,:,i),qin(3:nx,:,i),u(2:nx-1,:,i),nx-2,nz,1,f1)  !Ux1
           call flux2(qin(1:nx-2,:,i),qin(2:nx-1,:,i),u(1:nx-2,:,i),nx-2,nz,1,f2)  !Ux0
           call flux2(qin(2:nx-1,:,i),qin(2:nx-1,:,i+1),v(2:nx-1,:,i),nx-2,nz,1,f3)  !Vy1
           call flux2(qin(2:nx-1,:,i-1),qin(2:nx-1,:,i),v(2:nx-1,:,i-1),nx-2,nz,1,f4)  !Vy0
           call flux2(qin(2:nx-1,2:nz-1,i),qin(2:nx-1,3:nz,i),w(2:nx-1,2:nz-1,i),nx-2,nz-2,1,f5)
           call flux2(qin(2:nx-1,1:nz-2,i),qin(2:nx-1,2:nz-1,i),w(2:nx-1,1:nz-2,i),nx-2,nz-2,1,f6)
		   
           ! perform horizontal advection
           q(2:nx-1,:,i)=q(2:nx-1,:,i) - (f1(:,:,1)-f2(:,:,1)) - (f3(:,:,1)-f4(:,:,1))
           ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
           ! add fluxes to middle layers
           q(2:nx-1,2:nz-1,i)=q(2:nx-1,2:nz-1,i)-(f5(:,:,1)-f6(:,:,1))
           ! add fluxes to bottom layer
           q(2:nx-1,1,i)=q(2:nx-1,1,i)-f6(:,1,1)
           ! add fluxes to top layer
           q(2:nx-1,nz,i)=q(2:nx-1,nz,i)-(qin(2:nx-1,nz,i)*w(2:nx-1,nz,i)-f5(:,nz-2,1))
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine advect3d

! 	primary entry point, advect all scalars in domain
	subroutine advect(domain,options,dt)
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		real,intent(in)::dt
		
		real::dx
		integer::nx,nz,ny
		real,dimension(:,:,:),allocatable::U,V,W
			
		dx=domain%dx
		nx=size(domain%dz,1)
		nz=size(domain%dz,2)
		ny=size(domain%dz,3)
		
! 		calculate U,V,W normalized for dt/dx
		allocate(U(nx-1,nz,ny))
		U=domain%u(1:nx-1,:,:)*dt/dx
		allocate(V(nx,nz,ny-1))
		V=domain%v(:,:,1:ny-1)*dt/dx
		allocate(W(nx,nz,ny))
! 		note, even though dz!=dx, W is computed from the divergence in U/V so it is scaled by dx/dz already
		W=domain%w*dt/dx
! 		should probably be converting to mass (q*rho) before advecting, then back again... but testing showed minimal difference
		call advect3d(domain%th,   U,V,W,nx,nz,ny,0)
		call advect3d(domain%qv,   U,V,W,nx,nz,ny,0)
		call advect3d(domain%cloud,U,V,W,nx,nz,ny,0)
		call advect3d(domain%ice,  U,V,W,nx,nz,ny,0)
		call advect3d(domain%nice, U,V,W,nx,nz,ny,0)
		call advect3d(domain%qrain,U,V,W,nx,nz,ny,0)
		call advect3d(domain%nrain,U,V,W,nx,nz,ny,0)
		call advect3d(domain%qsnow,U,V,W,nx,nz,ny,0)
		call advect3d(domain%qgrau,U,V,W,nx,nz,ny,0)
		
		deallocate(U,V,W)
		
	end subroutine advect


    subroutine flux2d(l,r,U,ny,nz,f)
    !     Calculate the donor cell flux function
    !     l = left gridcell scalar 
    !     r = right gridcell scalar
    !     U = Courant number (u*dt/dx)
    !     
    !     If U is positive, return l*U if U is negative return r*U
    !     By using the mathematical form instead of the logical form, 
    !     we can run on the entire grid simultaneously, and avoid branches

    !   arguments
        real, dimension(1:ny,1:nz), intent(in) :: l,r,U
        real, dimension(1:ny,1:nz), intent(inout) :: f
        integer,intent(in) :: ny,nz
        !   internal parameter
!         integer ::  err,i!,j,Ny,Nz,Nx
        !   main code
        f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2

    end subroutine flux2d

    subroutine advect2d(q,u,w,nx,nz)
	    real,dimension(nx,nz), intent(inout) :: q
        real,dimension(nx,nz),intent(in) :: u
	    real,dimension(nx,nz), intent(in) :: w
        ! interal parameters
		integer,intent(in) :: nx,nz
        integer :: err,i
	    real,dimension (nx,nz) :: qin
        real, dimension(nx-1) :: f1
        real, dimension(nx) :: f2
        !$omp parallel shared(q,u,w) firstprivate(nx,nz) private(i,f1,f2)
        !$omp do
        do i=1,nz
			qin(:,i)=q(:,i)
		enddo
        !$omp end do
        !$omp do
        do i=1,nz-1
! 			write(*,*) i
           call flux2d(qin(1:nx-1,i),qin(2:nx,i),u(2:nx,i),nx-1,1,f1)
           call flux2d(qin(1:nx,i),qin(1:nx,i+1),w(1:nx,i),nx,1,f2)

           ! perform horizontal advection
		   q(2:nx-1,i)=q(2:nx-1,i)+(f1(1:nx-2)-f1(2:nx-1))
		   ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
           ! add fluxes to middle layers
		   q(:,i)=q(:,i)-f2
		   if (i<nz-1) then
			   q(:,i+1)=q(:,i+1)+f2
		   endif
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine advect2d

! 	subroutine advect2d_driver(domain,options,dt,dx)
! 		type(domain_type),intent(inout)::domain
! 		type(options_type),intent(in)::options
! 		real,intent(in)::dt,dx
! 		
! 		integer::nx,nz
! 		real,dimension(:,:),allocatable::U,W
! 		nx=size(domain%dz,1)
! 		nz=size(domain%dz,2)
! 		
! 		allocate(U(nx,nz))
! 		U=domain%u*dt/dx
! 		allocate(W(nx,nz))
! 		W=domain%w*dt/domain%dz
! 		call advect2d(domain%th,   U,W,nx,nz)
! 		call advect2d(domain%qv,   U,W,nx,nz)
! 		call advect2d(domain%cloud,U,W,nx,nz)
! 		call advect2d(domain%ice,  U,W,nx,nz)
! 		call advect2d(domain%nice, U,W,nx,nz)
! 		call advect2d(domain%qrain,U,W,nx,nz)
! 		call advect2d(domain%nrain,U,W,nx,nz)
! 		call advect2d(domain%qsnow,U,W,nx,nz)
! 		call advect2d(domain%qgrau,U,W,nx,nz)
! 		
! 	end subroutine advect2d_driver
	
end module advection
