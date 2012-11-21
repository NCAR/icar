module advect
    
    implicit none
    contains
    subroutine flux2(l,r,U,ny,nz,nx,f)
    !     Calculate the donor cell flux function
    !     l = left gridcell scalar 
    !     r = right gridcell scalar
    !     U = Courant number (u*dt/dx)
    !     
    !     If U is positive, return l*U if U is negative return r*U
    !     By using the mathematical form instead of the logical form, 
    !     we can run on the entire grid simultaneously, and avoid branches

    !   arguments
        real, dimension(1:ny,1:nz,1:nx), intent(in) :: l,r,U
        real, dimension(1:ny,1:nz,1:nx), intent(inout) :: f
        integer,intent(in) :: ny,nz,nx
        !   internal parameter
        integer ::  err,i!,j,Ny,Nz,Nx
        !   main code
        f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2

    end subroutine flux2

    subroutine advect3d(q,u,v,w,ny,nz,nx)
	    real,dimension(1:ny,1:nz,1:nx), intent(inout) :: q
	    real,dimension(1:ny,1:nz,1:nx), intent(in) :: w
        real,dimension(1:ny,1:nz,1:nx-1),intent(in) :: u
        real,dimension(1:ny-1,1:nz,1:nx),intent(in) :: v
		integer, intent(in) :: ny,nz,nx
        ! interal parameters
        integer :: err,i
        real, dimension(1:ny-2,1:nz,1) :: f1,f2,f3,f4
        real, dimension(1:ny-2,1:nz-2,1) ::f5,f6
        !$omp parallel shared(q,u,v,w) firstprivate(nx,ny,nz) private(i,f1,f2,f3,f4,f5,f6)
        !$omp do
        do i=2,nx-1
           call flux2(q(2:ny-1,:,i),q(2:ny-1,:,i+1),u(2:ny-1,:,i-1),ny-2,nz,1,f1)  !Ux1
           call flux2(q(2:ny-1,:,i-1),q(2:ny-1,:,i),u(2:ny-1,:,i),ny-2,nz,1,f2)  !Ux0
           call flux2(q(2:ny-1,:,i),q(3:ny,:,i),v(1:ny-2,:,i),ny-2,nz,1,f3)  !Vy1
           call flux2(q(1:ny-2,:,i),q(2:ny-1,:,i),v(2:ny-1,:,i),ny-2,nz,1,f4)  !Vy0
           call flux2(q(2:ny-1,2:nz-1,i),q(2:ny-1,3:nz,i),w(2:ny-1,2:nz-1,i),ny-2,nz-2,1,f5)
           call flux2(q(2:ny-1,1:nz-2,i),q(2:ny-1,2:nz-1,i),w(2:ny-1,1:nz-2,i),ny-2,nz-2,1,f6)

           ! perform horizontal advection
           q(2:ny-1,:,i)=q(2:ny-1,:,i) - (f1(:,:,1)-f2(:,:,1)) - (f3(:,:,1)-f4(:,:,1))
           ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
           ! add fluxes to middle layers
           q(2:ny-1,2:nz-1,i)=q(2:ny-1,2:nz-1,i)-(f5(:,:,1)-f6(:,:,1))
           ! add fluxes to bottom layer
           q(2:ny-1,1,i)=q(2:ny-1,1,i)-f6(:,1,1)
           ! add fluxes to top layer
           q(2:ny-1,nz,i)=q(2:ny-1,nz,i)-(q(2:ny-1,nz,i)*w(2:ny-1,nz,i)-f5(:,nz-2,1))
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine advect3d
	
end module advect
