program test_mpdata
    use adv_mpdata
    
    print*, "------------------------------"
    print*, "Testing U"
    print*, "------------------------------"
    call test_u()
    print*, "------------------------------"
    print*, ""
    print*, "------------------------------"
    print*, "Testing V"
    print*, "------------------------------"
    call test_v()
    print*, "------------------------------"
    print*, ""
    print*, "------------------------------"
    print*, "Testing w"
    print*, "------------------------------"
    call test_w()
    print*, "------------------------------"
    
contains
    subroutine test_v()
        implicit none
        
        integer, parameter :: nx=3, ny=100, nz=3
        real, parameter :: cfl=0.25
        real, allocatable, dimension(:,:,:) :: q
        real, allocatable, dimension(:,:,:) :: u
        integer :: i, loop, nloops
        
        allocate(q(nx,nz,ny))
        allocate(u(nx,nz,ny-1))
        
        nloops=10
        print*, "CFL=",cfl
        
        ! setup a sine curve for the initial conditions.  
        do i=1,ny
            q(:,:,i)=sin(i/real(ny-2) * 2*3.141592)+1
        end do
        u=cfl
        
        write(*, "(A,F6.4)") "Initial mean: ", sum(q(2,1,:))/ny
        write(*, "(A,10F6.3)") "Initial       ", q(2,1,::ny/10)
        do loop=1,nloops
            do i=1,(ny-2)/cfl
                call advect_v(q,u,nx,nz,ny)
                ! Wrap around boundary conditions
                q(:,:,1)=q(:,:,ny-1)
                q(:,:,ny)=q(:,:,2)
            end do
            write(*, "(A,I2,10F6.3)") "After loop: ", loop, q(2,1,::ny/10)
        end do
        write(*, "(A,F6.4)") "Final mean: ", sum(q(2,1,:))/ny
        
    end subroutine test_v


    subroutine test_u()
        implicit none
        
        integer, parameter :: nx=100, ny=3, nz=3
        real, parameter :: cfl=0.25
        real, allocatable, dimension(:,:,:) :: q
        real, allocatable, dimension(:,:,:) :: u
        integer :: i, loop, nloops
        
        allocate(q(nx,nz,ny))
        allocate(u(nx-1,nz,ny))
        
        nloops=10
        print*, "CFL=",cfl
        
        ! setup a sine curve for the initial conditions.  
        do i=1,nx
            q(i,:,:)=sin(i/real(nx-2) * 2*3.141592)+1
        end do
        u=cfl
        
        write(*, "(A,F6.4)") "Initial mean: ", sum(q(:,1,2))/nx
        write(*, "(A,10F6.3)") "Initial       ", q(::nx/10,1,2)
        do loop=1,nloops
            do i=1,(nx-2)/cfl
                call advect_u(q,u,nx,nz,ny)
                ! Wrap around boundary conditions
                q(1,:,:)=q(nx-1,:,:)
                q(nx,:,:)=q(2,:,:)
            end do
            write(*, "(A,I2,10F6.3)") "After loop: ", loop, q(::nx/10,1,2)
        end do
        write(*, "(A,F6.4)") "Final mean: ", sum(q(:,1,2))/nx
        
    end subroutine test_u
    
    subroutine test_w()
        implicit none
        
        integer, parameter :: nx=3, ny=3, nz=106
        real, parameter :: cfl=0.25
        real, allocatable, dimension(:,:,:) :: q
        real, allocatable, dimension(:,:,:) :: u
        integer :: i, loop, nloops
        
        allocate(q(nx,nz,ny))
        allocate(u(nx,nz,ny))
        
        nloops=10
        print*, "CFL=",cfl
        
        ! setup a sine curve for the initial conditions.  
        ! offset in space to account for large buffered boundaries for w test
        do i=1,nz
            q(:,i,:)=sin((i-3)/real(nz-7) * 2*3.141592)+1
        end do
        u=cfl
        
        write(*, "(A,F6.4)") "Initial mean: ", sum(q(2,4:nz-3,2))/(nz-6)
        write(*, "(A,11F6.3)") "Initial       ", q(2,4::nz/10,2)
        write(*, "(10F6.3)") q(2,:10,2)
        do loop=1,nloops
            do i=1,(nz-7)/cfl
                call advect_w(q,u,nx,nz,ny)
                ! Wrap around boundary conditions
                ! keep a large buffer because we have issues with boundaries in w
                ! w has special boundary conditions, so we are only examining internal
                ! advection with this test. 
                q(:,1,:)=q(:,nz-6,:)
                q(:,2,:)=q(:,nz-5,:)
                q(:,3,:)=q(:,nz-4,:)
                q(:,nz-3,:)=q(:,4,:)
                q(:,nz-2,:)=q(:,5,:)
                q(:,nz-1,:)=q(:,6,:)
                q(:,nz,  :)=q(:,7,:)
            end do
            write(*, "(A,I2,12F6.3)") "After loop: ", loop, q(2,4::nz/10,2), sum(q(2,4:nz-3,2))/(nz-6)
        end do
        write(*, "(A,F6.4)") "Final mean: ", sum(q(2,4:nz-3,2))/(nz-6)
        
    end subroutine test_w
    
end program