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
        
        integer, parameter :: nx=3, ny=20, nz=3
        real, allocatable, dimension(:,:,:) :: q
        real, allocatable, dimension(:,:,:) :: v
        integer :: i,j
        
        allocate(q(nx,nz,ny))
        allocate(v(nx,nz,ny-1))
        
        do i=1,ny
            q(:,:,i)=abs(i-3)
        end do
        v=0.25
        
        print*, "Testing individual steps CFL=",v(1,1,1)
        print*, "Initial ",q(2,1,::2)
        do i=1,10
            call advect_v(q,v,nx,nz,ny)
            print*, i, "  ",q(2,1,::2)
        end do
        print*, "Final Boundary ",q(1,1,::2)
        print*, "Final Internal ",q(2,1,::2)

        ! Reset domain
        do i=1,ny
            q(:,:,i)=abs(i-ny/2)
        end do
        ! Wrap around boundary conditions
        q(:,:,1)=q(:,:,ny-1)
        q(:,:,ny)=q(:,:,2)

        print*, "Testing complete cycles over domain CFL=",v(1,1,1)
        print*, "Initial ",q(2,1,::2)
        
        do j=1,10
            do i=1, nint(1/v(1,1,1)) * ny
                call advect_v(q,v,nx,nz,ny)
                ! Wrap around boundary conditions
                q(:,:,1)=q(:,:,ny-1)
                q(:,:,ny)=q(:,:,2)
            end do
            print*, j, q(2,1,::2)
        end do

        
    end subroutine test_v


    subroutine test_u()
        implicit none
        
        integer, parameter :: nx=10, ny=3, nz=3
        real, allocatable, dimension(:,:,:) :: q
        real, allocatable, dimension(:,:,:) :: u
        integer :: i
        
        allocate(q(nx,nz,ny))
        allocate(u(nx-1,nz,ny))
        
        do i=1,nx
            q(i,:,:)=abs(i-3)
        end do
        u=0.25
        
        print*, "Initial ",q(:,1,2)
        do i=1,10
            call advect_u(q,u,nx,nz,ny)
            print*, i, " ",q(:,1,2)
        end do
        print*, "Final Boundary ",q(:,1,1)
        print*, "Final Internal ",q(:,1,2)
        
    end subroutine test_u
    
    subroutine test_w()
        implicit none
        
        integer, parameter :: nx=3, ny=3, nz=10
        real, allocatable, dimension(:,:,:) :: q
        real, allocatable, dimension(:,:,:) :: u
        integer :: i
        
        allocate(q(nx,nz,ny))
        allocate(u(nx,nz,ny))
        
        do i=1,nz
            q(:,i,:)=abs(i-3)
        end do
        u=0.25
        
        print*, "Initial ",q(2,:,2)
        do i=1,10
            call advect_w(q,u,nx,nz,ny)
            print*, i, " ",q(2,:,2)
        end do
        print*, "Final Boundary ",q(1,:,1)
        print*, "Final Internal ",q(2,:,2)
        
    end subroutine test_w
    
end program