!>------------------------------------------------------------
!!  Test MPDATA advection code.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
program test_mpdata
    use adv_mpdata
    use data_structures

    integer, parameter :: kSTEP_FUNCTION=0
    integer, parameter :: kSINE_CURVE=1

    logical :: FCT
    integer :: initial_conditions

    FCT=.True.
!     FCT=.False.

!     initial_conditions=kSINE_CURVE
    initial_conditions=kSTEP_FUNCTION

! NOTE: would do well to update these to all be calls to advect3D, but they have served their purpose
!     print*, "------------------------------"
!     print*, "Testing U"
!     print*, "------------------------------"
!     call test_u(FCT,initial_conditions)
!     print*, "------------------------------"
!     print*, ""
!     print*, "------------------------------"
!     print*, "Testing V"
!     print*, "------------------------------"
!     call test_v(FCT,initial_conditions)
!     print*, "------------------------------"
!     print*, ""
!     print*, "------------------------------"
!     print*, "Testing w"
!     print*, "------------------------------"
!     call test_w(FCT,initial_conditions)
!     print*, "------------------------------"
    print*, ""
    print*, "------------------------------"
    print*, "Testing 3D advection"
    print*, "------------------------------"
    call test_3d(FCT,initial_conditions)
    print*, "------------------------------"

contains

    subroutine test_3d(FCT,initial_conditions)
        implicit none
        logical, intent(in) :: FCT
        integer, intent(in) :: initial_conditions

        type(options_type) :: options
        integer, parameter :: nx=3, ny=100, nz=3
        real, parameter :: cfl=0.25
        real, allocatable, dimension(:,:,:) :: q,rho,dz
        real, allocatable, dimension(:,:,:) :: u,v,w
        integer :: i, loop, nloops, err

        allocate(q(nx,nz,ny))
        allocate(rho(nx,nz,ny))
        allocate(dz(nx,nz,ny))
        allocate(u(nx-1,nz,ny))
        allocate(v(nx,nz,ny-1))
        allocate(w(nx,nz,ny))

        nloops=10
        print*, "CFL=",cfl
        options%adv_options%flux_corrected_transport=FCT
        options%adv_options%mpdata_order=2
        dz=1
        rho=1
        err=0

        if (initial_conditions==kSINE_CURVE) then
            ! setup a sine curve for the initial conditions.
            do i=1,ny
                q(:,:,i)=sin(i/real(ny-2) * 2*3.141592)+1
            end do
        elseif (initial_conditions==kSTEP_FUNCTION) then
            q=1
            do i=1,ny
                if (i>ny/2) then
                    q(:,:,i)=2
                endif
            end do
        endif

        v=cfl
        u=0
        w=0

        q(:,:,1)=q(:,:,ny-1)
        q(:,:,ny)=q(:,:,2)
        write(*, "(A,F6.4)") "Initial mean: ", sum(q(2,1,:))/ny
        write(*, "(A,12F6.3)") "Initial       ", q(2,1,::ny/10), maxval(q(2,1,:)), minval(q(2,1,:))
        do loop=1,nloops
            do i=1,(ny-2)/cfl
                call advect3d(q,u,v,w,rho,dz,nx,nz,ny, options,err)
                ! Wrap around boundary conditions
                q(:,:,1)=q(:,:,ny-1)
                q(:,:,ny)=q(:,:,2)
            end do
            write(*, "(A,I2,12F6.3)") "After loop: ", loop, q(2,1,::ny/10), maxval(q(2,1,:)), minval(q(2,1,:))
        end do
        write(*, "(A,3F7.4)") "Final mean: ", sum(q(2,1,:))/ny, maxval(q(2,1,:)), minval(q(2,1,:))

    end subroutine test_3d


!     subroutine test_v(FCT, initial_conditions)
!         implicit none
!         logical, intent(in) :: FCT
!         integer, intent(in) :: initial_conditions
!
!         integer, parameter :: nx=3, ny=100, nz=3
!         real, parameter :: cfl=0.25
!         real, allocatable, dimension(:,:,:) :: q
!         real, allocatable, dimension(:,:,:) :: u
!         integer :: i, loop, nloops
!
!         allocate(q(nx,nz,ny))
!         allocate(u(nx,nz,ny-1))
!
!         nloops=10
!         print*, "CFL=",cfl
!
!         if (initial_conditions==kSINE_CURVE) then
!             ! setup a sine curve for the initial conditions.
!             do i=1,ny
!                 q(:,:,i)=sin(i/real(ny-2) * 2*3.141592)+1
!             end do
!         elseif (initial_conditions==kSTEP_FUNCTION) then
!             q=1
!             do i=1,ny
!                 if (i>ny/2) then
!                     q(:,:,i)=2
!                 endif
!             end do
!         endif
!
!         u=cfl
!
!         q(:,:,1)=q(:,:,ny-1)
!         q(:,:,ny)=q(:,:,2)
!         write(*, "(A,F6.4)") "Initial mean: ", sum(q(2,1,:))/ny
!         write(*, "(A,12F6.3)") "Initial       ", q(2,1,::ny/10), maxval(q(2,1,:)), minval(q(2,1,:))
!         do loop=1,nloops
!             do i=1,(ny-2)/cfl
!                 call advect_v(q,u,nx,nz,ny, FCT)
!                 ! Wrap around boundary conditions
!                 q(:,:,1)=q(:,:,ny-1)
!                 q(:,:,ny)=q(:,:,2)
!             end do
!             write(*, "(A,I2,12F6.3)") "After loop: ", loop, q(2,1,::ny/10), maxval(q(2,1,:)), minval(q(2,1,:))
!         end do
!         write(*, "(A,3F7.4)") "Final mean: ", sum(q(2,1,:))/ny, maxval(q(2,1,:)), minval(q(2,1,:))
!
!     end subroutine test_v
!
!
!     subroutine test_u(FCT, initial_conditions)
!         implicit none
!         logical, intent(in) :: FCT
!         integer, intent(in) :: initial_conditions
!
!         integer, parameter :: nx=100, ny=3, nz=3
!         real, parameter :: cfl=0.25
!         real, allocatable, dimension(:,:,:) :: q
!         real, allocatable, dimension(:,:,:) :: u
!         integer :: i, loop, nloops
!
!         allocate(q(nx,nz,ny))
!         allocate(u(nx-1,nz,ny))
!
!         nloops=10
!         print*, "CFL=",cfl
!
!         if (initial_conditions==kSINE_CURVE) then
!             ! setup a sine curve for the initial conditions.
!             do i=1,nx
!                 q(i,:,:)=sin(i/real(nx-2) * 2*3.141592)+1
!             end do
!         elseif (initial_conditions==kSTEP_FUNCTION) then
!             q=0
!             do i=1,nx
!                 if (i>nx/2) then
!                     q(i,:,:)=1
!                 endif
!             end do
!         endif
!
!         u=cfl
!
!         write(*, "(A,F6.4)") "Initial mean: ", sum(q(:,1,2))/nx
!         write(*, "(A,12F6.3)") "Initial       ", q(::nx/10,1,2), maxval(q(:,1,2)), minval(q(:,1,2))
!         do loop=1,nloops
!             do i=1,(nx-2)/cfl
!                 call advect_u(q,u,nx,nz,ny, FCT)
!                 ! Wrap around boundary conditions
!                 q(1,:,:)=q(nx-1,:,:)
!                 q(nx,:,:)=q(2,:,:)
!             end do
!             write(*, "(A,I2,12F6.3)") "After loop: ", loop, q(::nx/10,1,2), maxval(q(:,1,2)), minval(q(:,1,2))
!         end do
!         write(*, "(A,3F7.4)") "Final mean: ", sum(q(:,1,2))/nx, maxval(q(:,1,2)), minval(q(:,1,2))
!
!     end subroutine test_u
!
!     subroutine test_w(FCT, initial_conditions)
!         implicit none
!         logical, intent(in) :: FCT
!         integer, intent(in) :: initial_conditions
!
!         integer, parameter :: nx=3, ny=3, nz=106
!         real, parameter :: cfl=0.25
!         real, allocatable, dimension(:,:,:) :: q
!         real, allocatable, dimension(:,:,:) :: u
!         integer :: i, loop, nloops
!
!         allocate(q(nx,nz,ny))
!         allocate(u(nx,nz,ny))
!
!         nloops=10
!         print*, "CFL=",cfl
!
!         if (initial_conditions==kSINE_CURVE) then
!             ! setup a sine curve for the initial conditions.
!             ! offset in space to account for large buffered boundaries for w test
!             do i=1,nz
!                 q(:,i,:)=sin((i-3)/real(nz-7) * 2*3.141592)+1
!             end do
!         elseif (initial_conditions==kSTEP_FUNCTION) then
!             q=1
!             do i=1,nz
!                 if (i>nz/2) then
!                     q(:,i,:)=2
!                 endif
!             end do
!         endif
!
!         u=cfl
!
!         q(:,1,:)=q(:,nz-6,:)
!         q(:,2,:)=q(:,nz-5,:)
!         q(:,3,:)=q(:,nz-4,:)
!         q(:,nz-3,:)=q(:,4,:)
!         q(:,nz-2,:)=q(:,5,:)
!         q(:,nz-1,:)=q(:,6,:)
!         q(:,nz,  :)=q(:,7,:)
!         write(*, "(A,F6.4)") "Initial mean: ", sum(q(2,4:nz-3,2))/(nz-6)
!         write(*, "(A,13F6.3)") "Initial       ", q(2,4::nz/10,2), maxval(q(2,:,2)), minval(q(2,:,2))
!         do loop=1,nloops
!             do i=1,(nz-7)/cfl
!                 call advect_w(q,u,nx,nz,ny, FCT)
!                 ! Wrap around boundary conditions
!                 ! keep a large buffer because we have issues with boundaries in w
!                 ! w has special boundary conditions, so we are only examining internal
!                 ! advection with this test.
!                 q(:,1,:)=q(:,nz-6,:)
!                 q(:,2,:)=q(:,nz-5,:)
!                 q(:,3,:)=q(:,nz-4,:)
!                 q(:,nz-3,:)=q(:,4,:)
!                 q(:,nz-2,:)=q(:,5,:)
!                 q(:,nz-1,:)=q(:,6,:)
!                 q(:,nz,  :)=q(:,7,:)
!             end do
!             write(*, "(A,I2,14F6.3)") "After loop: ", loop, q(2,4::nz/10,2), sum(q(2,4:nz-3,2))/(nz-6), maxval(q(2,:,2)), minval(q(2,:,2))
!         end do
!         write(*, "(A,3F7.4)") "Final mean: ", sum(q(2,4:nz-3,2))/(nz-6), maxval(q(2,:,2)), minval(q(2,:,2))
!
!
!     end subroutine test_w
!
end program
