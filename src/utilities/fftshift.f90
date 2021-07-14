!>------------------------------------------------------------
!! Supplies 1D and 2D FFT shift procedures ala matlab fftshift
!!
!! Uses a generic interface so that procedures can be called with
!! any variation of complex, real, 1D or 2D arrays
!! 2D can also be called with C_DOUBLE_COMPLEX variables
!!
!! fftshift swaps the left and right sides of an array.  For a
!!  2D array, the top and bottom halves are also swapped.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module fftshifter
    use fft
    implicit none
    !> Generic interface to 1D, 2D real, complex, and C-complex fftshift routines
    interface fftshift
        module procedure fftshift2cc, fftshift2c, fftshift2r, fftshift1c, fftshift1r
    end interface
    !> Generic interface to 1D, 2D real, complex, and C-complex inverse fftshift routines
    interface ifftshift
        module procedure ifftshift2cc_z, ifftshift2cc, ifftshift2c, ifftshift2r, ifftshift1c, ifftshift1r
    end interface
contains
    subroutine fftshift1c(fftimage)
        implicit none
        complex, intent(inout),dimension(:):: fftimage
        integer::nx
        complex,allocatable,dimension(:) :: tmp
        integer::i,ii

        nx=size(fftimage)
        allocate(tmp(nx))
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(ii) = fftimage(i)
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift1c
    subroutine fftshift1r(fftimage)
        implicit none
        real, intent(inout),dimension(:):: fftimage
        integer::nx
        real,allocatable,dimension(:) :: tmp
        integer::i,ii

        nx=size(fftimage)
        allocate(tmp(nx))
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(ii) = fftimage(i)
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift1r
    subroutine fftshift2cc_z(fftimage,fixed_axis)
        implicit none
        complex(C_DOUBLE_COMPLEX), intent(inout),dimension(:,:,:):: fftimage
        integer, intent(in) :: fixed_axis
        integer::nx,ny
        complex, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i,j,fixed_axis)
            enddo
        enddo
        fftimage(:,:,fixed_axis) = tmp
        deallocate(tmp)
    end subroutine fftshift2cc_z
    subroutine fftshift2cc(fftimage)
        implicit none
        complex(C_DOUBLE_COMPLEX), intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        complex, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i,j)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift2cc
    subroutine fftshift2c(fftimage)
        implicit none
        complex, intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        complex, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i,j)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift2c
    subroutine fftshift2r(fftimage)
        implicit none
        real, intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        real, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i,j)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift2r

    subroutine ifftshift1c(fftimage)
        implicit none
        complex, intent(inout),dimension(:):: fftimage
        integer::nx
        complex,allocatable,dimension(:) :: tmp
        integer::i,ii

        nx=size(fftimage)
        allocate(tmp(nx))
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(i) = fftimage(ii)
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift1c
    subroutine ifftshift1r(fftimage)
        implicit none
        real, intent(inout),dimension(:):: fftimage
        integer::nx
        real,allocatable,dimension(:) :: tmp
        integer::i,ii

        nx=size(fftimage)
        allocate(tmp(nx))
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(i) = fftimage(ii)
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift1r
    subroutine ifftshift2cc_z(fftimage, fixed_axis)
        implicit none
        complex(C_DOUBLE_COMPLEX), intent(inout),dimension(:,:,:):: fftimage
        integer, intent(in) :: fixed_axis
        integer::nx,ny
        complex, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(i,j) = fftimage(ii,jj,fixed_axis)
            enddo
        enddo
        fftimage(:,:,fixed_axis) = tmp
        deallocate(tmp)
    end subroutine ifftshift2cc_z
    subroutine ifftshift2cc(fftimage)
        implicit none
        complex(C_DOUBLE_COMPLEX), intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        complex, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(i,j) = fftimage(ii,jj)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift2cc
    subroutine ifftshift2c(fftimage)
        implicit none
        complex, intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        complex, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(i,j) = fftimage(ii,jj)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift2c
    subroutine ifftshift2r(fftimage)
        implicit none
        real, intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        real, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(i,j) = fftimage(ii,jj)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift2r


end module fftshifter
