!>------------------------------------------------------------
!!  Simple unit test to confirm that the fft_shift module works
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
program test_fftshift
    use fftshifter
    implicit none

    integer,parameter::n=5

    real,dimension(:),allocatable::test1r
    real,dimension(:,:),allocatable::test2r
    complex,dimension(:),allocatable::test1c
    complex,dimension(:,:),allocatable::test2c
    integer::i,j

    allocate(test1r(n))
    allocate(test2r(n,n))
    allocate(test1c(n))
    allocate(test2c(n,n))
    do i=1,n
        test1r(i)=i
        test1c(i)=i
        do j=1,n
            test2r(i,j)=i+j*n*10
            test2c(i,j)=i+j*n*10
        enddo
    enddo
    print*, "1D real array"
    print*, " Original:"
    print*, aint(test1r)
    call fftshift(test1r)
    print*, " Shifted:"
    print*, aint(test1r)
    call ifftshift(test1r)
    print*, " Reverted:"
    print*, aint(test1r)

    print*, "2D real array"
    print*, " Original:"
    print*, aint(test2r)
    call fftshift(test2r)
    print*, " Shifted:"
    print*, aint(test2r)
    call ifftshift(test2r)
    print*, " Reverted:"
    print*, aint(test2r)

    print*, "1D complex array"
    print*, " Original:"
    print*, aint(real(test1c))
    call fftshift(test1c)
    print*, " Shifted:"
    print*, aint(real(test1c))
    call ifftshift(test1c)
    print*, " Reverted:"
    print*, aint(real(test1c))

    print*, "2D complex array"
    print*, " Original:"
    print*, aint(real(test2c))
    call fftshift(test2c)
    print*, " Shifted:"
    print*, aint(real(test2c))
    call ifftshift(test2c)
    print*, " Reverted:"
    print*, aint(real(test2c))


end program test_fftshift
