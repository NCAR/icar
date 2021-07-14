!>------------------------------------------------------------
!! Simple test program to verify that fftw libraries are properly linked
!!
!!  @author
!!  Ethan Gutmann
!!
!!------------------------------------------------------------
program fftw_test
    use, intrinsic :: iso_c_binding
    use fft

    implicit none

    integer*8 :: n_elements, nx, ny, nz

    type(C_PTR) :: fftw_aligned_data
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: fftw_fdata

    nx = 100
    ny = 100
    nz = 20

    n_elements = nx * ny * nz
    print*, "Testing fftw_alloc_complex"
    fftw_aligned_data = fftw_alloc_complex(n_elements)
    call c_f_pointer(fftw_aligned_data, fftw_fdata, [nx,ny,nz])

    fftw_fdata(1,1,1) = 10
    ! print*, fftw_fdata(1,1,1)
    if (real(fftw_fdata(1,1,1)) == 10) then
        print*, "PASSED"
    endif

end program fftw_test
