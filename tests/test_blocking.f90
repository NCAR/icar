program test_blocking

    use linear_theory_winds
    use data_structures
    use icar_constants
    use mod_blocking,    only : compute_fft_topography, compute_blocked_flow_for_wind
    use io_routines,     only : io_write

    use iso_c_binding

    implicit none
    integer, parameter :: nx   = 200
    integer, parameter :: ny   = 200
    integer, parameter :: nz   = 15
    real    :: zmax = 1000
    real    :: dx   = 2000
    real    :: u    = 10
    real    :: v    = 0
    real, allocatable :: terrain(:,:)
    real :: z_layers(nz)
    type(linear_theory_type) :: lt_data
    integer :: i
    real, allocatable :: wfield(:,:,:)
    real, allocatable :: ufield(:,:,:)
    real, allocatable :: vfield(:,:,:)

    complex(C_DOUBLE_COMPLEX), allocatable  :: complex_terrain(:,:)
    complex(C_DOUBLE_COMPLEX), allocatable  :: fft_terrain(:,:)

    allocate(wfield(nx,nz,ny))
    wfield = 0
    allocate(ufield(nx,nz,ny))
    ufield = 0
    allocate(vfield(nx,nz,ny))
    vfield = 0
    allocate(complex_terrain(nx,ny))
    allocate(fft_terrain(nx,ny))

    call ideal_topography(nx,ny,zmax,terrain, 50)

    complex_terrain = terrain
    call compute_fft_topography(complex_terrain, fft_terrain)

    z_layers = [  50, 100, 300, 500, 500, 500, 500, 500, 500, 500, &
                 500, 500, 500, 500, 500]
    do i=2,nz
        z_layers(i) = z_layers(i-1) + z_layers(i)
    enddo

    call compute_blocked_flow_for_wind(u, v, z_layers, dx, fft_terrain, lt_data, ufield, vfield, debug=.False.)

    call io_write("ideal_terrain.nc","data",terrain)
    call io_write("ufield.nc", "data", ufield)
    call io_write("vfield.nc", "data", vfield)

    do i=1,nz
        wfield(2:nx-1, i, 2:ny-1) = ufield(3:nx,   i, 2:ny-1) - ufield(1:nx-2, i, 2:ny-1) &
                                  + vfield(2:nx-1, i, 3:ny  ) - vfield(2:nx-1, i, 1:ny-2)
    enddo
    call io_write("divergence.nc","data",wfield)

    wfield = wfield * (-1)

    do i=2,nz
        wfield(:,i,:) = wfield(:,i,:) + wfield(:,i-1,:)
    enddo
    call io_write("raww.nc","data",wfield)

    ! do i=1,nz
    !     wfield(2:nx-1, i, 2:ny-1) = wfield(2:nx-1, i, 2:ny-1)        &
    !                               + ufield(2:nx-1, i, 2:ny-1) * dzdx &
    !                               + vfield(2:nx-1, i, 2:ny-1) * dzdy
    ! enddo
    ! call io_write("w_full.nc","data",wfield)

contains

    subroutine ideal_topography(nx, ny, zmax, terrain, buffer_input)
        implicit none
        integer, intent(in) :: nx, ny
        real,    intent(in) :: zmax
        real,    intent(inout), allocatable :: terrain(:,:)
        integer, intent(in), optional :: buffer_input

        integer :: i,j
        integer :: buffer

        buffer = 0
        if (present(buffer_input)) buffer = buffer_input
        if (allocated(terrain)) deallocate(terrain)

        allocate(terrain(nx, ny))

        do i=buffer+1,nx-buffer
            do j=buffer+1,ny-buffer
                terrain(i,j) = (sin((i-buffer)/real(nx-(buffer*2)) * 2*pi - pi/2)/2+0.5) &
                                * (sin((j-buffer)/real(ny-(buffer*2)) * 2*pi - pi/2)/2+0.5) &
                                * zmax
            enddo
        enddo

    end subroutine ideal_topography

end program test_blocking
