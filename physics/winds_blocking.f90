module mod_blocking

    use data_structures
    use linear_theory_winds, only : linear_perturbation_at_height, initialize_linear_theory_data
    use mod_atm_utilities
    use array_utilities,     only : linear_space

    use fft
    use fftshifter

    use, intrinsic :: iso_c_binding

    implicit none
    private
    public add_blocked_flow

    ! although public, these are only exposed for testing purposes
    public compute_blocked_flow_for_wind, compute_fft_topography, compute_flow_for_wind

    type(linear_theory_type) :: lt_data
    real, allocatable, dimension(:,:,:,:,:) :: blocked_u_lut, blocked_v_lut
    real, allocatable, dimension(:)         :: dir_values, spd_values
    integer :: ndir, nspeed
    real :: dir_min, dir_max
    real :: spd_min, spd_max
    real :: Nsq = 1e-4
    real :: fraction_continued_divergence

    logical :: initialized = .False.

contains

    subroutine add_blocked_flow(domain, options)
        implicit none
        type(domain_type),  intent(inout) :: domain
        type(options_type), intent(in)    :: options

        if (.not.initialized) then
            call initialize_blocking(domain, options)
        endif

    end subroutine add_blocked_flow

    !>---------------------------------------------
    !> Initialize module level variables for blocking parameterization
    !>
    !> This includes generating the wind field LUT
    !>
    !----------------------------------------------
    subroutine initialize_blocking(domain, options)
        implicit none
        type(domain_type),  intent(inout) :: domain
        type(options_type), intent(in)    :: options

        integer :: nx, ny, nz
        real :: dx

        dx = domain%dx
        nx = size(domain%terrain,1)
        ny = size(domain%terrain,2)
        nz = size(domain%z, 2)

        fraction_continued_divergence = 0.05

        write(*,*) ""
        write(*,*) "Initializing Blocking Parameterization"
        call initialize_directions(options)
        call initialize_speeds(options)
        call initialize_linear_theory_data(lt_data, nx, ny, dx)

        allocate(blocked_u_lut(nx+1, nz, ny,   ndir, nspeed))
        allocate(blocked_v_lut(nx,   nz, ny+1, ndir, nspeed))

        call generate_blocked_flow_lut(domain%terrain, domain%z(1,:,1)-domain%terrain(1,1), dx)

        write(*,*) "Blocking Parameterization Initialized"
        initialized = .True.
    end subroutine initialize_blocking

    subroutine generate_blocked_flow_lut(terrain, layers, dx)
        implicit none
        real, intent(in) :: terrain(:,:)
        real, intent(in) :: layers(:)
        real, intent(in) :: dx

        integer :: nx, ny, nz, dir, speed
        real    :: u, v, direction, magnitude
        real, allocatable, dimension(:,:,:)     :: u_field, v_field
        complex(C_DOUBLE_COMPLEX), allocatable  :: fft_terrain(:,:)       !> Fourier transform of the terrain

        nx  = size(terrain,1)
        ny  = size(terrain,2)
        nz  = size(layers)

        allocate(u_field(nx, nz, ny))
        allocate(v_field(nx, nz, ny))

        call create_fft_terrain(terrain, fft_terrain)

        do dir = 1,ndir
            print*, dir, ndir
            do speed = 1,nspeed

                direction = dir_values(direction)
                magnitude = spd_values(speed)

                u = calc_u(direction, magnitude)
                v = calc_v(direction, magnitude)

                call compute_blocked_flow_for_wind(u, v, layers, dx, fft_terrain, lt_data, u_field, v_field)

                blocked_u_lut(2:nx,:, :,  dir, speed) = ( u_field(2:nx,:,:) + u_field(1:nx-1,:,:) ) /2
                blocked_u_lut( 1,  :, :,  dir, speed) = u_field(1,:,:)
                blocked_u_lut(nx+1,:, :,  dir, speed) = u_field(nx,:,:)

                blocked_v_lut( :,  :,2:ny,dir, speed) = ( v_field(:,:,2:ny) + v_field(:,:,1:ny-1) ) /2
                blocked_v_lut( :,  :, 1,  dir, speed) = v_field(:,:,1)
                blocked_v_lut( :,  :,ny+1,dir, speed) = v_field(:,:,ny)
            enddo
        enddo

    end subroutine generate_blocked_flow_lut

    subroutine compute_flow_for_wind(u, v, z_layers, dx, fft_terrain, lt_data, ufield, vfield)
        implicit none
        real,                       intent(in)    :: u, v              !> background u and v magnitudes             [m / s]
        real,                       intent(in)    :: z_layers(:)       !> vertical coordinate =height of each layer [m]
        real,                       intent(in)    :: dx                !> horizontal grid spacing                   [m]
        complex(C_DOUBLE_COMPLEX),  intent(in)    :: fft_terrain(:,:)  !> Fourier transform of the terrain
        type(linear_theory_type),   intent(inout) :: lt_data           !> linear theory data structure
        real,                       intent(inout) :: ufield(:,:,:)     !> output 3D U values for blocked flow       [m / s]
        real,                       intent(inout) :: vfield(:,:,:)     !> output 3D V values for blocked flow       [m / s]

        integer :: nx, ny, nz, i
        real :: z

        nx = size(fft_terrain,1)
        ny = size(fft_terrain,2)
        nz = size(z_layers)

        if (.not.allocated(lt_data%sig)) then
            call initialize_linear_theory_data(lt_data, nx, ny, dx)
        endif

        do i=1,nz
            z = z_layers(i)
            call linear_perturbation_at_height(u, v, Nsq, z, fft_terrain, lt_data)
            ufield(:,i,:) = real(real(lt_data%u_perturb))
            vfield(:,i,:) = real(real(lt_data%v_perturb))
        enddo

    end subroutine compute_flow_for_wind


    subroutine compute_blocked_flow_for_wind(u, v, z_layers, dx, fft_terrain, lt_data, ufield, vfield)
        implicit none
        real,                       intent(in)    :: u, v              !> background u and v magnitudes             [m / s]
        real,                       intent(in)    :: z_layers(:)       !> vertical coordinate =height of each layer [m]
        real,                       intent(in)    :: dx                !> horizontal grid spacing                   [m]
        complex(C_DOUBLE_COMPLEX),  intent(in)    :: fft_terrain(:,:)  !> Fourier transform of the terrain
        type(linear_theory_type),   intent(inout) :: lt_data           !> linear theory data structure
        real,                       intent(inout) :: ufield(:,:,:)     !> output 3D U values for blocked flow       [m / s]
        real,                       intent(inout) :: vfield(:,:,:)     !> output 3D V values for blocked flow       [m / s]


        real, allocatable :: wfield(:,:,:)
        integer :: nx, ny, nz, i, key_level
        real :: z

        nx = size(fft_terrain,1)
        ny = size(fft_terrain,2)
        nz = size(z_layers)

        allocate(wfield(nx,nz,ny))
        wfield = 0

        if (.not.allocated(lt_data%sig)) then
            call initialize_linear_theory_data(lt_data, nx, ny, dx)
        endif

        do i=1,nz
            z = z_layers(i)
            call linear_perturbation_at_height(u, v, Nsq, z, fft_terrain, lt_data)

            ufield(:,i,:) = real(real(lt_data%u_perturb))
            vfield(:,i,:) = real(real(lt_data%v_perturb))

            wfield(2:nx-1, i, 2:ny-1) = ufield(1:nx-2, i, 2:ny-1) - ufield(3:nx,   i, 2:ny-1) &
                                      + vfield(2:nx-1, i, 1:ny-2) - vfield(2:nx-1, i, 3:ny  )
            if (i>1) then
                wfield(:,i,:) = wfield(:,i,:) + wfield(:,i-1,:)
            endif
        enddo

        where(wfield > 0) wfield = 0

        key_level = find_maximum_downward_motion(wfield)

        if (key_level < nz) then
            do i=key_level+1,nz
                ufield(:,i,:) = ufield(:,key_level,:) * fraction_continued_divergence
                vfield(:,i,:) = vfield(:,key_level,:) * fraction_continued_divergence
            enddo
        endif


    end subroutine compute_blocked_flow_for_wind

    function find_maximum_downward_motion(wfield) result(max_level)
        implicit none
        real, intent(in) :: wfield(:,:,:)
        integer :: max_level

        integer :: nz, i
        real    :: minw, current_w

        nz = size(wfield,2)
        minw = 999999
        max_level = 1
        do i = 1, nz
            current_w = sum(wfield(:,i,:))
            if (current_w < minw) then
                max_level = i
                minw = current_w
            else
                ! if we have already found a maximum level, and we are now on the way down, then just return
                if (max_level /= 1) then
                    return
                endif
            endif
        end do

    end function find_maximum_downward_motion

    subroutine create_fft_terrain(terrain, fft_terrain)
        implicit none
        real,                                    intent(in)    :: terrain(:,:)
        complex(C_DOUBLE_COMPLEX), allocatable,  intent(inout) :: fft_terrain(:,:)       !> Fourier transform of the terrain

        integer :: nx, ny
        complex(C_DOUBLE_COMPLEX), allocatable :: complex_terrain(:,:)   !> Terrain as a complex data array

        nx = size(terrain,1)
        ny = size(terrain,2)
        allocate(complex_terrain(nx,ny))
        complex_terrain = terrain
        call compute_fft_topography(complex_terrain,fft_terrain)

    end subroutine create_fft_terrain

    subroutine compute_fft_topography(complex_terrain, fft_terrain)
        implicit none
        complex(C_DOUBLE_COMPLEX),               intent(inout) :: complex_terrain(:,:)   !> Terrain as a complex data array
        complex(C_DOUBLE_COMPLEX), allocatable,  intent(inout) :: fft_terrain(:,:)       !> Fourier transform of the terrain

        type(C_PTR) :: plan
        integer :: nx, ny

        nx = size(complex_terrain,1)
        ny = size(complex_terrain,2)

        if (allocated(fft_terrain)) then
            if ((size(fft_terrain,1) /= nx).or.(size(fft_terrain,2) /= ny)) then
                deallocate(fft_terrain)
            endif
        endif
        if (.not.allocated(fft_terrain)) then
            allocate(fft_terrain(nx,ny))
        endif

        plan = fftw_plan_dft_2d(ny, nx, complex_terrain, fft_terrain, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, complex_terrain, fft_terrain)
        call fftw_destroy_plan(plan)
        ! normalize FFT by N - grid cells
        fft_terrain = fft_terrain / (nx * ny)

        ! shift the grid cell quadrants
        call fftshift(fft_terrain)

    end subroutine compute_fft_topography

    subroutine initialize_directions(options)
        implicit none
        type(options_type), intent(in) :: options

        ndir    = options%lt_options%n_dir_values
        dir_min = options%lt_options%dirmin
        dir_max = options%lt_options%dirmax

        call linear_space(dir_values, dir_min, dir_max, ndir)

    end subroutine initialize_directions

    subroutine initialize_speeds(options)
        implicit none
        type(options_type), intent(in) :: options

        nspeed  = options%lt_options%n_spd_values
        spd_min = options%lt_options%spdmin
        spd_max = options%lt_options%spdmax

        call linear_space(spd_values, spd_min, spd_max, nspeed)

    end subroutine initialize_speeds

end module mod_blocking
