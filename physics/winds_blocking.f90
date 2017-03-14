module mod_blocking

    use data_structures
    use linear_theory_winds, only : linear_perturbation_at_height, initialize_linear_theory_data, linear_perturbation, add_buffer_topo
    use io_routines,         only : io_write, io_read, file_exists
    use mod_atm_utilities,   only : calc_u, calc_v, calc_direction, calc_speed, &
                                    calc_froude, calc_dry_stability, blocking_fraction
    use array_utilities,     only : linear_space, calc_weight, check_array_dims
    use string,              only : str

    use fft
    use fftshifter

    use, intrinsic :: iso_c_binding

    implicit none
    private
    public add_blocked_flow

    ! although public, these are only exposed for testing purposes
    public compute_blocked_flow_for_wind, compute_fft_topography, compute_flow_for_wind

    type(linear_theory_type) :: lt_data
    !$omp threadprivate(lt_data)
    real, allocatable, dimension(:,:,:,:,:) :: blocked_u_lut, blocked_v_lut
    real, allocatable, dimension(:,:,:)     :: u_perturbation, v_perturbation
    real, allocatable, dimension(:)         :: dir_values, spd_values

    ! note this could be another config option, but it *shouldnt* matter much
    ! this just sets the maximum size of the smoothing window used when buffer topography
    ! ultimately, buffered topography should probably be created exactly as it is in linear_winds...
    integer, parameter :: smooth_window = 5
    ! these are just short cuts for values from the options structure
    integer :: ndir, nspeed
    real :: dir_min, dir_max
    real :: spd_min, spd_max
    real :: Nsq = 1e-4
    real :: fraction_continued_divergence = 0.05
    real :: minimum_step = 100
    integer :: buffer = 0

    real    :: blocking_contribution  ! fractional contribution of flow blocking perturbation that is added [0-1]
    real    :: smooth_froude_distance ! distance (m) over which Froude number is smoothed
    integer :: nsmooth_gridcells      ! number of grid cells corresponding to smooth_froude_distance/dx
    integer :: n_smoothing_passes     ! number of times the smoothing window is applied
    real    :: max_froude             ! max froude no at which flow is only partially blocked above, no blocking

    logical :: initialized = .False.

contains

    subroutine add_blocked_flow(domain, options)
        implicit none
        type(domain_type),  intent(inout) :: domain
        type(options_type), intent(in)    :: options

        if (.not.initialized) then
            call initialize_blocking(domain, options)
        endif

        call update_froude_number(domain)
        print*, "adding blocked flow"
        call spatial_blocking(domain, options%lt_options%stability_window_size)

    end subroutine add_blocked_flow

    subroutine update_froude_number(domain)
        implicit none
        type(domain_type), intent(inout) :: domain

        integer :: i, j, k, nx, ny, nz
        real :: wind_speed, u, v, stability
        real :: z_top, z_bot, th_top, th_bot
        integer :: ymin, ymax, xmin, xmax
        real, allocatable :: temp_froude(:,:)


        nx = size(domain%terrain_blocking,1)
        ny = size(domain%terrain_blocking,2)
        nz = size(domain%z_layers)

        allocate(temp_froude(nx,ny))

        if (.not.domain%blocking_initialized) then
            call compute_terrain_blocking_heights(domain%terrain_blocking, domain%terrain)
            domain%blocking_initialized = .True.
        endif

        ! use a homogenous background temperature gradient and wind field from the boundaryies
        th_bot = (sum(domain%th(:, 1,1)) / nx + sum(domain%th(:, 1,ny)) / nx)  / 2
        th_top = (sum(domain%th(:,nz,1)) / nx + sum(domain%th(:,nz,ny)) / nx)  / 2

        u = (sum(domain%u(:,:,1)) / (nx*nz) + sum(domain%u(:,:,ny)) / (nx*nz)) / 2
        v = (sum(domain%v(:,:,1)) / (nx*nz) + sum(domain%v(:,:,ny)) / (nx*nz)) / 2
        wind_speed = sqrt(u**2 + v**2)

        z_bot  = domain%z(1, 1,1)
        z_top  = domain%z(1,nz,1)
        do j=1,ny
            do i=1,nx
                ! u = sum(domain%u(i:i+1,:,j    )) / (nz*2)
                ! v = sum(domain%v(i,    :,j:j+1)) / (nz*2)
                ! wind_speed = sqrt(u**2 + v**2)

                ! th_bot = domain%th(i,1,j)
                ! th_top = domain%th(i,nz,j)
                ! z_bot  = domain%z(i,1,j)
                ! z_top  = domain%z(i,nz,j)
                stability = calc_dry_stability(th_top, th_bot, z_top, z_bot) !sum(domain%nsquared(i,:,j)) / nz
                stability = sqrt(max(stability, 0.))

                temp_froude(i,j) = calc_froude(stability, domain%terrain_blocking(i,j), wind_speed)
            enddo
        enddo

        do k = 1,n_smoothing_passes
            do j=1,ny
                ymin = max(j-nsmooth_gridcells, 1)
                ymax = min(j+nsmooth_gridcells, ny)
                do i=1,nx
                    xmin = max(i-nsmooth_gridcells, 1)
                    xmax = min(i+nsmooth_gridcells, nx)

                    domain%froude(i,j) = sum(temp_froude(xmin:xmax,ymin:ymax)) / ((xmax-xmin+1) * (ymax-ymin+1))
                enddo
            enddo

            if (k/=n_smoothing_passes) then
                temp_froude = domain%froude
            endif
        enddo

    end subroutine update_froude_number

    !>----------------------------------------------------------
    !! Compute a spatially variable linear wind perturbation
    !! based off of look uptables computed in via setup
    !! for each grid point, find the closest LUT data in U and V space
    !! then bilinearly interpolate the nearest LUT values for that points linear wind field
    !!
    !!----------------------------------------------------------
    subroutine spatial_blocking(domain, winsz)
        implicit none
        type(domain_type),intent(inout)::domain
        integer, intent(in) :: winsz

        integer :: nx,nxu, ny,nyv, nz, i,j,k, smoothz
        integer :: uk, vi !store a separate value of i for v and of k for u to we can handle nx+1, ny+1
        integer :: step, dpos, spos, nexts, nextd
        integer :: north, south, east, west, top, bottom, n
        real :: u, v, froude
        real :: dweight, sweight, curspd, curdir, wind_first, wind_second

        nx=size(domain%lat,1)
        ny=size(domain%lat,2)
        nz=size(domain%u,2)
        nxu=size(domain%u,1)
        nyv=size(domain%v,3)

        ! $omp parallel firstprivate(nx,nxu,ny,nyv,nz, winsz), default(none), &
        ! $omp private(i,j,k,step, uk, vi, east, west, north, south, top, bottom), &
        ! $omp private(spos, dpos, nexts,nextd, n, smoothz, u, v), &
        ! $omp private(curspd, curdir, sweight, dweight, froude), &
        ! $omp shared(domain, spd_values, dir_values, blocked_u_lut, blocked_v_lut), &
        ! $omp shared(u_perturbation, v_perturbation), &
        ! $omp shared(ndir, nspeed)
        ! $omp do
        do k=1, nyv
            do j=1, nz
                do i=1, nxu

                    froude = domain%froude(min(i,nx), min(k,ny))
                    if (froude < max_froude) then

                        uk = min(k,ny)
                        vi = min(i,nx)

                        !   First find the bounds of the region to average over
                        west  = max(i - winsz, 1)
                        east  = min(i + winsz,nx)
                        bottom= max(j - winsz, 1)
                        top   = min(j + winsz,nz)
                        south = max(k - winsz, 1)
                        north = min(k + winsz,ny)

                        ! smooth the winds vertically first
                        u = sum(domain%u( i, bottom:top,uk)) / (top-bottom+1)
                        v = sum(domain%v(vi, bottom:top, k)) / (top-bottom+1)

                        n = (((east-west)+1) * ((north-south)+1))
                        n = n * ((top-bottom)+1)

                        ! Calculate the direction of the current grid cell wind
                        dpos = 1
                        curdir = calc_direction( u, v )
                        ! and find the corresponding position in the Look up Table
                        do step=1, ndir
                            if (curdir > dir_values(step)) then
                                dpos = step
                            endif
                        end do

                        ! Calculate the wind speed of the current grid cell
                        spos = 1
                        curspd = calc_speed( u, v )
                        ! and find the corresponding position in the Look up Table
                        do step=1, nspeed
                            if (curspd > spd_values(step)) then
                                spos = step
                            endif
                        end do

                        ! Calculate the weights and the "next" u/v position
                        ! "next" usually = pos+1 but for edge cases next = 1 or n
                        dweight = calc_weight(dir_values, dpos, nextd, curdir)
                        sweight = calc_weight(spd_values, spos, nexts, curspd)

                        if ((i==20).and.((k>15).or.(k<25))) then
                            print*, u, j, froude, blocking_fraction(froude)
                        endif
                        ! perform linear interpolation between LUT values
                        if (k<=ny) then
                            u_perturbation(i,j,k) =      sweight  * (dweight    * blocked_u_lut(i,j,k, dpos, spos)     &
                                                                  + (1-dweight) * blocked_u_lut(i,j,k, nextd,spos))    &
                                                   +  (1-sweight) * (dweight    * blocked_u_lut(i,j,k, dpos, nexts)    &
                                                                  + (1-dweight) * blocked_u_lut(i,j,k, nextd,nexts))

                            ! u_perturbation(i,j,k) = u_perturbation(i,j,k) * (1-linear_update_fraction) &
                            !             + linear_update_fraction * (sweight*wind_first + (1-sweight)*wind_second)
                            u_perturbation(i,j,k) = u_perturbation(i,j,k) * blocking_fraction(froude) * blocking_contribution
                            domain%u(i,j,k) = domain%u(i,j,k) + u_perturbation(i,j,k) !* linear_mask(min(nx,i),min(ny,k))
                        endif
                        if (i<=nx) then
                            v_perturbation(i,j,k) =      sweight  * (dweight    * blocked_v_lut(i,j,k, dpos, spos)     &
                                                                  + (1-dweight) * blocked_v_lut(i,j,k, nextd,spos))    &
                                                   +  (1-sweight) * (dweight    * blocked_v_lut(i,j,k, dpos, nexts)    &
                                                                  + (1-dweight) * blocked_v_lut(i,j,k, nextd,nexts))

                            ! v_perturbation(i,j,k) = v_perturbation(i,j,k) * (1-linear_update_fraction) &
                            !             + linear_update_fraction * (sweight*wind_first + (1-sweight)*wind_second)

                            v_perturbation(i,j,k) = v_perturbation(i,j,k) * blocking_fraction(froude) * blocking_contribution
                            domain%v(i,j,k) = domain%v(i,j,k) + v_perturbation(i,j,k) !* linear_mask(min(nx,i),min(ny,k))
                        endif
                    endif
                end do
            end do
        end do
        ! $omp end do
        ! $omp end parallel
    end subroutine spatial_blocking


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

        integer :: nx, ny, nz, err
        real :: dx

        dx = domain%dx
        nx = size(domain%terrain,1)
        ny = size(domain%terrain,2)
        nz = size(domain%z, 2)

        write(*,*) ""
        write(*,*) "Initializing Blocking Parameterization"
        fraction_continued_divergence = 0.05
        minimum_step = options%lt_options%minimum_layer_size
        buffer = options%lt_options%buffer

        blocking_contribution   = options%block_options%blocking_contribution
        smooth_froude_distance  = options%block_options%smooth_froude_distance
        nsmooth_gridcells       = nint(smooth_froude_distance / options%dx)
        n_smoothing_passes      = options%block_options%n_smoothing_passes
        max_froude              = options%block_options%block_fr_max

        call compute_terrain_blocking_heights(domain%terrain_blocking, domain%terrain)

        call initialize_directions(options)
        call initialize_speeds(options)
        ! call initialize_linear_theory_data(lt_data, nx+buffer*2, ny+buffer*2, dx)

        allocate(u_perturbation(nx+1,nz,ny  ))
        allocate(v_perturbation(nx,  nz,ny+1))

        err = 0
        if (file_exists("u_blocked_lut.nc")) then
            call io_read("u_blocked_lut.nc","data",blocked_u_lut)
            if (.not.check_array_dims(blocked_u_lut, nx+1, nz, ny,   ndir, nspeed)) then
                deallocate(blocked_u_lut)
                err = 1
            endif
        else
            err = 1
        endif

        if (allocated(blocked_u_lut)) then
            if (file_exists("v_blocked_lut.nc")) then
                call io_read("v_blocked_lut.nc","data",blocked_v_lut)
                if (.not.check_array_dims(blocked_v_lut, nx, nz, ny+1,   ndir, nspeed)) then
                    deallocate(blocked_v_lut)
                    deallocate(blocked_u_lut)
                    err = 1
                endif
            else
                err = 1
            endif
        endif

        if (err/=0) then
            allocate(blocked_u_lut(nx+1, nz, ny,   ndir, nspeed))
            allocate(blocked_v_lut(nx,   nz, ny+1, ndir, nspeed))

            call generate_blocked_flow_lut(domain%terrain, domain%z_interface_layers, dx, buffer)

            call io_write("u_blocked_lut.nc","data",blocked_u_lut)
            call io_write("v_blocked_lut.nc","data",blocked_v_lut)
        endif

        write(*,*) "Blocking Parameterization Initialized"

        domain%blocking_initialized = .True.
        initialized = .True.

    end subroutine initialize_blocking

    !>-----------------------------------------
    !> Compute a smoothed terrain varience field for use in Froude number calculation
    !>
    !------------------------------------------
    subroutine compute_terrain_blocking_heights(terrain_blocking, terrain)
        implicit none
        real, intent(inout) :: terrain_blocking(:,:)
        real, intent(in)    :: terrain(:,:)

        integer :: nx, ny, x, y
        integer :: xs,xe, ys,ye, n
        integer :: window_size, smooth_window
        real, allocatable :: temp_terrain(:,:)
        integer :: i

        window_size   = 5
        smooth_window = 2

        nx = size(terrain,1)
        ny = size(terrain,2)

        allocate(temp_terrain(nx,ny))

        ! first compute lightly smoothed terrain
        do y=1,ny
            ys = max( y - smooth_window, 1)
            ye = min( y + smooth_window, ny)
            do x=1,nx
                xs = max( x - smooth_window, 1)
                xe = min( x + smooth_window, nx)
                n = (xe-xs+1) * (ye-ys+1)
                terrain_blocking(x,y) = sum(terrain(xs:xe,ys:ye)) / n
            enddo
        enddo

        ! then compute the range of terrain (max-min) in a given window
        do y=1,ny
            ys = max( y - window_size, 1)
            ye = min( y + window_size, ny)
            do x=1,nx
                xs = max( x - window_size, 1)
                xe = min( x + window_size, nx)
                temp_terrain(x,y) = maxval(terrain_blocking(xs:xe,ys:ye)) - minval(terrain_blocking(xs:xe,ys:ye))
            enddo
        enddo
        ! call io_write("initial_terrain_delta.nc","data",temp_terrain)

        ! finally smooth that terrain delta field a few times as well
        smooth_window = 2
        do i=1,n_smoothing_passes
            do y=1,ny
                ys = max( y - smooth_window, 1)
                ye = min( y + smooth_window, ny)
                do x=1,nx
                    xs = max( x - smooth_window, 1)
                    xe = min( x + smooth_window, nx)
                    n = (xe-xs+1) * (ye-ys+1)
                    terrain_blocking(x,y) = sum(temp_terrain(xs:xe,ys:ye)) / n
                enddo
            enddo
            if (i /= n_smoothing_passes) then
                temp_terrain = terrain_blocking
            endif
        enddo
        ! call io_write("terrain_blocking.nc","data",terrain_blocking)

    end subroutine compute_terrain_blocking_heights

    subroutine generate_blocked_flow_lut(terrain, layers, dx, buf_in)
        implicit none
        real,    intent(in) :: terrain(:,:)
        real,    intent(in) :: layers(:)
        real,    intent(in) :: dx
        integer, intent(in), optional :: buf_in

        integer :: nx, ny, nz, dir, speed, buf, lut_iteration
        real    :: u, v, direction, magnitude
        real, allocatable, dimension(:,:,:)     :: u_field, v_field
        complex(C_DOUBLE_COMPLEX), allocatable  :: fft_terrain(:,:)       !> Fourier transform of the terrain


        nx  = size(terrain,1)
        ny  = size(terrain,2)
        nz  = size(layers)

        buf = buffer
        if (present(buf_in)) buf = buf_in

        allocate(u_field(nx+buf*2, nz, ny+buf*2))
        allocate(v_field(nx+buf*2, nz, ny+buf*2))

        call create_fft_terrain(terrain, fft_terrain, buf)

        !$omp parallel default(shared) &! though not listed, lt_data is defined as threadprivate at the module level
        !$omp private(lut_iteration, dir, speed, u, v, direction, magnitude, u_field, v_field)  &
        !$omp firstprivate(ndir, nspeed, dir_values, spd_values)                                &
        !$omp firstprivate(layers, dx, nx, ny, buf)
        !$omp do
        do lut_iteration = 0, ndir*nspeed-1
            ! make two loops into one to increase the granularity of parallelization
            dir     = lut_iteration / nspeed + 1
            speed   = mod(lut_iteration,nspeed) + 1

            direction = dir_values(dir)
            magnitude = spd_values(speed)

            u = calc_u(direction, magnitude)
            v = calc_v(direction, magnitude)

            call compute_blocked_flow_for_wind(u, v, layers, dx, fft_terrain, lt_data, u_field, v_field, debug=.False.)

            blocked_u_lut(2:nx,:, :,  dir, speed) = ( u_field(buf+2:buf+nx,:,buf:ny+buf) + u_field(buf+1:buf+nx-1,:,buf:ny+buf) ) /2
            blocked_u_lut( 1,  :, :,  dir, speed) =   u_field(buf+1, :,buf:ny+buf)
            blocked_u_lut(nx+1,:, :,  dir, speed) =   u_field(buf+nx,:,buf:ny+buf)

            blocked_v_lut( :,  :,2:ny,dir, speed) = ( v_field(buf:nx+buf,:,buf+2:buf+ny) + v_field(buf:nx+buf,:,buf+1:buf+ny-1) ) /2
            blocked_v_lut( :,  :, 1,  dir, speed) =   v_field(buf:nx+buf,:,buf+1)
            blocked_v_lut( :,  :,ny+1,dir, speed) =   v_field(buf:nx+buf,:,buf+ny)
        enddo
        !$omp end do
        !$omp end parallel

        ! oddly, it was getting non-zero values for 0 wind speed, though all other values looked OK
        ! should probably make something like this into a small function for general usage
        ! if (abs(spd_values(1)/max(abs(spd_values(1)),kSMALL_VALUE)) < kSMALL_VALUE) then
        !     blocked_u_lut(:,:,:,:,1) = 0
        !     blocked_v_lut(:,:,:,:,1) = 0
        ! endif


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


    subroutine compute_blocked_flow_for_wind(u, v, z_layers, dx, fft_terrain, lt_data, ufield, vfield, debug)
        implicit none
        real,                       intent(in)    :: u, v              !> background u and v magnitudes             [m / s]
        real,                       intent(in)    :: z_layers(:)       !> vertical coordinate =height of each layer [m]
        real,                       intent(in)    :: dx                !> horizontal grid spacing                   [m]
        complex(C_DOUBLE_COMPLEX),  intent(in)    :: fft_terrain(:,:)  !> Fourier transform of the terrain
        type(linear_theory_type),   intent(inout) :: lt_data           !> linear theory data structure
        real,                       intent(inout) :: ufield(:,:,:)     !> output 3D U values for blocked flow       [m / s]
        real,                       intent(inout) :: vfield(:,:,:)     !> output 3D V values for blocked flow       [m / s]
        logical,                    intent(in)    :: debug


        real, allocatable :: wfield(:,:,:)
        integer :: nx, ny, nz, i, key_level
        real :: z, z_bottom, z_top

        nx = size(fft_terrain,1)
        ny = size(fft_terrain,2)
        nz = size(z_layers) - 1

        allocate(wfield(nx,nz,ny))
        wfield = 0

        if (.not.allocated(lt_data%sig)) then
            call initialize_linear_theory_data(lt_data, nx, ny, dx)
        endif

        do i=1,nz
            z_bottom = z_layers(i)
            z_top    = z_layers(i+1)

            call linear_perturbation(u, v, Nsq, z_bottom, z_top, minimum_step, fft_terrain, lt_data)
            ! call linear_perturbation_at_height(u, v, Nsq, z, fft_terrain, lt_data)

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
        if (debug) then
            !$omp critical (print_lock)
            print*, key_level, u, v
            !$omp end critical (print_lock)
        endif
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

    subroutine create_fft_terrain(terrain, fft_terrain, buffer)
        implicit none
        real,                                    intent(in)    :: terrain(:,:)
        complex(C_DOUBLE_COMPLEX), allocatable,  intent(inout) :: fft_terrain(:,:)       !> Fourier transform of the terrain
        integer,                                 intent(in)    :: buffer                 !> number of grid cells to add as a buffer around the terrain

        integer :: nx, ny
        complex(C_DOUBLE_COMPLEX), allocatable :: complex_terrain(:,:)   !> Terrain as a complex data array

        nx = size(terrain,1)
        ny = size(terrain,2)

        call add_buffer_topo(terrain, complex_terrain, smooth_window, buffer)

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
