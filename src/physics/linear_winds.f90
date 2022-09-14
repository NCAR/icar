!>----------------------------------------------------------
!! This module provides the linear wind theory calculations
!!
!!  Code is based primarily off equations in Barstad and Gronas (2006)
!!  @see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
!!
!! The main entry point to the code is:
!!      linear_perturb(domain, options, vsmooth, reverse, useDensity)
!!
!! <pre>
!! Call tree graph :
!!  linear_perturb->[ setup_linwinds -> add_buffer_topo,
!!                    calc_domain_stability,
!!                    linear_winds -> various fft routines]
!!
!! High level routine descriptions / purpose
!!   calc_domain_stability    - calculates a mean Brunt Vaisala frequency over the domain
!!   linear_winds             - primary routine that calculates the linear wind perturbation
!!   add_buffer_topo          - generates a topo grid that with a surrounding buffer for the fft
!!   initialize_spatial_winds - generated the look up tables to use spatially varying linear winds
!!   setup_linwinds           - sets up module level variables and calls add_buffer_topo
!!   linear_perturb           - main entry point, calls setup on first entry for a given domain
!!
!! Inputs: domain, options, vsmooth, reverse, useDensity
!!      domain,options  = as defined in data_structures
!!      vsmooth         = number of vertical levels to smooth winds over
!!      reverse         = remove linear perturbation instead of adding it
!!      useDensity      = create a linear field that attempts to mitigate the
!!                          boussinesq approx that is embedded in the linear theory
!!                          so it advection can properly incorporate density.
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module linear_theory_winds
    use iso_c_binding
    use fft,                        only: fftw_execute_dft,                    &
                                          fftw_alloc_complex, fftw_free,       &
                                          fftw_plan_dft_2d, fftw_destroy_plan, &
                                          FFTW_FORWARD, FFTW_MEASURE,          &
                                          FFTW_BACKWARD, FFTW_ESTIMATE ! note fft module is defined in fftshift.f90
    use fftshifter,                 only: ifftshift, fftshift
    use data_structures
    use domain_interface,           only: domain_t
    !use io_routines,                only: io_read, io_write
    use string,                     only: str
    use grid_interface,             only: grid_t
    use linear_theory_lut_disk_io,  only: read_LUT, write_LUT
    use mod_atm_utilities,         only : calc_stability, calc_u, calc_v, &
                                          calc_speed, calc_direction,     &
                                          blocking_fraction, options_t
    use array_utilities,            only: smooth_array, calc_weight, &
                                          linear_space
    use icar_constants,             only: kMAX_FILE_LENGTH

    implicit none

    private
    public :: linear_perturb, linear_perturbation, linear_perturbation_at_height
    public :: add_buffer_topo, initialize_linear_theory_data, setup_linwinds

    logical :: module_initialized = .False.
    logical :: variable_N
    logical :: smooth_nsq
    logical :: using_blocked_flow
    real    :: N_squared

    logical :: debug_print = .false.

    interface linear_perturbation
        module procedure linear_perturbation_constz, linear_perturbation_varyingz
    end interface

    !! unfortunately these have to be allocated every call because we could be calling on both the high-res
    !! domain and the low res domain (to "remove" the linear winds)
    real,                       allocatable,    dimension(:,:)  :: k, l, kl, sig
    complex(C_DOUBLE_COMPLEX),  allocatable,    dimension(:,:)  :: denom, m, ineta, msq, mimag
    complex(C_DOUBLE_COMPLEX),  pointer,        dimension(:,:,:):: uhat, u_hat, vhat, v_hat
    type(C_PTR),                allocatable,    dimension(:)    :: uplans, vplans
    type(C_PTR)                 :: uh_aligned_data, u_h_aligned_data, vh_aligned_data, v_h_aligned_data
    logical                     :: data_allocated=.False. ! need a boolean because you cant test if a cptr is allocated?

    ! note this data structure holds a lot of important variables for the new linear_perturbation subroutine
    ! this is created at the module scope to prevent some apparent stack issues related to OpenMP and ifort
    type(linear_theory_type) :: lt_data_m
    !$omp threadprivate(lt_data_m)

    !! Linear wind look up table values and tables
    real, allocatable,          dimension(:)           :: dir_values, nsq_values, spd_values
    ! Look Up Tables for linear perturbation are nspd x n_dir_values x n_nsq_values x nx x nz x ny
    real, allocatable,          dimension(:,:,:,:,:,:) :: hi_u_LUT[:], hi_v_LUT[:] !, rev_u_LUT, rev_v_LUT
    ! real, pointer,              dimension(:,:,:,:,:,:) :: u_LUT, v_LUT
    real, allocatable,          dimension(:,:)         :: linear_mask, nsq_calibration

    ! store the linear perturbation so we can update it slightly each time step
    ! this permits the linear field to build up over time.
    real, allocatable, target,  dimension(:,:,:) :: hi_u_perturbation, hi_v_perturbation, lo_u_perturbation, lo_v_perturbation
    real, pointer,              dimension(:,:,:) :: u_perturbation, v_perturbation

    logical :: use_spatial_linear_fields
    logical :: use_linear_mask     ! use a spatial mask for the linear wind field
    logical :: use_nsq_calibration ! use a spatial mask to calibrate the nsquared (brunt vaisala frequency) field

    integer :: buffer, original_buffer ! number of grid cells to buffer around the domain MUST be >=1
    integer :: stability_window_size

    real :: max_stability ! limits on the calculated Brunt Vaisala Frequency
    real :: min_stability ! these may need to be a little narrower.
    real :: linear_contribution = 1.0 ! multiplier on uhat,vhat before adding to u,v
    ! controls the rate at which the linearfield updates (should be calculated as f(in_dt))
    ! new/current perturbation is multiplited by linear_update and added to (1-linear_update) * the previous combined perturbation
    real :: linear_update_fraction = 1.0

    real :: dirmax ! =2*pi
    real :: dirmin ! =0
    real :: spdmax ! =30
    real :: spdmin ! =0
    real :: nsqmax ! =log(max_stability)
    real :: nsqmin ! =log(min_stability)
    real :: minimum_layer_size

    integer :: n_dir_values=36
    integer :: n_nsq_values=10
    integer :: n_spd_values=10

    complex,parameter :: imaginary_number = (0,1)

    real, parameter :: SMALL_VALUE = 1e-15

contains

    !>----------------------------------------------------------
    !! Calculate a single brunt vaisala frequency effectively averaged across the entire domain
    !!
    !!----------------------------------------------------------
    ! pure function calc_domain_stability(domain) result(BV_freq)
    !     implicit none
    !     class(domain_t),intent(in)::domain
    !     real :: BV_freq
    !     real, allocatable, dimension(:) :: vertical_N
    !     integer :: i, nx,nz,ny, top_layer,bottom_layer
    !     real :: dz,dtheta,t_mean
    !
    !     nx=size(domain%th,1)
    !     nz=size(domain%th,2)
    !     ny=size(domain%th,3)
    !     allocate(vertical_N(nz))
    !
    !     if (variable_N) then
    !         do i=1,nz
    !             top_layer    = min(i+stability_window_size, nz)
    !             bottom_layer = max(i-stability_window_size,  1)
    !
    !             dz     = sum(domain%z (:,top_layer,:)-domain%z (:,bottom_layer,:))/(nx*ny) ! distance between top layer and bottom layer
    !             dtheta = sum(domain%th(:,top_layer,:)-domain%th(:,bottom_layer,:))/(nx*ny) ! average temperature change between layers
    !             t_mean = sum(domain%th(:,i,:))/(nx*ny) ! mean temperature in the current layer
    !
    !             ! BV frequency = gravity/theta * dtheta/dz
    !             vertical_N(i)=gravity/t_mean * dtheta/dz
    !         end do
    !         ! calculate the mean stability over the entire profile
    !         BV_freq=sum(vertical_N)/nz
    !         ! impose limits so the linear solution doesn't go crazy in really stable or unstable air
    !         BV_freq=min(max(BV_freq, min_stability), max_stability)
    !     else
    !         ! or just use the supplied BV frequency
    !         BV_freq=N_squared
    !     endif
    !
    !     deallocate(vertical_N)
    ! end function calc_domain_stability

    !>----------------------------------------------------------
    !! Compute linear wind perturbations given background U, V, Nsq
    !!
    !! see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
    !!
    !!----------------------------------------------------------
    subroutine linear_perturbation_at_height(U, V, Nsq, z, fourier_terrain, lt_data)
        use, intrinsic :: iso_c_binding
        implicit none
        real,                     intent(in)    :: U, V ! U and V components of background wind
        real,                     intent(in)    :: Nsq  ! Brunt-Vaisalla frequency (N squared)
        real,                     intent(in)    :: z    ! elevation at which to compute the linear solution
        complex(C_DOUBLE_COMPLEX),intent(in)    :: fourier_terrain(:,:) ! FFT(terrain)
        type(linear_theory_type), intent(inout) :: lt_data


        integer :: nx, ny ! store the size of the grid

        if ((U==0).and.(V==0)) then
            lt_data%u_perturb = 0
            lt_data%v_perturb = 0
            return
        endif

        nx = size(fourier_terrain, 1)
        ny = size(fourier_terrain, 2)

        lt_data%sig  = U * lt_data%k + V * lt_data%l
        ! prevent divide by zero
        where(lt_data%sig == 0) lt_data%sig = SMALL_VALUE

        lt_data%denom = lt_data%sig**2 ! -f**2 ! to add coriolis

        lt_data%msq   = Nsq / lt_data%denom * lt_data%kl
        lt_data%mimag = (0, 0) ! be sure to reset real and imaginary components
        lt_data%mimag = lt_data%mimag + (real(sqrt(-lt_data%msq)) * imaginary_number)

        lt_data%m = sqrt(lt_data%msq)                         ! vertical wave number, hydrostatic
        where(lt_data%sig < 0) lt_data%m = lt_data%m * (-1)   ! equivalent to m = m * sign(sig)
        where(real(lt_data%msq) < 0) lt_data%m = lt_data%mimag

        lt_data%ineta = imaginary_number * fourier_terrain * exp(imaginary_number * lt_data%m * z)
        !  what=sig*ineta

        ! with coriolis : [-/+] j * [l/k] * f
        !uhat = (0 - m) * ((sig * k) - (j * l * f)) * ineta / kl
        !vhat = (0 - m) * ((sig * l) + (j * k * f)) * ineta / kl
        ! removed coriolis term assumes the scale coriolis operates at is largely defined by the coarse model
        ! if using e.g. a sounding or otherwise spatially constant u/v, then coriolis should be defined.
        lt_data%ineta = lt_data%ineta / (lt_data%kl / ((0 - lt_data%m) * lt_data%sig))
        lt_data%uhat  = lt_data%k * lt_data%ineta
        lt_data%vhat  = lt_data%l * lt_data%ineta

        ! pull it back out of fourier space.
        ! The fftw transform scales by N so the previous Fzs/Nx/Ny provides the only normalization necessary (?)
        call ifftshift(lt_data%uhat)
        call ifftshift(lt_data%vhat)
        ! plans are only created once and re-used, this is a limit to further parallelization at the moment.
        call fftw_execute_dft(lt_data%uplan, lt_data%uhat, lt_data%u_perturb)
        call fftw_execute_dft(lt_data%vplan, lt_data%vhat, lt_data%v_perturb)

        ! returns u_perturb and v_perturb
    end subroutine linear_perturbation_at_height

    subroutine linear_perturbation_constz(U, V, Nsq, z_bottom, z_top, minimum_step, fourier_terrain, lt_data)
        use, intrinsic :: iso_c_binding
        implicit none
        real,                     intent(in)    :: U, V             ! U and V components of background wind
        real,                     intent(in)    :: Nsq              ! Brunt-Vaisalla frequency (N squared)
        real,                     intent(in)    :: z_top, z_bottom  ! elevation for the top and bottom bound of a layer
        real,                     intent(in)    :: minimum_step     ! minimum layer step size to compute LT for
        complex(C_DOUBLE_COMPLEX),intent(in)    :: fourier_terrain(:,:) ! FFT(terrain)
        type(linear_theory_type), intent(inout) :: lt_data          ! intermediate arrays needed for LT calc

        integer :: n_steps, i
        real    :: step_size, current_z

        ! Handle the trivial case quickly
        if ((U==0).and.(V==0)) then
            lt_data%u_perturb = 0
            lt_data%v_perturb = 0
            return
        endif

        n_steps = max(1,ceiling((z_top - z_bottom) / minimum_step))
        step_size = (z_top - z_bottom) / n_steps

        lt_data%u_accumulator = 0
        lt_data%v_accumulator = 0

        current_z = z_bottom + step_size / 2

        do i=1,n_steps
            call linear_perturbation_at_height(U,V,Nsq,current_z, fourier_terrain, lt_data)
            lt_data%u_accumulator = lt_data%u_accumulator + lt_data%u_perturb
            lt_data%v_accumulator = lt_data%v_accumulator + lt_data%v_perturb

            current_z = current_z + step_size
        enddo

        lt_data%u_perturb = lt_data%u_accumulator / n_steps
        lt_data%v_perturb = lt_data%v_accumulator / n_steps
    end subroutine linear_perturbation_constz


    subroutine linear_perturbation_varyingz(U, V, Nsq, z_bottom, z_top, minimum_step, fourier_terrain, lt_data)
        use, intrinsic :: iso_c_binding
        implicit none
        real,                     intent(in)    :: U, V             ! U and V components of background wind
        real,                     intent(in)    :: Nsq              ! Brunt-Vaisalla frequency (N squared)
        real,                     intent(in)    :: z_top(:,:), z_bottom(:,:)  ! elevation for the top and bottom bound of a layer
        real,                     intent(in)    :: minimum_step     ! minimum layer step size to compute LT for
        complex(C_DOUBLE_COMPLEX),intent(in)    :: fourier_terrain(:,:) ! FFT(terrain)
        type(linear_theory_type), intent(inout) :: lt_data          ! intermediate arrays needed for LT calc

        integer :: n_steps, i
        real    :: step_size, current_z, start_z, end_z
        real,   allocatable :: layer_count(:,:), layer_fraction(:,:)
        real,   allocatable :: internal_z_top(:,:), internal_z_bottom(:,:)

        ! Handle the trivial case quickly
        if ((U==0).and.(V==0)) then
            lt_data%u_perturb = 0
            lt_data%v_perturb = 0
            return
        endif

        start_z = minval(z_bottom)
        end_z = maxval(z_top)

        allocate(internal_z_top(size(lt_data%u_perturb,1),size(lt_data%u_perturb,2)), source=maxval(z_top))
        internal_z_top(buffer:buffer+size(z_top,1)-1, buffer:buffer+size(z_top,2)-1) = z_top

        allocate(internal_z_bottom(size(lt_data%u_perturb,1),size(lt_data%u_perturb,2)), source=minval(z_bottom))
        internal_z_bottom(buffer:buffer+size(z_bottom,1)-1, buffer:buffer+size(z_bottom,2)-1) = z_bottom

        allocate(layer_count(size(lt_data%u_perturb,1), size(lt_data%u_perturb,2)))
        layer_count = 0
        allocate(layer_fraction, source=layer_count)


        step_size = min(minimum_step, minval(z_top - z_bottom))
        if (step_size < 5) then
            if (this_image()==1) write(*,*) "WARNING: Very small vertical step size in linear winds, check layer thicknesses in output"
        endif

        current_z = start_z + step_size/2 ! we want the value in the middle of each theoretical layer

        ! sum up u/v perturbations over all layers evaluated
        lt_data%u_accumulator = 0
        lt_data%v_accumulator = 0

        do while (current_z < end_z)

            call linear_perturbation_at_height(U,V,Nsq,current_z, fourier_terrain, lt_data)

            layer_fraction = max(0.0,                                                                                       &
                                min(step_size/2, current_z - internal_z_bottom) + min(0.0, internal_z_top - current_z)      &
                              + min(step_size/2, internal_z_top - current_z)    + min(0.0, current_z - internal_z_bottom)   ) / step_size

            layer_count = layer_count + layer_fraction
            lt_data%u_accumulator = lt_data%u_accumulator + lt_data%u_perturb * layer_fraction
            lt_data%v_accumulator = lt_data%v_accumulator + lt_data%v_perturb * layer_fraction

            current_z = current_z + step_size
        enddo

        lt_data%u_perturb = lt_data%u_accumulator / layer_count
        lt_data%v_perturb = lt_data%v_accumulator / layer_count
    end subroutine linear_perturbation_varyingz

    !>----------------------------------------------------------
    !! Add a smoothed buffer around the edge of the terrain to prevent crazy wrap around effects
    !! in the FFT due to discontinuities between the left and right (or top and bottom) edges of the domain
    !!
    !!----------------------------------------------------------
    subroutine add_buffer_topo(terrain, buffer_topo, smooth_window, buffer, debug)
        implicit none
        real, dimension(:,:), intent(in) :: terrain
        complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:), intent(inout) :: buffer_topo
        integer, intent(in) :: smooth_window
        integer, intent(in) :: buffer
        logical, intent(in), optional :: debug
        real, dimension(:,:),allocatable :: real_terrain
        integer::nx,ny,i,j,pos, xs,xe,ys,ye, window
        real::weight

        nx=size(terrain,1)+buffer*2
        ny=size(terrain,2)+buffer*2
        allocate(buffer_topo(nx,ny))
        buffer_topo=minval(terrain)

        buffer_topo(1+buffer:nx-buffer,1+buffer:ny-buffer)=terrain
        do i=1,buffer
            weight=i/(real(buffer)*2)
            pos=buffer-i
            buffer_topo(pos+1,1+buffer:ny-buffer)  =terrain(1,1:ny-buffer*2)*(1-weight)+terrain(nx-buffer*2,1:ny-buffer*2)*   weight
            buffer_topo(nx-pos,1+buffer:ny-buffer) =terrain(1,1:ny-buffer*2)*(  weight)+terrain(nx-buffer*2,1:ny-buffer*2)*(1-weight)
        enddo
        do i=1,buffer
            weight=i/(real(buffer)*2)
            pos=buffer-i
            buffer_topo(:,pos+1)  =buffer_topo(:,buffer+1)*(1-weight) + buffer_topo(:,ny-buffer)*   weight
            buffer_topo(:,ny-pos) =buffer_topo(:,buffer+1)*(  weight) + buffer_topo(:,ny-buffer)*(1-weight)
        enddo

        ! smooth the outer most grid cells in all directions to minimize artifacts at the borders of real terrain
        ! smoothing effectively increases as it gets further from the real terrain border (window=min(j,smooth_window))
        if (smooth_window>0) then
            do j=1,buffer
                window=min(j,smooth_window)
                ! smooth the top and bottom borders
                do i=1,nx
                    xs=max(1, i-window)
                    xe=min(nx,i+window)

                    ys=max(1, buffer-j+1-window)
                    ye=min(ny,buffer-j+1+window)

                    buffer_topo(i,buffer-j+1)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))

                    ys=max(1, ny-(buffer-j)-window)
                    ye=min(ny,ny-(buffer-j)+window)

                    buffer_topo(i,ny-(buffer-j))=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))
                end do
                ! smooth the left and right borders
                do i=1,ny
                    xs=max(1, buffer-j+1-window)
                    xe=min(nx,buffer-j+1+window)
                    ys=max(1, i-window)
                    ye=min(ny,i+window)

                    buffer_topo(buffer-j+1,i)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))

                    xs=max(1, nx-(buffer-j)-window)
                    xe=min(nx,nx-(buffer-j)+window)

                    buffer_topo(nx-(buffer-j),i)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))
                end do
            end do
        endif

    end subroutine add_buffer_topo

    !>----------------------------------------------------------
    !! Allocate and initialize arrays in lt_data structure
    !!
    !! Initialize constant arrays in lt_data, e.g. k, l, kl wave number arrays, fftw plans
    !!
    !!----------------------------------------------------------
    subroutine initialize_linear_theory_data(lt_data, nx, ny, dx)
        implicit none
        type(linear_theory_type), intent(inout) :: lt_data
        integer, intent(in) :: nx, ny
        real,    intent(in) :: dx

        integer(C_SIZE_T) :: n_elements
        real :: gain, offset
        integer :: i


        allocate(lt_data%k(nx,ny))
        allocate(lt_data%l(nx,ny))
        allocate(lt_data%kl(nx,ny))
        allocate(lt_data%sig(nx,ny))
        allocate(lt_data%denom(nx,ny))
        allocate(lt_data%m(nx,ny))
        allocate(lt_data%msq(nx,ny))
        allocate(lt_data%mimag(nx,ny))
        allocate(lt_data%ineta(nx,ny))

        ! Compute 2D k and l wavenumber fields
        offset = pi / dx
        gain = 2 * offset / (nx-1)

        ! compute k wave numbers for the first row
        do i=1, nx
            lt_data%k(i,1) = ((i-1) * gain - offset)
        end do
        ! copy k wave numbers to all other rows
        do i=2, ny
            lt_data%k(:,i) = lt_data%k(:,1)
        enddo

        gain = 2 * offset / (ny-1)
        ! compute l wave numbers for the first column
        do i=1,ny
            lt_data%l(1,i) = ((i-1)*gain-offset)
        end do
        ! copy l wave numbers to all other columns
        do i=2,nx
            lt_data%l(i,:) = lt_data%l(1,:)
        enddo

        ! finally compute the kl combination array
        lt_data%kl = lt_data%k**2 + lt_data%l**2
        WHERE (lt_data%kl == 0.0) lt_data%kl = SMALL_VALUE

        ! using fftw_alloc routines to ensure better allignment for vectorization may not be threadsafe
        n_elements = nx * ny
        !$omp critical (fftw_lock)
        lt_data%uh_aligned_data = fftw_alloc_complex(n_elements)
        lt_data%up_aligned_data = fftw_alloc_complex(n_elements)
        lt_data%ua_aligned_data = fftw_alloc_complex(n_elements)
        lt_data%vh_aligned_data = fftw_alloc_complex(n_elements)
        lt_data%vp_aligned_data = fftw_alloc_complex(n_elements)
        lt_data%va_aligned_data = fftw_alloc_complex(n_elements)
        !$omp end critical (fftw_lock)

        call c_f_pointer(lt_data%uh_aligned_data,   lt_data%uhat,           [nx,ny])
        call c_f_pointer(lt_data%up_aligned_data,   lt_data%u_perturb,      [nx,ny])
        call c_f_pointer(lt_data%ua_aligned_data,   lt_data%u_accumulator,  [nx,ny])
        call c_f_pointer(lt_data%vh_aligned_data,   lt_data%vhat,           [nx,ny])
        call c_f_pointer(lt_data%vp_aligned_data,   lt_data%v_perturb,      [nx,ny])
        call c_f_pointer(lt_data%va_aligned_data,   lt_data%v_accumulator,  [nx,ny])

        ! note FFTW plan creation is not threadsafe
        !$omp critical (fftw_lock)
        lt_data%uplan = fftw_plan_dft_2d(ny,nx, lt_data%uhat, lt_data%u_perturb, FFTW_BACKWARD, FFTW_MEASURE) ! alternatives to MEASURE are PATIENT, or ESTIMATE
        lt_data%vplan = fftw_plan_dft_2d(ny,nx, lt_data%vhat, lt_data%v_perturb, FFTW_BACKWARD, FFTW_MEASURE) ! alternatives to MEASURE are PATIENT, or ESTIMATE
        !$omp end critical (fftw_lock)


    end subroutine initialize_linear_theory_data


    !>----------------------------------------------------------
    !! Deallocate lt_data arrays if allocated, includes fftw_* calls where necessary
    !!
    !!----------------------------------------------------------
    subroutine destroy_linear_theory_data(lt_data)
        implicit none
        type(linear_theory_type), intent(inout) :: lt_data

        if (allocated(lt_data%k))       deallocate(lt_data%k)
        if (allocated(lt_data%l))       deallocate(lt_data%l)
        if (allocated(lt_data%kl))      deallocate(lt_data%kl)
        if (allocated(lt_data%sig))     deallocate(lt_data%sig)
        if (allocated(lt_data%denom))   deallocate(lt_data%denom)
        if (allocated(lt_data%m))       deallocate(lt_data%m)
        if (allocated(lt_data%msq))     deallocate(lt_data%msq)
        if (allocated(lt_data%mimag))   deallocate(lt_data%mimag)
        if (allocated(lt_data%ineta))   deallocate(lt_data%ineta)

        !$omp critical (fftw_lock)
        call fftw_free(lt_data%uh_aligned_data)
        call fftw_free(lt_data%up_aligned_data)
        call fftw_free(lt_data%ua_aligned_data)
        call fftw_free(lt_data%vh_aligned_data)
        call fftw_free(lt_data%vp_aligned_data)
        call fftw_free(lt_data%va_aligned_data)
        !$omp end critical (fftw_lock)

        NULLIFY(lt_data%uhat)
        NULLIFY(lt_data%u_perturb)
        NULLIFY(lt_data%u_accumulator)
        NULLIFY(lt_data%vhat)
        NULLIFY(lt_data%v_perturb)
        NULLIFY(lt_data%v_accumulator)

        ! note FFTW plan creation is not threadsafe
        !$omp critical (fftw_lock)
        call fftw_destroy_plan(lt_data%uplan)
        call fftw_destroy_plan(lt_data%vplan)
        !$omp end critical (fftw_lock)
    end subroutine destroy_linear_theory_data


    subroutine setup_remote_grids(u_grids, v_grids, terrain, nz)
        implicit none
        type(grid_t), intent(inout), allocatable :: u_grids(:), v_grids(:)
        real,         intent(in)    :: terrain(:,:)
        integer,      intent(in)    :: nz

        integer :: nx, ny, i

        if (allocated(u_grids)) deallocate(u_grids)
        if (allocated(v_grids)) deallocate(v_grids)

        allocate(u_grids(num_images()))
        allocate(v_grids(num_images()))

        nx = size(terrain, 1)
        ny = size(terrain, 2)

        do i=1,num_images()
            call u_grids(i)%set_grid_dimensions(nx, ny, nz, nx_extra=1, for_image=i)
            call v_grids(i)%set_grid_dimensions(nx, ny, nz, ny_extra=1, for_image=i)
        enddo

    end subroutine setup_remote_grids

    subroutine copy_data_to_remote(wind, grids, LUT, i,j,k, z)
        implicit none
        real,           intent(in)  :: wind(:,:)
        type(grid_t),   intent(in)  :: grids(:)
        real,           intent(inout):: LUT(:,:,:,:,:,:)[*]
        integer,        intent(in)  :: i,j,k, z

        integer :: img

        do img = 1, num_images()
            associate(ims => grids(img)%ims, &
                      ime => grids(img)%ime, &
                      jms => grids(img)%jms, &
                      jme => grids(img)%jme  &
                )
            !$omp critical
            LUT(k,i,j, 1:ime-ims+1, z, 1:jme-jms+1)[img] = wind(ims:ime,jms:jme)
            !$omp end critical

            end associate
        enddo

    end subroutine copy_data_to_remote

    !>----------------------------------------------------------
    !! Compute look up tables for all combinations of U, V, and Nsq
    !!
    !!----------------------------------------------------------
    subroutine initialize_spatial_winds(domain,options,reverse)
        implicit none
        type(domain_t),  intent(inout)::domain
        type(options_t), intent(in) :: options
        logical, intent(in) :: reverse

        ! local variables used to calculate the LUT
        real :: u,v, layer_height, layer_height_bottom, layer_height_top
        integer :: nx,ny,nz, nxu,nyv, i,j,k,z,ik, error
        integer :: fftnx, fftny
        integer, dimension(3,2) :: LUT_dims
        integer :: loops_completed ! this is just used to measure progress in the LUT creation
        integer :: total_LUT_entries, ijk, start_pos, stop_pos
        integer :: ims, jms, this_n
        real, allocatable :: temporary_u(:,:), temporary_v(:,:)
        character(len=kMAX_FILE_LENGTH) :: LUT_file

        type(grid_t), allocatable :: u_grids(:), v_grids(:)

        ! append total number of images and the current image number to the LUT filename
        LUT_file = trim(options%lt_options%u_LUT_Filename) // "_" // trim(str(num_images())) // "_" // trim(str(this_image())) // ".nc"

        ims = lbound(domain%z%data_3d,1)
        jms = lbound(domain%z%data_3d,3)

        ! the domain to work over
        nz = size(domain%u%data_3d,  2)

        call setup_remote_grids(u_grids, v_grids, domain%global_terrain, nz)

        ! ensure these are at their required size for all images
        nx = maxval(v_grids%nx)
        ny = maxval(u_grids%ny)
        nxu = maxval(u_grids%nx)
        nyv = maxval(v_grids%ny)

        fftnx = size(domain%terrain_frequency, 1)
        fftny = size(domain%terrain_frequency, 2)
        ! note:
        ! buffer = (fftnx - nx)/2

        ! default assumes no errors in reading the LUT
        error = 0

        ! store to make it easy to check dim sizes in read_LUT
        LUT_dims(:,1) = [nxu,nz,ny]
        LUT_dims(:,2) = [nx,nz,nyv]

        total_LUT_entries = n_dir_values * n_spd_values * n_nsq_values

        start_pos = nint((real(this_image()-1) / num_images()) * total_LUT_entries)
        if (this_image()==num_images()) then
            stop_pos = total_LUT_entries - 1
        else
            stop_pos  = nint((real(this_image()) / num_images()) * total_LUT_entries) - 1
        endif

        ! create the array of spd, dir, and nsq values to create LUTs for
        ! generates the values for each look up table dimension
        ! generate table of wind directions to be used
        call linear_space(dir_values,dirmin,dirmax,n_dir_values)
        ! generate table of wind speeds to be used
        call linear_space(nsq_values,nsqmin,nsqmax,n_nsq_values)
        ! generate table of Brunt-Vaisalla frequencies (Nsq) to be used
        call linear_space(spd_values,spdmin,spdmax,n_spd_values)

        ! Allocate the (LARGE) look up tables for both U and V
        if (.not.options%lt_options%read_LUT) then
            allocate(hi_u_LUT(n_spd_values, n_dir_values, n_nsq_values, nxu, nz, ny)[*], source=0.0)
            allocate(hi_v_LUT(n_spd_values, n_dir_values, n_nsq_values, nx,  nz, nyv)[*], source=0.0)
            error=0
        else
            if (this_image()==1) write(*,*) "    Reading LUT from file: ", trim(LUT_file)
            error=1
            error = read_LUT(LUT_file, hi_u_LUT, hi_v_LUT, options%parameters%dz_levels(:nz), LUT_dims, options%lt_options)
            if (error/=0) then
                if (this_image()==1) write(*,*) "WARNING: LUT on disk does not match that specified in the namelist or does not exist."
                if (this_image()==1) write(*,*) "    LUT will be recreated"
                if (allocated(hi_u_LUT)) deallocate(hi_u_LUT)
                allocate(hi_u_LUT(n_spd_values, n_dir_values, n_nsq_values, nxu, nz, ny)[*], source=0.0)
                if (allocated(hi_v_LUT)) deallocate(hi_v_LUT)
                allocate(hi_v_LUT(n_spd_values, n_dir_values, n_nsq_values, nx,  nz, nyv)[*], source=0.0)
            endif
        endif

        ! if (options%parameters%debug) then
        if (this_image()==1) write(*,*) "Local Look up Table size:", 4*product(shape(hi_u_LUT))/real(2**20), "MB"
        if (this_image()==1) write(*,*) "Wind Speeds:",spd_values
        if (this_image()==1) write(*,*) "Directions:",360*dir_values/(2*pi)
        if (this_image()==1) write(*,*) "Stabilities:",exp(nsq_values)
        ! endif


        if (reverse.or.(.not.((options%lt_options%read_LUT).and.(error==0)))) then
            ! loop over combinations of U, V, and Nsq values
            loops_completed = 0
            if (this_image()==1) write(*,*) "    Initializing linear theory"
            ! $omp parallel default(shared) &
            ! $omp private(i,j,k,ik,ijk, z, u,v, layer_height, layer_height_bottom, layer_height_top, temporary_u, temporary_v) &
            ! $omp firstprivate(minimum_layer_size, n_dir_values, n_spd_values, n_nsq_values, nz,nx,ny,nxu,nyv,fftnx,fftny) &
            ! $omp firstprivate(start_pos, stop_pos, total_LUT_entries)
            ! $omp threadprivate(lt_data_m) declared at the top of the module

            ! initialization has to happen in each thread so each thread has its own copy
            ! lt_data_m is a threadprivate variable, within initialization, there are omp critical sections for fftw calls
            call initialize_linear_theory_data(lt_data_m, fftnx, fftny, domain%dx)
            allocate(temporary_u(fftnx - buffer*2+1, fftny - buffer*2  ), source=0.0)
            allocate(temporary_v(fftnx - buffer*2,   fftny - buffer*2+1), source=0.0)
            this_n = stop_pos-start_pos+1
            ! $omp do

            do ijk = start_pos, stop_pos
                ! loop over the combined ijk space to improve parallelization (more granular parallelization)
                ! because it is one combined loop, we have to calculate the i,k indicies from the combined ik variable

                ! do ik=0, n_dir_values*n_spd_values*n_nsq_values-1
                ik = ijk / n_nsq_values ! no +1 yet because this still needs to go through another div / mod iteration to compute i and k

                ! do j=1, n_nsq_values
                j = mod(ijk,n_nsq_values) + 1

                ! do i=1, n_spd_values
                i = ik/n_spd_values + 1
                ! do k=1, n_spd_values
                k = mod(ik,n_spd_values) + 1

                ! set the domain wide U and V values to the current u and v values
                u = calc_u( dir_values(i), spd_values(k) )
                v = calc_v( dir_values(i), spd_values(k) )

                debug_print = .False.
                if ((spd_values(k)==15).and.(abs(u-15) < 1).and.(j==15)) debug_print=.True.
                ! calculate the linear wind field for the current u and v values

                do z=1,nz
                    ! print the current status if this is being run "interactively"
                    if (options%parameters%interactive) then
                        !$omp critical (print_lock)
                        if (this_image()==1) write(*,"(f5.1,A)") loops_completed/real(nz*(stop_pos-start_pos+1))*100," %"
                        !$omp end critical (print_lock)
                    endif

                    if (options%parameters%space_varying_dz) then

                        ! call update_irregular_grid(u,v, nsq_values(j), z, domain, minimum_layer_size)
                        call linear_perturbation(u, v, exp(nsq_values(j)),                                                                      &
                                                 domain%global_z_interface(:,z,:) - domain%global_terrain,                                      &
                                                 domain%global_z_interface(:,z,:) - domain%global_terrain + domain%global_dz_interface(:,z,:),  &
                                                 minimum_layer_size, domain%terrain_frequency, lt_data_m)

                    else
                        layer_height = domain%z%data_3d(ims,z,jms) - domain%terrain%data_2d(ims,jms)
                        layer_height_bottom = layer_height - (options%parameters%dz_levels(z) / 2)
                        layer_height_top    = layer_height + (options%parameters%dz_levels(z) / 2)

                        if (z>1) then
                            if (layer_height_bottom /= sum(options%parameters%dz_levels(1:z-1))) then
                                print*, this_image(), layer_height_bottom - sum(options%parameters%dz_levels(1:z-1)), "layer_height_bottom = ", layer_height_bottom, "sum(dz) = ", sum(options%parameters%dz_levels(1:z-1))
                            endif
                            if (layer_height_top /= sum(options%parameters%dz_levels(1:z))) then
                                print*, this_image(), layer_height_top - sum(options%parameters%dz_levels(1:z)), "layer_height_top = ", layer_height_top, "sum(dz) = ", sum(options%parameters%dz_levels(1:z))
                            endif
                        endif

                        call linear_perturbation(u, v, exp(nsq_values(j)),                                  &
                                                 layer_height_bottom, layer_height_top, minimum_layer_size, &
                                                 domain%terrain_frequency, lt_data_m)
                    endif
                    ! need to handle stagger (nxu /= nx) and the buffer around edges of the domain
                    if (nxu /= nx) then
                        temporary_u(:,:) = real( real(                                              &
                                ( lt_data_m%u_perturb(buffer:fftnx-buffer,     1+buffer:fftny-buffer)           &
                                + lt_data_m%u_perturb(1+buffer:fftnx-buffer+1,     1+buffer:fftny-buffer)) )) / 2

                        temporary_v(:,:) = real( real(                                              &
                                ( lt_data_m%v_perturb(1+buffer:fftnx-buffer,     buffer:fftny-buffer)         &
                                + lt_data_m%v_perturb(1+buffer:fftnx-buffer,     1+buffer:fftny-buffer+1)) )) / 2

                        call copy_data_to_remote(temporary_u, u_grids, hi_u_LUT, i,j,k, z)
                        call copy_data_to_remote(temporary_v, v_grids, hi_v_LUT, i,j,k, z)

                    else
                        stop "ERROR: linear wind LUT creation not set up for non-staggered grids yet"
                        ! hi_u_LUT(k,i,j,:,z,:) = real( real(                                                    &
                        !         lt_data_m%u_perturb(1+buffer:fftnx-buffer,     1+buffer:fftny-buffer) ))
                        !
                        ! hi_v_LUT(k,i,j,:,z,:) = real( real(                                                    &
                        !         lt_data_m%v_perturb(1+buffer:fftnx-buffer,     1+buffer:fftny-buffer) ))
                    endif

                    !$omp critical (print_lock)
                    loops_completed = loops_completed+1
                    !$omp end critical (print_lock)

                    ! for now sync all has to be inside the z loop to conserve memory for large domains

                    !$omp critical (print_lock)
                    sync all
                    !$omp end critical (print_lock)
                enddo

            end do
            ! $omp end do

            ! If this image doesn't have as many steps to run as some other images might,
            ! then it has to run additional ghost step(s) so that it syncs with the others correctly
            do i=1, (total_LUT_entries/num_images()+1) - this_n
                ! the syncs should be outside of the z loop, but this is a little more forgiving with memory requirements
                do z=1,nz
                    !$omp critical
                    sync all
                    !$omp end critical
                enddo
            enddo

            ! memory needs to be freed so this structure can be used again when removing linear winds
            call destroy_linear_theory_data(lt_data_m)
            ! $omp end parallel

            sync all

            if (this_image()==1) write(*,*) "All images: 100 % Complete"
            if (this_image()==1) write(*,*) char(10),"--------  Linear wind look up table generation complete ---------"
        endif

        if ((options%lt_options%write_LUT).and.(.not.reverse)) then
            if ((options%lt_options%read_LUT) .and. (error == 0)) then
                if (this_image()==1) write(*,*) "    Not writing Linear Theory LUT to file because LUT was read from file"
            else
                if (this_image()==1) write(*,*) "    Writing Linear Theory LUT to file: ", trim(LUT_file)
                error = write_LUT(LUT_file, hi_u_LUT, hi_v_LUT, options%parameters%dz_levels(:nz), options%lt_options)
            endif
        endif

    end subroutine initialize_spatial_winds


    !>----------------------------------------------------------
    !! Compute a spatially variable linear wind perturbation
    !! based off of look uptables computed in via setup
    !! for each grid point, find the closest LUT data in U and V space
    !! then bilinearly interpolate the nearest LUT values for that points linear wind field
    !!
    !!----------------------------------------------------------
    subroutine spatial_winds(domain,reverse, vsmooth, winsz, update)
        implicit none
        type(domain_t), intent(inout):: domain
        logical,        intent(in)   :: reverse
        integer,        intent(in)   :: vsmooth
        integer,        intent(in)   :: winsz
        logical,        intent(in)   :: update

        logical :: externalN ! flag whether N is read from a field in the forcing dataset
        integer :: nx,nxu, ny,nyv, nz, i,j,k, smoothz
        integer :: uk, vi !store a separate value of i for v and of k for u to we can handle nx+1, ny+1
        integer :: step, dpos, npos, spos, nexts, nextd, nextn
        integer :: north, south, east, west, top, bottom, n
        real :: u, v
        real,allocatable :: u1d(:),v1d(:)
        integer :: ims, ime, jms, jme, ims_u, ime_u, jms_v, jme_v, kms, kme
        real :: dweight, nweight, sweight, curspd, curdir, curnsq, wind_first, wind_second
        real :: blocked, hydrometeors

        ! pointers to the u/v data to be updated so they can point to different places depending on the update flag
        real, pointer :: u3d(:,:,:), v3d(:,:,:), nsquared(:,:,:)


        if (update) then
            u3d => domain%u%meta_data%dqdt_3d
            v3d => domain%v%meta_data%dqdt_3d
        else
            u3d => domain%u%data_3d
            v3d => domain%v%data_3d
        endif
        nsquared => domain%nsquared%data_3d

        ims_u = lbound(u3d,1)
        ime_u = ubound(u3d,1)
        jms = lbound(u3d,3)
        jme = ubound(u3d,3)

        ims = lbound(v3d,1)
        ime = ubound(v3d,1)
        jms_v = lbound(v3d,3)
        jme_v = ubound(v3d,3)

        kms = lbound(u3d,2)
        kme = ubound(u3d,2)

        nx  = size(domain%latitude%data_2d,1)
        ny  = size(domain%latitude%data_2d,2)
        nz  = size(u3d,2)
        nxu = size(u3d,1)
        nyv = size(v3d,3)

        if (reverse) then
            ! u_LUT=>rev_u_LUT
            ! v_LUT=>rev_v_LUT
            u_perturbation=>lo_u_perturbation
            v_perturbation=>lo_v_perturbation
        else
            ! u_LUT=>hi_u_LUT
            ! v_LUT=>hi_v_LUT
            u_perturbation=>hi_u_perturbation
            v_perturbation=>hi_v_perturbation
        endif

        ! if (reverse) print*, "WARNING using fixed nsq for linear wind removal: 3e-6"
        ! !$omp parallel firstprivate(nx,nxu,ny,nyv,nz, kms, kme, reverse, vsmooth, winsz, using_blocked_flow), default(none), &
        ! !$omp private(i,j,k,step, uk, vi, east, west, north, south, top, bottom, u1d, v1d), &
        ! !$omp private(spos, dpos, npos, nexts,nextd, nextn,n, smoothz, u, v, blocked), &
        ! !$omp private(wind_first, wind_second, curspd, curdir, curnsq, sweight,dweight, nweight), &
        ! !$omp shared(domain, u3d,v3d, spd_values, dir_values, nsq_values, u_LUT, v_LUT, linear_mask), &
        ! !$omp shared(u_perturbation, v_perturbation, linear_update_fraction, linear_contribution, nsq_calibration), &
        ! !$omp shared(min_stability, max_stability, n_dir_values, n_spd_values, n_nsq_values, smooth_nsq)
        !
        ! !$omp do
        do k=1,ny

            do j=1,nz
                do i=1,nx

                    if (variable_N) then
                        ! look up vsmooth gridcells up to nz at the maximum
                        top = min(j+vsmooth, nz)
                        ! if (top-j)/=vsmooth, then look down enough layers to make the window vsmooth in size
                        bottom = max(1, j - (vsmooth - (top-j)))

                        hydrometeors = 0
                        if (associated(domain%cloud_water_mass%data_3d)) &
                            hydrometeors = hydrometeors + domain%cloud_water_mass%data_3d(i+ims-1,j+kms-1,k+jms-1)
                        if (associated(domain%cloud_ice_mass%data_3d)) &
                            hydrometeors = hydrometeors + domain%cloud_ice_mass%data_3d(i+ims-1,j+kms-1,k+jms-1)
                        if (associated(domain%rain_mass%data_3d)) &
                            hydrometeors = hydrometeors + domain%rain_mass%data_3d(i+ims-1,j+kms-1,k+jms-1)
                        if (associated(domain%snow_mass%data_3d)) &
                            hydrometeors = hydrometeors + domain%snow_mass%data_3d(i+ims-1,j+kms-1,k+jms-1)


                        if (.not.reverse) then
                            nsquared(i+ims-1,j+kms-1,k+jms-1) = calc_stability(domain%potential_temperature%data_3d(i+ims-1,bottom+kms-1,k+jms-1),                              &
                                                             domain%potential_temperature%data_3d(i+ims-1,top+kms-1,k+jms-1),                                                   &
                                                             domain%exner%data_3d(i+ims-1,bottom+kms-1,k+jms-1),domain%exner%data_3d(i+ims-1,top+kms-1,k+jms-1),                &
                                                             domain%z%data_3d(i+ims-1,bottom+kms-1,k+jms-1),  domain%z%data_3d(i+ims-1,top+kms-1,k+jms-1),                      &
                                                             domain%water_vapor%data_3d(i+ims-1,bottom+kms-1,k+jms-1), domain%water_vapor%data_3d(i+ims-1,top+kms-1,k+jms-1),   &
                                                             hydrometeors)

                            nsquared(i+ims-1,j+kms-1,k+jms-1) = max(min_stability, min(max_stability, &
                                                    nsquared(i+ims-1,j+kms-1,k+jms-1)))! * nsq_calibration(i,k)))
                        else
                            ! Low-res boundary condition variables will be in a different array format.  It should be
                            ! easy enough to call calc_stability after e.g. transposing z and y dimension, but some
                            ! e.g. pii will not be set in the forcing data, so this may need a little thought.
                            nsquared(i+ims-1,j+kms-1,k+jms-1) = 3e-6
                        endif
                    else
                        nsquared(i+ims-1,j+kms-1,k+jms-1) = N_squared
                    endif
                end do
                ! look up table is computed in log space
                nsquared(:,j+kms-1,k+jms-1) = log(nsquared(:,j+kms-1,k+jms-1))
            end do

            if (smooth_nsq) then
                do j=1,nz
                    ! compute window as above.
                    top = min(j+vsmooth,nz)
                    bottom = max(1, j - (vsmooth - (top-j)) )

                    do smoothz = bottom, j-1
                        nsquared(:,j+kms-1,k+jms-1) = nsquared(:,j+kms-1,k+jms-1) + nsquared(:,smoothz+kms-1,k+jms-1)
                    end do
                    do smoothz = j+1, top
                        nsquared(:,j+kms-1,k+jms-1) = nsquared(:,j+kms-1,k+jms-1) + nsquared(:,smoothz+kms-1,k+jms-1)
                    end do
                    nsquared(:,j+kms-1,k+jms-1) = nsquared(:,j+kms-1,k+jms-1)/(top-bottom+1)
                end do
            endif
        end do
        ! !$omp end do
        ! !$omp end parallel


        ! smooth array has it's own parallelization, so this probably can't go in a critical section
        if (smooth_nsq) then
            call smooth_array(nsquared, winsz, ydim=3)
        endif

        !$omp parallel firstprivate(ims, ime, jms, jme, ims_u, ime_u, jms_v, jme_v, kms, kme, nx,nxu,ny,nyv,nz), &
        !$omp firstprivate(reverse, vsmooth, winsz, using_blocked_flow), default(none), &
        !$omp private(i,j,k,step, uk, vi, east, west, north, south, top, bottom, u1d, v1d), &
        !$omp private(spos, dpos, npos, nexts,nextd, nextn,n, smoothz, u, v, blocked), &
        !$omp private(wind_first, wind_second, curspd, curdir, curnsq, sweight,dweight, nweight), &
        !$omp shared(domain, u3d,v3d, nsquared, spd_values, dir_values, nsq_values, hi_u_LUT, hi_v_LUT, linear_mask), &
        !$omp shared(u_perturbation, v_perturbation, linear_update_fraction, linear_contribution, nsq_calibration), &
        !$omp shared(min_stability, max_stability, n_dir_values, n_spd_values, n_nsq_values, smooth_nsq)
        allocate(u1d(nxu), v1d(nxu))
        !$omp do
        do k=1, nyv

            uk = min(k,ny)
            do i=1,nxu
                vi = min(i,nx)
                u1d(i)  = sum(u3d(i+ims_u-1,:,uk+jms-1)) / nz
                v1d(i)  = sum(v3d(vi+ims-1, :,k+jms_v-1)) / nz
            enddo

            do j=1, nz
                do i=1, nxu

                    ! if (using_blocked_flow) then
                    !     blocked = blocking_fraction(domain%froude(min(i,nx),min(k,ny)))
                    ! else
                        blocked = 0
                    ! endif
                    if (blocked<1) then

                        uk = min(k,ny)
                        vi = min(i,nx)

                        !   First find the bounds of the region to average over
                        west  = max(i - winsz, 1)
                        east  = min(i + winsz,nx)
                        bottom= max(j - winsz, 1)
                        top   = min(j + winsz,nz)
                        south = max(k - winsz, 1)
                        north = min(k + winsz,ny)

                        n = (((east-west)+1) * ((north-south)+1))
                        if (reverse) then
                            u = u3d(i+ims_u-1,j,uk+jms-1 ) ! u3d(i,j,uk)
                            v = v3d(vi+ims-1, j,k+jms_v-1) ! v3d(vi,j,k)
                            ! WARNING: see below for why this does not work (yet)
                            ! u = sum( domain%u(west:east,j,south:north) ) / n
                            ! v = sum( domain%v(west:east,j,south:north) ) / n
                        else
                            ! smooth the winds vertically first
                            u = u1d(i)
                            v = v1d(i)
                            ! WARNING for now domain%u,v are updated inside the loop, so the code below will end
                            ! up using u and v values that have had the linear theory applied already (not good)
                            ! eventually this should be pulled out and only the pertubration should be updated in the
                            ! parallel region
                            ! u = sum(domain%u( i, bottom:top,uk)) / (top-bottom+1)
                            ! v = sum(domain%v(vi, bottom:top, k)) / (top-bottom+1)
                        endif
                        n = n * ((top-bottom)+1)

                        ! Calculate the direction of the current grid cell wind
                        dpos = 1
                        curdir = calc_direction( u, v )
                        ! and find the corresponding position in the Look up Table
                        do step=1, n_dir_values
                            if (curdir > dir_values(step)) then
                                dpos = step
                            endif
                        end do

                        ! Calculate the wind speed of the current grid cell
                        spos = 1
                        curspd = calc_speed( u, v )
                        ! and find the corresponding position in the Look up Table
                        do step=1, n_spd_values
                            if (curspd > spd_values(step)) then
                                spos = step
                            endif
                        end do

                        ! Calculate the Brunt-Vaisalla frequency of the current grid cell
                        !   Then compute the mean in log space
                        ! curnsq = log(2.75e-5)
                        curnsq = sum(nsquared(vi+ims-1,bottom+kms-1:top+kms-1,uk+jms-1)) / (top - bottom + 1)
                        !   and find the corresponding position in the Look up Table
                        npos = 1
                        do step=1, n_nsq_values
                            if (curnsq > nsq_values(step)) then
                                npos = step
                            endif
                        end do

                        ! Calculate the weights and the "next" u/v position
                        ! "next" usually = pos+1 but for edge cases next = 1 or n
                        dweight = calc_weight(dir_values, dpos, nextd, curdir)
                        sweight = calc_weight(spd_values, spos, nexts, curspd)
                        nweight = calc_weight(nsq_values, npos, nextn, curnsq)

                        ! perform linear interpolation between LUT values
                        if (k<=ny) then
                            wind_first =      nweight  * (dweight * hi_u_LUT(spos, dpos,npos, i,j,k) + (1-dweight) * hi_u_LUT(spos, nextd,npos, i,j,k))   &
                                        +  (1-nweight) * (dweight * hi_u_LUT(spos, dpos,nextn,i,j,k) + (1-dweight) * hi_u_LUT(spos, nextd,nextn,i,j,k))

                            wind_second=      nweight  * (dweight * hi_u_LUT(nexts,dpos,npos, i,j,k) + (1-dweight) * hi_u_LUT(nexts,nextd,npos, i,j,k))   &
                                        +  (1-nweight) * (dweight * hi_u_LUT(nexts,dpos,nextn,i,j,k) + (1-dweight) * hi_u_LUT(nexts,nextd,nextn,i,j,k))

                            u_perturbation(i,j,k) = u_perturbation(i,j,k) * (1-linear_update_fraction) &
                                        + linear_update_fraction * (sweight*wind_first + (1-sweight)*wind_second)

                            ! if (reverse) then
                            !     u3d(i,j,k) = u3d(i,j,k) - u_perturbation(i,j,k) * linear_contribution
                            ! else
                                u3d(i+ims_u-1,j,k+jms-1) = u3d(i+ims_u-1,j,k+jms-1) + u_perturbation(i,j,k) *linear_contribution ! * linear_mask(min(nx,i),min(ny,k)) * (1-blocked)
                            ! endif
                        endif
                        if (i<=nx) then
                            wind_first =      nweight  * (dweight * hi_v_LUT(spos, dpos,npos, i,j,k) + (1-dweight) * hi_v_LUT(spos, nextd,npos, i,j,k))    &
                                        +  (1-nweight) * (dweight * hi_v_LUT(spos, dpos,nextn,i,j,k) + (1-dweight) * hi_v_LUT(spos, nextd,nextn,i,j,k))

                            wind_second=      nweight  * (dweight * hi_v_LUT(nexts,dpos,npos, i,j,k) + (1-dweight) * hi_v_LUT(nexts,nextd,npos, i,j,k))    &
                                        +  (1-nweight) * (dweight * hi_v_LUT(nexts,dpos,nextn,i,j,k) + (1-dweight) * hi_v_LUT(nexts,nextd,nextn,i,j,k))

                            v_perturbation(i,j,k) = v_perturbation(i,j,k) * (1-linear_update_fraction) &
                                        + linear_update_fraction * (sweight*wind_first + (1-sweight)*wind_second)

                            ! if (reverse) then
                            !     v3d(i,j,k) = v3d(i,j,k) - v_perturbation(i,j,k) * linear_contribution
                            ! else
                                ! for the high res domain, linear_mask should incorporate linear_contribution
                                v3d(i+ims-1,j,k+jms_v-1) = v3d(i+ims-1,j,k+jms_v-1) + v_perturbation(i,j,k) *linear_contribution! * linear_mask(min(nx,i),min(ny,k)) * (1-blocked)
                            ! endif
                        endif

                    endif
                end do
            end do
        end do
        !$omp end do

        deallocate(u1d, v1d)
        !$omp end parallel
        nsquared = exp(nsquared)

    end subroutine spatial_winds


    !>----------------------------------------------------------
    !! Copy options from the options data structure into module level variables
    !!
    !!----------------------------------------------------------
    subroutine set_module_options(options)
        implicit none
        type(options_t), intent(in) :: options

        original_buffer       = options%lt_options%buffer
        variable_N            = options%lt_options%variable_N
        smooth_nsq            = options%lt_options%smooth_nsq

        stability_window_size = options%lt_options%stability_window_size        ! window to average nsq over
        max_stability         = options%lt_options%max_stability                ! limits on the calculated Brunt Vaisala Frequency
        min_stability         = options%lt_options%min_stability                ! these may need to be a little narrower.
        linear_contribution   = options%lt_options%linear_contribution          ! multiplier on uhat,vhat before adding to u,v
                                                                                ! =fractional contribution of linear perturbation to wind field
        ! these are defined per call to permit different values for low res and high res domains
        ! N_squared              = options%lt_options%N_squared                 ! static Brunt Vaisala Frequency (N^2) to use
        ! rm_N_squared           = options%lt_options%rm_N_squared              ! static BV Frequency (N^2) to use in removing linear wind field
        ! rm_linear_contribution = options%lt_options%rm_linear_contribution    ! fractional contribution of linear perturbation to remove wind field
        ! remove_lowres_linear   = options%lt_options%remove_lowres_linear      ! remove the linear mountain wave from low res forcing model

        linear_update_fraction    = options%lt_options%linear_update_fraction   ! controls the rate at which the linearfield updates
                                                                                ! =fraction of linear perturbation to add each time step
        use_spatial_linear_fields = options%lt_options%spatial_linear_fields    ! use a spatially varying linear wind perturbation
        use_linear_mask           = options%lt_options%linear_mask              ! use a spatial mask for the linear wind field
        use_nsq_calibration       = options%lt_options%nsq_calibration          ! use a spatial mask to calibrate the nsquared (brunt vaisala frequency) field

        using_blocked_flow        = options%block_options%block_flow

        ! Look up table generation parameters, range for each parameter, and number of steps to cover that range
        dirmax = options%lt_options%dirmax
        dirmin = options%lt_options%dirmin
        spdmax = options%lt_options%spdmax
        spdmin = options%lt_options%spdmin
        nsqmax = options%lt_options%nsqmax
        nsqmin = options%lt_options%nsqmin
        n_dir_values = options%lt_options%n_dir_values
        n_nsq_values = options%lt_options%n_nsq_values
        n_spd_values = options%lt_options%n_spd_values
        minimum_layer_size = options%lt_options%minimum_layer_size

    end subroutine set_module_options

    !>----------------------------------------------------------
    !! Called from linear_perturb the first time perturb is called
    !! compute FFT(terrain), and dzdx,dzdy components
    !!
    !!----------------------------------------------------------
    subroutine setup_linwinds(domain, options, reverse, useDensity)
        implicit none
        type(domain_t),     intent(inout)  :: domain
        type(options_t),    intent(in)     :: options
        logical,            intent(in)     :: reverse, useDensity
        ! locals
        complex(C_DOUBLE_COMPLEX), allocatable  :: complex_terrain_firstpass(:,:)
        complex(C_DOUBLE_COMPLEX), allocatable  :: complex_terrain(:,:)
        type(C_PTR) :: plan
        integer     :: nx, ny, nz
        real        :: saved_linear_contribution ! temporary variable to save the linear contribution state so we can restore
        integer     :: save_buffer

        ! store module level variables so we don't have to pass options through everytime
        ! lots of these things probably need to be moved to the linearizable class so they
        ! can be separated for the domain and bc fields
        ! this is a little tricky, because we don't want to have to calculate the LUTs
        ! twice, once for domain and once for bc%next_domain
        call set_module_options(options)

        ! Create a buffer zone around the topography to smooth the edges
        buffer = original_buffer
        ! first create it including a 5 grid cell smoothing function
        call add_buffer_topo(domain%global_terrain, complex_terrain_firstpass, 5, buffer)
        buffer = 2
        ! then further add a small (~2) grid cell buffer where all cells have the same value
        call add_buffer_topo(real(real(complex_terrain_firstpass)), complex_terrain, 0, buffer, debug=options%parameters%debug)
        buffer = buffer + original_buffer

        nx = size(complex_terrain, 1)
        ny = size(complex_terrain, 2)

        if (this_image()==1) write(*,*) "Initializing linear winds"
        allocate(domain%terrain_frequency(nx,ny))

        ! calculate the fourier transform of the terrain for use in linear winds
        plan = fftw_plan_dft_2d(ny, nx, complex_terrain, domain%terrain_frequency, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, complex_terrain, domain%terrain_frequency)
        call fftw_destroy_plan(plan)
        ! normalize FFT by N - grid cells
        domain%terrain_frequency = domain%terrain_frequency / (nx * ny)
        ! shift the grid cell quadrants
        ! need to test what effect all of the related shifts actually have...
        call fftshift(domain%terrain_frequency)

        if (linear_contribution/=1) then
            if (this_image()==1) write(*,*) "  Using a fraction of the linear perturbation:",linear_contribution
        endif

        nx = size(domain%global_terrain, 1)
        nz = size(domain%u%data_3d,       2)
        ny = size(domain%global_terrain, 2)


        ! set up linear_mask variable
        if (.not.reverse) then
            if (allocated(linear_mask)) deallocate(linear_mask)
            allocate(linear_mask(nx,ny))
            linear_mask = linear_contribution

            ! if (use_linear_mask) then
            !     if (this_image()==1) write(*,*) "  Reading Linear Mask"
            !     if (this_image()==1) write(*,*) "    from file: " // trim(options%linear_mask_file)
            !     if (this_image()==1) write(*,*) "    with var: "  // trim(options%linear_mask_var)
            !     call io_read(options%linear_mask_file, options%linear_mask_var, domain%linear_mask)
            !
            !     linear_mask = domain%linear_mask * linear_contribution
            ! endif

            ! Stupidly simple adjustment to the linear wind field to account for using density.
            ! If we are using density in the advection calculations, modify the linear perturbation
            ! to get the vertical velocities closer to what they would be without density (boussinesq)
            ! need to check if this makes the most sense when close to the surface
            ! if (useDensity) then
                ! linear_mask = linear_mask * 2
            ! endif
        endif

        ! set up nsq_calibration variable
        if ( (.not.reverse) .and. (.not.allocated(nsq_calibration)) ) then
            allocate(nsq_calibration(nx,ny))
            nsq_calibration = 1

            ! if (use_nsq_calibration) then
            !     if (this_image()==1) write(*,*) "  Reading Linear Mask"
            !     if (this_image()==1) write(*,*) "    from file: " // trim(options%nsq_calibration_file)
            !     if (this_image()==1) write(*,*) "    with var: "  // trim(options%nsq_calibration_var)
            !     call io_read(options%nsq_calibration_file, options%nsq_calibration_var, domain%nsq_calibration)
            !     nsq_calibration = domain%nsq_calibration
            !
            !     where(nsq_calibration<1) nsq_calibration = 1 + 1/( (1-1/nsq_calibration)/100 )
            !     where(nsq_calibration>1) nsq_calibration = 1 + (nsq_calibration-1)/100
            ! endif
        endif

        ! allocate the fields that will hold the perturbation only so we can update it
        ! slowly and add the total to the domain%u,v
        if (reverse) then
            if (.not.allocated(lo_u_perturbation)) then
                allocate(lo_u_perturbation(nx+1,nz,ny))
                lo_u_perturbation=0
                allocate(lo_v_perturbation(nx,nz,ny+1))
                lo_v_perturbation=0
            endif
            u_perturbation=>lo_u_perturbation
            v_perturbation=>lo_v_perturbation
        else
            if (.not.allocated(hi_u_perturbation)) then
                allocate(hi_u_perturbation(nx+1,nz,ny))
                hi_u_perturbation=0
                allocate(hi_v_perturbation(nx,nz,ny+1))
                hi_v_perturbation=0
            endif
            u_perturbation=>hi_u_perturbation
            v_perturbation=>hi_v_perturbation
        endif

        if (use_spatial_linear_fields) then
            ! if    ((.not.allocated(hi_u_LUT)  .and. (.not.reverse)) &
            !  .or.  (.not.allocated(rev_u_LUT) .and. reverse)) then

                if (this_image()==1) write(*,*) "  Generating a spatially variable linear perturbation look up table"
                call initialize_spatial_winds(domain, options, reverse)

            ! endif
        endif

        module_initialized = .True.

    end subroutine setup_linwinds

    !>----------------------------------------------------------
    !! Initialize and/or apply linear wind solution.
    !!
    !! Called from ICAR to update the U and V wind fields based on linear theory (W is calculated to balance U/V)
    !!
    !!----------------------------------------------------------
    subroutine linear_perturb(domain,options,vsmooth,reverse,useDensity, update)
        implicit none
        type(domain_t),     intent(inout):: domain
        type(options_t),    intent(in)   :: options
        integer,            intent(in)   :: vsmooth
        logical,            intent(in),  optional :: reverse, useDensity, update

        logical :: rev, useD, updt
        ! logical, save :: debug=.True.
        real :: stability

        ! these probably need to get moved to options...
        rev = reverse
        if (present(reverse)) rev = .False.

        useD = .False.
        if (present(useDensity)) useD=useDensity

        updt = .False.
        if (present(update)) updt = update

        ! this is a little trickier, because it does have to be domain dependant... could at least be stored in the domain though...
        if (rev) then
            linear_contribution = options%lt_options%rm_linear_contribution
            N_squared = options%lt_options%rm_N_squared
        else
            linear_contribution = options%lt_options%linear_contribution
            N_squared = options%lt_options%N_squared
        endif

        ! if linear_perturb hasn't been called before we need to perform some setup actions.
        if (.not. module_initialized) then
            call setup_linwinds(domain, options, rev, useD)
        endif

        ! add the spatially variable linear field
        ! if we are reverseing the effects, that means we are in the low-res domain
        ! that domain does not have a spatial LUT calculated, so it can not be performed

        ! if (use_spatial_linear_fields)then
            call spatial_winds(domain,rev, vsmooth, stability_window_size, update=updt)
        ! else
        !     ! Nsq = squared Brunt Vaisalla frequency (1/s) typically from dry static stability
        !     stability = calc_domain_stability(domain)
        !     ! This should probably be called twice, once for dry, and once or moist regions
        !     call linear_winds(domain,stability,vsmooth,rev,useD,debug)
        ! endif
        ! debug=.False.

    end subroutine linear_perturb
end module linear_theory_winds
