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
    use fft ! note fft module is defined in fftshift.f90
    use fftshifter
    use data_structures
    use io_routines,               only : io_write2d, io_read2d
    use string,                    only : str
    use linear_theory_lut_disk_io, only : read_LUT, write_LUT
    implicit none
    private
    public::linear_perturb

    logical :: variable_N
    logical :: smooth_nsq
    real    :: N_squared

    !! unfortunately these have to be allocated every call because we could be calling on both the high-res
    !! domain and the low res domain (to "remove" the linear winds)
    real,                       allocatable,    dimension(:,:)  :: k, l, kl, sig
    complex(C_DOUBLE_COMPLEX),  allocatable,    dimension(:,:)  :: denom, m, ineta, msq, mimag
    complex(C_DOUBLE_COMPLEX),  pointer,        dimension(:,:,:):: uhat, u_hat, vhat, v_hat
    type(C_PTR),                allocatable,    dimension(:)    :: uplans, vplans
    type(C_PTR)              :: uh_aligned_data, u_h_aligned_data, vh_aligned_data, v_h_aligned_data
    logical                  :: data_allocated=.False. ! need a boolean because you cant test if a cptr is allocated?

    ! note this data structure holds a lot of important variables for the new linear_perturbation subroutine
    ! this is created at the module scope to prevent some apparent stack issues related to OpenMP and ifort
    type(linear_theory_type) :: lt_data_m
    !$omp threadprivate(lt_data_m)

    !! Linear wind look up table values and tables
    real, allocatable,          dimension(:)           :: dir_values, nsq_values, spd_values
    ! Look Up Tables for linear perturbation are nspd x n_dir_values x n_nsq_values x nx x nz x ny
    real, allocatable, target,  dimension(:,:,:,:,:,:) :: hi_u_LUT, hi_v_LUT, rev_u_LUT, rev_v_LUT
    real, pointer,              dimension(:,:,:,:,:,:) :: u_LUT, v_LUT
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

    integer :: n_dir_values=36
    integer :: n_nsq_values=10
    integer :: n_spd_values=10

    complex,parameter :: j= (0,1)

    real, parameter :: SMALL_VALUE = 1e-15

contains

    !>----------------------------------------------------------
    !! Calculate direction [0-2*pi) from u and v wind speeds
    !!
    !!----------------------------------------------------------
    pure function calc_direction(u,v) result(direction)
        implicit none
        real, intent(in) :: u,v
        real :: direction

        if (v<0) then
            direction = atan(u/v) + pi
        elseif (v==0) then
            if (u>0) then
                direction=pi/2.0
            else
                direction=pi*1.5
            endif
        else
            if (u>=0) then
                direction = atan(u/v)
            else
                direction = atan(u/v) + (2*pi)
            endif
        endif

    end function calc_direction

    !>----------------------------------------------------------
    !! Calculate the strength of the u wind field given a direction [0-2*pi] and magnitude
    !!
    !!----------------------------------------------------------
    pure function calc_speed(u, v) result(speed)
        implicit none
        real, intent(in) :: u,v
        real :: speed

        speed = sqrt(u**2 + v**2)
    end function calc_speed

    !>----------------------------------------------------------
    !! Calculate the strength of the u wind field given a direction [0-2*pi] and magnitude
    !!
    !!----------------------------------------------------------
    pure function calc_u(direction, magnitude) result(u)
        implicit none
        real, intent(in) :: direction, magnitude
        real :: u

        u = sin(direction) * magnitude
    end function calc_u

    !>----------------------------------------------------------
    !! Calculate the strength of the v wind field given a direction [0-2*pi] and magnitude
    !!
    !!----------------------------------------------------------
    pure function calc_v(direction, magnitude) result(v)
        implicit none
        real, intent(in) :: direction, magnitude
        real :: v

        v = cos(direction) * magnitude
    end function calc_v

    !>----------------------------------------------------------
    !! Calculate the saturated adiabatic lapse rate from a T/Moisture input
    !!
    !! return the moist / saturated adiabatic lapse rate for a given
    !! Temperature and mixing ratio (really MR could be calculated as f(T))
    !! from http://glossary.ametsoc.org/wiki/Saturation-adiabatic_lapse_rate
    !!
    !!----------------------------------------------------------
    pure function calc_sat_lapse_rate(T,mr) result(sat_lapse)
        implicit none
        real, intent(in) :: T,mr  ! inputs T in K and mr in kg/kg
        real :: L
        real :: sat_lapse

        L=LH_vaporization ! short cut for imported parameter
        sat_lapse = gravity*((1 + (L*mr) / (Rd*T))          &
                    / (cp + (L*L*mr*(Rd/Rw)) / (Rd*T*T) ))
    end function calc_sat_lapse_rate

    !>----------------------------------------------------------
    !! Calculate the moist brunt vaisala frequency (Nm^2)
    !! formula from Durran and Klemp, 1982 after Lalas and Einaudi 1974
    !!
    !!----------------------------------------------------------
    pure function calc_moist_stability(t_top, t_bot, z_top, z_bot, qv_top, qv_bot, qc) result(BV_freq)
        implicit none
        real, intent(in) :: t_top, t_bot, z_top, z_bot, qv_top, qv_bot, qc
        real :: t,qv, dz, sat_lapse
        real :: BV_freq

        t  = ( t_top +  t_bot)/2
        qv = (qv_top + qv_bot)/2
        dz = ( z_top - z_bot)
        sat_lapse = calc_sat_lapse_rate(t,qv)

        BV_freq = (gravity/t) * ((t_top-t_bot)/dz + sat_lapse) * &
                  (1 + (LH_vaporization*qv)/(Rd*t)) - (gravity/(1+qv+qc) * (qv_top-qv_bot)/dz)
    end function calc_moist_stability

    !>----------------------------------------------------------
    !! Calculate the dry brunt vaisala frequency (Nd^2)
    !!
    !!----------------------------------------------------------
    pure function calc_dry_stability(th_top, th_bot, z_top, z_bot) result(BV_freq)
        implicit none
        real, intent(in) :: th_top, th_bot, z_top, z_bot
        real :: BV_freq

        BV_freq = gravity * (log(th_top)-log(th_bot)) / (z_top - z_bot)
    end function calc_dry_stability

    !>----------------------------------------------------------
    !! Calculate either moist or dry brunt vaisala frequency and keep within min and max bounds
    !!
    !!----------------------------------------------------------
    pure function calc_stability(th_top, th_bot, pii_top, pii_bot, z_top, z_bot, qv_top, qv_bot, qc) result(BV_freq)
        implicit none
        real, intent(in) :: th_top, th_bot, pii_top, pii_bot, z_top, z_bot, qv_top, qv_bot, qc
        real :: BV_freq

        if (qc<1e-7) then
            if (variable_N) then
                BV_freq=calc_dry_stability(th_top, th_bot, z_top, z_bot)
            else
                BV_freq=N_squared
            endif
        else
            if (variable_N) then
                BV_freq=calc_moist_stability(th_top*pii_top, th_bot*pii_bot, z_top, z_bot, qv_top, qv_bot, qc)
            else
                BV_freq=N_squared/10.0 ! might be better as max(1e-7,N_squared-(1e-4))
            endif
        endif

        BV_freq = min(max(BV_freq,min_stability),max_stability)
    end function calc_stability

    !>----------------------------------------------------------
    !! Calculate a single brunt vaisala frequency effectively averaged across the entire domain
    !!
    !!----------------------------------------------------------
    pure function calc_domain_stability(domain) result(BV_freq)
        implicit none
        class(linearizable_type),intent(in)::domain
        real :: BV_freq
        real, allocatable, dimension(:) :: vertical_N
        integer :: i, nx,nz,ny, top_layer,bottom_layer
        real :: dz,dtheta,t_mean

        nx=size(domain%th,1)
        nz=size(domain%th,2)
        ny=size(domain%th,3)
        allocate(vertical_N(nz))

        if (variable_N) then
            do i=1,nz
                top_layer    = min(i+stability_window_size, nz)
                bottom_layer = max(i-stability_window_size,  1)

                dz     = sum(domain%z (:,top_layer,:)-domain%z (:,bottom_layer,:))/(nx*ny) ! distance between top layer and bottom layer
                dtheta = sum(domain%th(:,top_layer,:)-domain%th(:,bottom_layer,:))/(nx*ny) ! average temperature change between layers
                t_mean = sum(domain%th(:,i,:))/(nx*ny) ! mean temperature in the current layer

                ! BV frequency = gravity/theta * dtheta/dz
                vertical_N(i)=gravity/t_mean * dtheta/dz
            end do
            ! calculate the mean stability over the entire profile
            BV_freq=sum(vertical_N)/nz
            ! impose limits so the linear solution doesn't go crazy in really stable or unstable air
            BV_freq=min(max(BV_freq, min_stability), max_stability)
        else
            ! or just use the supplied BV frequency
            BV_freq=N_squared
        endif

        deallocate(vertical_N)
    end function calc_domain_stability

    !>----------------------------------------------------------
    !! Compute linear wind perturbations given background U, V, Nsq
    !!
    !! see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
    !!
    !!----------------------------------------------------------
    subroutine linear_perturbation(U, V, Nsq, z, fourier_terrain, lt_data)
        use, intrinsic :: iso_c_binding
        implicit none
        real,                     intent(in)    :: U, V ! U and V components of background wind
        real,                     intent(in)    :: Nsq  ! Brunt-Vaisalla frequency (N squared)
        real,                     intent(in)    :: z    ! elevation at which to compute the linear solution
        complex(C_DOUBLE_COMPLEX),intent(in)    :: fourier_terrain(:,:) ! FFT(terrain)
        type(linear_theory_type), intent(inout) :: lt_data


        integer :: nx, ny ! store the size of the grid

        nx = size(fourier_terrain, 1)
        ny = size(fourier_terrain, 2)

        lt_data%sig  = U * lt_data%k + V * lt_data%l
        ! prevent divide by zero
        where(lt_data%sig == 0) lt_data%sig = SMALL_VALUE

        lt_data%denom = lt_data%sig**2 ! -f**2 ! to add coriolis

        lt_data%msq   = Nsq / lt_data%denom * lt_data%kl
        lt_data%mimag = 0 + 0 * j ! be sure to reset real and imaginary components
        lt_data%mimag = lt_data%mimag + (real(sqrt(-lt_data%msq)) * j)

        lt_data%m = sqrt(lt_data%msq)                         ! vertical wave number, hydrostatic
        where(lt_data%sig < 0) lt_data%m = lt_data%m * (-1)   ! equivalent to m = m * sign(sig)
        where(real(lt_data%msq) < 0) lt_data%m = lt_data%mimag

        lt_data%ineta = j * fourier_terrain * exp(j * lt_data%m * z)
        !  what=sig*ineta

        ! with coriolis : [+/-]j*[l/k]*f
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
    end subroutine linear_perturbation

    !>----------------------------------------------------------
    !! Compute linear wind perturbations to U and V and add them back to the domain
    !!
    !! see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
    !!
    !!----------------------------------------------------------
    subroutine linear_winds(domain, Nsq, vsmooth, reverse, useDensity, debug)
        use, intrinsic :: iso_c_binding
        implicit none
        class(linearizable_type),intent(inout) :: domain
        real,    intent(in) :: Nsq      !> Input brunt-vaisalla frequency to use
        integer, intent(in) :: vsmooth  !> Number of vertical layers to smooth winds over
        logical, intent(in) :: reverse  !> reverse the linear solution (to "remove" linear wind effects)
        logical, intent(in) :: useDensity !> Specify that we are using density in the advection calculations
                                        !> if useDensity: use hueristics to modify the linear perturbation accordingly
        logical, intent(in) :: debug    !> Debug flag to print additional information

        ! Local variables
        real    :: gain, offset !used in setting up k and l arrays
        integer :: nx, ny, nz, i, z, realnx, realny, realnx_u, realny_v, bottom,top
        logical :: staggered, zaxis_is_third
        real    :: U, V         ! used to store the domain wide average values
        real,dimension(:),allocatable :: U_layers, V_layers, preU_layers, preV_layers
        integer(C_SIZE_T) :: n_elements

        nx=size(domain%fzs,1)
        nz=size(domain%u,2)
        ny=size(domain%fzs,2)
        allocate(U_layers(nz))
        allocate(V_layers(nz))
        allocate(preU_layers(nz))
        allocate(preV_layers(nz))

        realnx=size(domain%z,1)
        realnx_u=size(domain%u,1)
        realny=size(domain%z,3)
        realny_v=size(domain%v,3)
        ! this isn't perfect because if the number of vertical levels happen to equal forcing levels...
        zaxis_is_third = (size(domain%z,3)==nz)
        if (zaxis_is_third) then
            if (size(domain%z,2)==nz) then
                ! zaxis is not determined, test more rigorously
                zaxis_is_third = (minval(domain%z(:,:,2) - domain%z(:,:,1)) >  0) .and. &
                                 (minval(domain%z(:,2,:) - domain%z(:,1,:)) <= 0)
                if (debug) then
                    !$omp critical (print_lock)
                    write(*,*) "Warning: z axis in a domain is not determined:"
                    write(*,*) "Z shape = ", shape(domain%z)
                    write(*,*) "U shape = ", shape(domain%u)
                    if (zaxis_is_third) then
                        write(*,*) "Assuming zaxis is the third axis"
                    else
                        write(*,*) "Assuming zaxis is the second axis"
                    endif
                    !$omp end critical (print_lock)
                endif
            endif
        endif

        ! determine whether or not we are working on a staggered grid.
        staggered = (realnx/=realnx_u).and.(realny/=realny_v)

        do i=1,nz
            preU_layers(i)=sum(domain%u(:realnx_u,i,:))/(realnx_u * realny)
            preV_layers(i)=sum(domain%v(:,i,:realny_v))/(realnx * realny_v)
        enddo
        do i=1,nz
            bottom=i-vsmooth
            top=i+vsmooth
            if (bottom<1) then
                bottom=1
            endif
            if (top>nz) then
                top=nz
            endif
            U_layers(i)=sum(preU_layers(bottom:top))/(top-bottom+1)
            V_layers(i)=sum(preV_layers(bottom:top))/(top-bottom+1)
        enddo

        if (.not.data_allocated) then
            if (debug) then
                !$omp critical (print_lock)
                write(*,"(A,A)") char(13),"Allocating Linear Wind Data and FFTW plans"
                !$omp end critical (print_lock)
            endif
            ! using fftw_alloc routines to ensure better allignment for vectorization
            n_elements = nx*ny*nz
            uh_aligned_data = fftw_alloc_complex(n_elements)
            call c_f_pointer(uh_aligned_data, uhat, [nx,ny,nz])
            u_h_aligned_data = fftw_alloc_complex(n_elements)
            call c_f_pointer(u_h_aligned_data, u_hat, [nx,ny,nz])
            vh_aligned_data = fftw_alloc_complex(n_elements)
            call c_f_pointer(vh_aligned_data, vhat, [nx,ny,nz])
            v_h_aligned_data = fftw_alloc_complex(n_elements)
            call c_f_pointer(v_h_aligned_data, v_hat, [nx,ny,nz])
            data_allocated = .True.

            allocate(uplans(nz))
            allocate(vplans(nz))
            do i=1,nz
                uplans(i) = fftw_plan_dft_2d(ny,nx, uhat(:,:,i),u_hat(:,:,i), FFTW_BACKWARD, FFTW_MEASURE)!PATIENT)!FFTW_ESTIMATE)
                vplans(i) = fftw_plan_dft_2d(ny,nx, vhat(:,:,i),v_hat(:,:,i), FFTW_BACKWARD, FFTW_MEASURE)!PATIENT)!FFTW_ESTIMATE)
            enddo
            if (debug) then
                !$omp critical (print_lock)
                write(*,*) "Allocation Complete:",trim(str(nx))," ",trim(str(ny))," ",trim(str(nz))
                !$omp end critical (print_lock)
            endif
        endif
        if (reverse) then
            !$omp critical (print_lock)
            write(*,*) "ERROR: reversing linear winds not set up for parallel fftw computation yet"
            write(*,*) "also not set up for the new spatial wind LUTs"
            !$omp end critical (print_lock)
            stop
        endif

        !$omp parallel default(shared), private(k,l,kl,sig,denom,m,msq,mimag,ineta,z,gain,offset,i,U,V) &
        !$omp firstprivate(nx,ny,nz, zaxis_is_third, Nsq, realnx,realny, realnx_u, realny_v) &
        !$omp firstprivate(useDensity, debug, staggered, reverse, buffer)

        ! these should be stored in a separate data structure... and allocated/deallocated in a subroutine
        ! for now these are reallocated/deallocated everytime so we can use it for different sized domains (e.g. coarse and fine)
        ! maybe linear winds need to be embedded in an object instead of a module to avoid this problem...
        if (.not.allocated(k)) then
            allocate(k(nx,ny))
            allocate(l(nx,ny))
            allocate(kl(nx,ny))
            allocate(sig(nx,ny))
            allocate(denom(nx,ny))
            allocate(m(nx,ny))
            allocate(msq(nx,ny))
            allocate(mimag(nx,ny))
            allocate(ineta(nx,ny))

            ! Compute 2D k and l wavenumber fields (could be stored in options or domain or something)
            offset=pi/domain%dx
            gain=2*offset/(nx-1)
            ! don't use an implied do loop because PGF90 crashes with it and openmp...
            ! k(:,1) = (/((i*gain-offset),i=0,nx-1)/)
            do i=1,nx
                k(i,1) = ((i-1)*gain-offset)
            end do
            do i=2,ny
                k(:,i)=k(:,1)
            enddo

            gain=2*offset/(ny-1)
            ! don't use an implied do loop because PGF90 crashes with it and openmp...
            ! l(1,:) = (/((i*gain-offset),i=0,ny-1)/)
            do i=1,ny
                l(1,i) = ((i-1)*gain-offset)
            end do
            do i=2,nx
                l(i,:)=l(1,:)
            enddo

            ! finally compute the kl combination array
            kl = k**2+l**2
            WHERE (kl==0.0) kl=1e-15
        endif

        m=1

        !$omp do
        do z=1,nz
            U = U_layers(z)
            V = V_layers(z)
            if ((abs(U)+abs(V))>0.5) then
                sig  = U*k+V*l
                where(sig==0.0) sig=1e-15
                denom = sig**2 ! -f**2

                !   where(denom.eq.0) denom=1e-20
                ! # two possible non-hydrostatic versions
                ! not yet converted from python...
                ! # mimag=np.zeros((Ny,Nx)).astype('complex')
                ! # msq = (Nsq/denom * kl).astype('complex')                   # % vertical wave number, hydrostatic
                ! # msq = ((Nsq-sig**2)/denom * kl).astype('complex')          # % vertical wave number, hydrostatic
                ! # mimag.imag=(np.sqrt(-msq)).real
                ! # m=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)
                !
                ! mimag=np.zeros(Nx).astype('complex')
                ! mimag.imag=(np.sqrt(-msq)).real
                ! m=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)

                msq   = Nsq / denom * kl
                mimag = 0 + 0 * j ! be sure to reset real and imaginary components
                mimag = mimag + (real(sqrt(-msq)) * j)

                m = sqrt(msq)         ! # % vertical wave number, hydrostatic
                where(sig < 0) m = m * (-1)   ! equivalent to m=m*sign(sig)
                where(real(msq) < 0) m = mimag

                if (zaxis_is_third) then
                    ineta = j * domain%fzs * exp(j * m * &
                            sum(domain%z(:,:,z) - domain%z(:,:,1)) / (realnx*realny))
                else
                    ineta=j*domain%fzs*exp(j*m* &
                            sum(domain%z(:,z,:)-domain%z(:,1,:)) / (realnx*realny))
                endif
                !  what=sig*ineta

                ! with coriolis : [+/-]j*[l/k]*f
                !uhat = (0 - m) * ((sig * k) - (j * l * f)) * ineta / kl
                !vhat = (0 - m) * ((sig * l) + (j * k * f)) * ineta / kl
                ! removed coriolis term assumes the scale coriolis operates at is largely defined by the coarse model
                ! if using e.g. a sounding or otherwise spatially constant u/v, then coriolis should be defined.
                ineta = ineta / (kl / ((0 - m) * sig))
                uhat(:,:,z) = k * ineta
                vhat(:,:,z) = l * ineta

                ! pull it back out of fourier space.
                ! NOTE, the fftw transform inherently scales by N so the Fzs/Nx/Ny provides the only normalization necessary (I think)
                call ifftshift(uhat, fixed_axis=z)
                call ifftshift(vhat, fixed_axis=z)
                ! plans are only created once and re-used, this is a limit to further parallelization at the moment.
                call fftw_execute_dft(uplans(z), uhat(:,:,z),u_hat(:,:,z))
                call fftw_execute_dft(vplans(z), vhat(:,:,z),v_hat(:,:,z))

                ! The linear_mask field only applies to the high res grid (for now at least)
                ! NOTE: we should be able to do this without the loop, but ifort -O was giving the wrong answer...
                ! possible compiler bug version 12.1.x?
                if (.not.reverse) then
                    do i=1,realny
                        u_hat(buffer+1:buffer+realnx,buffer+i,z) = &
                            u_hat(buffer+1:buffer+realnx,buffer+i,z) * linear_mask(:,i)
                        v_hat(buffer+1:buffer+realnx,buffer+i,z) = &
                            v_hat(buffer+1:buffer+realnx,buffer+i,z) * linear_mask(:,i)
                    enddo
                endif

                ! u/vhat are first staggered to apply to u/v appropriately if on a staggered grid (realnx/=real_nx_u)
                ! when removing linear winds from forcing data, it may NOT be on a staggered grid
                ! NOTE: we should be able to do this without the loop, but ifort -O was giving the wrong answer...
                ! possible compiler bug version 12.1.x?
                if (staggered) then
                    do i=1,ny-1
                        u_hat(1:nx-1,i,z) = (u_hat(1:nx-1,i,z) + u_hat(2:nx,i,z)) /2
                        v_hat(:,i,z)      = (v_hat(:,i,z)      + v_hat(:,i+1,z))  /2
                    enddo
                    i=ny
                    u_hat(1:nx-1,i,z) = (u_hat(1:nx-1,i,z) + u_hat(2:nx,i,z)) /2
                endif

                ! If we are using density in the advection calculations, modify the linear perturbation
                ! to get the vertical velocities closer to what they would be without density (boussinesq)
                ! need to check if this makes the most sense when close to the surface
                ! if (useDensity) then
                !     if (debug) then
                !         !$omp critical (print_lock)
                !         write(*,*) "Using a density correction in linear winds"
                !         !$omp end critical (print_lock)
                !     endif
                !     u_hat(buffer:realnx_u+buffer,buffer:realny+buffer,z) = &
                !         2*real(u_hat(buffer:realnx_u+buffer,buffer:realny+buffer,z))! / domain%rho(1:realnx,z,1:realny)
                !     v_hat(buffer:realnx+buffer,buffer:realny_v+buffer,z) = &
                !         2*real(v_hat(buffer:realnx+buffer,buffer:realny_v+buffer,z))! / domain%rho(1:realnx,z,1:realny)
                ! endif

                ! If we are removing linear winds from a low res field, subtract u_hat v_hat instead
                ! real(real()) extracts real component of complex, then converts to a real data type (may not be necessary except for IO?)
                if (reverse) then
                    domain%u(:,z,:) = domain%u(:,z,:) - &
                        real(real( u_hat(1+buffer:realnx_u+buffer, 1+buffer:realny+buffer  ,z) ))*linear_contribution

                    domain%v(:,z,:) = domain%v(:,z,:) - &
                        real(real( v_hat(1+buffer:realnx+buffer,   1+buffer:realny_v+buffer,z) ))*linear_contribution
                else
                    ! note, linear_contribution component comes from the linear mask applied above on the mass grid
                    if (staggered) then
                        domain%u(2:realnx,z,:) = domain%u(2:realnx,z,:) + &
                            real(real(u_hat(1+buffer:realnx-1+buffer,1+buffer:realny+buffer  ,z) ))
                        domain%v(:,z,2:realny) = domain%v(:,z,2:realny) + &
                            real(real(v_hat(1+buffer:realnx+buffer,  1+buffer:realny-1+buffer,z) ))
                    else
                        domain%u(:,z,:) = domain%u(:,z,:) + &
                            real(real(u_hat(1+buffer:realnx_u+buffer,1+buffer:realny+buffer  ,z) ))
                        domain%v(:,z,:) = domain%v(:,z,:) + &
                            real(real(v_hat(1+buffer:realnx+buffer,  1+buffer:realny_v+buffer,z) ))
                    endif
                endif

                if (debug) then
                    if (z==1)then
                        !$omp critical (print_lock)
                        write(*,*) "Nsq = ", Nsq
                        write(*,*) "U=",U, "    V=",V
                        write(*,*) "realnx=",realnx, "; nx=",nx, "; buffer=",buffer
                        write(*,*) "realny=",realny, "; ny=",ny!, buffer
                        !$omp end critical (print_lock)
                    endif
                endif
            endif
        end do ! z-loop
        !$omp end do

        ! finally deallocate all temporary arrays that were created... chould be a datastructure and a subroutine...
        deallocate(k,l,kl,sig,denom,m,ineta,msq,mimag)

        !$omp end parallel

        ! these are subroutine scoped not module, they should be deallocated automatically anyway
        deallocate(U_layers,V_layers,preU_layers,preV_layers)
    end subroutine linear_winds


    !>----------------------------------------------------------
    !! Add a smoothed buffer around the edge of the terrain to prevent crazy wrap around effects
    !! in the FFT due to discontinuities between the left and right (or top and bottom) edges of the domain
    !!
    !!----------------------------------------------------------
    subroutine add_buffer_topo(terrain, buffer_topo, smooth_window, debug)
        implicit none
        real, dimension(:,:), intent(in) :: terrain
        complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:), intent(inout) :: buffer_topo
        integer, intent(in) :: smooth_window
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
        lt_data%vh_aligned_data = fftw_alloc_complex(n_elements)
        lt_data%vp_aligned_data = fftw_alloc_complex(n_elements)
        !$omp end critical (fftw_lock)

        call c_f_pointer(lt_data%uh_aligned_data,   lt_data%uhat,       [nx,ny])
        call c_f_pointer(lt_data%up_aligned_data,   lt_data%u_perturb,  [nx,ny])
        call c_f_pointer(lt_data%vh_aligned_data,   lt_data%vhat,       [nx,ny])
        call c_f_pointer(lt_data%vp_aligned_data,   lt_data%v_perturb,  [nx,ny])

        ! note FFTW plan creation is not threadsafe
        !$omp critical (fftw_lock)
        lt_data%uplan = fftw_plan_dft_2d(ny,nx, lt_data%uhat, lt_data%u_perturb, FFTW_BACKWARD, FFTW_MEASURE) ! alternatives to MEASURE are PATIENT, or ESTIMATE
        lt_data%vplan = fftw_plan_dft_2d(ny,nx, lt_data%vhat, lt_data%v_perturb, FFTW_BACKWARD, FFTW_MEASURE) ! alternatives to MEASURE are PATIENT, or ESTIMATE
        !$omp end critical (fftw_lock)


    end subroutine initialize_linear_theory_data


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
        call fftw_free(lt_data%vh_aligned_data)
        call fftw_free(lt_data%vp_aligned_data)
        !$omp end critical (fftw_lock)

        NULLIFY(lt_data%uhat)
        NULLIFY(lt_data%u_perturb)
        NULLIFY(lt_data%vhat)
        NULLIFY(lt_data%v_perturb)

        ! note FFTW plan creation is not threadsafe
        !$omp critical (fftw_lock)
        call fftw_destroy_plan(lt_data%uplan)
        call fftw_destroy_plan(lt_data%vplan)
        !$omp end critical (fftw_lock)
    end subroutine destroy_linear_theory_data

    !>----------------------------------------------------------
    !! Compute look up tables for all combinations of U, V, and Nsq
    !!
    !!----------------------------------------------------------
    subroutine initialize_spatial_winds(domain,options,reverse)
        implicit none
        class(linearizable_type),intent(inout)::domain
        type(options_type), intent(in) :: options
        logical, intent(in) :: reverse

        ! local variables used to calculate the LUT
        real :: u,v, layer_height
        integer :: nx,ny,nz, i,j,k,z, nxu,nyv, error
        integer :: fftnx, fftny
        integer, dimension(3,2) :: LUT_dims
        integer :: loops_completed ! this is just used to measure progress in the LUT creation

        ! the domain to work over
        nx = size(domain%lat,1)
        nz = size(domain%u,2)
        ny = size(domain%lat,2)

        nxu = size(domain%u,1)
        nyv = size(domain%v,3)

        fftnx = size(domain%fzs,1)
        fftny = size(domain%fzs,2)
        ! note:
        ! buffer = (fftnx - nx)/2

        ! default assumes no errors in reading the LUT
        error = 0

        ! store to make it easy to check dim sizes in read_LUT
        LUT_dims(:,1) = [nxu,nz,ny]
        LUT_dims(:,2) = [nx,nz,nyv]

        ! create the array of dir and nsq values to create LUTs for
        if (.not.allocated(dir_values)) then
            allocate( dir_values(n_dir_values) )
            allocate( nsq_values(n_nsq_values) )
            allocate( spd_values(n_spd_values) )
        endif
        ! generate the values for each look up table dimension
        ! generate table of wind directions to be used
        do i=1,n_dir_values
            dir_values(i) = (i - 1) / real(n_dir_values - 1) * (dirmax - dirmin) + dirmin
        enddo
        ! generate table of wind speeds to be used
        do i=1,n_spd_values
            spd_values(i) = (i - 1) / real(n_spd_values - 1) * (spdmax - spdmin) + spdmin
        enddo
        ! generate table of Brunt-Vaisalla frequencies (Nsq) to be used
        do i=1,n_nsq_values
            nsq_values(i) = (i - 1) / real(n_nsq_values - 1) * (nsqmax - nsqmin) + nsqmin
        enddo

        ! Allocate the (LARGE) look up tables for both U and V
        if (reverse) then
            ! This is the reverse LUT to remove the linear contribution from the low res field
            allocate(rev_u_LUT(n_spd_values,n_dir_values,n_nsq_values,nxu,nz,ny))
            allocate(rev_v_LUT(n_spd_values,n_dir_values,n_nsq_values,nx,nz,nyv))
            u_LUT=>rev_u_LUT
            v_LUT=>rev_v_LUT
        else
            ! this is the more common forward transform
            if (.not.options%lt_options%read_LUT) then
                allocate(hi_u_LUT(n_spd_values,n_dir_values,n_nsq_values,nxu,nz,ny))
                allocate(hi_v_LUT(n_spd_values,n_dir_values,n_nsq_values,nx,nz,nyv))
                error=0
            else
                print*, "Reading LUT from file: ", trim(options%lt_options%u_LUT_Filename)
                error = read_LUT(options%lt_options%u_LUT_Filename, hi_u_LUT, hi_v_LUT, options%dz_levels, LUT_dims, options%lt_options)
                if (error/=0) then
                    write(*,*) "WARNING: LUT on disk does not match that specified in the namelist or does not exist."
                    write(*,*) "    LUT will be recreated"
                    if (allocated(hi_u_LUT)) deallocate(hi_u_LUT)
                    allocate(hi_u_LUT(n_spd_values,n_dir_values,n_nsq_values,nxu,nz,ny))
                    if (allocated(hi_v_LUT)) deallocate(hi_v_LUT)
                    allocate(hi_v_LUT(n_spd_values,n_dir_values,n_nsq_values,nx,nz,nyv))
                endif
            endif
            u_LUT=>hi_u_LUT
            v_LUT=>hi_v_LUT

        endif

        if (options%debug) then
            write(*,*) "Wind Speeds:",spd_values
            write(*,*) "Directions:",360*dir_values/(2*pi)
            write(*,*) "Stabilities:",exp(nsq_values)
        endif

        if (reverse.or.(.not.((options%lt_options%read_LUT).and.(error==0)))) then
            ! loop over combinations of U, V, and Nsq values
            loops_completed = 0
            write(*,*) "Percent Completed:"
            !$omp parallel default(shared) &
            !$omp private(i,j,k,z, u,v, layer_height)
            ! $omp private(lt_data_m)

            ! initialization has to happen in each thread so each thread has its own copy
            ! lt_data_m is a threadprivate variable, within initialization, there are omp critical sections for fftw calls
            call initialize_linear_theory_data(lt_data_m, fftnx, fftny, domain%dx)
            !$omp do
            do i=1, n_dir_values
                ! set the domain wide U and V values to the current u and v values
                ! this could use u/v_perturbation, but those would need to be put in a linearizable structure...
                do k=1, n_spd_values
                    if (options%interactive) then
                        !$omp critical (print_lock)
                        write(*,"(A,f5.1,A$)") char(13), loops_completed/real(n_dir_values*n_spd_values)*100," %"
                        !$omp end critical (print_lock)
                    endif
                    do j=1, n_nsq_values
                        u = calc_u( dir_values(i), spd_values(k) )
                        v = calc_v( dir_values(i), spd_values(k) )

                        ! calculate the linear wind field for the current u and v values
                        do z=1,nz
                            if (reverse) then
                                layer_height = domain%z(1,1,z)-domain%terrain(1,1)
                            else
                                layer_height = domain%z(1,z,1)-domain%terrain(1,1)
                            endif
                            call linear_perturbation(u, v, exp(nsq_values(j)), layer_height, domain%fzs, lt_data_m)

                            ! need to handle stagger (nxu /= nx) and the buffer around edges of the domain
                            if (nxu /= nx) then
                                u_LUT(k,i,j,2:nx,z, :  ) = real( real(                                  &
                                        ( lt_data_m%u_perturb(1+buffer:nx+buffer-1,   1+buffer:ny+buffer) &
                                        + lt_data_m%u_perturb(2+buffer:nx+buffer,     1+buffer:ny+buffer)) )) / 2

                                v_LUT(k,i,j, :,  z,2:ny) = real( real(                                      &
                                        ( lt_data_m%v_perturb(1+buffer:nx+buffer,     1+buffer:ny+buffer-1)   &
                                        + lt_data_m%v_perturb(1+buffer:nx+buffer,     2+buffer:ny+buffer)) )) / 2
                            else
                                u_LUT(k,i,j,:,z,:) = real( real(                                  &
                                        lt_data_m%u_perturb(1+buffer:nx+buffer,     1+buffer:ny+buffer) ))
                                v_LUT(k,i,j,:,z,:) = real( real(                                  &
                                        lt_data_m%v_perturb(1+buffer:nx+buffer,     1+buffer:ny+buffer) ))
                            endif
                        enddo

                    end do
                    !$omp critical (print_lock)
                    loops_completed = loops_completed+1
                    !$omp end critical (print_lock)
                end do

            end do
            !$omp end do
            ! memory needs to be freed so this structure can be used again when removing linear winds
            call destroy_linear_theory_data(lt_data_m)
            !$omp end parallel
            write(*,"(A,f5.1,A$)") char(13), loops_completed/real(n_dir_values*n_spd_values)*100," %"
            write(*,*) char(10),"--------  Linear wind look up table generation complete ---------"
        endif

        if ((options%lt_options%write_LUT).and.(.not.reverse)) then
            if ((options%lt_options%read_LUT) .and. (error == 0)) then
                print*, "Not writing Linear Theory LUT to file because LUT was read from file"
            else
                print*, "Writing Linear Theory LUT to file: ", trim(options%lt_options%u_LUT_Filename)
                error = write_LUT(options%lt_options%u_LUT_Filename, hi_u_LUT, hi_v_LUT, options%dz_levels, options%lt_options)
            endif
        endif

    end subroutine initialize_spatial_winds


    !>----------------------------------------------------------
    !! Calculate the weights between the positions bestpos and nextpos
    !! based on the distance between match and indata(nextpos) (normalized by nextpos - bestpos)
    !! assumes indata is monotonically increasing,
    !! bestpos must be set prior to entry
    !! nextpos is calculated internally (either 1, bestpos+1, or n)
    !!
    !!----------------------------------------------------------
    function calc_weight(indata, bestpos, nextpos, match) result(weight)
        implicit none
        real :: weight
        real, dimension(:), intent(in) :: indata
        integer, intent(in) :: bestpos
        integer, intent(inout) :: nextpos
        real, intent(in) :: match

        integer :: n

        n=size(indata)

        if (match<indata(1)) then
            nextpos=1
            weight=1
        else
            if (bestpos==n) then
                nextpos=n
                weight=1
            else
                nextpos=bestpos+1
                weight=(indata(nextpos)-match) / (indata(nextpos) - indata(bestpos))
            endif
        endif

    end function

    !>----------------------------------------------------------
    !! Compute a spatially variable linear wind perturbation
    !! based off of look uptables computed in via setup
    !! for each grid point, find the closest LUT data in U and V space
    !! then bilinearly interpolate the nearest LUT values for that points linear wind field
    !!
    !!----------------------------------------------------------
    subroutine spatial_winds(domain,reverse, vsmooth, winsz)
        implicit none
        class(linearizable_type),intent(inout)::domain
        logical, intent(in) :: reverse
        integer, intent(in) :: vsmooth
        integer, intent(in) :: winsz

        integer :: nx,nxu, ny,nyv, nz, i,j,k, smoothz
        integer :: uk, vi !store a separate value of i for v and of k for u to we can handle nx+1, ny+1
        integer :: step, dpos, npos, spos, nexts, nextd, nextn
        integer :: north, south, east, west, top, bottom, n
        real :: u, v
        real :: dweight, nweight, sweight, curspd, curdir, curnsq, wind_first, wind_second

        nx=size(domain%lat,1)
        ny=size(domain%lat,2)
        nz=size(domain%u,2)
        nxu=size(domain%u,1)
        nyv=size(domain%v,3)

        if (reverse) then
            u_LUT=>rev_u_LUT
            v_LUT=>rev_v_LUT
            u_perturbation=>lo_u_perturbation
            v_perturbation=>lo_v_perturbation
        else
            u_LUT=>hi_u_LUT
            v_LUT=>hi_v_LUT
            u_perturbation=>hi_u_perturbation
            v_perturbation=>hi_v_perturbation
        endif

        if (reverse) print*, "WARNING using fixed nsq for linear wind removal: 3e-6"
        !$omp parallel firstprivate(nx,nxu,ny,nyv,nz, reverse, vsmooth, winsz), default(none), &
        !$omp private(i,j,k,step, uk, vi, east, west, north, south, top, bottom), &
        !$omp private(spos, dpos, npos, nexts,nextd, nextn,n, smoothz, u, v), &
        !$omp private(wind_first, wind_second, curspd, curdir, curnsq, sweight,dweight, nweight), &
        !$omp shared(domain, spd_values, dir_values, nsq_values, u_LUT, v_LUT, linear_mask), &
        !$omp shared(u_perturbation, v_perturbation, linear_update_fraction, linear_contribution, nsq_calibration), &
        !$omp shared(min_stability, max_stability, n_dir_values, n_spd_values, n_nsq_values, smooth_nsq)
        !$omp do
        do k=1,ny
            do j=1,nz
                do i=1,nx

                    ! look up vsmooth gridcells up to nz at the maximum
                    top = min(j+vsmooth, nz)
                    ! if (top-j)/=vsmooth, then look down enough layers to make the window vsmooth in size
                    bottom = max(1, j - (vsmooth - (top-j)) )

                    if (.not.reverse) then
                        domain%nsquared(i,j,k) = calc_stability(domain%th(i,bottom,k), domain%th(i,top,k),  &
                                                                domain%pii(i,bottom,k),domain%pii(i,top,k), &
                                                                domain%z(i,bottom,k),  domain%z(i,top,k),   &
                                                                domain%qv(i,bottom,k), domain%qv(i,top,k),  &
                                                                domain%cloud(i,j,k)+domain%ice(i,j,k)       &
                                                                +domain%qrain(i,j,k)+domain%qsnow(i,j,k))
                        domain%nsquared(i,j,k) = max(min(domain%nsquared(i,j,k) * nsq_calibration(i,k), max_stability), min_stability)
                    else
                        ! Low-res boundary condition variables will be in a different array format.  It should be
                        ! easy enough to call calc_stability after e.g. transposing z and y dimension, but some
                        ! e.g. pii will not be set in the forcing data, so this may need a little thought.
                        domain%nsquared(i,j,k) = 3e-6
                    endif
                end do
            end do

            if (smooth_nsq) then
                do j=1,nz
                    ! compute window as above.
                    top = min(j+vsmooth,nz)
                    bottom = max(1, j - (vsmooth - (top-j)) )

                    do smoothz = bottom, j-1
                        domain%nsquared(:,j,k) = domain%nsquared(:,j,k) + domain%nsquared(:,smoothz,k)
                    end do
                    do smoothz = j+1, top
                        domain%nsquared(:,j,k) = domain%nsquared(:,j,k) + domain%nsquared(:,smoothz,k)
                    end do
                    domain%nsquared(:,j,k) = domain%nsquared(:,j,k)/(top-bottom+1)
                end do
            endif
        end do
        !$omp end do
        !$omp barrier
        !$omp do
        do k=1, nyv
            do j=1, nz
                do i=1, nxu
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
                        u = sum( domain%u(west:east,j,south:north) ) / n
                        v = sum( domain%v(west:east,j,south:north) ) / n
                    else
                        u = domain%u(i,j,uk)
                        v = domain%v(vi,j,k)
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
                    curnsq = sum(log( domain%nsquared(west:east, bottom:top, south:north) ))
                    curnsq = curnsq / n
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
                        wind_first =      nweight  * (dweight * u_LUT(spos, dpos,npos, i,j,k) + (1-dweight) * u_LUT(spos, nextd,npos, i,j,k))   &
                                    +  (1-nweight) * (dweight * u_LUT(spos, dpos,nextn,i,j,k) + (1-dweight) * u_LUT(spos, nextd,nextn,i,j,k))

                        wind_second=      nweight  * (dweight * u_LUT(nexts,dpos,npos, i,j,k) + (1-dweight) * u_LUT(nexts,nextd,npos, i,j,k))   &
                                    +  (1-nweight) * (dweight * u_LUT(nexts,dpos,nextn,i,j,k) + (1-dweight) * u_LUT(nexts,nextd,nextn,i,j,k))

                        u_perturbation(i,j,k) = u_perturbation(i,j,k) * (1-linear_update_fraction) &
                                    + linear_update_fraction * (sweight*wind_first + (1-sweight)*wind_second)

                        if (reverse) then
                            domain%u(i,j,k) = domain%u(i,j,k) - u_perturbation(i,j,k) * linear_contribution
                        else
                            domain%u(i,j,k) = domain%u(i,j,k) + u_perturbation(i,j,k) * linear_mask(min(nx,i),min(ny,k))
                        endif
                    endif
                    if (i<=nx) then
                        wind_first =      nweight  * (dweight * v_LUT(spos, dpos,npos, i,j,k) + (1-dweight) * v_LUT(spos, nextd,npos, i,j,k))    &
                                    +  (1-nweight) * (dweight * v_LUT(spos, dpos,nextn,i,j,k) + (1-dweight) * v_LUT(spos, nextd,nextn,i,j,k))

                        wind_second=      nweight  * (dweight * v_LUT(nexts,dpos,npos, i,j,k) + (1-dweight) * v_LUT(nexts,nextd,npos, i,j,k))    &
                                    +  (1-nweight) * (dweight * v_LUT(nexts,dpos,nextn,i,j,k) + (1-dweight) * v_LUT(nexts,nextd,nextn,i,j,k))

                        v_perturbation(i,j,k) = v_perturbation(i,j,k) * (1-linear_update_fraction) &
                                    + linear_update_fraction * (sweight*wind_first + (1-sweight)*wind_second)

                        if (reverse) then
                            domain%v(i,j,k) = domain%v(i,j,k) - v_perturbation(i,j,k) * linear_contribution
                        else
                            ! for the high res domain, linear_mask should incorporate linear_contribution
                            domain%v(i,j,k) = domain%v(i,j,k) + v_perturbation(i,j,k) * linear_mask(min(nx,i),min(ny,k))
                        endif
                    endif
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine spatial_winds

    !>----------------------------------------------------------
    !! Copy options from the options data structure into module level variables
    !!
    !!----------------------------------------------------------
    subroutine set_module_options(options)
        implicit none
        type(options_type), intent(in) :: options

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

    end subroutine set_module_options

    !>----------------------------------------------------------
    !! Called from linear_perturb the first time perturb is called
    !! compute FFT(terrain), and dzdx,dzdy components
    !!
    !!----------------------------------------------------------
    subroutine setup_linwinds(domain, options, reverse, useDensity)
        implicit none
        class(linearizable_type),intent(inout)  :: domain
        type(options_type),      intent(in)     :: options
        logical,                 intent(in)     :: reverse, useDensity
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
        call add_buffer_topo(domain%terrain, complex_terrain_firstpass, 5)
        buffer = 2
        ! then further add a small (~2) grid cell buffer where all cells have the same value
        call add_buffer_topo(real(real(complex_terrain_firstpass)), complex_terrain, 0, debug=options%debug)
        buffer = buffer + original_buffer

        nx = size(complex_terrain, 1)
        ny = size(complex_terrain, 2)

        write(*,*) "Initializing linear winds : ",nx,ny
        allocate(domain%fzs(nx,ny))

        ! calculate the fourier transform of the terrain for use in linear winds
        plan = fftw_plan_dft_2d(ny, nx, complex_terrain, domain%fzs, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, complex_terrain, domain%fzs)
        call fftw_destroy_plan(plan)
        ! normalize FFT by N - grid cells
        domain%fzs = domain%fzs / (nx * ny)
        ! shift the grid cell quadrants
        ! need to test what effect all of the related shifts actually have...
        call fftshift(domain%fzs)

        if (linear_contribution/=1) then
            write(*,*) "Using a fraction of the linear perturbation:",linear_contribution
        endif

        nx = size(domain%terrain, 1)
        nz = size(domain%u,       2)
        ny = size(domain%terrain, 2)

        if (.not.allocated(domain%nsquared)) allocate(domain%nsquared(nx,nz,ny))

        ! set up linear_mask variable
        if (.not.reverse) then
            if (allocated(linear_mask)) deallocate(linear_mask)
            allocate(linear_mask(nx,ny))
            linear_mask = linear_contribution

            if (use_linear_mask) then
                write(*,*) "Reading Linear Mask"
                write(*,*) "  from file: " // trim(options%linear_mask_file)
                write(*,*) "  with var: "  // trim(options%linear_mask_var)
                call io_read2d(options%linear_mask_file, options%linear_mask_var, domain%linear_mask)

                linear_mask = domain%linear_mask * linear_contribution
            endif

            ! Stupidly simple adjustment to the linear wind field to account for using density.
            ! If we are using density in the advection calculations, modify the linear perturbation
            ! to get the vertical velocities closer to what they would be without density (boussinesq)
            ! need to check if this makes the most sense when close to the surface
            if (useDensity) then
                linear_mask = linear_mask * 2
            endif
        endif

        ! set up nsq_calibration variable
        if ( (.not.reverse) .and. (.not.allocated(nsq_calibration)) ) then
            allocate(nsq_calibration(nx,ny))
            nsq_calibration = 1

            if (use_nsq_calibration) then
                write(*,*) "Reading Linear Mask"
                write(*,*) "  from file: " // trim(options%nsq_calibration_file)
                write(*,*) "  with var: "  // trim(options%nsq_calibration_var)
                call io_read2d(options%nsq_calibration_file, options%nsq_calibration_var, domain%nsq_calibration)
                nsq_calibration = domain%nsq_calibration

                where(nsq_calibration<1) nsq_calibration = 1 + 1/( (1-1/nsq_calibration)/100 )
                where(nsq_calibration>1) nsq_calibration = 1 + (nsq_calibration-1)/100
            endif
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
            if    ((.not.allocated(hi_u_LUT)  .and. (.not.reverse)) &
             .or.  (.not.allocated(rev_u_LUT) .and. reverse)) then

                write(*,*) "Generating a spatially variable linear perturbation look up table"
                call initialize_spatial_winds(domain, options, reverse)

            endif
        endif


    end subroutine setup_linwinds

    !>----------------------------------------------------------
    !! Initialize and/or apply linear wind solution.
    !!
    !! Called from ICAR to update the U and V wind fields based on linear theory (W is calculated to balance U/V)
    !!
    !!----------------------------------------------------------
    subroutine linear_perturb(domain,options,vsmooth,reverse,useDensity)
        implicit none
        class(linearizable_type),intent(inout)::domain
        type(options_type), intent(in) :: options
        integer, intent(in) :: vsmooth
        logical, intent(in), optional :: reverse,useDensity
        logical :: rev, useD
        logical, save :: debug=.True.
        real::stability

        ! these probably need to get moved to options...
        if (present(reverse)) then
            rev=reverse
        else
            rev=.False.
        endif

        if (present(useDensity)) then
            useD=useDensity
        else
            useD=.False.
        endif

        ! this is a little trickier, because it does have to be domain dependant... could at least be stored in the domain though...
        if (rev) then
            linear_contribution = options%lt_options%rm_linear_contribution
            N_squared = options%lt_options%rm_N_squared
        else
            linear_contribution = options%lt_options%linear_contribution
            N_squared = options%lt_options%N_squared
        endif

        ! if linear_perturb hasn't been called before we need to perform some setup actions.
        if (.not.allocated(domain%fzs)) then
            call setup_linwinds(domain, options, rev, useD)
        endif

        ! add the spatially variable linear field
        ! if we are reverseing the effects, that means we are in the low-res domain
        ! that domain does not have a spatial LUT calculated, so it can not be performed
        if (use_spatial_linear_fields)then
            call spatial_winds(domain,rev, vsmooth, stability_window_size)
        else
            ! Nsq = squared Brunt Vaisalla frequency (1/s) typically from dry static stability
            stability = calc_domain_stability(domain)
            ! This should probably be called twice, once for dry, and once or moist regions
            call linear_winds(domain,stability,vsmooth,reverse,useDensity,debug)
        endif
        debug=.False.

    end subroutine linear_perturb
end module linear_theory_winds
