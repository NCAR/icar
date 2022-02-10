!> ----------------------------------------------------------------------------
!!  Main time stepping module.
!!  Calculates a stable time step (dt) and loops over physics calls
!!  Also updates boundaries every time step.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module time_step
    use data_structures             ! *_type  types and kCONSTANTS
    use icar_constants,             only : Rd
    use microphysics,               only : mp
    use advection,                  only : advect
    use mod_atm_utilities,          only : exner_function
    use convection,                 only : convect
    use land_surface,               only : lsm
    use planetary_boundary_layer,   only : pbl
    use radiation,                  only : rad

    use domain_interface,           only : domain_t
    use options_interface,          only : options_t
    use debug_module,               only : domain_check
    use string,                     only : str

    implicit none
    private
    public :: step

    !  temporary variables used to update the w_real state. This should move to mod_atm_utilities
    real, allocatable :: lastw(:,:)
    real, allocatable :: currw(:,:)
    real, allocatable :: uw(:,:)
    real, allocatable :: vw(:,:)

contains


    !>------------------------------------------------------------
    !! Update model diagnostic fields
    !!
    !! Calculates most model diagnostic fields such as Psfc, 10m height winds and ustar
    !!
    !! @param domain    Model domain data structure to be updated
    !! @param options   Model options (not used at present)
    !!
    !!------------------------------------------------------------
    subroutine diagnostic_update(domain, options)
        implicit none
        type(domain_t),  intent(inout)   :: domain
        type(options_t), intent(in)      :: options
        integer :: z

        logical :: use_delta_terrain

        associate(ims => domain%ims, ime => domain%ime,                         &
                  jms => domain%jms, jme => domain%jme,                         &
                  kms => domain%kms, kme => domain%kme,                         &
                  exner                 => domain%exner%data_3d,                &
                  pressure              => domain%pressure%data_3d,             &
                  pressure_i            => domain%pressure_interface%data_3d,   &
                  dz_interface          => domain%dz_interface%data_3d,         &
                  psfc                  => domain%surface_pressure%data_2d,     &
                  density               => domain%density%data_3d,              &
                  temperature           => domain%temperature%data_3d,          &
                  temperature_i         => domain%temperature_interface%data_3d,&
                  u                     => domain%u%data_3d,                    &
                  v                     => domain%v%data_3d,                    &
                  w                     => domain%w%data_3d,                    &
                  w_real                => domain%w_real%data_3d,               &
                  u_mass                => domain%u_mass%data_3d,               &
                  v_mass                => domain%v_mass%data_3d,               &
                  potential_temperature => domain%potential_temperature%data_3d )


        exner = exner_function(pressure)

        ! domain%p_inter=domain%p
        ! call update_pressure(domain%p_inter, domain%z, domain%z_inter, domain%t)
        pressure_i(:,kms+1:kme, :) = (pressure(:,kms:kme-1, :) + pressure(:,kms+1:kme, :)) / 2
        pressure_i(:, kms, :) = pressure(:, kms, :) + (pressure(:, kms, :) - pressure(:, kms+1, :)) / 2
        ! this isn't correct, we should be using update_pressure or similar to solve this
        ! domain%ptop = 2*domain%p(:,nz,:) - domain%p(:,nz-1,:)
        if (associated(domain%surface_pressure%data_2d)) then
            psfc = pressure_i(:, kms, :)
        endif

        temperature = potential_temperature * exner
        temperature_i(:,kms+1:kme, :) = (temperature(:,kms:kme-1, :) + temperature(:,kms+1:kme, :)) / 2
        temperature_i(:, kms, :) = temperature(:, kms, :) + (temperature(:, kms, :) - temperature(:, kms+1, :)) / 2

        if (associated(domain%density%data_3d)) then
            density =  pressure / &
                        (Rd * temperature) ! kg/m^3
        endif
        if (associated(domain%u_mass%data_3d)) then
            u_mass = (u(ims+1:ime+1,:,:) + u(ims:ime,:,:)) / 2
        endif
        if (associated(domain%v_mass%data_3d)) then
            v_mass = (v(:,:,jms+1:jme+1) + v(:,:,jms:jme)) / 2
        endif


    ! NOTE: all code below is not implemented in ICAR 2.0 yet
    ! it is left as a reminder of what needs to be done, and example when the time comes
    !
    !     ! update mut
    !
    !     domain%p_inter=domain%p
    !     call update_pressure(domain%p_inter, domain%z, domain%z_inter, domain%t)
    !     domain%psfc = domain%p_inter(:,1,:)
    !     ! technically this isn't correct, we should be using update_pressure or similar to solve this
    !     domain%ptop = 2*domain%p(:,nz,:) - domain%p(:,nz-1,:)
    !
    !     ! dry mass in the gridcell is equivalent to the difference in pressure from top to bottom
    !     domain%mut(:,1:nz-1,:) = domain%p_inter(:,1:nz-1,:) - domain%p_inter(:,2:nz,:)
    !     domain%mut(:,nz,:) = domain%p_inter(:,nz,:) - domain%ptop
    !
        if (.not.allocated(lastw)) then
            allocate( lastw( ims+1:ime-1, jms+1:jme-1))
            allocate( currw( ims+1:ime-1, jms+1:jme-1))
            allocate(    uw( ims+1:ime,   jms+1:jme-1))
            allocate(    vw( ims+1:ime-1, jms+1:jme  ))
        endif

        ! temporary constant
        if (associated(domain%roughness_z0%data_2d)) then
            ! use log-law of the wall to convert from first model level to surface
            currw = karman / log((domain%z%data_3d(ims+1:ime-1,kms,jms+1:jme-1) - domain%terrain%data_2d(ims+1:ime-1,jms+1:jme-1)) / domain%roughness_z0%data_2d(ims+1:ime-1,jms+1:jme-1))
            ! use log-law of the wall to convert from surface to 10m height
            lastw = log(10.0 / domain%roughness_z0%data_2d(ims+1:ime-1,jms+1:jme-1)) / karman
        endif

        if (associated(domain%u_10m%data_2d)) then
            domain%ustar        (ims+1:ime-1,jms+1:jme-1) = u_mass      (ims+1:ime-1,kms,jms+1:jme-1) * currw
            domain%u_10m%data_2d(ims+1:ime-1,jms+1:jme-1) = domain%ustar(ims+1:ime-1,jms+1:jme-1)     * lastw
            domain%ustar        (ims+1:ime-1,jms+1:jme-1) = v_mass      (ims+1:ime-1,kms,jms+1:jme-1) * currw
            domain%v_10m%data_2d(ims+1:ime-1,jms+1:jme-1) = domain%ustar(ims+1:ime-1,jms+1:jme-1)     * lastw
        endif

        if (allocated(domain%ustar)) then
            ! now calculate master ustar based on U and V combined in quadrature
            domain%ustar(ims+1:ime-1,jms+1:jme-1) = sqrt(u_mass(ims+1:ime-1,kms,jms+1:jme-1)**2 + v_mass(ims+1:ime-1,kms,jms+1:jme-1)**2) * currw
        endif

        ! finally, calculate the real vertical motions (including U*dzdx + V*dzdy)
        lastw = 0
        do z = kms, kme

            ! ! if(options%parameters%use_terrain_difference) then
            !                 ! compute the U * dz/dx component of vertical motion
            !     uw    = u(ims+1:ime,   z, jms+1:jme-1) * domain%delta_dzdx(:,z,jms+1:jme-1)
            !     ! compute the V * dz/dy component of vertical motion
            !     vw    = v(ims+1:ime-1, z, jms+1:jme  ) * domain%delta_dzdy(ims+1:ime-1,z,:)
            ! else
                ! compute the U * dz/dx component of vertical motion
                uw    = u(ims+1:ime,   z, jms+1:jme-1) * domain%dzdx(ims+1:ime,z,jms+1:jme-1)
                ! compute the V * dz/dy component of vertical motion
                vw    = v(ims+1:ime-1, z, jms+1:jme  ) * domain%dzdy(ims+1:ime-1,z,jms+1:jme)
            ! endif
            ! ! convert the W grid relative motion to m/s
            ! currw = w(ims+1:ime-1, z, jms+1:jme-1) * dz_interface(ims+1:ime-1, z, jms+1:jme-1) / domain%dx

            ! the W grid relative motion
            currw = w(ims+1:ime-1, z, jms+1:jme-1)

            ! if (options%physics%convection>0) then
            !     currw = currw + domain%w_cu(2:nx-1,z,2:ny-1) * domain%dz_inter(2:nx-1,z,2:ny-1) / domain%dx
            ! endif

            ! compute the real vertical velocity of air by combining the different components onto the mass grid
            ! includes vertical interpolation between w_z-1/2 and w_z+1/2
            w_real(ims+1:ime-1, z, jms+1:jme-1) = (uw(ims+1:ime-1,:) + uw(ims+2:ime,:))*0.5 &
                                                 +(vw(:,jms+1:jme-1) + vw(:,jms+2:jme))*0.5 &
                                                 +domain%jacobian(ims+1:ime-1,z,jms+1:jme-1)*(lastw + currw) * 0.5
            lastw = currw ! could avoid this memcopy cost using pointers or a single manual loop unroll
        end do

        end associate

    end subroutine diagnostic_update


    !>------------------------------------------------------------
    !!  Calculate the maximum stable time step given some CFL criteria
    !!
    !!  For each grid cell, find the mean of the wind speeds from each
    !!  direction * sqrt(3) for the 3D advection CFL limited time step
    !!  Also find the maximum wind speed anywhere in the domain to check
    !!  against a 1D advection limit.
    !!
    !! @param dx  [ scalar ]        horizontal grid cell width  [m]
    !! @param u   [nx+1 x nz x ny]  east west wind speeds       [m/s]
    !! @param v   [nx x nz x ny+1]  North South wind speed      [m/s]
    !! @param w   [nx x nz x ny]    vertical wind speed         [m/s]
    !! @param CFL [ scalar ]        CFL limit to use (e.g. 1.0)
    !! @return dt [ scalar ]        Maximum stable time step    [s]
    !!
    !!------------------------------------------------------------
    function compute_dt(dx, u, v, w, rho, dz, CFL, cfl_strictness, use_density) result(dt)
        real,       intent(in)                   :: dx
        real,       intent(in), dimension(:,:,:) :: u, v, w, rho
        real,       intent(in), dimension(:)     :: dz
        real,       intent(in)                   :: CFL
        integer,    intent(in)                   :: cfl_strictness
        logical,    intent(in)                   :: use_density
        ! output value
        real :: dt
        ! locals
        real :: three_d_cfl = 0.577350269 ! = sqrt(3)/3
        integer :: i, j, k, nx, nz, ny, zoffset
        real :: maxwind3d, maxwind1d, current_wind, sqrt3

        sqrt3 = sqrt(3.0) * 1.001 ! with a safety factor

        maxwind1d = 0
        maxwind3d = 0

        nx = size(w,1)
        nz = size(w,2)
        ny = size(w,3)

        if (cfl_strictness==1) then
            ! to ensure we are stable for 1D advection:
            if (use_density) then
                !maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                !maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
            else
                maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                maxwind1d = max( maxwind1d, maxval(abs(w)))
            endif

            maxwind3d = maxwind1d * sqrt3
        else if (cfl_strictness==5) then

            if (use_density) then
                !maxwind1d = maxval(abs(u(2:,:,:) / (rho*dz*dx) )) &
                !          + maxval(abs(v(:,:,2:) / (rho*dz*dx) )) &
                !          + maxval(abs(w(:,:, :) / (rho*dz*dx) ))
            else
                maxwind3d = maxval(abs(u)) + maxval(abs(v)) + maxval(abs(w))
            endif

        else
            ! to ensure we are stable for 3D advection we'll use the average "max" wind speed
            ! but that average has to be divided by sqrt(3) for stability in 3 dimensional advection
            do j=1,ny
                do k=1,nz
                    if (k==1) then
                        zoffset = 0
                    else
                        zoffset = -1
                    endif

                    do i=1,nx
                        ! just compute the sum of the wind speeds, but take the max over the two
                        ! faces of the grid cell (e.g. east and west sides)
                        ! this will be divided by 3 later by three_d_cfl
                        if (use_density) then
                            !current_wind = (max(abs(u(i,k,j)), abs(u(i+1,k,j))) &
                            !              + max(abs(v(i,k,j)), abs(v(i,k,j+1))) &
                            !              + max(abs(w(i,k,j)), abs(w(i,k+zoffset,j))) ) &
                            !              / (rho(i,k,j) * dz(i,k,j) * dx)
                        else
                            current_wind = max(abs(u(i,k,j)), abs(u(i+1,k,j))) / dx &
                                          +max(abs(v(i,k,j)), abs(v(i,k,j+1))) / dx &
                                          +max(abs(w(i,k,j)), abs(w(i,k+zoffset,j))) / dz(k)
                        endif
                        maxwind3d = max(maxwind3d, current_wind)
                    ENDDO
                ENDDO
            ENDDO

            if (cfl_strictness==2) then
                ! effectively divides by 3 to take the mean and multiplies by the sqrt(3) for the 3D advection limit
                maxwind3d = maxwind3d * three_d_cfl

                ! to ensure we are stable for 1D advection:
                if (use_density) then
                    !maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                    !maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
                else
                    maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                    maxwind1d = max( maxwind1d, maxval(abs(w)))
                endif
                ! also insure stability for 1D advection
                maxwind3d = max(maxwind1d,maxwind3d)

            ! else if (cfl_strictness==3) then
            !   leave maxwind3d as the sum of the max winds
            ! This should be the default, does it need to be multiplied by sqrt(3)?
            elseif (cfl_strictness==4) then
                maxwind3d = maxwind3d * sqrt3
            endif

        endif

        !TESTING: Do we need to multiply maxwind3d by sqrt3 as the comment above suggests?
        ! maxwind3d = maxwind3d * sqrt3

        dt = CFL / maxwind3d

        ! If we have too small a time step throw an error
        ! something is probably wrong in the physics or input data
        if (dt<1e-1) then
            write(*,*) "dt   = ", dt
            write(*,*) "Umax = ", maxval(abs(u))
            write(*,*) "Vmax = ", maxval(abs(v))
            write(*,*) "Wmax = ", maxval(abs(w))
            stop "ERROR time step too small"
        endif

    end function compute_dt


    !>------------------------------------------------------------
    !!  Prints progress to the terminal if requested
    !!
    !! @param current_time  the current state of the model time
    !! @param end_time      the end of the current full time step (when step will return)
    !! @param time_step     length of full time step to be integrated by step
    !! @param dt            numerical timestep to print for information
    !!
    !!------------------------------------------------------------
    subroutine print_progress(current_time, end_time, time_step, dt, last_time)
        implicit none
        type(Time_type),    intent(in)    :: current_time,    end_time
        type(time_delta_t), intent(in)    :: time_step,       dt
        real,               intent(inout) :: last_time

        type(time_delta_t) :: progress_dt
        real :: time_percent

        ! first compute the current time until reaching the end
        progress_dt  = (end_time - current_time)

        ! convert that to a percentage of the total time required
        time_percent = 100 - progress_dt%seconds() / time_step%seconds()  * 100

        ! finally if it has been at least 5% of the time since the last time we printed output, print output
        if (time_percent > (last_time + 5.0)) then
            last_time = last_time + 5.0
            ! this used to just use the nice $ (or advance="NO") trick, but at least with some mpi implementations, it buffers this output until it crashes
            write(*,"(A,f5.1,A,A)") char(13), max(0.0, time_percent)," %  dt=",trim(dt%as_string())
        endif

    end subroutine print_progress

    !>------------------------------------------------------------
    !! Update the numerical timestep to use
    !!
    !! @param dt            numerical timestep to use
    !! @param options       set options for how to update the time step
    !! @param domain        the full domain structure (need winds to compute dt)
    !! @param end_time      the end of the current full time step (when step will return)
    !!
    !!------------------------------------------------------------
    subroutine update_dt(dt, options, domain, end_time)
        implicit none
        type(time_delta_t), intent(inout) :: dt
        type(options_t),    intent(in)    :: options
        type(domain_t),     intent(in)    :: domain
        type(Time_type),    intent(in)    :: end_time

        double precision :: seconds

        ! compute internal timestep dt to maintain stability
        ! courant condition for 3D advection. Note that w is normalized by dx/dz

        ! Note this needs to be performed when advect_density is enabled
        ! if (options%parameters%advect_density) then
            ! call dt%set(seconds=compute_dt(domain%dx, domain%ur, domain%vr, domain%wr, domain%rho, domain%dz_inter,&
            !                 options%parameters%cfl_reduction_factor, cfl_strictness=options%parameters%cfl_strictness,                   &
            !                 use_density=.True.))

        ! else
        ! compute the dt to meet the CFL criteria specified given dx, u, v, w, dz
        associate(dx         => domain%dx,                              &
                  u          => domain%u%data_3d,                       &
                  v          => domain%v%data_3d,                       &
                  w          => domain%w%data_3d,                       &
                  density    => domain%density%data_3d,                 &
                  dz         => options%parameters%dz_levels,           &
                  reduction  => options%parameters%cfl_reduction_factor,&
                  strictness => options%parameters%cfl_strictness       &
            )

            seconds = compute_dt(dx, u, v, w, density, dz, reduction, &
                                 cfl_strictness=strictness, use_density=.false.)

        end associate
        ! endif

        ! perform a reduction across all images to find the minimum time step required
#ifndef __INTEL_COMPILER
        call co_min(seconds)
#endif
#ifdef __INTEL_COMPILER
        seconds = domain%dx / 100
#endif

        ! set an upper bound on dt to keep microphysics and convection stable (?)
        ! store this back in the dt time_delta data structure
        call dt%set(seconds=min(seconds,120.0D0))

    end subroutine update_dt


    !>------------------------------------------------------------
    !!  Step forward one IO time step.
    !!
    !!  Calculated the internal model time step to satisfy the CFL criteria,
    !!  then updates all forcing update increments for that dt and loops through
    !!  time calling physics modules.
    !!  Also checks to see if it is time to write a model output file.
    !!
    !! @param domain    domain data structure containing model state
    !! @param options   model options structure
    !! @param bc        model boundary conditions data structure
    !! @param next_output   Next time to write an output file (in "model_time")
    !!
    !!------------------------------------------------------------
    subroutine step(domain, end_time, options)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(Time_type),    intent(in)      :: end_time
        type(options_t),    intent(in)      :: options

        real :: last_print_time
        type(time_delta_t) :: dt, time_step_size
        integer :: ims, jms, kms
        integer :: test_x, test_y, test_z, test_image

        test_x = 87
        test_z = 3
        test_y = 198
        test_image = 33

        ims = domain%grid%ims
        jms = domain%grid%jms
        kms = domain%grid%kms

        last_print_time = 0.0
        time_step_size = end_time - domain%model_time

        ! now just loop over internal timesteps computing all physics in order (operator splitting...)
        do while (domain%model_time < end_time)

            call update_dt(dt, options, domain, end_time)

            ! Make sure we don't over step the forcing or output period
            if ((domain%model_time + dt) > end_time) then
                dt = end_time - domain%model_time
            endif

            ! ensure internal model consistency
            call diagnostic_update(domain, options)


            ! if an interactive run was requested than print status updates everytime at least 5% of the progress has been made
            if (options%parameters%interactive .and. (this_image()==1)) then
                call print_progress(domain%model_time, end_time, time_step_size, dt, last_print_time)
            endif

            ! this if is to avoid round off errors causing an additional physics call that won't really do anything
            if (dt%seconds() > 1e-3) then

                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" init", fix=.True.)

                ! first process the halo section of the domain (currently hard coded at 1 should come from domain?)
                call rad(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" rad(domain", fix=.True.)

                call lsm(domain, options, real(dt%seconds()))!, halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" lsm")

                call pbl(domain, options, real(dt%seconds()))!, halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" pbl")

                call convect(domain, options, real(dt%seconds()))!, halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" convect")

                call mp(domain, options, real(dt%seconds()), halo=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" mp_halo", fix=.True.)

                call domain%halo_send()
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%halo_send", fix=.True.)

                ! call rad(domain, options, real(dt%seconds()), subset=1)
                ! call lsm(domain, options, real(dt%seconds()))!, subset=1)
                ! call pbl(domain, options, real(dt%seconds()))!, subset=1)
                ! call convect(domain, options, real(dt%seconds()), subset=1)

                call mp(domain, options, real(dt%seconds()), subset=1)
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" mp(domain", fix=.True.)

                call domain%halo_retrieve()
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%halo_retrieve", fix=.True.)

                call advect(domain, options, real(dt%seconds()))
                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" advect(domain", fix=.True.)


                ! ! apply/update boundary conditions including internal wind and pressure changes.
                call domain%apply_forcing(dt)

                !If we are in the last 2 updates of a time step and a variable drops below 0, we have probably over-shot a value of 0. Force back to 0
                if ((end_time%seconds() - domain%model_time%seconds()) < (dt%seconds()*2)) then
                    call domain%enforce_limits()
                endif


                if (options%parameters%debug) call domain_check(domain, "img: "//trim(str(this_image()))//" domain%apply_forcing", fix=.True.)

            endif

            ! step model_time forward
            domain%model_time = domain%model_time + dt

        enddo

    end subroutine step
end module time_step
