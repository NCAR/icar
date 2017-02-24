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
    use microphysics,               only : mp
    use convection,                 only : convect
    use land_surface,               only : lsm
    use wind,                       only : balance_uvw
    use advection,                  only : advect
    use output,                     only : write_domain
    use planetary_boundary_layer,   only : pbl
    use radiation,                  only : rad
    use boundary_conditions,        only : update_pressure
    use debug_module,               only : domain_check
    implicit none
    private
    public :: step

    real, dimension(:,:), allocatable :: lastw, currw, uw, vw !> temporaries used to compute diagnostic w_real field
contains
    
    !>------------------------------------------------------------
    !!  Update the edges of curdata by adding dXdt
    !!
    !!  Apply dXdt to the boundaries of a given data array. 
    !!  For these variables, dXdt is a small array just storing the boundary increments
    !!
    !! @param curdata   data array to be updated
    !! @param dXdt      Array containing increments to be applied to the boundaries. 
    !!                  (nz x max(nx,ny) x 4)
    !!
    !!------------------------------------------------------------
    subroutine boundary_update(curdata,dXdt)
        implicit none
        real,dimension(:,:,:), intent(inout) :: curdata
        real,dimension(:,:,:), intent(in) :: dXdt
        integer::nx,nz,ny,i
        
        nx=size(curdata,1)
        nz=size(curdata,2)
        ny=size(curdata,3)

        do i=1,nz
            curdata(1, i,2:ny-1) = curdata(1, i, 2:ny-1) + dXdt(i,2:ny-1,1)
            curdata(nx,i,2:ny-1) = curdata(nx,i, 2:ny-1) + dXdt(i,2:ny-1,2)
            curdata(:, i,1 )     = curdata(:, i, 1)      + dXdt(i,1:nx,  3)
            curdata(:, i,ny)     = curdata(:, i, ny)     + dXdt(i,1:nx,  4)
        enddo
        ! correct possible rounding errors, primarily an issue of clouds...
        where(curdata(1, :,: )<0) curdata(1, :,: ) = 0
        where(curdata(nx,:,: )<0) curdata(nx,:,: ) = 0
        where(curdata(:, :,1 )<0) curdata(:, :,1 ) = 0
        where(curdata(:, :,ny)<0) curdata(:, :,ny) = 0
        
    end subroutine boundary_update
    
    !>------------------------------------------------------------
    !!  Updated all fields in domain using the respective dXdt variables
    !!
    !!  
    !!
    !! @param domain    full domain data structure to be updated
    !! @param bc        boundary conditions (containing dXdt increments)
    !! @param options   options structure specifies which variables need to be updated
    !!
    !!------------------------------------------------------------
    subroutine forcing_update(domain, bc, options, dt)
        implicit none
        type(domain_type),  intent(inout) :: domain
        type(bc_type),      intent(inout) :: bc
        type(options_type), intent(in)    :: options
        real,               intent(in)    :: dt
        integer :: j,ny
        
        ny = size(domain%p, 3)
        
        !$omp parallel firstprivate(ny, dt) &
        !$omp private(j) &
        !$omp shared(bc, domain)
        !$omp do schedule(static)
        do j=1,ny
            domain%u(:,:,j) = domain%u(:,:,j) + (bc%du_dt(:,:,j) * dt)
            domain%v(:,:,j) = domain%v(:,:,j) + (bc%dv_dt(:,:,j) * dt)
            if (.not.options%ideal) then
                domain%p(:,:,j)   =  domain%p(:,:,j) + (bc%dp_dt(:,:,j) * dt)
                ! update the exner function and model density while we are at it should this be in diagnostic_update(?)
                domain%pii(:,:,j) = (domain%p(:,:,j) / 100000.0)**(Rd/cp)
                domain%rho(:,:,j) =  domain%p(:,:,j) / (Rd * domain%th(:,:,j) * domain%pii(:,:,j)) ! kg/m^3
            endif
            
            ! these only get updated if we are using the fluxes derived from the forcing model
            if (options%physics%landsurface==kLSM_BASIC) then
                domain%sensible_heat(:,j)  = domain%sensible_heat(:,j)+ (bc%dsh_dt(:,j) * dt)
                domain%latent_heat(:,j)    = domain%latent_heat(:,j)  + (bc%dlh_dt(:,j) * dt)
            endif
            
            ! these only get updated if we are using the fluxes (and PBL height) derived from the forcing model
            if (options%physics%boundarylayer==kPBL_BASIC) then
                domain%pbl_height(:,j) = domain%pbl_height(:,j)+ (bc%dpblh_dt(:,j) * dt)
            endif
            
            ! these only get updated if we are using the fluxes derived from the forcing model
            if (options%physics%radiation==kRA_BASIC) then
                domain%swdown(:,j)  = domain%swdown(:,j)  + (bc%dsw_dt(:,j) * dt)
                domain%lwdown(:,j)  = domain%lwdown(:,j)  + (bc%dlw_dt(:,j) * dt)
            endif
            domain%sst(:,j)  = domain%sst(:,j)  + (bc%dsst_dt(:,j) * dt)

            if (trim(options%rain_var)/="") then
                domain%rain(:,j) = domain%rain(:,j) + (bc%drain_dt(:,j) * dt)
            endif

        enddo
        !$omp end do
        !$omp end parallel
        ! v has one more y element than others
        ny = ny + 1
        domain%v(:,:,ny) = domain%v(:,:,ny) + (bc%dv_dt(:,:,ny) * dt)
        ! dXdt for qv,qc,th are only applied to the boundarys
        if (.not.options%ideal) then
            call boundary_update(domain%th, bc%dth_dt * dt)
            call boundary_update(domain%qv, bc%dqv_dt * dt)
            call boundary_update(domain%cloud, bc%dqc_dt * dt)
        endif
        
        ! because density changes with each time step, u/v/w have to be rebalanced as well. 
        ! could avoid this by assuming density doesn't change... but would need to keep an "old" density around
        ! also, convection can modify u and v so it needs to rebalanced
        call balance_uvw(domain,options)
        
        call diagnostic_update(domain,options)
    end subroutine forcing_update
    
    !>------------------------------------------------------------
    !! Update model diagnostic fields
    !!
    !! Calculates most model diagnostic fields such as Psfc, 10m height winds and ustar
    !!
    !! @param domain    Model domain data structure to be updated
    !! @param options   Model options (not used at present)
    !!
    !!------------------------------------------------------------
    subroutine diagnostic_update(domain,options)
        implicit none
        type(domain_type),  intent(inout)   :: domain
        type(options_type), intent(in)      :: options
        
        integer :: nx,ny,nz, y, z
        
        nx=size(domain%p,1)
        nz=size(domain%p,2)
        ny=size(domain%p,3)
        
        ! update p_inter, psfc, ptop, Um, Vm, mut
        domain%Um = 0.5*(domain%u(1:nx-1,:,:)+domain%u(2:nx,:,:))
        domain%Vm = 0.5*(domain%v(:,:,1:ny-1)+domain%v(:,:,2:ny))
        domain%t  = domain%th * domain%pii
        
        domain%p_inter=domain%p
        call update_pressure(domain%p_inter, domain%z, domain%z_inter, domain%t)
        domain%psfc = domain%p_inter(:,1,:)
        ! technically this isn't correct, we should be using update_pressure or similar to solve this
        domain%ptop = 2*domain%p(:,nz,:) - domain%p(:,nz-1,:)
        
        ! dry mass in the gridcell is equivalent to the difference in pressure from top to bottom
        domain%mut(:,1:nz-1,:) = domain%p_inter(:,1:nz-1,:) - domain%p_inter(:,2:nz,:)
        domain%mut(:,nz,:) = domain%p_inter(:,nz,:) - domain%ptop
        
        if (.not.allocated(lastw)) then
            allocate(lastw(nx-2,ny-2))
            allocate(currw(nx-2,ny-2))
            allocate(uw(nx-1,ny-2))
            allocate(vw(nx-2,ny-1))
        endif
        
        ! temporary constant
        ! use log-law of the wall to convert from first model level to surface
        currw = karman / log((domain%z(2:nx-1,1,2:ny-1)-domain%terrain(2:nx-1,2:ny-1)) / domain%znt(2:nx-1,2:ny-1))
        ! use log-law of the wall to convert from surface to 10m height
        lastw = log(10.0 / domain%znt(2:nx-1,2:ny-1)) / karman
        domain%ustar(2:nx-1,2:ny-1) = domain%Um   (2:nx-1,1,2:ny-1) * currw
        domain%u10  (2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1)   * lastw
        domain%ustar(2:nx-1,2:ny-1) = domain%Vm   (2:nx-1,1,2:ny-1) * currw
        domain%v10  (2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1)   * lastw
        
        ! now calculate master ustar based on U and V combined in quadrature
        domain%ustar(2:nx-1,2:ny-1) = sqrt(domain%Um(2:nx-1,1,2:ny-1)**2 + domain%Vm(2:nx-1,1,2:ny-1)**2) * currw
        
        ! finally, calculate the real vertical motions (including U*dzdx + V*dzdy)
        lastw=0
        do z=1,nz
            ! compute the U * dz/dx component of vertical motion
            uw    = domain%u(2:nx,  z,2:ny-1) * domain%dzdx(:,2:ny-1)
            ! compute the V * dz/dy component of vertical motion
            vw    = domain%v(2:nx-1,z,2:ny  ) * domain%dzdy(2:nx-1,:)
            ! convert the W grid relative motion to m/s
            currw = domain%w(2:nx-1,z,2:ny-1) * domain%dz_inter(2:nx-1,z,2:ny-1) / domain%dx
            
            ! compute the real vertical velocity of air by combining the different components onto the mass grid
            domain%w_real(2:nx-1,z,2:ny-1) = (uw(1:nx-2,:) + uw(2:nx-1,:))*0.5 &
                                            +(vw(:,1:ny-2) + vw(:,2:ny-1))*0.5 &
                                            +(lastw + currw) * 0.5
                                            
            lastw=currw ! could avoid this memcopy cost using pointers or a single manual loop unroll
        end do
        
    end subroutine diagnostic_update


    !>------------------------------------------------------------
    !!  Uses the dt from step() to convert the forcing increments to per/time step increments
    !!
    !!  Divides dXdt variables by the length of the input timestep so that it can be multiplied by the internal
    !!  timestep to get the each update value. 
    !!
    !! @param bc        Boundary conditions structure containing dXdt variables
    !! @param dt        Length of time between input forcing data for which the dXdt updates were calculated
    !! @param options   Model options structure specifies which variables need to be updated
    !!
    !!------------------------------------------------------------
    subroutine apply_dt(bc, dt, options)
        implicit none
        type(bc_type),      intent(inout)   :: bc
        real,               intent(in)      :: dt
        type(options_type), intent(in)      :: options
        ! internal variables
        integer :: j, ny, nx
        
        ny=size(bc%du_dt,3)
        nx=size(bc%du_dt,1)

        bc%du_dt  = bc%du_dt  / dt
        bc%dv_dt  = bc%dv_dt  / dt
        bc%dp_dt  = bc%dp_dt  / dt
        bc%dth_dt = bc%dth_dt / dt
        bc%dqv_dt = bc%dqv_dt / dt
        bc%dqc_dt = bc%dqc_dt / dt
        
        if (trim(options%rain_var)/="") then
            bc%drain_dt = bc%drain_dt / dt
        endif

        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%landsurface==kLSM_BASIC) then
            bc%dsh_dt   = bc%dsh_dt   / dt
            bc%dlh_dt   = bc%dlh_dt   / dt
        endif
        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%boundarylayer==kPBL_BASIC) then
            bc%dpblh_dt = bc%dpblh_dt / dt
        endif
        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%radiation==kRA_BASIC) then
            bc%dsw_dt   = bc%dsw_dt   / dt
            bc%dlw_dt   = bc%dlw_dt   / dt
        endif
        bc%dsst_dt   = bc%dsst_dt   / dt
            
    end subroutine apply_dt

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
        real,       intent(in), dimension(:,:,:) :: u, v, w, rho, dz
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
        
        nx = size(rho,1)
        nz = size(rho,2)
        ny = size(rho,3)

        if (cfl_strictness==1) then
            ! to ensure we are stable for 1D advection:
            if (use_density) then
                maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
            else
                maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                maxwind1d = max( maxwind1d, maxval(abs(w)))
            endif
                        
            maxwind3d = maxwind1d * sqrt3
        else if (cfl_strictness==5) then
            
            if (use_density) then
                maxwind1d = maxval(abs(u(2:,:,:) / (rho*dz*dx) )) &
                          + maxval(abs(v(:,:,2:) / (rho*dz*dx) )) &
                          + maxval(abs(w(:,:, :) / (rho*dz*dx) ))
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
                            current_wind = (max(abs(u(i,k,j)), abs(u(i+1,k,j))) &
                                          + max(abs(v(i,k,j)), abs(v(i,k,j+1))) & 
                                          + max(abs(w(i,k,j)), abs(w(i,k+zoffset,j))) ) &
                                          / (rho(i,k,j) * dz(i,k,j) * dx)
                        else
                            current_wind = max(abs(u(i,k,j)), abs(u(i+1,k,j))) &
                                          +max(abs(v(i,k,j)), abs(v(i,k,j+1))) & 
                                          +max(abs(w(i,k,j)), abs(w(i,k+zoffset,j)))
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
                    maxwind1d = max( maxval(abs(u(2:,:,:) / (rho*dz*dx) )), maxval(abs(v(:,:,2:) / (rho*dz*dx))) )
                    maxwind1d = max( maxwind1d, maxval(abs(w/(rho*dz*dx))) )
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
        
        dt = CFL * dx / maxwind3d
        
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
    !! @param model_time    Current internal model time (seconds since start of run)
    !! @param next_output   Next time to write an output file (in "model_time")
    !!
    !!------------------------------------------------------------
    subroutine step(domain,options,bc,model_time,next_output)
        implicit none
        type(domain_type),  intent(inout)   :: domain
        type(bc_type),      intent(inout)   :: bc
        type(options_type), intent(in)      :: options
        real*8,             intent(inout)   :: model_time, next_output
        
        real*8  :: end_time
        real    :: dt, CFL, cfl_reduction_factor
        
        CFL = 1.0
        cfl_reduction_factor = options%cfl_reduction_factor
        CFL = CFL * cfl_reduction_factor
        
        ! calculate the number of timesteps
        end_time = model_time + options%in_dt
        
        ! Make the boundary condition dXdt values into units of [X]/s
        call apply_dt(bc, options%in_dt, options)
        
        ! ensure internal model consistency (should only need to be called here when the model starts...)
        ! e.g. for potential temperature and temperature
        call diagnostic_update(domain,options)
        
        ! now just loop over internal timesteps computing all physics in order (operator splitting...)
        do while (model_time < end_time)
            
            ! compute internal timestep dt to maintain stability
            ! courant condition for 3D advection. Note that w is normalized by dx/dz
            ! pick the minimum dt from the begining or the end of the current timestep
            if (options%advect_density) then
                dt = compute_dt(domain%dx, domain%ur, domain%vr, domain%wr, domain%rho, domain%dz_inter, &
                                CFL, cfl_strictness=options%cfl_strictness, use_density=options%advect_density)
            else
                dt = compute_dt(domain%dx, domain%u, domain%v, domain%w, domain%rho, domain%dz_inter, &
                                CFL, cfl_strictness=options%cfl_strictness, use_density=options%advect_density)
            endif
            ! set an upper bound on dt to keep microphysics and convection stable (?) not sure what time is required here. 
            dt = min(dt,120.0) !better min=180?
            if (options%interactive) then
                write(*,"(A,f5.1,A,f5.1,A$)") char(13), 100-max(0.0,(end_time-model_time-dt)/options%in_dt*100)," %  dt=",dt,"s  "
            endif
            ! Make sure we don't over step the forcing period
            if ((model_time + dt) > end_time) then
                dt = end_time - model_time
            endif
            
            ! this if is to avoid round off errors causing an additional physics call that won't really do anything
            if (dt > 1e-5) then
                if (options%debug) call domain_check(domain,"Time step loop start")
                
                call advect(domain,options,dt)
                if (options%debug) call domain_check(domain,"After advection")
                
                call mp(domain,options,dt, model_time)
                if (options%debug) call domain_check(domain,"After microphysics")
                
                call rad(domain,options,model_time/86400.0+50000, dt)
                if (options%debug) call domain_check(domain,"After radiation")
                
                call lsm(domain,options,dt,model_time)
                if (options%debug) call domain_check(domain,"After LSM")
                
                call pbl(domain,options,dt)
                if (options%debug) call domain_check(domain,"After PBL")
                
                call convect(domain,options,dt)
                if (options%debug) call domain_check(domain,"After Convection")
                
                ! apply/update boundary conditions including internal wind and pressure changes. 
                call forcing_update(domain, bc, options, dt)
                if (options%debug) call domain_check(domain,"After Forcing update")

            endif
            
            ! step model_time forward
            model_time = model_time + dt
            domain%model_time = model_time
            
            if ((abs(model_time-next_output)<1e-1).or.(model_time>next_output)) then
                call write_domain(domain,options,nint((model_time-options%time_zero)/options%out_dt))
                next_output = next_output + options%out_dt
            endif
            
        enddo
        
    end subroutine step
end module time_step
