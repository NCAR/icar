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
    use data_structures     ! *_type  types
    use microphysics,               only : mp
    use convection,                 only : convect
    use land_surface,               only : lsm
    use wind,                       only : balance_uvw
    use advection,                  only : advect
    use output,                     only : write_domain
    use planetary_boundary_layer,   only : pbl
    use radiation,                  only : rad
    use boundary_conditions,        only : update_pressure
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
            curdata(1,i,2:ny-1)  = curdata(1,i,2:ny-1)  + dXdt(i,2:ny-1,1)
            curdata(nx,i,2:ny-1) = curdata(nx,i,2:ny-1) + dXdt(i,2:ny-1,2)
            curdata(:,i,1)  = curdata(:,i,1)  + dXdt(i,1:nx,3)
            curdata(:,i,ny) = curdata(:,i,ny) + dXdt(i,1:nx,4)
        enddo
        ! correct possible rounding errors, primarily an issue of clouds...
        where(curdata(1,:,:)<0) curdata(1,:,:)=0
        where(curdata(nx,:,:)<0) curdata(nx,:,:)=0
        where(curdata(:,:,1)<0) curdata(:,:,1)=0
        where(curdata(:,:,ny)<0) curdata(:,:,ny)=0
        
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
    subroutine forcing_update(domain,bc,options)
        implicit none
        type(domain_type),intent(inout)::domain
        type(bc_type),intent(inout)::bc
        type(options_type),intent(in)::options
        integer::j,ny
        ny=size(domain%p,3)
        !$omp parallel firstprivate(ny) &
        !$omp private(j) &
        !$omp shared(bc,domain)
        !$omp do schedule(static)
        do j=1,ny
            domain%u(:,:,j)=domain%u(:,:,j)+bc%du_dt(:,:,j)
            domain%v(:,:,j)=domain%v(:,:,j)+bc%dv_dt(:,:,j)
            if (.not.options%ideal) then
                domain%p(:,:,j)=domain%p(:,:,j)+bc%dp_dt(:,:,j)
                domain%pii(:,:,j)=(domain%p(:,:,j)/100000.0)**(Rd/cp)
                domain%rho(:,:,j)=domain%p(:,:,j)/(Rd*domain%th(:,:,j)*domain%pii(:,:,j)) ! kg/m^3
            endif
            
            ! these only get updated if we are using the fluxes derived from the forcing model
            if (options%physics%landsurface==kLSM_BASIC) then
                domain%sensible_heat(:,j)  = domain%sensible_heat(:,j)+ bc%dsh_dt(:,j)
                domain%latent_heat(:,j)    = domain%latent_heat(:,j)  + bc%dlh_dt(:,j)
            endif
            
            ! these only get updated if we are using the fluxes (and PBL height) derived from the forcing model
            if (options%physics%boundarylayer==kPBL_BASIC) then
                domain%pbl_height(:,j) = domain%pbl_height(:,j)+ bc%dpblh_dt(:,j)
            endif
            
            ! these only get updated if we are using the fluxes derived from the forcing model
            if (options%physics%radiation==kRA_BASIC) then
                domain%swdown(:,j)  = domain%swdown(:,j)  + bc%dsw_dt(:,j)
                domain%lwdown(:,j)  = domain%lwdown(:,j)  + bc%dlw_dt(:,j)
            endif
            domain%sst(:,j)  = domain%sst(:,j)  + bc%dsst_dt(:,j)
        enddo
        !$omp end do
        !$omp end parallel
        ! v has one more y element than others
        ny=ny+1
        domain%v(:,:,ny)=domain%v(:,:,ny)+bc%dv_dt(:,:,ny)
        ! dXdt for qv,qc,th are only applied to the boundarys
        if (.not.options%ideal) then
            call boundary_update(domain%th,bc%dth_dt)
            call boundary_update(domain%qv,bc%dqv_dt)
            call boundary_update(domain%cloud,bc%dqc_dt)
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
        type(domain_type), intent(inout) :: domain
        type(options_type), intent(in) :: options
        
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
        domain%ustar(2:nx-1,2:ny-1) = domain%Um(2:nx-1,1,2:ny-1) * currw
        domain%u10(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1) * lastw
        domain%ustar(2:nx-1,2:ny-1) = domain%Vm(2:nx-1,1,2:ny-1) * currw
        domain%v10(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1) * lastw
        
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
    !!  Divides dXdt variables by n timesteps so that after adding it N times we will be at the
    !!  correct final value. 
    !!
    !! @param bc        Boundary conditions structure containing dXdt variables
    !! @param nsteps    Number of time steps calculated in step() function
    !! @param options   Model options structure specifies which variables need to be updated
    !!
    !!------------------------------------------------------------
    subroutine apply_dt(bc,nsteps, options)
        implicit none
        type(bc_type), intent(inout) :: bc
        integer,intent(in)::nsteps
        type(options_type), intent(in) :: options
        integer::j, ny, nx
        ny=size(bc%du_dt,3)
        nx=size(bc%du_dt,1)

        bc%du_dt  = bc%du_dt  / nsteps
        bc%dv_dt  = bc%dv_dt  / nsteps
        bc%dp_dt  = bc%dp_dt  / nsteps
        bc%dth_dt = bc%dth_dt / nsteps
        bc%dqv_dt = bc%dqv_dt / nsteps
        bc%dqc_dt = bc%dqc_dt / nsteps

        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%landsurface==kLSM_BASIC) then
            bc%dsh_dt   = bc%dsh_dt   / nsteps
            bc%dlh_dt   = bc%dlh_dt   / nsteps
        endif
        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%boundarylayer==kPBL_BASIC) then
            bc%dpblh_dt = bc%dpblh_dt / nsteps
        endif
        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%radiation==kRA_BASIC) then
            bc%dsw_dt   = bc%dsw_dt   / nsteps
            bc%dlw_dt   = bc%dlw_dt   / nsteps
        endif
        bc%dsst_dt   = bc%dsst_dt   / nsteps

        ! parallel version didn't seem to work
!         !$omp parallel firstprivate(ny,nsteps) &
!         !$omp private(j) &
!         !$omp shared(bc)
!         !$omp do schedule(static)
!         do j=1,ny
!             bc%du_dt(:,:,j)  = bc%du_dt(:,:,j)  / nsteps
!             bc%dv_dt(:,:,j)  = bc%dv_dt(:,:,j)  / nsteps
!             bc%dp_dt(:,:,j)  = bc%dp_dt(:,:,j)  / nsteps
!             bc%dth_dt(:,j,:) = bc%dth_dt(:,j,:) / nsteps
!             bc%dqv_dt(:,j,:) = bc%dqv_dt(:,j,:) / nsteps
!             bc%dqc_dt(:,j,:) = bc%dqc_dt(:,j,:) / nsteps
!
!             ! these only get updated if we are using the fluxes derived from the forcing model
!             if (options%physics%landsurface==kLSM_BASIC) then
!                 bc%dsh_dt(:,j)   = bc%dsh_dt(:,j)   / nsteps
!                 bc%dlh_dt(:,j)   = bc%dlh_dt(:,j)   / nsteps
!             endif
!             ! these only get updated if we are using the fluxes derived from the forcing model
!             if (options%physics%boundarylayer==kPBL_BASIC) then
!                 bc%dpblh_dt(:,j) = bc%dpblh_dt(:,j) / nsteps
!             endif
!             ! these only get updated if we are using the fluxes derived from the forcing model
!             if (options%physics%radiation==kRA_BASIC) then
!                 bc%dsw_dt(:,j)   = bc%dsw_dt(:,j)   / nsteps
!                 bc%dlw_dt(:,j)   = bc%dlw_dt(:,j)   / nsteps
!             endif
!             bc%dsw_dt(:,j)   = bc%dsw_dt(:,j)       / nsteps
!             
!         end do
!         !$omp end do
!         !$omp end parallel
!         if (nx>ny) then
!             do j=ny+1,nx
!                 bc%dth_dt(:,j,:) = bc%dth_dt(:,j,:) / nsteps
!                 bc%dqv_dt(:,j,:) = bc%dqv_dt(:,j,:) / nsteps
!                 bc%dqc_dt(:,j,:) = bc%dqc_dt(:,j,:) / nsteps
!             end do
!         endif
            
    end subroutine apply_dt
    
    
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
        type(domain_type),intent(inout)::domain
        type(bc_type),intent(inout)::bc
        type(options_type),intent(in)::options
        real*8,intent(inout)::model_time,next_output
        real*8::end_time
        integer::i,ntimesteps,tenp
        real::dt,dtnext
        
!       compute internal timestep dt to maintain stability
!       courant condition for 3D advection. Note that w is normalized by dx/dz
        dt = domain%dx/(maxval(abs(domain%u)) &
                      + maxval(abs(domain%v)) &
                      + maxval(abs(domain%w)))
!       pick the minimum dt from the begining or the end of the current timestep
        dtnext = domain%dx/(maxval(abs(bc%next_domain%u)) &
                          + maxval(abs(bc%next_domain%v)) &
                          + maxval(abs(bc%next_domain%w)))
        
        dt=min(dt,dtnext)
!       set an upper bound on dt to keep microphysics and convection stable (?) not sure what time is required here. 
        dt=min(dt,120.0) !better min=180?
!       if we have too small a time step just throw an error
        if (dt<1e-1) then
            write(*,*) "dt=",dt
            write(*,*) "Umax",maxval(abs(domain%u)),maxval(abs(bc%next_domain%u))
            write(*,*) "Vmax",maxval(abs(domain%v)),maxval(abs(bc%next_domain%v))
            write(*,*) "Wmax",maxval(abs(domain%w)),maxval(abs(bc%next_domain%w))
            call write_domain(domain,options,99998)
            call write_domain(bc%next_domain,options,99999,'timestep_error_file.nc')
            stop "ERROR time step too small"
        endif
        
!       make dt an integer fraction of the full timestep
        dt=options%in_dt/ceiling(options%in_dt/dt)
!       calculate the number of timesteps
        ntimesteps=nint(options%in_dt/dt)
        end_time=model_time+options%in_dt
        
!       adjust the boundary condition dXdt values for the number of time steps
        call apply_dt(bc,ntimesteps,options)
        write(*,*) "    dt=",dt, "nsteps=",ntimesteps
        
!       ensure internal model consistency (should only need to be called here when the model starts...)
!       e.g. for potential temperature and temperature
        call diagnostic_update(domain,options)
        
!       now just loop over internal timesteps computing all physics in order (operator splitting...)
        do i=1,ntimesteps
            model_time=model_time+dt
            if (dt>1e-5) then
                call advect(domain,options,dt)
                call mp(domain,options,dt)
                call rad(domain,options,model_time/86400.0+50000, dt)
                call lsm(domain,options,dt,model_time)
                call pbl(domain,options,dt)
                call convect(domain,options,dt)

    !           apply/update boundary conditions including internal wind and pressure changes. 
                call forcing_update(domain,bc,options)
            
    !           step model time forward
                domain%model_time=model_time
                if ((abs(model_time-next_output)<1e-1).or.(model_time>next_output)) then
                    call write_domain(domain,options,nint((model_time-options%time_zero)/options%out_dt))
                    next_output=next_output+options%out_dt
                endif
    !           in case out_dt and in_dt arent even multiples of each other.  Make sure we don't over step
                if ((model_time+dt)>end_time) then
                    dt=end_time-model_time
                endif
            endif
        enddo
        
    end subroutine step
end module time_step
