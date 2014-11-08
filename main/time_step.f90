module time_step
	use data_structures     ! *_type  types
	use microphysics,                only : mp
	use convection,                  only : convect
	use land_surface,                only : lsm
	use wind,                        only : balance_uvw
	use advection,                   only : advect
	use output,                      only : write_domain
	use planetary_boundary_layer,    only : pbl
	use radiation, 					 only : rad
	implicit none
	private
	public :: step
	
contains
	
! 	update just the edges of curdata by adding dXdt
	subroutine boundary_update(curdata,dXdt)
		implicit none
		real,dimension(:,:,:), intent(inout) :: curdata
		real,dimension(:,:,:), intent(in) :: dXdt
		integer::nx,nz,ny,i
		
		nx=size(curdata,1)
		nz=size(curdata,2)
		ny=size(curdata,3)

		do i=1,nz
			curdata(1,i,:) =curdata(1,i,:) +dXdt(i,1:ny,1)
			curdata(nx,i,:)=curdata(nx,i,:)+dXdt(i,1:ny,2)
			curdata(:,i,1) =curdata(:,i,1) +dXdt(i,1:nx,3)
			curdata(:,i,ny)=curdata(:,i,ny)+dXdt(i,1:nx,4)
		enddo
		! correct possible rounding errors, primarily an issue of clouds...
		where(curdata(1,:,:)<0) curdata(1,:,:)=0
		where(curdata(nx,:,:)<0) curdata(nx,:,:)=0
		where(curdata(:,:,1)<0) curdata(:,:,1)=0
		where(curdata(:,:,ny)<0) curdata(:,:,ny)=0
		
	end subroutine boundary_update
	
! 	updated all fields in domain using the respective dXdt variables
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
	! 		domain%w(:,:,j)=domain%w(:,:,j)+bc%dw_dt(:,:,j) ! this is now recalculated every time step (as are ur,vr,wr)
			domain%p(:,:,j)=domain%p(:,:,j)+bc%dp_dt(:,:,j)
			domain%pii(:,:,j)=(domain%p(:,:,j)/100000.0)**(R/cp)
	        domain%rho(:,:,j)=domain%p(:,:,j)/(R*domain%th(:,:,j)*domain%pii(:,:,j)) ! kg/m^3
		enddo
        !$omp end do
        !$omp end parallel
		! v has one more y element than others
		ny=ny+1
		domain%v(:,:,ny)=domain%v(:,:,ny)+bc%dv_dt(:,:,ny)
		! dXdt for qv,qc,th are only applied to the boundarys
		call boundary_update(domain%th,bc%dth_dt)
		call boundary_update(domain%qv,bc%dqv_dt)
		call boundary_update(domain%cloud,bc%dqc_dt)
		
		! because density changes with each time step, u/v/w have to be rebalanced as well. 
		! could avoid this by assuming density doesn't change... but would need to keep an "old" density around
		call balance_uvw(domain,options)
	end subroutine forcing_update		


! 	Divides dXdt variables by n timesteps so that after adding it N times we will be at the
! 	correct final value. 
	subroutine apply_dt(bc,nsteps)
		implicit none
		type(bc_type), intent(inout) :: bc
		integer,intent(in)::nsteps
		
		bc%du_dt  =bc%du_dt/nsteps
		bc%dv_dt  =bc%dv_dt/nsteps
! 		bc%dw_dt  =bc%dw_dt/nsteps
		bc%dp_dt  =bc%dp_dt/nsteps
! 		bc%drho_dt=bc%drho_dt/nsteps
		bc%dth_dt =bc%dth_dt/nsteps
		bc%dqv_dt =bc%dqv_dt/nsteps
		bc%dqc_dt =bc%dqc_dt/nsteps
	end subroutine apply_dt
	
	
!	Step forward one IO time step. 
	subroutine step(domain,options,bc,model_time,next_output)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		real*8,intent(inout)::model_time,next_output
		integer::i,ntimesteps,tenp
		real::dt,dtnext,end_time
		
! 		compute internal timestep dt to maintain stability
! 		courant condition for 3D advection. Note that w is normalized by dx/dz
		dt=(domain%dx/max(max(maxval(abs(domain%u)),maxval(abs(domain%v))),maxval(abs(domain%w)))/(3.0))
! 		pick the minimum dt from the begining or the end of the current timestep
		dtnext=(domain%dx/max(max(maxval(abs(bc%next_domain%u)), &
								  maxval(abs(bc%next_domain%v))), &
								  maxval(abs(bc%next_domain%w)))/(3.0))
		dt=min(dt,dtnext)
! 		set an upper bound on dt to keep microphysics and convection stable (?) not sure what time is required here. 
		dt=min(dt,120.0) !better min=180?
! 		if we have too small a time step just throw an error
		if (dt<1e-1) then
			write(*,*) "dt=",dt
			write(*,*) "Umax",maxval(abs(domain%u)),maxval(abs(bc%next_domain%u))
			write(*,*) "Vmax",maxval(abs(domain%v)),maxval(abs(bc%next_domain%v))
			write(*,*) "Wmax",maxval(abs(domain%w)),maxval(abs(bc%next_domain%w))
			call write_domain(domain,options,99998)
			call write_domain(bc%next_domain,options,99999)
			stop "ERROR time step too small"
		endif
		
! 		make dt an integer fraction of the full timestep
		dt=options%out_dt/ceiling(options%out_dt/dt)
! 		calculate the number of timesteps
		ntimesteps=nint(options%in_dt/dt)
		end_time=model_time+options%in_dt
		
! 		adjust the boundary condition dXdt values for the number of time steps
		call apply_dt(bc,ntimesteps)
		write(*,*) "    dt=",dt, "nsteps=",ntimesteps
! 		now just loop over internal timesteps computing all physics in order (operator splitting...)
		do i=1,ntimesteps
			call lsm(domain,options,dt,model_time)
			call pbl(domain,options,dt)
			call advect(domain,options,dt)
			call mp(domain,options,dt)
			call convect(domain,options,dt)
			call rad(domain,options,model_time/86400.0+50000, dt)

! 			apply/update boundary conditions including internal wind and pressure changes. 
			if (.not.options%ideal) then
				call forcing_update(domain,bc,options)
			endif
			
! 			step model time forward
			model_time=model_time+dt
			if ((abs(model_time-next_output)<1e-1).or.(model_time>next_output)) then
				call write_domain(domain,options,nint((model_time-options%time_zero)/options%out_dt))
				next_output=next_output+options%out_dt
			endif
! 			in case out_dt and in_dt arent even multiples of each other.  Make sure we don't over step
			if ((model_time+dt)>end_time) then
				dt=end_time-model_time
			endif
		enddo
		
	end subroutine step
end module time_step
