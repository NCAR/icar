module time_step
	use data_structures     ! *_type  types
	use microphysics        ! mp
	use wind                ! update_winds
	use advection           ! advect
	
	implicit none
	private
	public :: step
	
contains
	
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
	end subroutine boundary_update
	
	
	subroutine forcing_update(domain,bc)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		
		domain%u=domain%u+bc%dudt
		domain%v=domain%v+bc%dvdt
		domain%w=domain%w+bc%dwdt
		domain%p=domain%p+bc%dpdt
! 		dXdt for qv,qc,th are only applied to the boundarys
		call boundary_update(domain%th,bc%dthdt)
		call boundary_update(domain%qv,bc%dqvdt)
		call boundary_update(domain%cloud,bc%dqcdt)
	end subroutine forcing_update		


	subroutine apply_dt(bc,nsteps)
		implicit none
		type(bc_type), intent(inout) :: bc
		integer,intent(in)::nsteps
		
		bc%dudt  =bc%dudt/nsteps
		bc%dvdt  =bc%dvdt/nsteps
		bc%dwdt  =bc%dwdt/nsteps
		bc%dpdt  =bc%dpdt/nsteps
		bc%dthdt =bc%dthdt/nsteps
		bc%dqvdt =bc%dqvdt/nsteps
		bc%dqcdt =bc%dqcdt/nsteps
	end subroutine apply_dt
	
	
	subroutine step(domain,options,bc)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer::i,ntimesteps,tenp
		real::dt,dtnext
		
		! courant condition for 3D advection... could make 3 x 1D to maximize dt? esp. w/linear wind speedups...
		dt=(options%dx/max(max(maxval(abs(domain%u)),maxval(abs(domain%v))),maxval(abs(domain%w)))/3.0)
! 		pick the minimum dt from the begining or the end of the current timestep
		dtnext=(options%dx/max(max(maxval(abs(bc%next_domain%u)), &
										maxval(abs(bc%next_domain%v))), &
										maxval(abs(bc%next_domain%w)))/3.0)

		dt=min(dt,dtnext)
! 		make dt an integer fraction of the full timestep
		dt=min(dt,60.0)
		dt=options%io_dt/ceiling(options%io_dt/dt)
! 		calculate the number of timesteps
		ntimesteps=options%io_dt/dt
		
		call apply_dt(bc,ntimesteps)
		write(*,*) dt,ntimesteps
		tenp=10
		do i=1,ntimesteps
			if (i>=(ntimesteps*tenp/100.0)) then
				write(*,*) nint((100.0*i)/ntimesteps), "%"
				tenp=tenp+10
			endif
			call advect(domain,options,dt,options%dx)
			call mp(domain,options,dt)
	! 		call lsm(domain,options,dt)
	! 		call pbl(domain,options,dt)
	! 		call radiation(domain,options,dt)
			
			call forcing_update(domain,bc)
		enddo
	end subroutine step
end module time_step