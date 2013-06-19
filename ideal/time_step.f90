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
		integer::i,ntimesteps,nx,ny,nz
		real::dt,dtnext
		real,dimension(:,:,:),allocatable::rho,pii
	    real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
	    real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
		
		! courant condition for 3D advection... could make 3 x 1D to maximize dt? esp. w/linear wind speedups...
		dt=(options%dx/max(max(maxval(abs(domain%u)),maxval(abs(domain%v))),maxval(abs(domain%w)))/3.0)
! 		pick the minimum dt from the begining or the end of the current timestep
		dtnext=(options%dx/max(max(maxval(abs(bc%next_domain%u)), &
										maxval(abs(bc%next_domain%v))), &
										maxval(abs(bc%next_domain%w)))/3.0)
		write(*,*) maxval(bc%next_domain%w), maxval(bc%next_domain%v), maxval(bc%next_domain%u)
		write(*,*) minval(bc%next_domain%w), minval(bc%next_domain%v), minval(bc%next_domain%u)
		write(*,*) maxval(domain%w), maxval(domain%v), maxval(domain%u)
		dt=min(dt,dtnext)
! 		make dt an integer fraction of the full timestep
		dt=min(dt,60.0)
		dt=options%io_dt/ceiling(options%io_dt/dt)
! 		calcualte the number of timesteps
		ntimesteps=options%io_dt/dt
		
		call apply_dt(bc,ntimesteps)
		nx=size(domain%p,1)
		nz=size(domain%p,2)
		ny=size(domain%p,3)
! 		allocate(rho(nx,nz,ny),pii(nx,nz,ny))
		write(*,*) dt,ntimesteps
		do i=1,ntimesteps
! 			pii=(100000.0/domain%p)**(R/cp)
! 			rho=0.622*domain%p/(R*(domain%th/pii)*(domain%qv+0.622))
! 			write(*,*) minval(rho),maxval(rho), minval(domain%th), maxval(domain%th), minval(domain%cloud), maxval(domain%cloud), minval(domain%qv), maxval(domain%qv)
! 			write(*,*) "pre-adv", i,ntimesteps
			call advect(domain,options,dt,options%dx)
! 			pii=(100000.0/domain%p)**(R/cp)
! 			rho=0.622*domain%p/(R*(domain%th/pii)*(domain%qv+0.622))
! 			write(*,*) minval(rho),maxval(rho), minval(domain%th), maxval(domain%th), minval(domain%cloud), maxval(domain%cloud), minval(domain%qv), maxval(domain%qv)
! 			write(*,*) "pre-MP", i,ntimesteps
			call mp(domain,options,dt)
	! 		call lsm(domain,options,dt)
	! 		call pbl(domain,options,dt)
	! 		call radiation(domain,options,dt)
			write(*,*) i,ntimesteps
			
			call forcing_update(domain,bc)
		enddo
! 		deallocate(rho,pii)
	end subroutine step
end module time_step