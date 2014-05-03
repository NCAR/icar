module wind
	use linear_theory_winds
	use data_structures
	implicit none
	private
	public::update_winds,balance_uvw,init_winds
	real, parameter::deg2rad=0.017453293 !2*pi/360
contains

! Forces u,v, and w fields to balance
!       du/dx+dv/dy = dw/dz
! Starts by setting w out of the ground=0 then works through layers
	subroutine balance_uvw(domain,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		real, allocatable,dimension(:,:) ::du,dv,divergence,rhou,rhov,rhow
		integer ::nx,ny,nz,i
		nx=size(domain%w,1)
		nz=size(domain%w,2)
		ny=size(domain%w,3)
		
		! these could me module level to prevent lots of allocation/deallocation/reallocations
		if (options%advect_density) then
			allocate(rhou(nx-1,ny-2))
			allocate(rhov(nx-2,ny-1))
			allocate(rhow(nx-2,ny-2))
		endif
		allocate(du(nx-2,ny-2))
		allocate(dv(nx-2,ny-2))
		allocate(divergence(nx-2,ny-2))
		
! If this becomes a bottle neck in the code it could be parallelized over y  
! 		loop over domain levels
		do i=1,nz
			! if we are incorporating density into the advection equation (as we should be default now)
			if (options%advect_density) then
				! first calculate density on the u and v grid points
				rhou=(domain%rho(1:nx-1,i,2:ny-1) + domain%rho(2:nx,i,2:ny-1))/2
				rhov=(domain%rho(2:nx-1,i,1:ny-1) + domain%rho(2:nx-1,i,2:ny))/2
				! for most grid points interpolate between the current grid and the one above it
				if (i<nz) then
					rhow=(domain%rho(2:nx-1,i,2:ny-1) + domain%rho(2:nx-1,i+1,2:ny-1))/2
				else
					! for the top grid cell extrapolate upwards based on the current grid and the one below it
					rhow=2*domain%rho(2:nx-1,i,2:ny-1)-domain%rho(2:nx-1,i-1,2:ny-1)
				endif
				! calculate horizontal divergence based on the wind * density field
				domain%vr(2:nx-1,i,1:ny-1)=rhov*domain%v(2:nx-1,i,1:ny-1) * domain%dz(1,i,1)*domain%dx
				dv=domain%vr(2:nx-1,i,2:ny-1) - domain%vr(2:nx-1,i,1:ny-2)
				
				domain%ur(1:nx-1,i,2:ny-1)=rhou*domain%u(1:nx-1,i,2:ny-1) * domain%dz(1,i,1)*domain%dx
				du=domain%ur(2:nx-1,i,2:ny-1) - domain%ur(1:nx-2,i,2:ny-1)
				
				divergence=du+dv
				if (i==1) then
					! if this is the first model level start from 0 at the ground
					domain%wr(2:nx-1,i,2:ny-1)=0-divergence
					domain%w(2:nx-1,i,2:ny-1)=0-divergence/rhow/(domain%dx**2)
				else
					! else calculate w as a change from w at the level below
					domain%wr(2:nx-1,i,2:ny-1)=domain%wr(2:nx-1,i-1,2:ny-1)-divergence
					domain%w(2:nx-1,i,2:ny-1)=domain%w(2:nx-1,i-1,2:ny-1)-divergence/rhow/(domain%dx**2)
				endif
			else
				! calculate horizontal divergence
				dv=domain%v(2:nx-1,i,2:ny-1) - domain%v(2:nx-1,i,1:ny-2)
				du=domain%u(2:nx-1,i,2:ny-1) - domain%u(1:nx-2,i,2:ny-1)
				divergence=du+dv
				if (i==1) then
					! if this is the first model level start from 0 at the ground
					domain%w(2:nx-1,i,2:ny-1)=0-divergence
				else
					! else calculate w as a change from w at the level below
					domain%w(2:nx-1,i,2:ny-1)=domain%w(2:nx-1,i-1,2:ny-1)-divergence
				endif
			endif
		enddo
		!NOTE w is scaled by dx/dz because it is balancing divergence over dx with a flow through dz
! 		could rescale below but for now this is left in the advection code (the only place it matters for now)
!       this makes it easier to play with varying options for including varying dz and rho in advection. 
! 		domain%w(:,:nz-1,:)=domain%w(:,:nz-1,:)/domain%dx * (domain%dz(:,1:nz-1,:)+domain%dz(:,2:nz,:))/2
! 		domain%w(:,nz,:)=domain%w(:,nz,:)/domain%dx * domain%dz(:,nz,:)
		
		deallocate(du,dv,divergence)
		if (options%advect_density) then
			deallocate(rhou,rhov,rhow)
		endif
	end subroutine balance_uvw
	
	subroutine make_winds_grid_relative(domain)
		type(domain_type), intent(inout) :: domain
		real,dimension(:),allocatable :: u,v
		integer::nx,nz,ny,k,j
		
		nx=size(domain%p,1)
		nz=size(domain%p,2)
		ny=size(domain%p,3)
		
		allocate(u(nx))
		allocate(v(nx))
		!assumes u and v come in on a staggered grid with one additional grid point in x/y for u/v respectively
		domain%u(:nx,:,:) = (domain%u(:nx,:,:)+domain%u(2:,:,:))/2
		domain%v(:,:,:ny) = (domain%v(:,:,:ny)+domain%v(:,:,2:))/2
		do j=1,ny
			do k=1,nz
				u=domain%u(:nx,k,j)*domain%costheta(:nx,j) + domain%v(:nx,k,j)*domain%sintheta(:nx,j)
				v=domain%v(:nx,k,j)*domain%costheta(:nx,j) + domain%u(:nx,k,j)*domain%sintheta(:nx,j)
				domain%u(:nx,k,j)=u
				domain%v(:nx,k,j)=v
			enddo
		enddo
		deallocate(u,v)
		domain%u(1:nx-1,:,:) = (domain%u(:nx-1,:,:)+domain%u(2:nx,:,:))/2
		domain%v(:,:,1:ny-1) = (domain%v(:,:,:ny-1)+domain%v(:,:,2:ny))/2
		
	end subroutine make_winds_grid_relative

! 	rotate winds from real space back to terrain following grid (approximately)
!   assumes a simple slope transform in u and v independantly
	subroutine rotate_wind_field(domain)
        implicit none
        class(linearizable_type),intent(inout)::domain
		integer :: nx,ny,nz,i,j
	
		nx=size(domain%u,1)
		nz=size(domain%u,2)
		ny=size(domain%u,3)
		do j=1,ny
			do i=1,nz
				domain%u(1:nx-2,i,j)=domain%u(1:nx-2,i,j)*domain%dzdx(:,j)
				if (j<ny) then
					domain%v(:,i,j)=domain%v(:,i,j)*domain%dzdy(:,j)
				endif
			end do
		end do
	
	end subroutine rotate_wind_field

	! apply wind field physics and adjustments
	! this will call the linear wind module if necessary, otherwise it just updates for 
	subroutine update_winds(domain,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		real,allocatable,dimension(:,:,:)::temparray
		integer::nx,ny,nz,i,j
		
		!note, this should only be called once per time step, DO NOT use update winds internal to the time step
		print*, "rotating winds into the model grid"
		call make_winds_grid_relative(domain)
		
! 		linear winds
		if (options%physics%windtype==1) then
			call linear_perturb(domain)
		endif
! 		else assumes even flow over the mountains
		call rotate_wind_field(domain)
		call balance_uvw(domain,options)
		
	end subroutine update_winds
	
	!setup initial fields (i.e. grid relative rotation fields)
	subroutine init_winds(domain,options)
		type(domain_type), intent(inout) :: domain
		type(options_type), intent(in) :: options
		integer:: i,j,nx,ny,starti,endi
		real::dist,dlat,dlon
		
		call allocate_winds(domain)
		
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
		do j=1,ny
			do i=1,nx
				! in case we are in the first or last grid, reset boundaries
				if (i==1) then
					starti=i
					endi=i+1
				elseif (i==nx) then
					starti=i-1
					endi=i
				else
					starti=i-1
					endi=i+1
				endif
				
				! change in latitute
				dlat=domain%lat(endi,j)-domain%lat(starti,j)
				! change in longitude
				dlon=(domain%lon(endi,j)-domain%lon(starti,j))*cos(deg2rad*domain%lat(i,j))
				! distance between two points
				dist=sqrt(dlat**2+dlon**2)
				! sin/cos of angles for use in rotating fields later
				domain%sintheta(i,j)=(-1)*dlat/dist
				domain%costheta(i,j)=dlon/dist
			enddo
		enddo
		
! 		dzdx/y used in rotating windfield back to terrain following grid in a simple fashion
		allocate(domain%dzdx(nx-1,ny))
		allocate(domain%dzdy(nx,ny-1))
		domain%dzdx=sqrt((domain%terrain(2:nx,:)-domain%terrain(1:nx-1,:))**2+domain%dx**2)/domain%dx
		domain%dzdy=sqrt((domain%terrain(:,2:ny)-domain%terrain(:,1:ny-1))**2+domain%dx**2)/domain%dx
		
	end subroutine init_winds
	
! allocate memory
	subroutine allocate_winds(domain)
		type(domain_type), intent(inout) :: domain
		integer::nx,ny
		
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
		
		if (.not.allocated(domain%sintheta)) then
			allocate(domain%sintheta(nx,ny))
		endif
		if (.not.allocated(domain%costheta)) then
			allocate(domain%costheta(nx,ny))
		endif
	end subroutine allocate_winds

! deallocate memory
	subroutine finalize_winds(domain)
		type(domain_type), intent(inout) :: domain
		
		if (allocated(domain%sintheta)) then
			deallocate(domain%sintheta)
		endif
		if (allocated(domain%costheta)) then
			deallocate(domain%costheta)
		endif
	
	end subroutine finalize_winds
end module wind