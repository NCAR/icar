module wind
	use linear_theory_winds
	use data_structures
	implicit none
	private
	public::update_winds
contains

	subroutine balance_uvw(domain,options)
! Forces u,v, and w fields to balance
!       du/dx+dv/dy = dw/dz
! Starts by setting w out of the ground=0 then works through layers
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		real, allocatable,dimension(:,:) ::du,dv,divergence,rhou,rhov,rhow
		integer ::nx,ny,nz,i
		nx=size(domain%w,1)
		nz=size(domain%w,2)
		ny=size(domain%w,3)
		
		if (options%advect_density) then
			allocate(rhou(nx-1,ny-2))
			allocate(rhov(nx-2,ny-1))
			allocate(rhow(nx-2,ny-2))
		endif
		allocate(du(nx-2,ny-2))
		allocate(dv(nx-2,ny-2))
		allocate(divergence(nx-2,ny-2))
		
! 		loop over domain levels
		do i=1,nz
			if (options%advect_density) then
				rhou=(domain%rho(1:nx-1,i,2:ny-1) + domain%rho(2:nx,i,2:ny-1))/2
				rhov=(domain%rho(2:nx-1,i,1:ny-1) + domain%rho(2:nx-1,i,2:ny))/2
				if (i<nz) then
					rhow=(domain%rho(2:nx-1,i,2:ny-1) + domain%rho(2:nx-1,i+1,2:ny-1))/2
				else
					rhow=domain%rho(2:nx-1,i,2:ny-1)
				endif
! 			calculate horizontal divergence
				dv=rhov(:,2:ny-1)*domain%v(2:nx-1,i,2:ny-1) - rhov(:,1:ny-2)*domain%v(2:nx-1,i,1:ny-2)
				du=rhou(2:nx-1,:)*domain%u(2:nx-1,i,2:ny-1) - rhou(1:nx-2,:)*domain%u(1:nx-2,i,2:ny-1)
				divergence=du+dv
				if (i==1) then
					! if this is the first model level start from 0 at the ground
					domain%w(2:nx-1,i,2:ny-1)=0-divergence/rhow
				else
					! else calculate w as a change from w at the level below
					domain%w(2:nx-1,i,2:ny-1)=domain%w(2:nx-1,i-1,2:ny-1)-divergence/rhow
				endif
			else
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
		
		deallocate(du,dv,divergence,rhou,rhov,rhow)
	end subroutine balance_uvw
	
! 	apply wind field physics and adjustments
	subroutine update_winds(domain,options)
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		real,allocatable,dimension(:,:,:)::temparray
		integer::nx,ny,nz,i,j
		
! 		linear winds
		if (options%physics%windtype==1) then
			call linear_perturb(domain)
			if (options%ideal) then
				domain%v=0
			endif
! 		assumes even flow over the mountains
		else
			nx=size(domain%u,1)
			nz=size(domain%u,2)
			ny=size(domain%u,3)
! I'm not sure the temparray is necessary. I was getting a segfault (sometimes) before changing this... 
!  
! 	U and V are now staggered on input so now we just offset by 1
! 
			allocate(temparray(nx,nz,ny))
			temparray=domain%u
			domain%u(1:nx-1,:,:)=temparray(2:nx,:,:)
			deallocate(temparray)
			
			ny=ny+1
			allocate(temparray(nx-1,nz,ny))
			temparray=domain%v
			domain%v(:,:,1:ny-1)=temparray(:,:,2:ny)
			deallocate(temparray)
		endif
		
		call balance_uvw(domain,options)
		
	end subroutine update_winds
end module wind