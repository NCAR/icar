module wind
	use linear_theory_winds
	use data_structures
	implicit none
	contains

	subroutine balance_uvw(domain)
! Forces u,v, and w fields to balance
!       du/dx+dv/dy = dw/dz
! Starts by setting w out of the ground=0 then works through layers
		type(domain_type),intent(inout)::domain
		real, allocatable,dimension(:,:) ::du,dv,divergence
		integer ::nx,ny,nz,i
		ny=size(domain%w,1)
		nz=size(domain%w,2)
		nx=size(domain%w,3)
		
		allocate(du(ny-2,nx-2))
		allocate(dv(ny-2,nx-2))
		allocate(divergence(ny-2,nx-2))
		
! 		loop over domain levels
		do i=1,nz
! 			calculate horizontal divergence
			du=domain%u(2:ny-1,i,2:nx-1) - domain%u(2:ny-1,i,1:nx-2)
			dv=domain%v(2:ny-1,i,2:nx-1) - domain%v(1:ny-2,i,2:nx-1)
			divergence=du+dv
			if (i==1) then
! 				if this is the first model level start from 0 at the ground
				domain%w(2:ny-1,i,2:nx-1)=0-divergence
			else
! 				else calculate w as a change from w at the level below
				domain%w(2:ny-1,i,2:nx-1)=domain%w(2:ny-1,i-1,2:nx-1)-divergence
			endif
		enddo
		
		deallocate(du,dv,divergence)
	end subroutine balance_uvw
	
	subroutine update_winds(domain,options)
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		
! 		linear winds
		if (options%physics%windtype==1) then
			call linear_perturb(domain)
		endif
		
		call balance_uvw(domain)
		
	end subroutine update_winds
end module wind