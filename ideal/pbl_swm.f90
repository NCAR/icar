! simple PBL diffusion package for the Simple Weather Model
!  Local-K diffusion type PBL as in Louis (1979) as documented in Hong and Pan (1996) = HP96
!  Hong and Pan used this for their free atmosphere diffusion, but noted differences
!  used in the "current operational model" notably the asymptotic length scale lambda 
!
! HP96 = Hong,S.-Y. and H.-L. Pan (1996) Monthly Weather Review v127 p2322
! 		 Nonlocal Boundary Layer Vertical Diffusion in a Medium Range Forecast Model
! 
! Implemented with K,shear,stability... on half levels
! 	rho on half levels for f=k*rho*dq/dz*dt
! 	rho on full levels for q=q+f/rho
!   q,U,V on full levels
module pbl_simple
	use data_structures
	private
	public :: simple_pbl, finalize_simple_pbl, init_simple_pbl
	
! 	note _m indicates module level variables
! 	these variables are declared as module level variables so that they do not need to be allocated
!   deallocated, and re-allocated all the time. 
	real, allocatable, dimension(:,:,:), target ::virt_pot_temp_zgradient_m
	real, allocatable, dimension(:,:,:) :: rig_m,shear_m
	real, allocatable, dimension(:,:,:) :: stability_m,prandtl_m,l_m,K_m,Kq_m
	integer :: nx,nz,ny
	
!	limits on Pr noted in HP96 page 2325 below eqn 13
	real, parameter :: pr_upper_limit = 4.0 !Prandtl number for stability
	real, parameter :: pr_lower_limit = 0.25 !Prandtl number for stability
	real, parameter :: asymp_length_scale = 1/250.0 !m from Hong and Pan (1996)
	real, parameter :: kappa =0.4 !von Karman constant
	! note, they actually use 30m because they only use this for free-atmosphere mixing
	! but they note that 250m is used in the operational model for the full PBL mixing
	
	
contains
	subroutine simple_pbl(domain,dt)
		type(domain_type), intent(inout) :: domain
		real, intent(in)::dt
		integer :: j,k
		
!       OpenMP parallelization		
!c      omp parallel shared(domain,l_m,k_m,stability_m,shear_m) &
!c      omp firstprivate(nz,ny,kappa,asymp_length_scale) private(k,j)
!c      omp do
		do j=1,ny
			call calc_shear(domain,j)
			call calc_virt_pot_temp_zgradient(domain,j)
			call calc_richardson_gradient(domain,j)
			call calc_pbl_stability_function(j)
			
			do k=1,nz
				l_m(:,k,j) = 1 / (1/(kappa*(domain%z(:,k,j)-domain%terrain(:,j))) + asymp_length_scale)
			enddo
! 			diffusion for momentum... can I ignore this term? 
!           k = l**2 * stability * shear * dt/dz
			K_m(:,:,j) = l_m(:,:,j)**2 * stability_m(:,:,j) * shear_m(:,:,j) * &
						 dt/((domain%dz(:,2:,j)+domain%dz(:,:nz-1,j))/2)
! 			diffusion for scalars
			Kq_m(:,:,j)=K_m(:,:,j)/prandtl_m(:,:,j)
			
			call pbl_diffusion(domain,j)
		enddo
!      omp end do
!      omp end parallel
	end subroutine simple_pbl
	
	subroutine pbl_diffusion(domain,j)
		type(domain_type), intent(inout) :: domain
		integer,intent(in) :: j
		integer::i
		real,dimension(nx,nz-1):: fluxes,rhomean,dz
		
		! If this works, dz could be moved into the domain
		! and calculated once at the begining of the simulation
		rhomean=(domain%rho(:,1:nz-1,j)+domain%rho(:,2:,j))/2
		dz=(domain%dz(:,1:nz-1,j)+domain%dz(:,2:,j))/2
! 		qmeans=(domain%qv(:,1:nz-1,j)+domain%qv(:,2:nz,j))/2
		
		! note Kq_m already has dt/dz embedded in it
		! diffusion fluxes within the PBL
		! q = q + (k dq/qz)/dz *dt
		
		! First water vapor
		fluxes=Kq_m(:,:,j)*rhomean*(domain%qv(:,:nz-1,j)-domain%qv(:,2:,j))
		! first layer assumes no flow through the surface, that comes from the LSM
		domain%qv(:,1,j) = domain%qv(:,1,j) + &
								fluxes(:,1) / (domain%dz(:,1,j)*domain%rho(:,1,j))
		! middle layers (no change for top layer assuming flux in = flux out)
		domain%qv(:,2:nz-1,j) = domain%qv(:,2:nz-1,j) + &
		   					   (fluxes(:,:nz-2)+fluxes(:,2:nz-1)) / (domain%dz(:,2:nz-1,j)*domain%rho(:,2:nz-1,j))

		! ditto for potential temperature
   		fluxes=Kq_m(:,:,j)*rhomean*(domain%th(:,:nz-1,j)-domain%th(:,2:,j))
   		domain%th(:,1,j) = domain%th(:,1,j) + &
   								fluxes(:,1) / (domain%dz(:,1,j)*domain%rho(:,1,j))
   		domain%th(:,2:nz-1,j) = domain%th(:,2:nz-1,j) + &
   		   					   (fluxes(:,:nz-2)+fluxes(:,2:nz-1)) / (domain%dz(:,2:nz-1,j)*domain%rho(:,2:nz-1,j))
	
		! and cloud water
		fluxes=Kq_m(:,:,j)*rhomean*(domain%cloud(:,:nz-1,j)-domain%cloud(:,2:,j))
		domain%cloud(:,1,j) = domain%cloud(:,1,j) + &
								fluxes(:,1) / (domain%dz(:,1,j)*domain%rho(:,1,j))
		domain%cloud(:,2:nz-1,j) = domain%cloud(:,2:nz-1,j) + &
		   					   (fluxes(:,:nz-2)+fluxes(:,2:nz-1)) / (domain%dz(:,2:nz-1,j)*domain%rho(:,2:nz-1,j))
        !
		! and cloud ice
		fluxes=Kq_m(:,:,j)*rhomean*(domain%ice(:,:nz-1,j)-domain%ice(:,2:,j))
		domain%ice(:,1,j) = domain%ice(:,1,j) + &
								fluxes(:,1) / (domain%dz(:,1,j)*domain%rho(:,1,j))
		domain%ice(:,2:nz-1,j) = domain%ice(:,2:nz-1,j) + &
		   					   (fluxes(:,:nz-2)+fluxes(:,2:nz-1)) / (domain%dz(:,2:nz-1,j)*domain%rho(:,2:nz-1,j))
		!
		! and snow
		fluxes=Kq_m(:,:,j)*rhomean*(domain%qsnow(:,:nz-1,j)-domain%qsnow(:,2:,j))
		domain%qsnow(:,1,j) = domain%qsnow(:,1,j) + &
							fluxes(:,1) / (domain%dz(:,1,j)*domain%rho(:,1,j))
		domain%qsnow(:,2:nz-1,j) = domain%qsnow(:,2:nz-1,j) + &
						   (fluxes(:,:nz-2)+fluxes(:,2:nz-1)) / (domain%dz(:,2:nz-1,j)*domain%rho(:,2:nz-1,j))
		!

	end subroutine pbl_diffusion
	
	subroutine calc_shear(domain,j)
		type(domain_type), intent(in) :: domain
		integer,intent(in) :: j
		real,dimension(nx,nz)::centered_winds
		integer::last_wind,k
		
		centered_winds(:,:)=sqrt( ((domain%u(1:nx-1,:,j)+domain%u(2:,:,j))/2)**2 &
			   				 +((domain%v(:,:,j)     +domain%v(:,:,j+1))/2)**2)
		
		shear_m(:,:,j)=abs(centered_winds(:,2:)-centered_winds(:,:nz-1))  &
		              /((domain%dz(:,:nz-1,j)+domain%dz(:,2:,j))*0.5)
	end subroutine calc_shear

! 	calculate the vertical gradient in virtual potential temperature
	subroutine calc_virt_pot_temp_zgradient(domain,j)
		type(domain_type), intent(in) :: domain
		integer,intent(in)::j
		real,dimension(size(domain%p,1))::vth_row,last_vth_row
		real,pointer::vth(:,:,:)
		integer::k
		
! 		set up vth as a pointer to the module variable to make reading this subroutine a little easier
		vth=>virt_pot_temp_zgradient_m
		
! 		first calculate the virtual potential temperature
! 		vth=th*(1+0.61*qv-(qc+qi+qr+qs))
		vth(:,:,j)=domain%th(:,:,j)* &
				(1+0.61*domain%qv(:,:,j) &
				 -(domain%cloud(:,:,j)+domain%ice(:,:,j)+domain%qrain(:,:,j)+domain%qsnow(:,:,j)))
		do k=1,nz-1
			vth(:,k,j)=(vth(:,k+1,j)-vth(:,k,j))/((domain%dz(:,k,j)+domain%dz(:,k+1,j))*0.5)
! 			
! 			OLDER code that will calculate virtual potential temperature gradient on full levels
! 			! for the first model level save this level for the next one and calculate the gradient
! 			! based only on the bottom two levels
! 			if (k==1) then
! 				last_vth_row=vth(:,k,j)
! 				vth(:,k,j)=(vth(:,k+1,j)-vth(:,k,j))/domain%dz(:,k,j)
! 			elseif (k==nz) then
! 				! for the top model level calculate the gradient based only on the top two levels
! 				vth(:,k,j)=(vth(:,k,j)-last_vth_row)/domain%dz(:,k,j)
! 			else
! 				! for all other levels, first save this level temporarily
! 				vth_row=vth(:,k,j)
! 				! then calculate the gradient across this level
! 				vth(:,k,j)=(vth(:,k+1,j)-last_vth_row)/ &
! 					((domain%dz(:,k-1,j)+domain%dz(:,k+1,j))*0.5+domain%dz(:,k,j))
! 				! then store the saved level in the "last" level for the next iteration
! 				last_vth_row=vth_row
! 			endif
		enddo
	end subroutine calc_virt_pot_temp_zgradient

	subroutine calc_pbl_stability_function(j)
		integer,intent(in)::j
		integer::i,k
! 		HP96 eqn 13
		stability_m(:,:,j) = exp(-8.5*rig_m(:,:,j)) + 0.15 / (rig_m(:,:,j)+3)
! 		HP96 eqn 13 continued
		prandtl_m(:,:,j) = 1.5+3.08*rig_m(:,:,j)
! 		limits on Pr noted in Hong and Pan (?)
		do k=1,size(prandtl_m,2)
			do i=1,size(prandtl_m,1)
				if (prandtl_m(i,k,j)>pr_upper_limit) then
					prandtl_m(i,k,j)=pr_upper_limit
				elseif (prandtl_m(i,k,j)<pr_lower_limit) then
					prandtl_m(i,k,j)=pr_lower_limit
				endif
			enddo
		enddo
! 		alternatively... which is faster? 
! 		where(prandtl_m(:,:,j)>4) prandtl_m(:,:,j)=4
! 		where(prandtl_m(:,:,j)<0.25) prandtl_m(:,:,j)=0.25
	end subroutine calc_pbl_stability_function
	
	subroutine calc_richardson_gradient(domain,j)
! 		calculate the local gradient richardson number as in eqn. between 11 and 12 in HP96
		type(domain_type), intent(inout) :: domain
		integer,intent(in)::j
		real,dimension(nx,nz)::temperature
		temperature=domain%th(:,:,j)*domain%pii(:,:,j)
! 		might be slightly better to interpolate theta to half levels, then recalc p and pii at half levels
		temperature(:,:nz-1)= (temperature(:,:nz-1)+temperature(:,2:))*0.5
		rig_m(:,:,j) =  g/temperature(:,:nz-1)  & 
					   * virt_pot_temp_zgradient_m(:,:nz-1,j) * 1/(shear_m(:,:,j)**2)
	end subroutine calc_richardson_gradient
	

! memory allocation and deallocation	
! Can/Should also add parameter definition from options%pbl (or something like that)
	subroutine init_simple_pbl(domain,options)
		type(domain_type), intent(in) :: domain
		type(options_type),intent(in) :: options
		nx=size(domain%p,1)
		nz=size(domain%p,2)
		ny=size(domain%p,3)
		
		allocate(virt_pot_temp_zgradient_m(nx,nz,ny))
		allocate(rig_m(nx,nz-1,ny))
		allocate(stability_m(nx,nz-1,ny))
		allocate(shear_m(nx,nz-1,ny))
		allocate(prandtl_m(nx,nz-1,ny))
		allocate(K_m(nx,nz-1,ny))
		allocate(Kq_m(nx,nz-1,ny))
		allocate(l_m(nx,nz-1,ny))
	end subroutine init_simple_pbl
	
	subroutine finalize_simple_pbl()
		if (allocated(virt_pot_temp_zgradient_m)) then
			deallocate(virt_pot_temp_zgradient_m)
		endif
		if (allocated(rig_m)) then
			deallocate(rig_m)
		endif
		if (allocated(stability_m)) then
			deallocate(stability_m)
		endif
		if (allocated(shear_m)) then
			deallocate(shear_m)
		endif
		if (allocated(prandtl_m)) then
			deallocate(prandtl_m)
		endif
		if (allocated(K_m)) then
			deallocate(K_m)
		endif
		if (allocated(Kq_m)) then
			deallocate(Kq_m)
		endif
		if (allocated(l_m)) then
			deallocate(l_m)
		endif
	end subroutine finalize_simple_pbl
end module pbl_simple