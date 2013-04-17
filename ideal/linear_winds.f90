module fft
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
end module fft

module linear_theory_winds
    use fft
    use data_structures
    implicit none
contains
	
    real function calc_stability(domain)
        implicit none
        type(domain_type),intent(in)::domain
        calc_stability=6.37e-5
!         real, parameter :: R  = 287.0
!         real, parameter :: Rv = 461.0
!         real, parameter :: cp = 1004.0
!         real, parameter :: L   = 2.5e6
!         real, parameter :: g  = 9.81
!         real, parameter :: ratio = 18.015/28.964
!         real, parameter :: t0 = 273.15
! 
!         p0=domain%.p(:,1,:)
!         pii=1.0/((100000.0/domain%p)**(R/cp))
!         T2m=domain%th(:,1,:)*pii(:,1,:)
!     
!         es = 611.21*exp(17.502*(T2m-t0)/(T2m-32.19))
!         qs0 = ratio * es/(p0-es)
! 
!         cap_gamma = -(g * (1.+(L*qs0)/(R*T2m)) / (cp + (L**2 * qs0*ratio) / (R*T2m**2)))
!         env_gamma = np.mean(np.diff(weather.th*pii,axis=0)/np.diff(base.hgt3d,axis=0),axis=0)
!     
!         dry_gamma=sum(env_gamma-cap_gamma)/size(env_gamma)
!         ndsq=(g/(sum(T2m)/size(T2m)))*(dry_gamma)
!         ndsq=max(min(1e-5,ndsq),1e-8)
!         calc_stability=ndsq
    end function calc_stability
    
    subroutine linear_winds(domain,Ndsq)
!         # see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
!         # -------------------------------------------------------------------------------
!         # Ndsq  = 0.003**2      # dry BV freq sq. was 0.01**2 initially, Idar suggested 0.005**2
!         # 1E-8 keeps it relatively stable, no wild oscillations in u,v
!         # but 1E-8 doesn't damp vertical winds with height very fast
!         # could/should be calculated from the real atm profile with limits
!         f  = 9.37e-5           # rad/s Coriolis frequency for 40deg north
!         # ---------------------------------------------------------------------------------
		use, intrinsic :: iso_c_binding
        implicit none
        type(domain_type),intent(inout)::domain
        real, intent(in)::Ndsq
        complex,parameter :: j= (0,1)
        real,allocatable,dimension(:,:)::k,l,kl,sig
        complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:)::denom,what,w_hat,uhat,u_hat,vhat,v_hat,m,ineta
        real::gain,offset !used in setting up k and l arrays
        integer::nx,ny,nz,i,midpoint
        real::U,V,z
		real,parameter::pi=3.1415927
        type(C_PTR) :: plan
        
        ny=size(domain%p,1)
        nz=size(domain%p,2)
        nx=size(domain%p,3)
        
        allocate(k(ny,nx))
        allocate(l(ny,nx))
        allocate(kl(ny,nx))
        allocate(sig(ny,nx))
        allocate(denom(ny,nx))
        allocate(what(ny,nx))
        allocate(w_hat(ny,nx))
        allocate(uhat(ny,nx))
        allocate(u_hat(ny,nx))
        allocate(vhat(ny,nx))
        allocate(v_hat(ny,nx))
        allocate(m(ny,nx))
        allocate(ineta(ny,nx))
        
!         # % Compute 2D k and l wavenumber fields (could be stored in options or domain or something)
        offset=pi/domain%dx
        gain=2*offset/(nx-1)
		k(1,:) = (/((i*gain-offset),i=0,nx-1)/)
		! foolishly inefficient ifftshift should just fix the array
		! creation above
		if (mod(nx,2)==1) then
			! odd
			k(2,:)=k(1,:)
			nx=nx-1
			do i=1,nx/2
				k(1,i)=k(2,i+nx/2)
				k(1,i+nx/2+1)=k(2,i)
			enddo
			k(1,nx/2+1)=k(2,nx+1)
			nx=nx+1
		else
			! even
			do i=1,nx/2
				gain=k(1,i)
				k(1,i)=k(1,i+nx/2)
				k(1,i+nx/2)=gain
			enddo
		endif
		
        gain=2*offset/(ny-1)
		l(:,1) = (/((i*gain-offset),i=0,ny-1)/)
		! foolishly inefficient ifftshift should just fix the array
		! creation above
		if (mod(ny,2)==1) then
			! odd
			l(:,2)=l(:,1)
			ny=ny-1
			do i=1,ny/2
				l(i,1)=l(i+ny/2,2)
				l(i+ny/2+1,1)=l(i,2)
			enddo
			l(ny/2+1,1)=l(ny+1,2)
			ny=ny+1
!  note, this is a standard fftshift, not an ifftshift
! 			ny=ny-1
! 			do i=1,ny/2
! 				l(1,i+1)=l(2,i+ny/2+1)
! 				l(1,i+ny/2+1)=l(2,i)
! 			enddo
! 			l(1,1)=l(2,ny/2+1)
! 			ny=ny+1
		else
			! even
			do i=1,ny/2
				gain=l(1,i)
				l(1,i)=l(1,i+ny/2)
				l(1,i+ny/2)=gain
			enddo
		endif

        do i=2,ny
            k(i,:)=k(1,:)
        end do
        do i=2,nx
            l(:,i)=l(:,1)
        end do
        
        kl = k**2+l**2
        WHERE (kl==0.0) kl=1e-20
		
        do z=1,nz
            U=sum(domain%u(:,z,:))/(nx*ny)
            V=sum(domain%v(:,z,:))/(nx*ny)
            sig  = U*k+V*l
            denom = sig**2!-f**2
            WHERE (sig==0.0) sig=1e-20
!         # mimag=np.zeros((Ny,Nx)).astype('complex')
! 		  # two possible non-hydrostatic versions
!         # msq = (Ndsq/denom * kl).astype('complex')          # % vertical wave number, hydrostatic
!         # msq = ((Ndsq-sig**2)/denom * kl).astype('complex')          # % vertical wave number, hydrostatic
!         # mimag.imag=(np.sqrt(-msq)).real
!         # m=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)
! 
            m = sqrt(((Ndsq-sig**2)/denom * kl))         ! # % vertical wave number, hydrostatic
    		
            ineta=j*domain%fzs*exp(j*m*sum(domain%z(:,z,:)-domain%z(:,1,:)+domain%dz(:,z,:)/2)/(nx*ny))
!             what=sig*ineta

!           # removed coriolis term
            ineta=ineta/(kl/((0-m)*sig))
            uhat=k*ineta
            vhat=l*ineta
!           # with coriolis : 
!           u_hat = -m*(sig*k-i*l*f)*ineta/kl
!           v_hat = -m*(sig*l+i*k*f)*ineta/kl
!           pull it back out of fourier space. 
! 			plan = fftw_plan_dft_2d(nx,ny, what,w_hat, FFTW_BACKWARD,FFTW_ESTIMATE)
! 			call fftw_execute_dft(plan, what,w_hat)
! 			call fftw_destroy_plan(plan)
			
            plan = fftw_plan_dft_2d(nx,ny, uhat,u_hat, FFTW_BACKWARD,FFTW_ESTIMATE)
            call fftw_execute_dft(plan, uhat,u_hat)
            call fftw_destroy_plan(plan)
			
            plan = fftw_plan_dft_2d(nx,ny, vhat,v_hat, FFTW_BACKWARD,FFTW_ESTIMATE)
            call fftw_execute_dft(plan, vhat,v_hat)
            call fftw_destroy_plan(plan)
			
!             domain%w(:,z,:)=domain%w(:,z,:)+real(w_hat)!*(ny*nx)
            domain%u(:,z,1:nx-1)=domain%u(:,z,1:nx-1)+real(u_hat(:,1:nx-1)+u_hat(:,2:nx))/2
            domain%v(1:ny-1,z,:)=domain%v(1:ny-1,z,:)+real(v_hat(1:ny-1,:)+v_hat(2:ny,:))/2
		end do
		deallocate(k,l,kl,sig,denom,uhat,u_hat,vhat,v_hat,m,ineta)
    end subroutine linear_winds
    
    
    subroutine setup_linwinds(domain)
        implicit none
        type(domain_type),intent(inout)::domain
		complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:)::complex_terrain
        type(C_PTR) :: plan
        integer::nx,ny
        
        ny=size(domain%p,1)
        nx=size(domain%p,3)

        write(*,*) "Fzs setup"
        allocate(domain%fzs(ny,nx))
        allocate(complex_terrain(ny,nx))
		complex_terrain=domain%terrain
! 		calculate the fourier transform of the terrain for use in linear winds
        plan = fftw_plan_dft_2d(nx,ny, complex_terrain,domain%fzs, FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan, complex_terrain,domain%fzs)
        call fftw_destroy_plan(plan)
! 		normalize FFT by N - grid cells
		domain%fzs=domain%fzs/(nx*ny)

		allocate(domain%dzdx(ny,nx-1))
		allocate(domain%dzdy(ny-1,nx))
		domain%dzdx=sqrt((domain%terrain(:,2:nx)-domain%terrain(:,1:nx-1))**2+domain%dx**2)/domain%dx
		domain%dzdy=sqrt((domain%terrain(2:ny,:)-domain%terrain(1:ny-1,:))**2+domain%dx**2)/domain%dx
		
		deallocate(complex_terrain)
        
    end subroutine
	
	subroutine stagger_winds(domain)
        implicit none
        type(domain_type),intent(inout)::domain
		integer :: nx,ny
		
		ny=size(domain%p,1)
		nx=size(domain%p,3)
		
		domain%u(:,:,1:nx-1)=(domain%u(:,:,1:nx-1)+domain%u(:,:,2:nx))/2
		domain%v(1:ny-1,:,:)=(domain%v(1:ny-1,:,:)+domain%v(2:ny,:,:))/2
		
	end subroutine stagger_winds
	
	subroutine rotate_wind_field(domain)
        implicit none
        type(domain_type),intent(inout)::domain
		integer :: nx,ny,nz,i
		
		ny=size(domain%p,1)
		nz=size(domain%p,2)
		nx=size(domain%p,3)
		do i=1,nz
			domain%u(:,i,1:nx-1)=domain%u(:,i,1:nx-1)*domain%dzdx
			domain%v(1:ny-1,i,:)=domain%v(1:ny-1,i,:)*domain%dzdy
		end do
		
	end subroutine rotate_wind_field
	
    subroutine linear_perturb(domain)
!         Called from the simple weather model to update the U,V,W wind fields based on linear theory
!         Ndsq = squared Brunt Vaisalla frequency (1/s) typically from dry static stability
        implicit none
        type(domain_type),intent(inout)::domain
		real::stability
        
        stability=calc_stability(domain)
        if (.not.allocated(domain%fzs)) then
            call setup_linwinds(domain)
        endif
		
! 		domain%dz=domain%dz*3
		
		call stagger_winds(domain)
		call linear_winds(domain,stability)
		call rotate_wind_field(domain)
        
        
    end subroutine linear_perturb
end module linear_theory_winds