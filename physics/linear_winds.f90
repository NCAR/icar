!>----------------------------------------------------------
!!
!! This module provides the linear wind theory calculations
!! 
!! The main entry point to the code is:
!! 	 	linear_perturb(domain, options, vsmooth, reverse, useDensity)
!!
!! Call tree graph :
!!	linear_perturb->[ setup_linwinds -> add_buffer_topo,
!! 					  calc_stability,
!!					  linear_winds -> various fft routines]
!! 
!! High level routine descriptions / purpose
!!   calc_stability		- calculates a mean Brunt Vaisala frequency over the domain
!!   linear_winds		- primary routine that calculates the linear wind perturbation
!!	 add_buffer_topo	- generates a topo grid that with a surrounding buffer for the fft
!! 	 setup_linwinds		- sets up module level variables and calls add_buffer_topo
!! 	 linear_perturb		- main entry point, calls setup on first entry for a given domain
!! 
!! Inputs: domain, options, vsmooth, reverse, useDensity
!! 		domain,options	= as defined in data_structures
!! 		vsmooth			= number of vertical levels to smooth winds over
!! 		reverse 		= remove linear perturbation instead of adding it
!!		useDensity		= create a linear field that attempts to mitigate the 
!! 							boussinesq approx that is embedded in the linear theory
!!							so it advection can properly incorporate density.
!!
!! Author : Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module linear_theory_winds
    use fft ! note fft module is defined in fftshift.f90
	use fftshifter
    use data_structures
	use io_routines, 		only : io_write2d
	use output, 			only : write_domain
    implicit none
	private
	public::linear_perturb
	
	logical :: variable_N
	real :: N_squared
	!! unfortunately these have to be allocated every call because we could be calling on both the high-res
	!! domain and the low res domain (to "remove" the linear winds)
    real,allocatable,dimension(:,:)::k,l,kl,sig
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:)::denom,uhat,u_hat,vhat,v_hat,m,ineta,msq,mimag
	integer,parameter::buffer=50 ! number of grid cells to buffer around the domain
	integer,parameter::stability_window_size=2
	
	real, parameter :: max_stability= 5e-4 ! limits on the calculated Brunt Vaisala Frequency
	real, parameter :: min_stability= 5e-7 ! these may need to be a little narrower. 
	
contains
	
    pure function calc_stability(domain) result(BV_freq)
        implicit none
        class(linearizable_type),intent(in)::domain
		real :: BV_freq
		real, allocatable, dimension(:) :: vertical_N
		integer :: i, nx,nz,ny, top_layer,bottom_layer
		real :: dz,dtheta,t_mean
		
		nx=size(domain%th,1)
		nz=size(domain%th,2)
		ny=size(domain%th,3)
		allocate(vertical_N(nz))
		
		if (variable_N) then
			do i=1,nz
				top_layer    = max(i+stability_window_size, nz)
				bottom_layer = min(i-stability_window_size,  i)
				
				dz     = sum(domain%z (:,top_layer,:)-domain%z (:,bottom_layer,:))/(nx*ny) ! distance between top layer and bottom layer
				dtheta = sum(domain%th(:,top_layer,:)-domain%th(:,bottom_layer,:))/(nx*ny) ! average temperature change between layers
				t_mean = sum(domain%th(:,i,:))/(nx*ny) ! mean temperature in the current layer
				
				! BV frequency = gravity/theta * dtheta/dz
				vertical_N(i)=gravity/t_mean * dtheta/dz
			end do
			! calculate the mean stability over the entire profile
			BV_freq=sum(vertical_N)/nz
			! impose limits so the linear solution doesn't go crazy in really stable or unstable air
			BV_freq=max(min(BV_freq, min_stability), max_stability)
		else
			! or just use the supplied BV frequency
			BV_freq=N_squared
		endif
		
		deallocate(vertical_N)
		! below are other possible calculations for Nd using the "dry" environmental lapse rate
		! real, parameter :: R  = 287.0
		! real, parameter :: Rv = 461.0
		! real, parameter :: cp = 1004.0
		! real, parameter :: L   = 2.5e6
		! real, parameter :: g  = 9.81
		! real, parameter :: ratio = 18.015/28.964
		! real, parameter :: t0 = 273.15
		! 
		! p0=domain%.p(:,1,:)
		! pii=1.0/((100000.0/domain%p)**(R/cp))
		! T2m=domain%th(:,1,:)*pii(:,1,:)
		!     
		!  es = 611.21*exp(17.502*(T2m-t0)/(T2m-32.19))
		!  qs0 = ratio * es/(p0-es)
		! 
		! cap_gamma = -(g * (1.+(L*qs0)/(R*T2m)) / (cp + (L**2 * qs0*ratio) / (R*T2m**2)))
		!! STILL NEEDS TO BE CONVERTED FROM PYTHON
		!! env_gamma = np.mean(np.diff(weather.th*pii,axis=0)/np.diff(base.hgt3d,axis=0),axis=0)
		! 
		! dry_gamma=sum(env_gamma-cap_gamma)/size(env_gamma)
		! ndsq=(g/(sum(T2m)/size(T2m)))*(dry_gamma)
		! ndsq=max(min(1e-4,ndsq),1e-8)
		! calc_stability=ndsq
    end function calc_stability

	! Compute linear wind perturbations to U and V and add them back to the domain
	!
	! see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
	! -------------------------------------------------------------------------------
	! Ndsq  = 0.003**2      # dry BV freq sq. was 0.01**2 initially, Idar suggested 0.005**2
	! 1E-8 keeps it relatively stable, no wild oscillations in u,v
	! but 1E-8 doesn't damp vertical winds with height very fast
	! could/should be calculated from the real atm profile with limits
	! f  = 9.37e-5           # rad/s Coriolis frequency for 40deg north
	! ---------------------------------------------------------------------------------
    subroutine linear_winds(domain,Ndsq,vsmooth,reverse_flag,useDensity,debug)
		use, intrinsic :: iso_c_binding
        implicit none
        class(linearizable_type),intent(inout)::domain
        real, intent(in)::Ndsq
		integer, intent(in)::vsmooth ! number of layers to smooth winds over in the vertical
		logical, intent(in), optional :: reverse_flag,useDensity,debug
		logical::reverse
        complex,parameter :: j= (0,1)
        real::gain,offset !used in setting up k and l arrays
        integer::nx,ny,nz,i,midpoint,z,realnx,realny,x,y,bottom,top
        real::U,V
		real,dimension(:),allocatable::U_layers,V_layers,preU_layers,preV_layers
		real,parameter::pi=3.1415927
        type(C_PTR) :: plan
        
		if (.not.present(reverse_flag)) then
			reverse=.False.
		else
			reverse=reverse_flag
		endif
		
		nx=size(domain%fzs,1)
		nz=size(domain%u,2)
		ny=size(domain%fzs,2)
		allocate(U_layers(nz))
		allocate(V_layers(nz))
		allocate(preU_layers(nz))
		allocate(preV_layers(nz))
		
		realnx=size(domain%z,1)
		realny=size(domain%z,3)
		
		do i=1,nz
			preU_layers(i)=sum(domain%u(:(realnx-1),i,:))/((realnx-1)*realny)
			preV_layers(i)=sum(domain%v(:,i,:(realny-1)))/(realnx*(realny-1))
		enddo
		do i=1,nz
			bottom=i-vsmooth
			top=i+vsmooth
			if (bottom<1) then
				bottom=1
			endif
			if (top>nz) then
				top=nz
			endif
			U_layers(i)=sum(preU_layers(bottom:top))/(top-bottom+1)
			V_layers(i)=sum(preV_layers(bottom:top))/(top-bottom+1)
		enddo
		
        
		! these should be stored in a separate data structure... and allocated/deallocated in a subroutine
		! for now these are reallocated/deallocated everytime so we can use it for different sized domains (e.g. coarse and fine)
		! maybe linear winds need to be embedded in an object instead of a module to avoid this problem...
		if (.not.allocated(k)) then
	        allocate(k(nx,ny))
	        allocate(l(nx,ny))
	        allocate(kl(nx,ny))
	        allocate(sig(nx,ny))
	        allocate(denom(nx,ny))
	        allocate(uhat(nx,ny))
	        allocate(u_hat(nx,ny))
	        allocate(vhat(nx,ny))
	        allocate(v_hat(nx,ny))
	        allocate(m(nx,ny))
	        allocate(msq(nx,ny))
	        allocate(mimag(nx,ny))
	        allocate(ineta(nx,ny))
        
			! Compute 2D k and l wavenumber fields (could be stored in options or domain or something)
	        offset=pi/domain%dx
	        gain=2*offset/(nx-1)
			k(:,1) = (/((i*gain-offset),i=0,nx-1)/)
			do i=2,ny
				k(:,i)=k(:,1)
			enddo
			
	        gain=2*offset/(ny-1)
			l(1,:) = (/((i*gain-offset),i=0,ny-1)/)
			do i=2,nx
				l(i,:)=l(1,:)
			enddo
			
			! finally compute the kl combination array
	        kl = k**2+l**2
	        WHERE (kl==0.0) kl=1e-15
		endif

		
		! process each array independantly
		! note fftw_plan_... may not be threadsafe so this can't be OMP parallelized without pulling that outside of the loop. 
		! to parallelize, compute uhat,vhat in parallel
		! then compute the plans serially
		! then perform ffts etc in parallel
		! finally destroy plans serially
		m=1
        do z=1,nz
            U=U_layers(z)
            V=V_layers(z)
			if ((abs(U)+abs(V))>0.5) then
	            sig  = U*k+V*l
	            where(sig==0.0) sig=1e-15
	            denom = sig**2!-f**2
			
				! 	where(denom.eq.0) denom=1e-20
				! # two possible non-hydrostatic versions
				! not yet converted from python...
				! # mimag=np.zeros((Ny,Nx)).astype('complex')
				! # msq = (Ndsq/denom * kl).astype('complex')          			# % vertical wave number, hydrostatic
				! # msq = ((Ndsq-sig**2)/denom * kl).astype('complex')          # % vertical wave number, hydrostatic
				! # mimag.imag=(np.sqrt(-msq)).real
				! # m=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)
				! 
				! mimag=np.zeros(Nx).astype('complex')
				! mimag.imag=(np.sqrt(-msq)).real
				! m=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)

				msq = Ndsq/denom * kl
				mimag=0 ! be sure to reset real and imaginary components
			    mimag=mimag+(real(sqrt(-msq))*j)
				! m=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)
	            m = sqrt(msq)         ! # % vertical wave number, hydrostatic
				where(sig<0) m=m*(-1)   ! equivilant to m=m*sign(sig)
	    		where(real(msq)<0) m=mimag
			
	            ineta=j*domain%fzs*exp(j*m* &
						sum(domain%z(:,z,:)-domain%z(:,1,:)) &
						/(realnx*realny))
				!  what=sig*ineta

				!# with coriolis : [+/-]j*[l/k]*f
				!uhat = (0-m)*(sig*k-j*l*f)*ineta/kl
				!vhat = (0-m)*(sig*l+j*k*f)*ineta/kl
				!# removed coriolis term
	            ineta=ineta/(kl/((0-m)*sig))
	            uhat=k*ineta
	            vhat=l*ineta

				! pull it back out of fourier space. 
				! NOTE, the fftw transform inherently scales by N so the Fzs/Nx/Ny provides the only normalization necessary (I think)

				! it should be possible to store the plan and execute it everytime rather than recreating it everytime, doesn't matter too much 
				call ifftshift(uhat)
				call ifftshift(vhat)
			
	            plan = fftw_plan_dft_2d(ny,nx, uhat,u_hat, FFTW_BACKWARD,FFTW_ESTIMATE)
	            call fftw_execute_dft(plan, uhat,u_hat)
	            call fftw_destroy_plan(plan)
			
	            plan = fftw_plan_dft_2d(ny,nx, vhat,v_hat, FFTW_BACKWARD,FFTW_ESTIMATE)
	            call fftw_execute_dft(plan, vhat,v_hat)
	            call fftw_destroy_plan(plan)
			
				! u/vhat are first staggered to apply to u/v appropriately...
				! NOTE: we should be able to do this without the loop, but ifort -O was giving the wrong answer... possible compiler bug version 12.1.x?
				do i=1,realny
					u_hat(1+buffer:nx-buffer-1,i+buffer) = (u_hat(1+buffer:nx-buffer-1,i+buffer)+u_hat(2+buffer:nx-buffer,i+buffer))/2
					v_hat(1+buffer:nx-buffer,i+buffer)   = (v_hat(1+buffer:nx-buffer,i+buffer)+v_hat(1+buffer:nx-buffer,i+buffer+1))/2
				enddo
				if (present(useDensity)) then
					! if we are using density in the advection calculations, modify the linear perturbation
					! to get the vertical velocities closer to what they would be without density (boussinesq)
					if (useDensity) then
						write(*,*), "Using a density correction in linear winds"
						u_hat(1+buffer:nx-buffer,1+buffer:ny-buffer) = &
							2*real(u_hat(1+buffer:nx-buffer,1+buffer:ny-buffer))! / domain%rho(1:realnx,z,1:realny)
						v_hat(1+buffer:nx-buffer,1+buffer:ny-buffer) = &
							2*real(v_hat(1+buffer:nx-buffer,1+buffer:ny-buffer))! / domain%rho(1:realnx,z,1:realny)
					endif
				endif
				! if we are removing linear winds from a low res field, subtract u_hat v_hat instead
				! real(real()) extracts real component of complex, then converts to a real data type (may not be necessary except for IO?)
				if (reverse) then
		            domain%u(1:realnx-1,z,:)=domain%u(1:realnx-1,z,:) - &
						real(real(u_hat(1+buffer:nx-buffer-1,1+buffer:realny+buffer  ) ))
		            domain%v(:,z,1:realny-1)=domain%v(:,z,1:realny-1) - &
						real(real(v_hat(1+buffer:nx-buffer,  1+buffer:realny+buffer-1) ))
				else
		            domain%u(1:realnx-1,z,:)=domain%u(1:realnx-1,z,:) + &
						real(real(u_hat(1+buffer:nx-buffer-1,1+buffer:realny+buffer  ) ))
		            domain%v(:,z,1:realny-1)=domain%v(:,z,1:realny-1) + &
						real(real(v_hat(1+buffer:nx-buffer,  1+buffer:realny+buffer-1) ))
				endif
			
				if (present(debug).and.(z==1))then
					if (debug) then
						write(*,*) "Ndsq = ", Ndsq
						write(*,*) "U=",U, "    V=",V
						write(*,*) "realnx=",realnx, "; nx=",nx, "; buffer=",buffer
						write(*,*) "realny=",realny, "; ny=",ny!, buffer
						write(*,*) "Writing internal linear wind data"
						call io_write2d("u_hat_sub2.nc","data",real(real(u_hat(1+buffer:nx-buffer-1,1+buffer:realny+buffer))) )
						call io_write2d("v_hat_sub2.nc","data",real(real(v_hat(1+buffer:nx-buffer,1+buffer:realny+buffer-1))) )
					endif
				endif
			endif
		end do ! z-loop
		
		
		! finally deallocate all temporary arrays that were created... chould be a datastructure and a subroutine...
		deallocate(k,l,kl,sig,denom,uhat,u_hat,vhat,v_hat,m,ineta,msq,mimag)
		deallocate(U_layers,V_layers,preU_layers,preV_layers)
    end subroutine linear_winds
    
    
	subroutine add_buffer_topo(terrain,buffer_topo)
		! add a smoothed buffer around the edge of the terrain to prevent crazy wrap around effects
		! in the FFT due to discontinuities between the left and right (top and bottom) edges of the domain
        implicit none
		real, dimension(:,:), intent(in) :: terrain
		complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:), intent(out):: buffer_topo
		real, dimension(:,:),allocatable :: real_terrain
		integer::nx,ny,i,pos
		real::weight
		nx=size(terrain,1)+buffer*2
		ny=size(terrain,2)+buffer*2
		allocate(buffer_topo(nx,ny))
		buffer_topo=minval(terrain)
		
		buffer_topo(1+buffer:nx-buffer,1+buffer:ny-buffer)=terrain
		do i=1,buffer
			weight=i/(real(buffer)*2)
			pos=buffer-i
			buffer_topo(pos+1,1+buffer:ny-buffer)  =terrain(1,1:ny-buffer*2)*(1-weight)+terrain(nx-buffer*2,1:ny-buffer*2)*   weight
			buffer_topo(nx-pos,1+buffer:ny-buffer) =terrain(1,1:ny-buffer*2)*(  weight)+terrain(nx-buffer*2,1:ny-buffer*2)*(1-weight)
		enddo
		do i=1,buffer
			weight=i/(real(buffer)*2)
			pos=buffer-i
			buffer_topo(:,pos+1)  =buffer_topo(:,buffer+1)*(1-weight) + buffer_topo(:,ny-buffer)*   weight
			buffer_topo(:,ny-pos) =buffer_topo(:,buffer+1)*(  weight) + buffer_topo(:,ny-buffer)*(1-weight)
		enddo
		allocate(real_terrain(nx,ny))
		real_terrain=buffer_topo
		call io_write2d("complex_terrain.nc","data",real_terrain)
		deallocate(real_terrain)
		
	end subroutine add_buffer_topo
	
	! called from linear_perturb the first time perturb is called
	! compute FFT(terrain), and dzdx,dzdy components
    subroutine setup_linwinds(domain,options)
        implicit none
        class(linearizable_type),intent(inout)::domain
		type(options_type),intent(in) :: options
		complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:)::complex_terrain
        type(C_PTR) :: plan
        integer::nx,ny
        
		! store module level variables so we don't have to pass options through everytime
		N_squared=options%N_squared
		variable_N=options%variable_N
		
		call add_buffer_topo(domain%terrain,complex_terrain)
        nx=size(complex_terrain,1)
        ny=size(complex_terrain,2)

        write(*,*) "Initializing linear winds : ",nx,ny
        allocate(domain%fzs(nx,ny))
		
		! calculate the fourier transform of the terrain for use in linear winds
        plan = fftw_plan_dft_2d(ny,nx, complex_terrain,domain%fzs, FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan, complex_terrain,domain%fzs)
        call fftw_destroy_plan(plan)
		! normalize FFT by N - grid cells
		domain%fzs=domain%fzs/(nx*ny)
		call fftshift(domain%fzs)
		
		! cleanup temporary array
		deallocate(complex_terrain)
        
    end subroutine
	
	! Primary entry point!
	! Called from ICAR to update the U,V,W wind fields based on linear theory
    subroutine linear_perturb(domain,options,vsmooth,reverse,useDensity)
        implicit none
        class(linearizable_type),intent(inout)::domain
		type(options_type), intent(in) :: options
		integer, intent(in) :: vsmooth
		logical, intent(in), optional :: reverse,useDensity
		logical, save :: debug=.True.
		real::stability
        
		! if linear_perturb hasn't been called before we need to perform some setup actions. 
        if (.not.allocated(domain%fzs)) then
            call setup_linwinds(domain,options)
        endif
		
		! Ndsq = squared Brunt Vaisalla frequency (1/s) typically from dry static stability
        stability=calc_stability(domain)
		
		call linear_winds(domain,stability,vsmooth,reverse,useDensity,debug)
		debug=.False.
        
    end subroutine linear_perturb
end module linear_theory_winds