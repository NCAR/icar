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
!!   initialize_spatial_winds - generated the look up tables to use spatially varying linear winds
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
	use io_routines, 		only : io_write3d, io_write2d, io_read2d
	use output, 			only : write_domain
	use string, 			only : str
    implicit none
	private
	public::linear_perturb
	
	logical :: variable_N
	real :: N_squared
	!! unfortunately these have to be allocated every call because we could be calling on both the high-res
	!! domain and the low res domain (to "remove" the linear winds)
    real,allocatable,dimension(:,:)::k,l,kl,sig
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:)::denom,uhat,u_hat,vhat,v_hat,m,ineta,msq,mimag
	integer::buffer=25 ! number of grid cells to buffer around the domain MUST be >=1
	integer,parameter::stability_window_size=2
	
	real, parameter :: max_stability= 5e-4 ! limits on the calculated Brunt Vaisala Frequency
	real, parameter :: min_stability= 5e-7 ! these may need to be a little narrower. 
	real :: linear_contribution = 1.0 ! multiplier on uhat,vhat before adding to u,v
	
	real, allocatable, dimension(:) :: u_values, v_values
	real, allocatable, target, dimension(:,:,:,:,:) :: hi_u_LUT, hi_v_LUT, rev_u_LUT, rev_v_LUT
	real, pointer, dimension(:,:,:,:,:) :: u_LUT, v_LUT
	real, allocatable, dimension(:,:) :: linear_mask
	
	real, parameter :: umax=30
	real, parameter :: umin=-20
	real, parameter :: vmax=20
	real, parameter :: vmin=-20
	
	integer, parameter :: n_U_values=30
	integer, parameter :: n_V_values=30
	
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
    subroutine linear_winds(domain,Ndsq,vsmooth,reverse,useDensity,debug,fixedU,fixedV)
		use, intrinsic :: iso_c_binding
        implicit none
        class(linearizable_type),intent(inout)::domain
        real, intent(in)::Ndsq
		integer, intent(in)::vsmooth ! number of layers to smooth winds over in the vertical
		logical, intent(in) :: reverse,useDensity
		logical, intent(in), optional ::debug
		real, intent(in), optional :: fixedU, fixedV !instead of computing U and V from domain, just use these
        complex,parameter :: j= (0,1)
        real::gain,offset !used in setting up k and l arrays
        integer::nx,ny,nz,i,midpoint,z,realnx,realny,realnx_u,realny_v,x,y,bottom,top
		logical :: staggered
        real::U,V
		real,dimension(:),allocatable::U_layers,V_layers,preU_layers,preV_layers
		real,parameter::pi=3.1415927
        type(C_PTR) :: plan
        
		nx=size(domain%fzs,1)
		nz=size(domain%u,2)
		ny=size(domain%fzs,2)
		allocate(U_layers(nz))
		allocate(V_layers(nz))
		allocate(preU_layers(nz))
		allocate(preV_layers(nz))
		
		realnx=size(domain%z,1)
		realnx_u=size(domain%u,1)
		realny=size(domain%z,3)
		realny_v=size(domain%v,3)
		staggered = (realnx/=realnx_u).and.(realny/=realny_v)
		
		if (present(fixedV)) then
			if (.not.present(fixedU)) then
				write(*,*) "Fixed U and V fields must both be supplied"
				stop
			endif
			V_layers=fixedV
			U_layers=fixedV
		else
			do i=1,nz
				preU_layers(i)=sum(domain%u(:(realnx_u-1),i,:))/((realnx_u-1)*realny)
				preV_layers(i)=sum(domain%v(:,i,:(realny_v-1)))/(realnx*(realny_v-1))
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
		endif
		
        
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
				mimag=0+0*j ! be sure to reset real and imaginary components
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
			
				! u/vhat are first staggered to apply to u/v appropriately if on a staggered grid (real_nx/=real_nx_u)
				! when removing linear winds from forcing data, it may NOT be on a staggered grid
				! NOTE: we should be able to do this without the loop, but ifort -O was giving the wrong answer... possible compiler bug version 12.1.x?
				if (staggered) then
					! note, buffer must be AT LEAST 1
					do i=0,realny_v
						u_hat(buffer:realnx_u+buffer,i+buffer) = &
							(u_hat(buffer:realnx_u+buffer,i+buffer) + u_hat(1+buffer:realnx_u+buffer+1,i+buffer))/2
						v_hat(1+buffer:realnx+buffer,i+buffer) = &
							(v_hat(1+buffer:realnx+buffer,i+buffer) + v_hat(1+buffer:realnx+buffer,i+buffer+1))/2
					enddo
				endif
				if (useDensity) then
					! if we are using density in the advection calculations, modify the linear perturbation
					! to get the vertical velocities closer to what they would be without density (boussinesq)
					! need to check if this makes the most sense when close to the surface
					if (debug) then
						write(*,*), "Using a density correction in linear winds"
					endif
					u_hat(buffer:realnx_u+buffer,buffer:realny+buffer) = &
						2*real(u_hat(buffer:realnx_u+buffer,buffer:realny+buffer))! / domain%rho(1:realnx,z,1:realny)
					v_hat(buffer:realnx+buffer,buffer:realny_v+buffer) = &
						2*real(v_hat(buffer:realnx+buffer,buffer:realny_v+buffer))! / domain%rho(1:realnx,z,1:realny)
				endif
				! if we are removing linear winds from a low res field, subtract u_hat v_hat instead
				! real(real()) extracts real component of complex, then converts to a real data type (may not be necessary except for IO?)
				if (reverse) then
					if (staggered) then
						! Forcing data need to go back to the correct position in the grid (e.g. 2:nx-1 not 1:nx-2)
			            domain%u(:,z,:)=domain%u(:,z,:) - &
							real(real( u_hat(1+buffer:realnx+buffer+1, 1+buffer:realny+buffer  ) ))*linear_contribution
							
			            domain%v(:,z,:)=domain%v(:,z,:) - &
							real(real( v_hat(1+buffer:realnx+buffer,   1+buffer:realny+buffer+1) ))*linear_contribution
					else
						! if we are not on a staggered grid, just get winds from inside the buffered area
			            domain%u(:,z,:)=domain%u(:,z,:) - &
							real(real( u_hat(1+buffer:realnx_u+buffer, 1+buffer:realny+buffer  ) ))*linear_contribution
							
			            domain%v(:,z,:)=domain%v(:,z,:) - &
							real(real( v_hat(1+buffer:realnx+buffer,   1+buffer:realny_v+buffer) ))*linear_contribution
					endif
				else
					! the internal domain is staggered by definition
		            domain%u(:,z,:)=domain%u(:,z,:) + &
						real(real(u_hat(1+buffer:realnx_u+buffer,1+buffer:realny+buffer  ) ))*linear_contribution
		            domain%v(:,z,:)=domain%v(:,z,:) + &
						real(real(v_hat(1+buffer:realnx+buffer,  1+buffer:realny_v+buffer) ))*linear_contribution
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
						call io_write2d("u_hat_full.nc","data",real(real(u_hat)) )
						call io_write2d("v_hat_full.nc","data",real(real(v_hat)) )
					endif
				endif
			endif
		end do ! z-loop
		
		
		! finally deallocate all temporary arrays that were created... chould be a datastructure and a subroutine...
		deallocate(k,l,kl,sig,denom,uhat,u_hat,vhat,v_hat,m,ineta,msq,mimag)
		deallocate(U_layers,V_layers,preU_layers,preV_layers)
    end subroutine linear_winds
    
    
	subroutine add_buffer_topo(terrain,buffer_topo,smooth_window)
		! add a smoothed buffer around the edge of the terrain to prevent crazy wrap around effects
		! in the FFT due to discontinuities between the left and right (or top and bottom) edges of the domain
        implicit none
		real, dimension(:,:), intent(in) :: terrain
		complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:), intent(inout) :: buffer_topo
		integer, intent(in) :: smooth_window
		real, dimension(:,:),allocatable :: real_terrain
		integer::nx,ny,i,j,pos, xs,xe,ys,ye, window
		real::weight
		nx=size(terrain,1)+buffer*2
		ny=size(terrain,2)+buffer*2
		allocate(buffer_topo(nx,ny))
		buffer_topo=minval(terrain)
		
		print*, buffer, nx, ny
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
		
		if (smooth_window>0) then
			do j=1,buffer
				window=min(j,smooth_window)
				do i=1,nx
					xs=max(1, i-window)
					xe=min(nx,i+window)
					
					ys=max(1, buffer-j+1-window)
					ye=min(ny,buffer-j+1+window)
				
					buffer_topo(i,buffer-j+1)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))

					ys=max(1, ny-(buffer-j)-window)
					ye=min(ny,ny-(buffer-j)+window)
				
					buffer_topo(i,ny-(buffer-j))=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))
				end do
				do i=1,ny
					xs=max(1, buffer-j+1-window)
					xe=min(nx,buffer-j+1+window)
					ys=max(1, i-window)
					ye=min(ny,i+window)
				
					buffer_topo(buffer-j+1,i)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))

					xs=max(1, nx-(buffer-j)-window)
					xe=min(nx,nx-(buffer-j)+window)
				
					buffer_topo(nx-(buffer-j),i)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))
				end do
			end do
		endif
		allocate(real_terrain(nx,ny))
		real_terrain=buffer_topo
		call io_write2d("complex_terrain.nc","data",real_terrain)
		deallocate(real_terrain)
		
	end subroutine add_buffer_topo
	
	subroutine initialize_spatial_winds(domain,options,reverse,useDensity)
		! compute look up tables for all combinations of  N different U speeds and N different V speeds
		implicit none
		type(linearizable_type), intent(inout) :: domain
		type(options_type), intent(in) :: options
		logical, intent(in) :: reverse,useDensity
		real, allocatable, dimension(:,:,:) :: savedU, savedV
		real :: u,v
		integer :: nx,ny,nz,i,j,k
		logical :: debug
		
		! the domain to work over
		nx=size(domain%u,1)-1
		nz=size(domain%u,2)
		ny=size(domain%u,3)
		
		! save the old U and V values so we can restore them
		allocate(savedU(nx+1,nz,ny))
		allocate(savedV(nx,nz,ny+1))
		savedU=domain%u
		savedV=domain%v
		
		! create the array of U and V values to create LUTs for
		if (.not.allocated(u_values)) then
			allocate(u_values(n_U_values))
			allocate(v_values(n_V_values))
		endif
		do i=1,n_U_values
			u_values(i)=(i-1)/real(n_U_values-1) * (umax-umin) + umin
		enddo
		do i=1,n_V_values
			v_values(i)=(i-1)/real(n_V_values-1) * (vmax-vmin) + vmin
		enddo
		
		! allocate the (LARGE) look up tables for both U and V
		if (reverse) then
			allocate(rev_u_LUT(n_U_values,n_V_values,nx+1,nz,ny))
			allocate(rev_v_LUT(n_U_values,n_V_values,nx,nz,ny+1))
			u_LUT=>rev_u_LUT
			v_LUT=>rev_v_LUT
		else
			allocate(hi_u_LUT(n_U_values,n_V_values,nx+1,nz,ny))
			allocate(hi_v_LUT(n_U_values,n_V_values,nx,nz,ny+1))
			u_LUT=>hi_u_LUT
			v_LUT=>hi_v_LUT
		endif
		
		! loop over combinations of U and V values
		write(*,*) "Percent Completed:"
		debug=.True.
! 		call io_write2d("internal_linear_mask.nc","data",linear_mask)
		do i=1,n_U_values
			write(*,*) i/real(n_U_values)*100," %"
			do j=1,n_V_values
				
				! set the domain wide U and V values to the current u and v values
				domain%u=u_values(i)
				domain%v=v_values(j)
				
				! calculate the linear wind field for the current u and v values
				call linear_winds(domain,N_squared, 0, reverse,useDensity,debug=debug)!,fixedU=u_values(i),fixedV=v_values(j))
				debug=.False.
				do k=1,nz
					u_LUT(i,j,:,k,:)=(domain%u(:,k,:)-u_values(i)) * linear_mask(:,:ny)
					v_LUT(i,j,:,k,:)=(domain%v(:,k,:)-v_values(j)) * linear_mask(:nx,:)
				end do
! 				call io_write3d("v_lut"//trim(str(i))//"_"//trim(str(j))//".nc","data",domain%v)
! 				call io_write3d("u_lut"//trim(str(i))//"_"//trim(str(j))//".nc","data",domain%u)
			end do
		end do
		domain%u=savedU
		domain%v=savedV
		
		deallocate(savedU,savedV)
		
	end subroutine initialize_spatial_winds
	
	function calc_weight(indata, bestpos, nextpos, match) result(weight)
		! simply calculate the weights between the positions bestpos and nextpos
		! based on the distance between match and indata(nextpos) (normatlized by nextpos - bestpos)
		! assumes indata is monotonically increasing, 
		! bestpos must be set prior to entry
		! nextpos is calculated internally (either bestpos+1 or n)
		implicit none
		real :: weight
		real, dimension(:), intent(in) :: indata
		integer, intent(in) :: bestpos
		integer, intent(inout) :: nextpos
		real, intent(in) :: match
		
		integer :: n
		
		n=size(indata)
		
		if (match<indata(1)) then
			nextpos=1
			weight=1
		else
			if (match>indata(n)) then
				nextpos=n
				weight=1
			else
				nextpos=bestpos+1
				weight=(indata(nextpos)-match) / (indata(nextpos) - indata(bestpos))
			endif
		endif
	
	end function
	
	subroutine spatial_winds(domain,reverse)
		! compute a spatially variable linear wind perturbation
		! based off of look uptables computed in via setup
		! for each grid point, find the closest LUT data in U and V space
		! then bilinearly interpolate the nearest LUT values for that points linear wind field
		implicit none
		type(linearizable_type), intent(inout) :: domain
		logical, intent(in) :: reverse
		integer :: nx,ny,nz,i,j,k
		integer :: uk, vi !store a separate value of i for v and of k for u to we can handle nx+1, ny+1
		integer :: step, upos, vpos, nextu, nextv
		real :: uweight, vweight
	
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
		nz=size(domain%u,2)
		
		if (reverse) then
			u_LUT=>rev_u_LUT
			v_LUT=>rev_v_LUT
		else
			u_LUT=>hi_u_LUT
			v_LUT=>hi_v_LUT
		endif
		
		
		!$omp parallel firstprivate(nx,ny,nz), &
		!$omp private(i,j,k,step, upos, vpos, nextu, nextv, uweight, vweight), &
		!$omp shared(domain, u_values, v_values, u_LUT, v_LUT)
		!$omp do
		do k=1,ny+1
			do j=1,nz
				do i=1,nx+1
					if (k<=ny) then
						uk=k
					else
						uk=k-1
					endif
					if (i<=nx) then
						vi=i
					else
						vi=i-1
					endif
					
					upos=1
					do step=1,n_U_values
						if (domain%u(i,j,uk)>u_values(step)) then
							upos=step
						endif
					end do
					vpos=1
					do step=1,n_V_values
						if (domain%v(vi,j,k)>v_values(step)) then
							vpos=step
						endif
					end do
					
					! calculate the weights and the "next" u/v position
					! "next" usually = pos+1 but for edge cases next = 1 or n
					uweight=calc_weight(u_values, upos,nextu,domain%u(i,j,uk))
					vweight=calc_weight(v_values, vpos,nextv,domain%v(vi,j,k))
					
					if (k<=ny) then
					! perform linear interpolation between LUT values
						domain%u(i,j,k)=domain%u(i,j,k) &
									+ (   vweight  * (uweight * u_LUT(upos,vpos,i,j,k)  + (1-uweight) * u_LUT(nextu,vpos,i,j,k))  &
									+  (1-vweight) * (uweight * u_LUT(upos,nextv,i,j,k) + (1-uweight) * u_LUT(nextu,nextv,i,j,k)) )
					endif
					if (i<=nx) then
						domain%v(i,j,k)=domain%v(i,j,k) &
									+ (   vweight  * (uweight * v_LUT(upos,vpos,i,j,k)  + (1-uweight) * v_LUT(nextu,vpos,i,j,k))  &
									+  (1-vweight) * (uweight * v_LUT(upos,nextv,i,j,k) + (1-uweight) * v_LUT(nextu,nextv,i,j,k)) )
					endif
					
				end do
			end do
		end do
		!$omp end do
		!$omp end parallel
		
	end subroutine spatial_winds
	
	! called from linear_perturb the first time perturb is called
	! compute FFT(terrain), and dzdx,dzdy components
    subroutine setup_linwinds(domain,options,reverse,useDensity)
        implicit none
        class(linearizable_type),intent(inout)::domain
		type(options_type),intent(in) :: options
		logical, intent(in) :: reverse,useDensity
		complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:)::complex_terrain_firstpass
		complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:)::complex_terrain
        type(C_PTR) :: plan
        integer::nx,ny, save_buffer
		
		! store module level variables so we don't have to pass options through everytime
		variable_N=options%variable_N
		
! 		call add_buffer_topo(domain%terrain,complex_terrain,5)
		call add_buffer_topo(domain%terrain,complex_terrain_firstpass,5)
		save_buffer=buffer
		buffer=2
		call add_buffer_topo(real(real(complex_terrain_firstpass)),complex_terrain,0)
		buffer=buffer+save_buffer
		
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
		if (allocated(complex_terrain_firstpass)) then
			deallocate(complex_terrain_firstpass)
		endif
		
		if (linear_contribution/=1) then
			write(*,*) "Using a fraction of the linear perturbation:",linear_contribution
		endif
		
		if (options%spatial_linear_fields) then
			if ((.not.allocated(hi_u_LUT) .and. (.not.reverse)) .or. ((.not.allocated(rev_u_LUT)) .and. reverse)) then
				
				nx=size(domain%terrain,1)+1
				ny=size(domain%terrain,2)+1
				allocate(linear_mask(nx,ny))
				linear_mask=1.0
				if (.not.reverse) then
		! 			if (options%linear_mask) then
					if (.True.) then
						! temporary hack to read in a mask from a hard coded filename
						! THIS MUST CHANGE
						print*, "WARNING READING HARD CODED linmask.nc FILE"
						call io_read2d("linmask.nc","data",domain%linear_mask)
						linear_mask(:nx-1,:ny-1)=1-(1-domain%linear_mask)*0.8
						linear_mask(nx,:ny-1)=linear_mask(nx-1,:ny-1)
						linear_mask(:nx-1,ny)=linear_mask(:nx-1,ny-1)
						linear_mask(nx,ny)=linear_mask(nx-1,ny-1)
					endif
				endif
			
				write(*,*) "Generating a spatially variable linear perturbation look up table"
				call initialize_spatial_winds(domain,options,reverse,useDensity)
				deallocate(linear_mask)
			else
				write(*,*) "Skipping spatial wind field for presumed domain repeat"
			endif
		endif
		
        
    end subroutine setup_linwinds
	
	! Primary entry point!
	! Called from ICAR to update the U,V,W wind fields based on linear theory
    subroutine linear_perturb(domain,options,vsmooth,reverse,useDensity)
        implicit none
        class(linearizable_type),intent(inout)::domain
		type(options_type), intent(in) :: options
		integer, intent(in) :: vsmooth
		logical, intent(in), optional :: reverse,useDensity
		logical :: rev, useD
		logical, save :: debug=.True.
		real::stability

		if (present(reverse)) then
			rev=reverse
		else
			rev=.False.
		endif
		if (present(useDensity)) then
			useD=useDensity
		else
			useD=.False.
		endif

		if (rev) then
			linear_contribution=options%rm_linear_contribution
			N_squared=options%rm_N_squared
		else
			linear_contribution=options%linear_contribution
			N_squared=options%N_squared
		endif
        
		! if linear_perturb hasn't been called before we need to perform some setup actions. 
        if (.not.allocated(domain%fzs)) then
            call setup_linwinds(domain,options,rev,useD)
        endif
		
		! add the spatially variable linear field
		! if we are reverseing the effects, that means we are in the low-res domain
		! that domain does not have a spatial LUT calculated, so it can not be performed
		if (options%spatial_linear_fields)then
			call spatial_winds(domain,rev)
		else
			! Ndsq = squared Brunt Vaisalla frequency (1/s) typically from dry static stability
	        stability=calc_stability(domain)
			! This should probably be called twice, once for dry, and once or moist regions
			call linear_winds(domain,stability,vsmooth,reverse,useDensity,debug)
		endif
		debug=.False.
        
    end subroutine linear_perturb
end module linear_theory_winds