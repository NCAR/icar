module swim_step
    use advect
    use module_mp_thompson
	use module_mp_simple
    implicit none
	real, parameter :: debug=0
    real, parameter :: LH_vaporization=2260000.0 ! J/kg
    real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
    real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
    real, parameter :: g=9.81 ! gravity m/s^2
    contains
!   just a pass through routine to call the microphysics initialization (warning it takes a while)
    subroutine init(physics)
		integer,intent(in)::physics
		if (physics.EQ.1) then
			call thompson_init
		else
			write(*,*) "No Initialization necessary for mp_simple"
		endif
    end subroutine init

! add the latent and sensible heat fluxes to the temprature and humidity variables in the first model level
    subroutine land_surface_fluxes(th,qv,p,dz,pii, dt, &
                                    sensible_heat, latent_heat, pblh, &
                                    ids,ide,kds,kde,jds,jde)
        implicit none
!       input potential temperature (K) and humidity mixing ratio (kg/kg)
        real, dimension(ids:ide,kds:kde,jds:jde),intent(inout) :: th,qv
!       input pressure (Pa), layer thickness (m), and potential temperature conversion function
        real, dimension(ids:ide,kds:kde,jds:jde),intent(inout) :: p,dz,pii
!       input land surface fluxes (W/m^2)
        real, dimension(ids:ide,jds:jde),intent(inout) :: sensible_heat, latent_heat
		integer, dimension(ids:ide,jds:jde),intent(inout) :: pblh
!       input time step (seconds)
        real, intent(in) :: dt
!       array subscripts, shouldn't be necessary in fortran90, but something is getting tripped up
        integer, intent(in) :: ids,ide, jds,jde, kds,kde
!       internal loop counter
        integer :: i,j
!       temporary variables internal to the loop and private to each thread
		real, dimension(ids:ide,kds:kde) :: rho
        real, dimension(ids+1:ide-1) :: temperature,dTemp,lhdQV
		real,dimension(jds+1:jde-1) :: curmax
		real, dimension(kds:kde)::means,f1,f2,rhomean
!       not really needed anymore variables for tracking statistics and modifying fluxes
        real :: factor,mxdT,mxdqv,temp,diffusionrate,coolingrate,sbconstant
        temp=0.
        mxdT=0.
        mxdqv=0.
        factor=1.0 ! factor to scale land surface fluxes for tuning purposes
		! diffusionrate*dq/dz must be < 1 (or even 0.5?) for model stability
		diffusionrate=dt/(1.0*60.0*60.0)*20 ! *20 = rate "constant" this could/should be based on stability?
		sbconstant=5.67e-8
		coolingrate=1.5*(dt/(60.0*60.0*24.0)) *sbconstant / 300.0 !1.5K/day radiative cooling rate (300 = W/m^2 at 270K)
		if (debug.eq.1) then
			write(*,*) "in land surface fluxes"
		endif
!       OpenMP commands
        !$omp parallel firstprivate(jds,jde,ids,ide,kds,kde,diffusionrate,coolingrate) &
		!$omp private(i,temperature,dTemp,lhdQV,rho,rhomean,means,f1,f2) &
        !$omp shared(pii,p,th,qv,dz,sensible_heat,latent_heat,pblh,curmax)
        !$omp do
        do i=jds+1,jde-1
! 			write(*,*) 0
!           loop over the last (slowest) dimension with the first (fastest) dimension processed internally
!           compute the real temperature from the potential temperature
            temperature=th(ids+1:ide-1,1,i)*pii(ids+1:ide-1,1,i) ! K
! 			write(*,*) 1
!           compute air density
            rho(ids+1:ide-1,kds:kde)=p(ids+1:ide-1,kds:kde,i)/(R*th(ids+1:ide-1,kds:kde,i)*pii(ids+1:ide-1,kds:kde,i)) ! kg/m^3
! 			write(*,*) 2
!           convert sensible heat flux to a temperature delta term
            ! J/s/m^2 * s / J/(kg*K) => kg*K/m^2 ... /((kg/m^3) * m) => K
            dTemp=(sensible_heat(ids+1:ide-1,i)*dt/cp) / (rho(ids+1:ide-1,1)*dz(ids+1:ide-1,1,i)) / factor! K
! 			write(*,*) 3
! 			curmax(i)=MAXVAL(dTemp)
!           add temperature delta and convert back to potential temperature
			temperature=temperature+dTemp
            th(ids+1:ide-1,1,i)=temperature/pii(ids+1:ide-1,1,i)
! 			write(*,*) 4
!           convert latent heat flux to a mixing ratio tendancy term
            ! J/s/m^2 * s / J/kg => kg/m^2 ... / (kg/m^3 * m) => kg/kg
            lhdQV=(latent_heat(ids+1:ide-1,i)/LH_vaporization*dt) / (rho(ids+1:ide-1,1)*dz(ids+1:ide-1,1,i)) / factor
! 			write(*,*) i
!           and add back to the mixing ratio
            qv(ids+1:ide-1,1,i)=qv(ids+1:ide-1,1,i)+lhdQV

!			stupid diffusion within the PBL	I should at least speed up/clean up the math at some point...
! 			ideally could use an implicit solution?
			do j=ids+1,ide-1
! 				write(*,*) i,j
				rhomean(kds:kde-1)=(rho(j,kds:kde-1)+rho(j,kds+1:kde))/2
				rhomean(kde)=rho(j,kde)
				
				means(kds:kde-1)=(qv(j,kds:kde-1,i)+qv(j,kds+1:kde,i))/2
! 				diffusion fluxes within the PBL
				f1(1:pblh(j,i)-1)=(qv(j,1:pblh(j,i)-1,i)-means(1:pblh(j,i)-1))*diffusionrate
				f2(1:pblh(j,i)-1)=(qv(j,2:pblh(j,i),i)-means(1:pblh(j,i)-1))*diffusionrate
! 				diffusion fluxes above the PBL
				f1(pblh(j,i):kde-1)=(qv(j,pblh(j,i):kde-1,i)-means(pblh(j,i):kde-1))*diffusionrate/5
				f2(pblh(j,i):kde-1)=(qv(j,pblh(j,i)+1:kde,i)-means(pblh(j,i):kde-1))*diffusionrate/5
! 				add fluxes to grid cells
				qv(j,kds:kde-1,i)=qv(j,kds:kde-1,i)-f1*rhomean(kds:kde-1)/rho(j,kds:kde-1)
				qv(j,kds+1:kde,i)=qv(j,kds+1:kde,i)-f2*rhomean(kds+1:kde)/rho(j,kds+1:kde)

! 				same process for potential temperature
				means(kds:kde-1)=(th(j,kds:kde-1,i)+th(j,kds+1:kde,i))/2
				f1(1:pblh(j,i)-1)=(th(j,1:pblh(j,i)-1,i)-means(1:pblh(j,i)-1))*diffusionrate
				f2(1:pblh(j,i)-1)=(th(j,2:pblh(j,i),i)-means(1:pblh(j,i)-1))*diffusionrate
				f1(pblh(j,i):kde-1)=(th(j,pblh(j,i):kde-1,i)-means(pblh(j,i):kde-1))*diffusionrate/5
				f2(pblh(j,i):kde-1)=(th(j,pblh(j,i)+1:kde,i)-means(pblh(j,i):kde-1))*diffusionrate/5
				th(j,kds:kde-1,i)=th(j,kds:kde-1,i)-f1*rhomean(kds:kde-1)/rho(j,kds:kde-1)
				th(j,kds+1:kde,i)=th(j,kds+1:kde,i)-f2*rhomean(kds+1:kde)/rho(j,kds+1:kde)
			enddo
! 			stupid radiative cooling th=th-coolingrate*(T^4)
			th(ids+1:ide-1,kds:kde,i)=th(ids+1:ide-1,kds:kde,i)-(((th(ids+1:ide-1,kds:kde,i)*pii(ids+1:ide-1,kds:kde,i))**4)*coolingrate)
        enddo
        !$omp end do
        !$omp end parallel
! 		write(*,*) MAXVAL(curmax),MAXVAL(pblh)
        
    end subroutine land_surface_fluxes

    SUBROUTINE timestep(nsteps,u,v,w,sensible_heat, latent_heat, pblh, &
                            qv, qc, qr, qi, qs, qg, ni, nr, &
                            th, pii, p, dz, dt_in, itimestep, &
                            dth,dqv,du,dv,dw,dp, &
                            RAINNC, RAINNCV, &
                            SNOWNC, SNOWNCV, &
                            GRAUPELNC, GRAUPELNCV, &
                            SR, physics,&
                            ids,ide, jds,jde, kds,kde, &             ! domain dims
                            ims,ime, jms,jme, kms,kme, &             ! memory dims
                            its,ite, jts,jte, kts,kte)               ! tile dims
        implicit none
        integer,intent(in) :: nsteps
        real,dimension(ids:ide,kds:kde,jds:jde-1),intent(inout) :: u,du
        real,dimension(ids:ide-1,kds:kde,jds:jde),intent(inout) :: v,dv
        real,dimension(ids:ide,jds:jde),intent(inout) :: sensible_heat, latent_heat
		integer, dimension(ids:ide,jds:jde),intent(inout) :: pblh
        real,dimension(ids:ide,kds:kde,jds:jde),intent(inout) :: w,pii,p,dz
        real,dimension(ids:ide,kds:kde,jds:jde),intent(in) :: dqv,dth,dw,dp
        real,dimension(ids:ide,kds:kde,jds:jde),intent(inout) :: qv,qc,qr,qi,qs,qg,ni,nr,th
        real,dimension(ids:ide,jds:jde),intent(inout) :: RAINNC,RAINNCV,SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV,SR
        integer, intent(in) :: ids,ide, jds,jde, kds,kde
        integer, intent(in) :: ims,ime, jms,jme, kms,kme
        integer, intent(in) :: its,ite, jts,jte, kts,kte
        real, intent(in) :: dt_in
        integer, intent(in) :: itimestep,physics
        real,dimension(ids:ide,kds:kde,jds:jde) :: temperature,rho!,dvol
        real ::dt_m
        integer :: t,i
        integer :: adv_steps,nx,ny,nz
		real ::dx2

		dx2=4000.0*4000.0
        ny=ide
        nz=kde
        nx=jde
        if (dt_in.LT.30.0) then
            adv_steps=floor(30/dt_in)
            dt_m=dt_in*adv_steps
        else
            adv_steps=1
            dt_m=dt_in
        endif
! 		dt_m=dt_in
! 		adv_steps=1
		if (debug.eq.1) then
			write(*,*) "in swim_step"
		endif
!         !$omp parallel private(i) firstprivate(jds,jde,dx2) & 
!         !$omp shared(dvol,dz)
!         !$omp do 
! 		do i=jds,jde
! 			dvol(:,:,i)=dx2*dz(:,:,i)
! 		enddo
!         !$omp end do
!         !$omp end parallel
		
        do t=1,nsteps
			
            if (mod(t,adv_steps).EQ.0) then
!               call the microphysics routine every adv_steps time steps
				if (physics.EQ.1) then
	                call mp_gt_driver(qv, qc, qr, qi, qs, qg, ni, nr, &
	                                    th, pii, p, dz, dt_m, itimestep, &
	                                    RAINNC, RAINNCV, &
	                                    SNOWNC, SNOWNCV, &
	                                    GRAUPELNC, GRAUPELNCV, &
	                                    SR, &
	                                    ids,ide, jds,jde, kds,kde, &             ! domain dims
	                                    ims,ime, jms,jme, kms,kme, &             ! memory dims
	                                    its,ite, jts,jte, kts,kte)               ! tile dims
				else
	                call mp_simple_driver(p,th,pii,qv, qc, qr, qs,RAINNC,SNOWNC, &
	                                    dt_m,dz, jde,ide,kde)
				endif				
            endif
			if (debug.eq.1) then
				write(*,*) "done microphysics"
			endif

            !$omp parallel private(i) firstprivate(jds,jde) & 
            !$omp shared(qv,qc,qr,qi,qs,qg,th,p,rho,temperature)
            !$omp do 
            do i=jds,jde
    !           compute the real temperature from the potential temperature
                temperature(:,:,i)=th(:,:,i)*pii(:,:,i) ! K should convert to J? =T *rho*cp*dvol??
    !           compute air density
                rho(:,:,i)=p(:,:,i)/(R*temperature(:,:,i))
! 				th(:,:,i)=th(:,:,i)*rho(:,:,i)
                qv(:,:,i)=qv(:,:,i)*rho(:,:,i)
                qc(:,:,i)=qc(:,:,i)*rho(:,:,i)
                qr(:,:,i)=qr(:,:,i)*rho(:,:,i)
                qs(:,:,i)=qs(:,:,i)*rho(:,:,i)
				if (physics.EQ.1) then
	                qi(:,:,i)=qi(:,:,i)*rho(:,:,i)
	                qg(:,:,i)=qg(:,:,i)*rho(:,:,i)
				endif
            enddo
            !$omp end do
            !$omp end parallel
! Run advection on each scalar independantly 
! (could all be done in their own thread? but advection is parallelized internally)
! 			write(*,*) "------------------------"
! 			write(*,*) th(100,13,100:104)
! 			write(*,*) th(101,13,100:104)
! 			write(*,*) th(102,13,100:104)
! 			write(*,*) th(103,13,100:104)
! 			write(*,*) "------------------------"
! 			write(*,*) th(100,12,100:104)
! 			write(*,*) th(101,12,100:104)
! 			write(*,*) th(102,12,100:104)
! 			write(*,*) th(103,12,100:104)
! 			write(*,*) "------------------------"
! 			write(*,*) th(100,11,100:104)
! 			write(*,*) th(101,11,100:104)
! 			write(*,*) th(102,11,100:104)
! 			write(*,*) th(103,11,100:104)
			
            call advect3d(qv,u,v,w,ny,nz,nx,0)
            call advect3d(th,u,v,w,ny,nz,nx,0)
            call advect3d(qc,u,v,w,ny,nz,nx,0)
            call advect3d(qr,u,v,w,ny,nz,nx,0)
            call advect3d(qs,u,v,w,ny,nz,nx,0)
			if (physics.EQ.1) then
	            call advect3d(qi,u,v,w,ny,nz,nx,0)
	            call advect3d(qg,u,v,w,ny,nz,nx,0)
	            call advect3d(nr,u,v,w,ny,nz,nx,0)
	            call advect3d(ni,u,v,w,ny,nz,nx,0)
			endif
			if (debug.eq.1) then
				write(*,*) "done advection"
			endif
! 			write(*,*) "----------q_post--------------"
! 			write(*,*) th(100,12,100:104)
! 			write(*,*) th(101,12,100:104)
! 			write(*,*) th(102,12,100:104)
! 			write(*,*) th(103,12,100:104)
! 			write(*,*) "------------------------"
! 			write(*,*) "----------U--------------"
! 			write(*,*) u(100,12,99:104)/dt_in*4000.0
! 			write(*,*) u(101,12,99:104)/dt_in*4000.0
! 			write(*,*) u(102,12,99:104)/dt_in*4000.0
! 			write(*,*) u(103,12,99:104)/dt_in*4000.0
! 			write(*,*) "----------V--------------"
! 			write(*,*) v(99,12,100:104)/dt_in*4000.0
! 			write(*,*) v(100,12,100:104)/dt_in*4000.0
! 			write(*,*) v(101,12,100:104)/dt_in*4000.0
! 			write(*,*) v(102,12,100:104)/dt_in*4000.0
! 			write(*,*) v(103,12,100:104)/dt_in*4000.0
! 			write(*,*) "----------W15--------------"
! 			write(*,*) w(100,12,100:104)/dt_in*4000.0
! 			write(*,*) w(101,12,100:104)/dt_in*4000.0
! 			write(*,*) w(102,12,100:104)/dt_in*4000.0
! 			write(*,*) w(103,12,100:104)/dt_in*4000.0
! 			write(*,*) "----------W14--------------"
! 			write(*,*) w(100,11,100:104)/dt_in*4000.0
! 			write(*,*) w(101,11,100:104)/dt_in*4000.0
! 			write(*,*) w(102,11,100:104)/dt_in*4000.0
! 			write(*,*) w(103,11,100:104)/dt_in*4000.0
! 			write(*,*) "------------------------"
! 			write(*,*) "------------------------"
			
            !$omp parallel firstprivate(jds,jde) private(i) &
            !$omp shared(pii,th,qv,qc,qr,qi,qs,qg,rho,p,u,v,w,dth,dqv,dp,du,dv,dw)
            !$omp do
            do i=jds,jde
! 				th(:,:,i)=th(:,:,i)/rho(:,:,i)
                qv(:,:,i)=qv(:,:,i)/rho(:,:,i)
                qc(:,:,i)=qc(:,:,i)/rho(:,:,i)
                qr(:,:,i)=qr(:,:,i)/rho(:,:,i)
                qs(:,:,i)=qs(:,:,i)/rho(:,:,i)
				if (physics.EQ.1) then
					qi(:,:,i)=qi(:,:,i)/rho(:,:,i)
	                qg(:,:,i)=qg(:,:,i)/rho(:,:,i)
				endif
!               add large scale time tendency terms to the u,v,w and recompute pii with the new pressure
                p(:,:,i)=p(:,:,i)+dp(:,:,i)
                pii(:,:,i)=1.0/((100000.0/p(:,:,i))**(R/cp))
                th(:,:,i)=th(:,:,i)+dth(:,:,i)
                qv(:,:,i)=qv(:,:,i)+dqv(:,:,i)
                v(:,:,i)=v(:,:,i)+dv(:,:,i)
                w(:,:,i)=w(:,:,i)+dw(:,:,i)
                if (i.LT.jde) then
                    u(:,:,i)=u(:,:,i)+du(:,:,i)
                endif
            enddo
            !$omp end do
            !$omp end parallel
! 			if (debug.eq.1) then
! 				write(*,*) "calling land surface"
! 			endif
! 			
!             call land_surface_fluxes(th,qv,p,dz,pii,dt_in,&
!                                     sensible_heat, latent_heat, pblh, &
!                                     ids,ide,kds,kde,jds,jde)
        enddo
    end subroutine timestep
end module swim_step
