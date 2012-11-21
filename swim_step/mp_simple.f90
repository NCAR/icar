module module_mp_simple
	implicit none
    real, parameter :: LH_vapor=2260000.0 ! J/kg
	real, parameter :: LH_liquid=334000.0 ! J/kg
    real, parameter :: heat_capacity = 1021.0 ! air heat capacity J/kg/K
!     real, parameter :: mp_R=287.058 ! J/(kg K) specific gas constant for air
!     real, parameter :: mp_g=9.81 ! gravity m/s^2
	!           compute air density
! 	            rho=p/(R*temperature) ! kg/m^3
	!           convert sensible heat flux to a temperature delta term
	            ! J/s/m^2 * s / J/(kg*K) => kg*K/m^2 ... /((kg/m^3) * m) => K
! 	            dTemp=(heat_flux*dt/cp) / (rho*dz) ! K
	
	contains
	real function sat_mr(t,p)
	! Calculate the saturated mixing ratio at a temperature (K), pressure (Pa)
		implicit none
		real,intent(in) :: t,p
		real :: e_s,mr_s
		
		! e_s = 6.112*exp(17.67*(t-273.15)/(t-29.65))
		!from : http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf
		e_s = 611.0*10.0**(7.5*(t-273.15)/(t-35.45))
		!from : http://www.srh.noaa.gov/images/epz/wxcalc/mixingRatio.pdf
		sat_mr=0.62197*e_s/(p-e_s) !(kg/kg)
	end function sat_mr
	
	subroutine cloud_conversion(p,t,qv,qc,qvsat,dt)
		implicit none
		real,intent(inout)::t,qv,qc,qvsat
		real,intent(in)::dt,p
		real :: vapor2temp,excess,deltat
		vapor2temp=LH_vapor/heat_capacity!*(p/(R*t))*dV)
		! write(*,*) LH_vapor/heat_capacity!*dV,p/(R*t),dV

		qvsat=sat_mr(t,p)
! 		write(*,*) "qvsat,qv",qvsat,qv
		if (qv>qvsat) then
			excess=qv-qvsat
			deltat=excess*vapor2temp
! 			write(*,*) "vapor-excess,dTemp",excess,excess*vapor2temp
			qvsat=sat_mr(t+deltat*0.25,p)
			excess=qv-qvsat
			t=t+(excess*vapor2temp)
			qv=qv-excess
			qc=qc+excess
! 			write(*,*) "vapor-excess,dTemp",excess,excess*vapor2temp
		else if (qc>1E-9) then
			excess=qvsat-qv
			if (excess<qc) then
				deltat=excess*vapor2temp
				qvsat=sat_mr(t-deltat*0.25,p)
				excess=qvsat-qv
				if (excess<qc) then
					t=t-(excess*vapor2temp)
					qv=qv+excess
					qc=qc-excess
! 					write(*,*) "cloud-excess,dTemp",excess,excess*vapor2temp
				else
					qv=qv+qc
					t=t-(qc*vapor2temp)
					qc=1E-15
				endif
			else
				qv=qv+qc
				t=t-(qc*vapor2temp)
				qc=1E-15
			endif
		endif
			
	end subroutine

	subroutine cloud2hydrometeor(qc,q,conversion)
		implicit none
		real,intent(inout) :: qc,q
		real,intent(in) :: conversion
		real::delta
	
		delta=qc*conversion
		if (delta<qc) then
			qc=qc-delta
			q=q+delta
		else
			q=q+qc
			qc=1E-15
		endif
	end subroutine
	
	subroutine phase_change(p,t,q1,qmax,q2,Lheat,change_rate)
		implicit none
		real, intent(inout)::t,q1,q2
		real,intent(in) :: p,qmax,change_rate,Lheat
		real :: mass2temp,delta
		mass2temp=Lheat/heat_capacity!*(p/(R*t)*dV))
		
		if (q1>3E-5) then
			delta=q1*change_rate
			!make sure we don't over saturate the air
			if (delta>(qmax-q2)) then
				delta=qmax-q2
			endif
			q1=q1-delta
			q2=q2+delta
		else
			delta=q1
			if (delta>(qmax-q2)) then
				delta=qmax-q2
			endif
			q2=q2+delta
			q1=1E-15
		endif
		t=t+delta*mass2temp
	
	end subroutine
	
	subroutine mp_conversions(p,t,qv,qc,qr,qs,dt)
		implicit none
		real, intent(inout) :: p,t,qv,qc,qr,qs
		real,intent(in)::dt
		real :: qvsat,rain_evap_time,snow_evap_time,cloud2rain_time,&
				cloud2snow_time,snow_melt_time,&
				L_evap,L_subl,L_melt
		
		L_melt=LH_liquid !kJ/kg
		L_evap=LH_vapor !kJ/kg
		L_subl=L_melt+L_evap
		!arbitrary calibratable timescales
		rain_evap_time=100.0 !seconds
		snow_evap_time=200.0 !seconds
		snow_melt_time=100.0 !seconds
		cloud2rain_time=500.0!seconds
		cloud2snow_time=500.0!seconds
		
		!convert cloud water to and from water vapor
		call cloud_conversion(p,t,qv,qc,qvsat,dt)
		! if there are no species to process we will just return
		if ((qc+qr+qs) >1E-9) then
			if (qc>1E-9) then
				if (t>273.15) then
					! convert cloud water to rain drops
					call cloud2hydrometeor(qc,qr,dt/cloud2rain_time)
! 					if (qs>1E-9) then
! 						! it is above freezing, so start melting any snow if present
! 						call phase_change(p,t,qs,100.,qr,L_melt,dt/snow_melt_time)
! ! 						write(*,*) "Snow Melt",t
! 					endif
				else
					! convert cloud water to snow flakes
					call cloud2hydrometeor(qc,qs,dt/cloud2snow_time)
					
				endif
			endif
			! if unsaturated, evaporate any existing snow and rain
! 			if (qv<qvsat) then
! 				if (qr>1E-9) then
! 					! evaporate rain
! 					call phase_change(p,t,qr,qvsat,qv,L_evap,dt/rain_evap_time)
! ! 					write(*,*) "evap rain"
! 				endif
! 				if (qs>1E-9) then
! 					! sublimate snow
! 					call phase_change(p,t,qs,qvsat,qv,L_subl,dt/snow_evap_time)
! ! 					write(*,*) "evap snow"
! 				endif
! 			endif
		endif
	end subroutine
	
	real function sediment(q,v,n)
		implicit none
		real,intent(inout),dimension(n)::q
		real,intent(in),dimension(n)::v
		integer,intent(in) :: n
		real :: flux
		integer :: i
		
		sediment=v(1)*q(1)
		if (sediment>q(1)) then
			sediment=q(1)
			q(1)=1E-15
		else
			q(1)=q(1)-sediment
		endif
		
		do i=1,n-1
			flux=v(i+1)*q(i+1)
			if (flux>q(i+1)) then
				! note it would be better to force v<=1 by sub_stepping in time
				! v should be <1 (courant condition) but it is not currently enforced 
				! that ~speed*dt < min(dz)
				q(i)=q(i)+q(i+1)
				q(i+1)=1E-15
			else
				q(i)=q(i)+flux
				q(i+1)=q(i+1)-flux
			endif
		enddo
	
	end function

	subroutine mp_simple(p,t,qv,qc,qr,qs,rain,snow,dt,dz,nz)
		implicit none
		real,intent(inout),dimension(nz)::p,t,qv,qc,qr,qs
		real,intent(inout)::rain,snow
		real,intent(in),dimension(nz)::dz
		real,intent(in)::dt
		integer,intent(in)::nz
		real::fall_rate
		integer::i
		
! 		fall_rate=2+(t-273.15)/5
		fall_rate=1.5
		do i=1,nz
			call mp_conversions(p(i),t(i),qv(i),qc(i),qr(i),qs(i),dt)
! 			if (fall_rate(i)<1) then
! 				fall_rate(i)=1.0 !m/s
! 			else if (fall_rate(i)>10) then
! 				fall_rate(i)=10.0 !m/s
! 			endif
		enddo
		rain=rain+sediment(qr,dt/dz*10.,nz)*dz(1)
		snow=snow+sediment(qs,dt/dz*fall_rate,nz)*dz(1)
		
	end subroutine mp_simple


! call mp_simple_driver(p,th,pii,qv, qc, qr, qs,RAINNC,SNOWNC, &
!                     dt_m,dz, jde,ide,kde)
	subroutine mp_simple_driver(p,th,pii,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz)
		implicit none
		real,intent(inout),dimension(ny,nz,nx)::p,th,pii,qv,qc,qs,qr,dz
		real,intent(inout),dimension(ny,nx)::rain,snow
		real,intent(in)::dt
		integer,intent(in)::nx,ny,nz
		real,dimension(ny,nz,nx)::t
		integer::i,j
		!$omp parallel private(i,j),&
		!$omp shared(p,t,th,pii,qv,qc,qs,qr,rain,snow,dz),&
		!$omp firstprivate(dt,nx,ny,nz)
		!$omp do
                do j=2,nx-1
		        do i=2,ny-1
				t(i,:,j)=th(i,:,j)*pii(i,:,j)
				call mp_simple(p(i,:,j),t(i,:,j),qv(i,:,j),&
							qc(i,:,j),qr(i,:,j),qs(i,:,j),&
							rain(i,j),snow(i,j),&
							dt,dz(i,:,j),nz)
				th(i,:,j)=t(i,:,j)/pii(i,:,j)
			enddo
		enddo
		!$omp end do
		!$omp end parallel
	end subroutine mp_simple_driver
end module

!!!! "unit" Test Code (sort of)
! program main
! 	use module_mp_simple
! 	real,dimension(5) ::p,t,qv,qc,qr,qs,dz
! 	real::rain,snow,dt,dx2
! 	integer::nz,i
! 	nz=5
! 	p(:)=80000.0
! 	t(:)=270.0
! 	qv(:)=0.005
! 	qc(:)=0.
! 	qr(:)=0.
! 	qs(:)=0.
! 	dz(:)=200.0
! 	dx2=2000*2000.0
! 	dt=20.0
! 	rain=0
! 	snow=0
! 	
! 	do i=1,10
! 		call mp_simple(p,t,qv,qc,qr,qs,rain,snow,dt,dx2,dz,nz)
! ! 		write(*,*) "p=",p
! 		write(*,*) t(2),qv(2),qc(2),qr(2),qs(2),rain,snow
! ! 		write(*,*) "t=",t
! ! 		write(*,*) "qv=",qv
! ! 		write(*,*) "qc=",qc
! ! 		write(*,*) "qr=",qr
! ! 		write(*,*) "qs=",qs
! ! 		write(*,*) "rain=",rain
! ! 		write(*,*) "snow=",snow
! 	end do
! end program
	
