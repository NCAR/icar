!>----------------------------------------------------------
!!
!! Very simple microphysics code modeled after the microphysics
!! used in Smith and Barstad '04 (the linear model). 
!! 
!! Clouds (solid and liquid) form and evaporate instantly
!! Clouds convert to rain (or snow) with a time constant tau [s]
!!   typically rain_tau ~500s snow_tau ~2000s
!! Fall speeds in SB04 assume 10m/s (rain) and 1.5m/s (snow). 
!! In SB04 this is treated by means of a time constant (=height/speed). 
!! Here it is modeled explicitly with these fall speeds
!!
!! The entry point to the code is mp_simple_driver. 
!!
!! <pre>
!! Call tree graph :
!! mp_simple_driver->mp_simple->
!!   [mp_conversions->
!!       [cloud_conversion->sat_mr,
!!        cloud2hydrometeor,
!!        phase_change],
!!   sediment]
!! 
!! High level routine descriptions / purpose
!!   mp_simple_driver    - loops over X,Y grid cells, calls mp_simple on columns
!!   mp_simple           - calls mp_conversions for all z, then calls sediment
!!   mp_conversions      - handles all microphysics conversions (e.g. vapor->cloud->rain[->snow->vapor],...)
!!   cloud_conversion    - uses sat_mr to calculate sub or supersaturation and changes vapor and cloud water to balance 
!!                          (adjusts T for Latent Heating)
!!   cloud2hydrometeor   - converts cloud water to rain or snow (T dependant)
!!   phase_change        - handles rain and snow evaporation (can take a time constant)
!!   sediment            - advects falling rain and snow due to gravity (not wind)
!!
!! Driver inputs: p,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz
!!   p   = pressure                      - 3D - input  - Pa     - (nx,nz,ny)
!!   th  = potential temperature         - 3D - in/out - K      - (nx,nz,ny)
!!   pii = exner function                - 3D - input  - []     - (nx,nz,ny)
!!   rho = air density                   - 3D - input  - kg/m^3 - (nx,nz,ny)
!!   qv  = specific humidity             - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qc  = cloud water content           - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qr  = rain water content            - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qs  = snow water content            - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   rain = accumulated rain             - 2D - output - mm     - (nx,ny)
!!   snow = accumulated snow             - 2D - output - mm     - (nx,ny)
!!   dt = time step                      - 0D - input  - seconds    - scalar
!!   nx = number of ew grid cells        - 0D - input  - n      - scalar
!!   ny = number of ns grid cells        - 0D - input  - n      - scalar
!!   nz = number of vertical grid cells  - 0D - input  - n      - scalar
!! </pre>
!!
!!  Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_mp_simple
    implicit none
    private
    public::mp_simple_driver
    
    real, parameter :: LH_vapor=2.26E6 ! J/kg
    real, parameter :: dLHvdt = 2400! APPROXIMATE increase in latent heat with a decrease in temperature (below 373.15K)
                                    ! derived from various curves plotted online
                                    ! not sure why this isn't 4186.0 (J/kg/K) (i.e. specific heat of water)
                                    ! NOTE, ice=2110 (J/kg/K), steam=2080 (J/kg/K)
    real, parameter :: LH_liquid=3.34E5 ! J/kg
    real, parameter :: heat_capacity = 1006.0 ! air heat capacity J/kg/K
    real, parameter :: SMALL_VALUE = 1E-15
    real, parameter :: SMALL_PRESSURE = 1000. ! in Pa this is actually pretty small... 
!     real, parameter :: mp_R=287.058 ! J/(kg K) specific gas constant for air
!     real, parameter :: mp_g=9.81 ! gravity m/s^2

! arbitrary calibratable timescales default values as used in the linear model
! these should be pulled out to a parameter file for calibration purposes 
! but this is approximately how they are implemented in SB04
    real,parameter :: snow_formation_time_const=1/2000.0 ! [1/s]
    real,parameter :: rain_formation_time_const=1/500.0  ! [1/s]
    real,parameter :: freezing_threshold=273.15          ! [K]
    real,parameter :: snow_fall_rate=1.5                 ! [m/s]   for a water vapor scale height of 3750m corresponds to tau_f = 2500
    real,parameter :: rain_fall_rate=10.0                ! [m/s]   for a water vapor scale height of 3750m corresponds to tau_f = 375
    real,parameter :: snow_cloud_init=0.0001            ! [kg/kg] cloud ice content before snow will start to form
    real,parameter :: rain_cloud_init=0.0001            ! [kg/kg] cloud water content before rain will start to form
    
    
!   these are recalculated every call because they are a function of dt
!   conversion "time" = exp( -1 * time_constant * dt)
    real :: rain_evap=0.999
    real :: snow_evap=0.999
    real :: snow_melt=0.999
    real :: cloud2rain=0.999
    real :: cloud2snow=0.999
    !$omp threadprivate(cloud2rain,cloud2snow,snow_melt,snow_evap,rain_evap)
    
    contains
    real function sat_mr(t,p)
    ! Calculate the saturated mixing ratio at a temperature (K), pressure (Pa)
        implicit none
        real,intent(in) :: t,p
        real :: e_s,mr_s,a,b

!       from http://www.dtic.mil/dtic/tr/fulltext/u2/778316.pdf
!           Lowe, P.R. and J.M. Ficke., 1974: THE COMPUTATION OF SATURATION VAPOR PRESSURE 
!               Environmental Prediction Research Facility, Technical Paper No. 4-74
!       which references:
!           Murray, F. W., 1967: On the computation of saturation vapor pressure. 
!               Journal of Applied Meteorology, Vol. 6, pp. 203-204.
!       Also notes a 6th order polynomial and look up table as viable options. 
        if (t<freezing_threshold) then
            a=21.8745584
            b=7.66
        else
            a=17.2693882
            b=35.86
        endif
        e_s = 610.78* exp(a*(t-273.16)/(t-b)) !(Pa)

!       alternate formulations
!       Polynomial:
!       e_s = ao + t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t))))) a0-6 defined separately for water and ice
!       e_s = 611.2*exp(17.67*(t-273.15)/(t-29.65)) ! (Pa)
        !from : http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf
!       e_s = 611.0*10.0**(7.5*(t-273.15)/(t-35.45))
        
        
!       if ((p-e_s)<0) then
!           print*, "e_s more than p!"
!           print*, p,e_s,t,p-e_s
!       endif
!       e_s=min(e_s,p-SMALL_PRESSURE)
        if ((p-e_s)<=0) then
            e_s=p*0.99999
        endif
        !from : http://www.srh.noaa.gov/images/epz/wxcalc/mixingRatio.pdf
        sat_mr=0.6219907*e_s/(p-e_s) !(kg/kg)
    end function sat_mr
    

    subroutine cloud_conversion(p,t,qv,qc,qvsat,dt)
        implicit none
        real,intent(inout)::t,qv,qc,qvsat
        real,intent(in)::dt,p
        real :: vapor2temp,excess,deltat,pre_qc,pre_qv,pre_t,lastqv
        integer :: iteration
        real :: maxerr
        
        maxerr=1e-6
        iteration=0
        lastqv=qv+maxerr*2
        vapor2temp=(LH_vapor+(373.15-t)*dLHvdt)/heat_capacity
        pre_qc=qc !DEBUG
        pre_qv=qv !DEBUG
        pre_t=t !DEBUG
        excess=0
        
        do while ((abs(lastqv-qv)>maxerr).and.(iteration<5))
            iteration=iteration+1
            lastqv=qv
            ! calculate the saturating mixing ratio
            qvsat=sat_mr(t,p)
            ! if saturated create clouds
            if (qv>qvsat) then
                excess=(qv-qvsat)*0.95
                ! temperature change if all vapor is converted
                t=t+(excess*vapor2temp)
                qv=qv-excess
                qc=qc+excess
            ! if unsaturated and clouds exist, evaporate clouds
            else if (qc>0) then
                excess=(qvsat-qv)*0.95
                if (excess<qc) then
                    t=t-(excess*vapor2temp)
                    qv=qv+excess
                    qc=qc-excess
                else
                    qv=qv+qc
                    t=t-(qc*vapor2temp)
                    excess=qc
                    qc=0.
                endif
                excess=excess*(-1) !DEBUG
            endif
        enddo
        if (iteration==5) then
!           print*, iteration
!           print*, pre_qc,pre_qv,pre_t,p
!           print*, qc, qv,t, qvsat
!           print*, sat_mr(pre_t,p),sat_mr(t,p), sat_mr(pre_t,p)-sat_mr(t,p)
            qv=sat_mr(pre_t,p)
            t=pre_t
            qc=pre_qc
        endif
!       if (t<170) then
!           t=170
!       else if (t>350) then
!           t=350
!       endif
        
        qc=max(qc,0.)
        if ((t>350).or.(qc>0.01).or.(qvsat>1)) then
            deltat=excess*vapor2temp
            print*, "mp_simple: data out of bounds"
            print*, "iter=",iteration
            print*, "preqc=",pre_qc,"preqv=",pre_qv,"pret=",pre_t,"p=",p
            print*, "qc=",qc, "qv=",qv,"t=",t, "qvsat=",qvsat
            print*, "mrs_pret=",sat_mr(pre_t,p),"mrs_t=",sat_mr(t,p), "mrs_delta=",sat_mr(pre_t,p)-sat_mr(t,p)

            print*, "qv=",qv, "qvs=",qvsat, "qc=",qc, "preqc=",pre_qc, "t=",t, "excess=",excess, "vapor2temp=",vapor2temp
            print*, "t-deltat",t-deltat, "p=",p, "mrs_t-dT=",sat_mr(t-deltat,p)
        endif
        
    end subroutine


    subroutine cloud2hydrometeor(qc,q,conversion,qcmin)
        implicit none
        real,intent(inout) :: qc,q
        real,intent(in) :: conversion, qcmin
        real::delta
        
        if (qc > qcmin) then
            delta=qc-(qc*conversion)
        else
            delta=0
        endif
        
        if (delta<qc) then
            qc=qc-delta
            q=q+delta
        else
            q=q+qc
            qc=0.
        endif
        qc=max(qc,0.)
    end subroutine
    
    subroutine phase_change(p,t,q1,qmax,q2,Lheat,change_rate)
        ! convert from q1 (e.g. rain) to q2 (e.g. vapor)
        ! qmax = the maximum value for q2 (e.g. qv_sat)
        ! p = pressure (not used)
        ! t = temperature
        ! Lheat = latent heat associated with the phase change
        ! change_rate = fraction of q1 that can change in a timestep
        implicit none
        real, intent(inout)::t,q1,q2
        real,intent(in) :: p,qmax,change_rate,Lheat
        real :: mass2temp,delta
        mass2temp=Lheat/heat_capacity!*(p/(R*t)*dV))
        
        delta=(qmax-q2)*change_rate
        if (delta>q1) delta=q1
        ! hopefully we don't over shoot saturation (use a 1% buffer)
        if (delta>((qmax-q2)*0.99)) then
            delta=(qmax-q2)*0.99
        endif
        
        q1=q1-delta
        if (q1<0) then 
            if ((q1+SMALL_VALUE)<0) then
                q1=0
            else
                print*, "phase_change"
                print*, q1,q2,delta,qmax,change_rate
                stop
            endif
        endif
        q2=q2+delta
        t=t+delta*mass2temp
        
    end subroutine
    
    subroutine mp_conversions(p,t,qv,qc,qr,qs,dt)
        implicit none
        real, intent(inout) :: p,t,qv,qc,qr,qs
        real,intent(in)::dt
        real :: qvsat,L_evap,L_subl,L_melt
        
        qvsat=1
        L_melt=-1*LH_liquid  ! J/kg (should change with temperature)
        L_evap=-1*(LH_vapor+(373.15-t)*dLHvdt)   ! J/kg
        L_subl=L_melt+L_evap ! J/kg
        !convert cloud water to and from water vapor
        call cloud_conversion(p,t,qv,qc,qvsat,dt)
        ! if there are no species to process we will just return
        if ((qc+qr+qs) >SMALL_VALUE) then
            if (qc>SMALL_VALUE) then
                if (t>freezing_threshold) then
                    ! convert cloud water to rain drops
                    call cloud2hydrometeor(qc,qr,cloud2rain,rain_cloud_init)
!                   if (qs>SMALL_VALUE) then
!                       ! it is above freezing, so start melting any snow if present
!                       call phase_change(p,t,qs,100.,qr,L_melt,dt/snow_melt_time)
! !                         write(*,*) "Snow Melt",t
!                   endif
                else
                    ! convert cloud water to snow flakes
                    call cloud2hydrometeor(qc,qs,cloud2snow,snow_cloud_init)
                    
                endif
            endif
            ! if unsaturated, evaporate any existing snow and rain
            if (qv<qvsat) then
                if (qr>SMALL_VALUE) then
                    ! evaporate rain
                    call phase_change(p,t,qr,qvsat,qv,L_evap,cloud2rain)
                endif
                if (qs>SMALL_VALUE) then
                    ! sublimate snow
                    call phase_change(p,t,qs,qvsat,qv,L_subl,cloud2snow)
                endif
            endif
        endif
    end subroutine
    
    real function sediment(q,v,rho,dz,n)
        implicit none
        real,intent(inout),dimension(n)::q
        real,intent(in),dimension(n)::v,rho,dz
        integer,intent(in) :: n
        real,dimension(n) :: flux
        integer :: i
        
!       calculate the mass of material falling out of the bottom model level
        sediment=v(1)*q(1)*rho(1) ![m] * [kg/kg] * [kg/m^3] = [kg/m^2]
!       remove that from the bottom model layer. 
        q(1)=q(1)-(sediment/dz(1)/rho(1)) ! [kg/m^2] / [m] / [kg/m^3] = [kg/kg]
        do i=1,n-1
            flux(i)=v(i+1)*q(i+1)*rho(i+1)
        enddo
        do i=1,n-1
            q(i)=q(i)+flux(i)/(rho(i)*dz(i))
            q(i+1)=q(i+1)-flux(i)/(rho(i+1)*dz(i+1))
        enddo
    
    end function

    subroutine mp_simple(p,t,rho,qv,qc,qr,qs,rain,snow,dt,dz,nz,debug)
        implicit none
        real,intent(inout),dimension(nz)::p,t,rho,qv,qc,qr,qs
        real,intent(inout)::rain,snow
        real,intent(in),dimension(nz)::dz
        real,intent(in)::dt
        integer,intent(in)::nz
        logical,intent(in)::debug
        
        real,dimension(nz)::fall_rate
        real::cfl,snowfall, qvsat
        integer::i,cfl_step
        real :: L_evap,L_subl,L_melt
        
        qvsat=1
        L_melt=-1*LH_liquid  ! J/kg (should change with temperature)
        
        if ((debug).and.(dt<1)) print*, "internal dt=",dt
        do i=1,nz
            ! convert specific humidity to mixing ratio
!           qv(i)=qv(i)/(1-qv(i))
            call mp_conversions(p(i),t(i),qv(i),qc(i),qr(i),qs(i),dt)
        enddo

        ! SEDIMENTATION for rain
        if (maxval(qr)>SMALL_VALUE) then
            fall_rate=rain_fall_rate
            cfl=ceiling(maxval(dt/dz*fall_rate))
            fall_rate=dt*fall_rate/cfl
            ! substepping to satisfy CFL criteria
            do cfl_step=1,nint(cfl)
                rain=rain+sediment(qr,fall_rate,rho,dz,nz)
                ! allow any rain that reached a unsaturated layer to evaporate
                do i=1,nz
                    L_evap=-1*(LH_vapor+(373.15-t(i))*dLHvdt)   ! J/kg
                    qvsat=sat_mr(t(i),p(i))
                    if (qv(i)<qvsat) then
                        if (qr(i)>SMALL_VALUE) then
                            ! evaporate rain
                            call phase_change(p(i),t(i),qr(i),qvsat,qv(i),L_evap,cloud2rain)
                        endif
                    endif
                enddo
            enddo
        endif
        
        ! SEDIMENTATION for snow
        if (maxval(qs)>SMALL_VALUE) then
            fall_rate=snow_fall_rate
            cfl=ceiling(maxval(dt/dz*fall_rate))
            fall_rate=dt*fall_rate/cfl
            ! substepping to satisfy CFL criteria
            do cfl_step=1,nint(cfl)
                snowfall=sediment(qs,fall_rate,rho,dz,nz)
                snow=snow+snowfall
                rain=rain+snowfall
                ! allow any snow that reached a unsaturated layer to sublimate
                do i=1,nz
                    L_evap=-1*(LH_vapor+(373.15-t(i))*dLHvdt)   ! J/kg
                    L_subl=L_melt+L_evap ! J/kg
                    qvsat=sat_mr(t(i),p(i))
                    if (qv(i)<qvsat) then
                        if (qs(i)>SMALL_VALUE) then
                            ! sublimate snow
                            call phase_change(p(i),t(i),qs(i),qvsat,qv(i),L_subl,cloud2snow)
                        endif
                    endif
                enddo
            enddo
        endif
!         do i=1,nz
!           ! convert mixing ratio back to specific humidity
!           qv(i)=qv(i)/(qv(i)+1)
!       enddo
    end subroutine mp_simple


    subroutine mp_simple_driver(p,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz)
        implicit none
        real,intent(inout),dimension(nx,nz,ny)::p,th,pii,rho,qv,qc,qs,qr,dz
        real,intent(inout),dimension(nx,ny)::rain,snow
        real,intent(in)::dt
        integer,intent(in)::nx,ny,nz
        real,dimension(:),allocatable::t
        integer::i,j
        
!       calculate these once for every call because they are only a function of dt
        cloud2snow=exp(-1.0*snow_formation_time_const*dt)
        cloud2rain=exp(-1.0*rain_formation_time_const*dt)
        !$omp parallel private(i,j,t), &
        !$omp copyin(cloud2rain,cloud2snow,snow_melt,snow_evap,rain_evap),&
        !$omp shared(p,th,pii,qv,qc,qs,qr,rain,snow,dz),&
        !$omp firstprivate(dt,nx,ny,nz)
        allocate(t(nz))
        if (dt<1e-10) print*, dt ! for some reason this seems to fix an issue with ifort -O(anything) not copying dt in to the parallel region(?)
        !$omp do
        do j=2,ny-1
            do i=2,nx-1
                t=th(i,:,j)*pii(i,:,j)
                call mp_simple(p(i,:,j),t,rho(i,:,j),qv(i,:,j),&
                            qc(i,:,j),qr(i,:,j),qs(i,:,j),&
                            rain(i,j),snow(i,j),&
                            dt,dz(i,:,j),nz,((i==(nx/2+20)).and.(j==2)))
                th(i,:,j)=t/pii(i,:,j)
                
            enddo
!             where(qs(:,:,j)<0) qs(:,:,j)=0
!             where(qr(:,:,j)<0) qr(:,:,j)=0
!             where(qc(:,:,j)<0) qc(:,:,j)=0
!             where(qv(:,:,j)<0) qv(:,:,j)=0
        enddo
        !$omp end do
        deallocate(t)
        !$omp end parallel
    end subroutine mp_simple_driver
end module

!!!! Old Test Code
! program main
!   use module_mp_simple
!   real,dimension(5) ::p,t,qv,qc,qr,qs,dz
!   real::rain,snow,dt,dx2
!   integer::nz,i
!   nz=5
!   p(:)=80000.0
!   t(:)=270.0
!   qv(:)=0.005
!   qc(:)=0.
!   qr(:)=0.
!   qs(:)=0.
!   dz(:)=200.0
!   dx2=2000*2000.0
!   dt=20.0
!   rain=0
!   snow=0
!   cloud2snow=exp(-1.0*snow_formation_time_const*dt)
!   cloud2rain=exp(-1.0*rain_formation_time_const*dt)
!   
!   do i=1,10
!       call mp_simple(p,t,qv,qc,qr,qs,rain,snow,dt,dx2,dz,nz,.True.)
! !         write(*,*) "p=",p
!       write(*,*) t(2),qv(2),qc(2),qr(2),qs(2),rain,snow
! !         write(*,*) "t=",t
! !         write(*,*) "qv=",qv
! !         write(*,*) "qc=",qc
! !         write(*,*) "qr=",qr
! !         write(*,*) "qs=",qs
! !         write(*,*) "rain=",rain
! !         write(*,*) "snow=",snow
!   end do
! end program
    
