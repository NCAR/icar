!>----------------------------------------------------------
!! Very simple microphysics code modeled after the microphysics of SB04
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
!!                          (adjusts temperature for Latent Heating)
!!   cloud2hydrometeor   - converts cloud water to rain or snow (temperature dependant)
!!   phase_change        - handles rain and snow evaporation (can take a time constant)
!!   sediment            - advects falling rain and snow due to gravity (not wind)
!!   sat_mr              - calculate saturated water vapor mixing ratio
!!
!! Driver inputs: pressure,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz
!!   pressure   = pressure                      - 3D - input  - Pa     - (nx,nz,ny)
!!   th         = potential temperature         - 3D - in/out - K      - (nx,nz,ny)
!!   pii        = exner function                - 3D - input  - []     - (nx,nz,ny)
!!   rho        = air density                   - 3D - input  - kg/m^3 - (nx,nz,ny)
!!   qv         = specific humidity             - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qc         = cloud water content           - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qr         = rain water content            - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qs         = snow water content            - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   rain       = accumulated rain              - 2D - output - mm     - (nx,ny)
!!   snow       = accumulated snow              - 2D - output - mm     - (nx,ny)
!!   dt         = time step                     - 0D - input  - sec.   - scalar
!!   nx         = number of ew grid cells       - 0D - input  - n      - scalar
!!   ny         = number of ns grid cells       - 0D - input  - n      - scalar
!!   nz         = number of vertical grid cells - 0D - input  - n      - scalar
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_mp_simple
    use options_interface, only : options_t
    use data_structures

    implicit none
    private
    public :: mp_simple_driver, mp_simple_var_request

    ! parameters should come from atm_utilties or icar_constants modules
    real, parameter :: LH_vapor = 2.26E6 ! J/kg
    real, parameter :: dLHvdt   = 2400! APPROXIMATE increase in latent heat with a decrease in temperature (below 373.15K)
                                    ! derived from various curves plotted online
                                    ! not sure why this isn't 4186.0 (J/kg/K) (i.e. specific heat of water)
                                    ! NOTE, ice = 2110 (J/kg/K), steam = 2080 (J/kg/K)
    real, parameter :: LH_liquid = 3.34E5 ! J/kg
    real, parameter :: heat_capacity = 1006.0 ! air heat capacity J/kg/K
    real, parameter :: SMALL_VALUE = 1E-30
    real, parameter :: SMALL_PRESSURE = 1000. ! in Pa this is actually pretty small...
!     real, parameter :: mp_R = 287.058 ! J/(kg K) specific gas constant for air
!     real, parameter :: mp_g = 9.81 ! gravity m/s^2

! arbitrary calibratable timescales default values as used in the linear model
! these should be pulled out to a parameter file for calibration purposes
! but this is approximately how they are implemented in SB04
    real, parameter :: snow_evap_time_const      = 1 / 2000.0   ! [1/s]
    real, parameter :: rain_evap_time_const      = 1 / 500.0    ! [1/s]
    real, parameter :: snow_formation_time_const = 1 / 2000.0   ! [1/s]
    real, parameter :: rain_formation_time_const = 1 / 500.0    ! [1/s]
    real, parameter :: freezing_threshold        = 273.15       ! [K]
    real, parameter :: snow_fall_rate            = 1.5          ! [m/s]   for a water vapor scale height of 3750m corresponds to tau_f = 2500
    real, parameter :: rain_fall_rate            = 10.0         ! [m/s]   for a water vapor scale height of 3750m corresponds to tau_f = 375
    real, parameter :: snow_cloud_init           = 0.0001       ! [kg/kg] cloud ice content before snow will start to form
    real, parameter :: rain_cloud_init           = 0.0001       ! [kg/kg] cloud water content before rain will start to form

!   these are recalculated every call because they are a function of dt
!   conversion "time" = exp( -1 * time_constant * dt)
    real :: rain_evap  = 0.999
    real :: snow_evap  = 0.999
    real :: snow_melt  = 0.999
    real :: cloud2rain = 0.999
    real :: cloud2snow = 0.999
    !$omp threadprivate(cloud2rain,cloud2snow,snow_melt,snow_evap,rain_evap)

contains

    !> ----------------------------------------------
    !! Communicate to the master process requesting the variables requred to be allocated, used for restart files, and advected
    !!
    !! ----------------------------------------------
    subroutine mp_simple_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options


        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density,      &
                      kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,  kVARS%snow_in_air,  &
                      kVARS%precipitation, kVARS%snowfall,              kVARS%dz])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water,   &
                       kVARS%rain_in_air,           kVARS%snow_in_air   ] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature,    kVARS%water_vapor,  &
                        kVARS%cloud_water,  kVARS%rain_in_air,              kVARS%snow_in_air,  &
                        kVARS%precipitation,kVARS%snowfall,                 kVARS%dz] )

    end subroutine mp_simple_var_request




    !>----------------------------------------------------------
    !!  Calculate the saturated mixing ratio for a given temperature and pressure
    !!
    !!  If temperature > 0C: returns the saturated mixing ratio with respect to liquid
    !!  If temperature < 0C: returns the saturated mixing ratio with respect to ice
    !!
    !!  @param temperature  Air Temperature [K]
    !!  @param pressure     Air Pressure [Pa]
    !!  @retval sat_mr      Saturated water vapor mixing ratio [kg/kg]
    !!
    !!  @see http://www.dtic.mil/dtic/tr/fulltext/u2/778316.pdf
    !!   Lowe, P.R. and J.M. Ficke., 1974: The Computation of Saturation Vapor Pressure
    !!   Environmental Prediction Research Facility, Technical Paper No. 4-74
    !!
    !!----------------------------------------------------------
    real function sat_mr(temperature, pressure)
    ! Calculate the saturated mixing ratio at a temperature (K), pressure (Pa)
        implicit none
        real,intent(in) :: temperature,pressure
        real :: e_s,mr_s,a,b

        ! from http://www.dtic.mil/dtic/tr/fulltext/u2/778316.pdf
        !   Lowe, P.R. and J.M. Ficke., 1974: THE COMPUTATION OF SATURATION VAPOR PRESSURE
        !       Environmental Prediction Research Facility, Technical Paper No. 4-74
        ! which references:
        !   Murray, F. W., 1967: On the computation of saturation vapor pressure.
        !       Journal of Applied Meteorology, Vol. 6, pp. 203-204.
        ! Also notes a 6th order polynomial and look up table as viable options.
        if (temperature < freezing_threshold) then
            a = 21.8745584
            b = 7.66
        else
            a = 17.2693882
            b = 35.86
        endif

        e_s = 610.78 * exp(a * (temperature - 273.16) / (temperature - b)) !(Pa)

        ! alternate formulations
        ! Polynomial:
        ! e_s = ao + t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t))))) a0-6 defined separately for water and ice
        ! e_s = 611.2*exp(17.67*(t-273.15)/(t-29.65)) ! (Pa)
        ! from : http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf
        ! e_s = 611.0*10.0**(7.5*(t-273.15)/(t-35.45))


        if ((pressure - e_s) <= 0) then
            e_s = pressure * 0.99999
        endif
        ! from : http://www.srh.noaa.gov/images/epz/wxcalc/mixingRatio.pdf
        sat_mr = 0.6219907 * e_s / (pressure - e_s) !(kg/kg)
    end function sat_mr


    !>----------------------------------------------------------
    !!  Convert cloud water to vapor and back
    !!
    !!  Iterates to find the conversion between liquid and vapor including temperature feedbacks
    !!
    !!  @param pressure     air pressure             [Pa]
    !!  @param temperature  air temperature          [K]
    !!  @param qv           water vapor mixing ratio [kg/kg]
    !!  @param qc           cloud water mixing ratio [kg/kg]
    !!  @param qvsat        saturated mixing ratio   [kg/kg]
    !!  @param dt           time step                [sec.]
    !!
    !!----------------------------------------------------------
    subroutine cloud_conversion(pressure,temperature,qv,qc,qvsat,dt)
        implicit none
        real,intent(inout)  :: temperature,qv,qc,qvsat
        real,intent(in)     :: dt,pressure
        real    :: vapor2temp,excess,deltat,pre_qc,pre_qv,pre_t,lastqv
        integer :: iteration
        real    :: maxerr
        real :: qv_iters(15), t_iters(15)

        maxerr = 1e-4
        iteration = 0
        lastqv = qv + maxerr * 2
        vapor2temp = (LH_vapor + (373.15 - temperature) * dLHvdt) / heat_capacity
        pre_qc = qc !DEBUG
        pre_qv = qv !DEBUG
        pre_t = temperature !DEBUG
        excess = 0

        ! this should be written to use a slightly smarter iteration scheme.
        do while ((abs(lastqv-qv)>maxerr).and.(iteration<15))
            iteration = iteration+1
            lastqv = qv
            ! calculate the saturating mixing ratio
            qvsat = sat_mr(temperature,pressure)
            ! if saturated create clouds
            if (qv>qvsat) then
                excess = (qv-qvsat)*0.5
                ! temperature change if all vapor is converted
                temperature = temperature+(excess*vapor2temp)
                qv = qv-excess
                qc = qc+excess
            ! if unsaturated and clouds exist, evaporate clouds
            else if (qc>0) then
                excess = (qvsat-qv)*0.5
                if (excess<qc) then
                    temperature = temperature-(excess*vapor2temp)
                    qv = qv+excess
                    qc = qc-excess
                else
                    qv = qv+qc
                    temperature = temperature-(qc*vapor2temp)
                    excess = qc
                    qc = 0.
                endif
                excess = excess*(-1) !DEBUG
            endif
            ! qv_iters(iteration) = qv
            ! t_iters(iteration) = temperature
        enddo
        ! print*, iteration
        if (iteration==15) then
            ! !$omp critical (print_lock)
            ! print*, "qv_simple failed to converge", pre_t, t_iters(1:6)!pre_t, temperature, pressure, qv, lastqv, pre_qv, qvsat, abs(lastqv-qv)
            ! !$omp end critical (print_lock)
            qv = sat_mr(pre_t,pressure)
            temperature = pre_t
            qc = pre_qc
        endif

        ! prevents precision limits from crashing later code
        qc = max(qc,0.)

        if ((temperature > 350).or.(qc > 0.01).or.(qvsat > 1)) then
            deltat = excess*vapor2temp
            !$omp critical (print_lock)
            print*, ""
            print*, ""
            print*, "mp_simple ERROR: "
            if (temperature > 350 ) print*, "Temperature out of bounds", temperature
            if (qc          > 0.01) print*, "Qc out of bounds", qc
            if (qvsat       > 1   ) print*, "saturated qv out of bounds", qvsat
            print*, "iter=",iteration
            print*, ""
            print*, "   preqc=",pre_qc,"       preqv=",pre_qv,"        pret=",pre_t
            print*, "   qc=",qc, "      qv=",qv,"       temperature=",temperature
            print*, "   pressure=",pressure, "        qvsat=",qvsat
            print*, "   mrs_pret=",sat_mr(pre_t,pressure),"     mrs_t=",sat_mr(temperature,pressure), "     mrs_delta=",sat_mr(pre_t,pressure)-sat_mr(temperature,pressure)
            print*, "   excess=",excess, "      vapor2temp=",vapor2temp
            print*, "   temperature-deltat=",temperature-deltat, "      mrs_(t-dT)=",sat_mr(temperature-deltat,pressure)
            !$omp end critical (print_lock)
        endif

    end subroutine


    !>----------------------------------------------------------
    !!  Convert cloud water or ice to rain or snow
    !!
    !!  Use a time constant to calculate the convertion between
    !!  cloud and hydrometeor and enforce reasonable bounds
    !!
    !!  @param qc           Cloud water (or ice) mixing ratio       [kg/kg]
    !!  @param q            Rain (or snow) mixing ratio             [kg/kg]
    !!  @param conversion   time constant for conversion (*dt)      []
    !!  @param qcmin        minimum cloud content before conversion [kg/kg]
    !!
    !!----------------------------------------------------------
    subroutine cloud2hydrometeor(qc,q,conversion,qcmin)
        implicit none
        real,intent(inout) :: qc,q
        real,intent(in) :: conversion, qcmin
        real::delta

        if (qc > qcmin) then
            delta = qc-(qc*conversion)
        else
            delta = 0
        endif

        if (delta<qc) then
            qc = qc-delta
            q = q+delta
        else
            q = q+qc
            qc = 0.
        endif
        qc = max(qc,0.)
    end subroutine

    !>----------------------------------------------------------
    !!  Change the "phase" of a hydrometeor (e.g. evaporate rain)
    !!
    !!  Written to apply generically to any conversion.
    !!  Convert from q1 (e.g. rain) to q2 (e.g. vapor) using Lheat to affect temperature
    !!  and change_rate to control the conversion time scale
    !!
    !! @param pressure    air pressure (not used)                       [Pa]
    !! @param temperature air temperature                               [K]
    !! @param q1          hydrometeor to convert from                   [kg/kg]
    !! @param qmax        the maximum value for q2 (e.g. qv_sat)        [kg/kg]
    !! @param q2          hydrometeor to convert to                     [kg/kg]
    !! @param Lheat       latent heat associated with the phase change  [J/kg]
    !! @param change_rate fraction of q1 that can change in a timestep  []
    !!
    !!----------------------------------------------------------
    subroutine phase_change(pressure,temperature,q1,qmax,q2,Lheat,change_rate)
        implicit none
        real, intent(inout) :: temperature, q1, q2
        real, intent(in)    :: pressure, qmax, change_rate, Lheat
        real :: mass2temp, delta

        mass2temp = Lheat/heat_capacity!*(pressure/(R*temperature)*dV))

        delta = (qmax-q2)*change_rate
        if (delta>q1) delta = q1
        ! hopefully we don't over shoot saturation (use a 1% buffer)
        if (delta>((qmax-q2)*0.99)) then
            delta = (qmax-q2)*0.99
        endif

        q1 = q1-delta
        ! bounds checking
        if (q1<0) then
            if ((q1+SMALL_VALUE)<0) then
                q1 = 0
            else
                print*, "phase_change"
                print*, q1,q2,delta,qmax,change_rate
                stop
            endif
        endif
        q2 = q2+delta
        temperature = temperature+delta*mass2temp

    end subroutine

    !>----------------------------------------------------------
    !!  Compute microphysical conversions
    !!
    !!  Convert cloud water to and from vapor
    !!  Convert cloud water to rain, and cloud ice to snow
    !!  Cloud water is assumed to be ice if temperature is less than 0C
    !!  Conversions also include latent heat feedbacks to temperature
    !!
    !!  @param pressure     air pressure             [Pa]
    !!  @param temperature  air temperature          [K]
    !!  @param qv           water vapor mixing ratio [kg/kg]
    !!  @param qc           cloud water mixing ratio [kg/kg]
    !!  @param qr           rain mixing ratio        [kg/kg]
    !!  @param qs           snow mixing ratio        [kg/kg]
    !!  @param dt           time step                [sec.]
    !!
    !!----------------------------------------------------------
    subroutine mp_conversions(pressure, temperature, qv, qc, qr, qs, dt)
        implicit none
        real, intent(inout) :: pressure, temperature, qv, qc, qr, qs
        real, intent(in)    :: dt
        real :: qvsat, L_evap, L_subl, L_melt

        L_melt = -1 * LH_liquid  ! J/kg (should change with temperature)
        L_evap = -1 * (LH_vapor + (373.15 - temperature) * dLHvdt)   ! J/kg
        L_subl = L_melt + L_evap ! J/kg
        ! convert cloud water to and from water vapor (also initializes qvsat)
        call cloud_conversion(pressure, temperature, qv, qc, qvsat, dt)
        ! if there are no species to process we will just return
        if ((qc+qr+qs) > SMALL_VALUE) then
            if (qc > SMALL_VALUE) then
                if (temperature > freezing_threshold) then
                    ! convert cloud water to rain drops
                    call cloud2hydrometeor(qc, qr, cloud2rain, rain_cloud_init)
                    if (qs > SMALL_VALUE) then
                        ! it is above freezing, so start melting any snow if present
                        call phase_change(pressure, temperature, qs, 100., qr, L_melt, cloud2rain)
                    endif
                else
                    ! convert cloud water to snow flakes
                    call cloud2hydrometeor(qc, qs, cloud2snow, snow_cloud_init)

                endif
            endif
            ! if unsaturated, evaporate any existing snow and rain
            if (qv < qvsat) then
                if (qr > SMALL_VALUE) then
                    ! evaporate rain
                    call phase_change(pressure, temperature, qr, qvsat, qv, L_evap, cloud2rain/2)
                endif
                if (qs > SMALL_VALUE) then
                    ! sublimate snow
                    call phase_change(pressure, temperature, qs, qvsat, qv, L_subl, cloud2snow/2)
                endif
            endif
        endif
    end subroutine

    !>----------------------------------------------------------
    !!  Compute sedimentation in a column
    !!
    !!  Takes a mixing ratio of some species along with the fall velocities,
    !!  air densities and layer thicknesses.
    !!
    !!  @param  q   1D microphysical species to sediment (snow, rain, etc) [kg/kg]
    !!  @param  v   1D vertical velocity [m/s]
    !!  @param  rho 1D air density [kg/m^3]
    !!  @param  dz  1D layer thickness [m]
    !!  @param  n   0D number of layers
    !!  @param  its 0D bottom layer to process
    !!  @param  ite 0D top layer to process
    !!
    !!----------------------------------------------------------
    real function sediment(q, v, rho, dz, kms, kme, kts, kte)
        implicit none
        real,   intent(inout):: q   (kms:kme)
        real,   intent(in)   :: v   (kms:kme)
        real,   intent(in)   :: rho (kms:kme)
        real,   intent(in)   :: dz  (kms:kme)
        integer,intent(in)   :: kms, kme, kts, kte
        real    :: flux(kms:kme)
        integer :: i

        ! calculate the mass of material falling out of the bottom model level
        sediment = v(kts) * q(kts) * rho(kts) ![m] * [kg/kg] * [kg/m^3] = [kg/m^2]
        ! remove that from the bottom model layer.
        q(kts) = q(kts) - (sediment / dz(kts) / rho(kts)) ! [kg/m^2] / [m] / [kg/m^3] = [kg/kg]
        do i = kts, min(kte, kme-1)
            flux(i) = v(i+1) * q(i+1) * rho(i+1)
        enddo
        do i = kts, min(kte, kme-1)
            q(i) = q(i) + flux(i) / (rho(i) * dz(i))
            q(i+1) = q(i+1) - flux(i) / (rho(i+1) * dz(i+1))
        enddo

    end function

    !>----------------------------------------------------------
    !!  Basic microphysics code for a column of air
    !!
    !!  Call microphysical conversion routines on each layer and sedimentation for the column
    !!
    !!  @param pressure   = pressure                      - 1D - input  - Pa     - (nz)
    !!  @param temperature= air temperature               - 1D - in/out - K      - (nz)
    !!  @param rho        = air density                   - 1D - input  - kg/m^3 - (nz)
    !!  @param qv         = specific humidity             - 1D - in/out - kg/kg  - (nz)
    !!  @param qc         = cloud water content           - 1D - in/out - kg/kg  - (nz)
    !!  @param qr         = rain water content            - 1D - in/out - kg/kg  - (nz)
    !!  @param qs         = snow water content            - 1D - in/out - kg/kg  - (nz)
    !!  @param rain       = accumulated rain              - 0D - output - mm     - ()
    !!  @param snow       = accumulated snow              - 0D - output - mm     - ()
    !!  @param dt         = time step                     - 0D - input  - sec.   - scalar
    !!  @param dz         = vertical thickness of layers  - 1D - input  - m      - (nz)
    !!  @param kms, kme   = start end of z array memory   - 0D - input  - n      - scalar
    !!  @param kts, kte   = start end of z tile to process- 0D - input  - n      - scalar
    !!
    !!----------------------------------------------------------
    subroutine mp_simple(pressure, temperature, rho, qv, qc, qr, qs, rain, snow, dt, dz, kms, kme, kts, kte)
        implicit none
        real,   intent(inout)   :: pressure     (kms:kme)
        real,   intent(inout)   :: temperature  (kms:kme)
        real,   intent(inout)   :: rho          (kms:kme)
        real,   intent(inout)   :: qv           (kms:kme)
        real,   intent(inout)   :: qc           (kms:kme)
        real,   intent(inout)   :: qr           (kms:kme)
        real,   intent(inout)   :: qs           (kms:kme)
        real,   intent(inout)   :: rain, snow
        real,   intent(in)      :: dz           (kms:kme)
        real,   intent(in)      :: dt
        integer,intent(in)      :: kms, kme, kts, kte

        real    :: fall_rate(kms:kme)
        real    :: cfl, snowfall, qvsat
        integer :: i, cfl_step
        real    :: L_evap, L_subl, L_melt

        L_melt = -1 * LH_liquid  ! J/kg (should change with temperature)

        do i = kts, kte
            call mp_conversions(pressure(i), temperature(i), qv(i), qc(i), qr(i), qs(i), dt)
        enddo

        ! SEDIMENTATION for rain
        if (maxval(qr)>SMALL_VALUE) then
            fall_rate = rain_fall_rate

            ! check the CFL criteria for the sedimentation routine
            cfl = ceiling(maxval( dt / dz * fall_rate))
            ! update the fall_rate to be grid cell and CFL time step relative
            fall_rate = dt * fall_rate / cfl
            ! substepping to satisfy CFL criteria
            do cfl_step = 1, nint(cfl)
                rain = rain + sediment(qr, fall_rate, rho, dz, kms, kme, kts, kte)
                ! allow any rain that reached an unsaturated layer to evaporate
                do i = kts, kte
                    L_evap = -1 * (LH_vapor + (373.15 - temperature(i)) * dLHvdt)   ! J/kg
                    qvsat  = sat_mr(temperature(i), pressure(i))
                    ! if the air is not saturated, try evaporating rain
                    if (qv(i) < qvsat) then
                        if (qr(i) > SMALL_VALUE) then
                            ! evaporate rain
                            call phase_change(pressure(i), temperature(i), qr(i), qvsat, qv(i), L_evap, cloud2rain/(2*nint(cfl)))
                        endif
                    endif
                enddo
            enddo

        endif

        ! SEDIMENTATION for snow
        if (maxval(qs) > SMALL_VALUE) then
            fall_rate = snow_fall_rate

            ! check the CFL criteria for the sedimentation routine
            cfl = ceiling(maxval( dt / dz * fall_rate))
            ! update the fall_rate to be grid cell and CFL time step relative
            fall_rate = dt * fall_rate / cfl
            ! substepping to satisfy CFL criteria
            do cfl_step = 1, nint(cfl)
                snowfall = sediment(qs, fall_rate, rho, dz, kms, kme, kts, kte)
                snow = snow + snowfall
                rain = rain + snowfall

                ! allow any snow that reached a unsaturated layer to sublimate
                do i = kts, kte
                    L_evap = -1 * (LH_vapor + (373.15 - temperature(i)) * dLHvdt)   ! J/kg
                    L_subl = L_melt + L_evap ! J/kg
                    qvsat  = sat_mr(temperature(i), pressure(i))
                    ! if not saturated, sublimate the snow
                    ! if tracking snow and rain separately, should probably think about
                    ! condensing water on the snow and melting some into rain...
                    if (qv(i) < qvsat) then
                        if (qs(i) > SMALL_VALUE) then
                            ! sublimate snow
                            call phase_change(pressure(i), temperature(i), qs(i),                  &
                                              qvsat, qv(i), L_subl, cloud2snow/(2*nint(cfl)))
                        endif
                    endif
                enddo
            enddo
        endif

    end subroutine mp_simple


    !>----------------------------------------------------------
    !!  Driver code to control simple microphysics
    !!
    !!  Handles horizontal spatial dimensions and parallelization in space,
    !!  Calculates time constants for the current dt,
    !!  Calculates real temperature from potential temperature
    !!
    !!  @param pressure   = pressure                      - 3D - input  - Pa     - (nx,nz,ny)
    !!  @param th         = potential temperature         - 3D - in/out - K      - (nx,nz,ny)
    !!  @param pii        = exner function                - 3D - input  - []     - (nx,nz,ny)
    !!  @param rho        = air density                   - 3D - input  - kg/m^3 - (nx,nz,ny)
    !!  @param qv         = specific humidity             - 3D - in/out - kg/kg  - (nx,nz,ny)
    !!  @param qc         = cloud water content           - 3D - in/out - kg/kg  - (nx,nz,ny)
    !!  @param qr         = rain water content            - 3D - in/out - kg/kg  - (nx,nz,ny)
    !!  @param qs         = snow water content            - 3D - in/out - kg/kg  - (nx,nz,ny)
    !!  @param rain       = accumulated rain              - 2D - output - mm     - (nx,ny)
    !!  @param snow       = accumulated snow              - 2D - output - mm     - (nx,ny)
    !!  @param dt         = time step                     - 0D - input  - sec.   - scalar
    !!  @param ims, ime   = start end of x array memory   - 0D - input  - n      - scalar
    !!  @param jms, jme   = start end of y array memory   - 0D - input  - n      - scalar
    !!  @param kms, kme   = start end of z array memory   - 0D - input  - n      - scalar
    !!  @param its, ite   = start end of x tile to process- 0D - input  - n      - scalar
    !!  @param jts, jte   = start end of y tile to process- 0D - input  - n      - scalar
    !!  @param kts, kte   = start end of z tile to process- 0D - input  - n      - scalar
    !!
    !!----------------------------------------------------------
    subroutine mp_simple_driver(pressure,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,    &
                                ims, ime, jms, jme, kms, kme, &
                                its, ite, jts, jte, kts, kte)
        implicit none
        real, intent(inout) :: pressure (ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: th       (ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: pii      (ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: rho      (ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: qv       (ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: qc       (ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: qs       (ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: qr       (ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: rain     (ims:ime,jms:jme)
        real, intent(inout) :: snow     (ims:ime,jms:jme)
        real, intent(in)    :: dt
        real, intent(in)    :: dz       (ims:ime,kms:kme,jms:jme)
        integer,intent(in)  :: ims, ime, jms, jme, kms, kme
        integer,intent(in)  :: its, ite, jts, jte, kts, kte

        ! local variables
        real,allocatable,dimension(:) :: temperature
        integer :: i,j

!       calculate these once for every call because they are only a function of dt
        cloud2snow = exp(-1.0*snow_formation_time_const*dt)
        cloud2rain = exp(-1.0*rain_formation_time_const*dt)
        snow_evap  = exp(-1.0*snow_evap_time_const*dt)
        rain_evap  = exp(-1.0*rain_evap_time_const*dt)

        !$omp parallel private(i,j,temperature), &
        !$omp copyin(cloud2rain,cloud2snow,snow_melt,snow_evap,rain_evap),&
        !$omp shared(pressure,th,pii,qv,qc,qs,qr,rain,snow,dz),&
        !$omp firstprivate(dt,ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        allocate(temperature(kts:kte))
        if (dt<1e-10) print*, "dt=",dt ! for some reason this seems to fix an issue with ifort -O(anything) not copying dt in to the parallel region(?)
        !$omp do
        do j = jts,jte
            do i = its, ite
                temperature = th(i,:,j) * pii(i,:,j)
                ! should probably test out explicit 1D temporaries to see which is faster
                call mp_simple(pressure(i,:,j), temperature, rho(i,:,j), qv(i,:,j), &
                               qc(i,:,j), qr(i,:,j), qs(i,:,j),                     &
                               rain(i,j), snow(i,j),                                &
                               dt, dz(i,:,j), kms, kme, kts, kte)
                th(i,:,j) = temperature / pii(i,:,j)

            enddo
        enddo
        !$omp end do
        deallocate(temperature)
        !$omp end parallel
    end subroutine mp_simple_driver
end module
