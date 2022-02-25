!>----------------------------------------------------------
!! This module provides basic atmospheric utility functions.
!!
!!  Utilities exist to convert u, v into speed and direction (and vice versa)
!!  Compute the dry and moist lapse rates and Brunt Vaisalla stabilities
!!  Compute the terrain blocking froude number (U / (terrain_height * brunt_vaisalla_frequency ))
!!  Compute the fraction of flow that is blocked (linear interpolation between froude_max and froude_min)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module mod_atm_utilities

    use icar_constants,           only : pi, gravity, Rd, Rw, cp, LH_vaporization
    ! use data_structures
    use options_interface,        only : options_t

    implicit none

    real,     private :: N_squared  = 1e-5
    logical,  private :: variable_N = .True.
    real,     private :: max_froude, min_froude, froude_gain

contains


    !>----------------------------------------------------------
    !! Compute a 3D height field given a surface (or sea level) pressure
    !! and 3D temperature, humidity and pressures.
    !!
    !! Input temperature is real temperature in Kelvin  [K]
    !! Input humidity is mixing ratio                   [kg/kg]
    !! Pressures (input and output) are in Pascals      [Pa]
    !! Change in height is in meters                    [m]
    !!
    !!----------------------------------------------------------
    subroutine compute_3d_z(p, ps, z, t, qv, zs)
        implicit none
        real, intent(in),            dimension(:,:,:)   :: p
        real, intent(in),            dimension(:,:)     :: ps  ! surface (or sea level) pressure
        real, intent(inout),         dimension(:,:,:)   :: z   ! height above sea level for each atmospheric level
        real, intent(in),            dimension(:,:,:)   :: t   ! air temperature (real) [K]
        real, intent(in),            dimension(:,:,:)   :: qv  ! water vapor mixing ratio
        real, intent(in),  optional, dimension(:,:)     :: zs  ! if present, this is the height above z that the first level is computed for.

        integer :: i

        i=1
        call compute_z_offset(z(:,i,:), p(:,i,:) / ps, t(:,i,:), qv(:,i,:))

        if (present(zs)) then
            z(:,i,:) = zs - z(:,i,:)
        else
            z(:,i,:) = 0 - z(:,i,:)
        endif

        do i=2, size(z, 2)
            call compute_z_offset(z(:,i,:), p(:,i,:) / p(:,i-1,:), (t(:,i,:)+t(:,i-1,:))/2, (qv(:,i,:)+qv(:,i-1,:))/2)
            z(:,i,:) = z(:,i-1,:) - z(:,i,:)
        enddo

    end subroutine


    !>----------------------------------------------------------
    !! Compute the change in height for a change in pressure,
    !! and the temperature and humidity in between them
    !!
    !! Input temperature is real temperature in Kelvin  [K]
    !! Input humidity is mixing ratio                   [kg/kg]
    !! Pressures (input and output) are in Pascals      [Pa]
    !! Change in height is in meters                    [m]
    !!
    !!----------------------------------------------------------
    subroutine compute_z_offset(dz_out, p_ratio, t, qv)
        implicit none
        real, intent(inout),    dimension(:,:)   :: dz_out  ! change in height caused by change in pressure [m]
        real, intent(in),       dimension(:,:)   :: p_ratio ! ratio of pressure_in and pressure_out [Pa/Pa]
        real, intent(in),       dimension(:,:)   :: t       ! temperature in layer between p_out and p [K]
        real, intent(in),       dimension(:,:)   :: qv      ! water vapor in layer between p_out and p [kg/kg]

        dz_out = &
            Rd / gravity * ( t * ( 1 + 0.608 * qv ) ) * &
            LOG ( p_ratio )

        ! WRF formulate to compute height from pressure
        ! z(k) = z(k-1) - &
        !     R_d / g * 0.5 * ( t(k) * ( 1 + 0.608 * qv(k) ) +   &
        !                     t(k-1) * ( 1 + 0.608 * qv(k-1) ) ) * &
        !     LOG ( p(k) / p(k-1) )

    end subroutine compute_z_offset




    !>----------------------------------------------------------
    !! Compute a 3D pressure field given a surface (or sea level) pressure
    !! and 3D temperature, humidity and their corresponding heights.
    !!
    !! Input temperature is real temperature in Kelvin  [K]
    !! Input humidity is mixing ratio                   [kg/kg]
    !! Pressures (input and output) are in Pascals      [Pa]
    !! Change in height is in meters                    [m]
    !!
    !!----------------------------------------------------------
    subroutine compute_3d_p(p, ps, z, t, qv, zs)
        implicit none
        real, intent(inout)        , dimension(:,:,:)   :: p
        real, intent(in)           , dimension(:,:)     :: ps  ! surface (or sea level) pressure
        real, intent(in)           , dimension(:,:,:)   :: z   ! height above sea level for each atmospheric level
        real, intent(in)           , dimension(:,:,:)   :: t   ! air temperature (real) [K]
        real, intent(in)           , dimension(:,:,:)   :: qv  ! water vapor mixing ratio
        real, intent(in),  optional, dimension(:,:)     :: zs  ! if present, this is the height above z that the first level is computed for.

        integer :: i

        i=1
        if (present(zs)) then
            call compute_p_offset(p(:,i,:), ps, z(:,i,:)-zs, t(:,i,:), qv(:,i,:))
        else
            call compute_p_offset(p(:,i,:), ps, z(:,i,:), t(:,i,:), qv(:,i,:))
        endif

        do i=2, size(p,2)
            call compute_p_offset(p(:,i,:), p(:,i-1,:), z(:,i,:)-z(:,i-1,:), (t(:,i,:)+t(:,i-1,:))/2, (qv(:,i,:)+qv(:,i-1,:))/2)
        enddo

    end subroutine

    !>----------------------------------------------------------
    !! Compute the pressure of level p_out based on the pressure dz meters below,
    !! and the temperature and humidity in between them
    !!
    !! Input temperature is real temperature in Kelvin  [K]
    !! Input humidity is mixing ratio                   [kg/kg]
    !! Pressures (input and output) are in Pascals      [Pa]
    !! Change in height is in meters                    [m]
    !!
    !!----------------------------------------------------------
    subroutine compute_p_offset(p_out, p, dz, t, qv)
        implicit none
        real, intent(inout),    dimension(:,:)   :: p_out   ! output as p+dz [Pa]
        real, intent(in),       dimension(:,:)   :: p       ! input pressure dz distance below the output pressure [Pa]
        real, intent(in),       dimension(:,:)   :: dz      ! height to raise p to get p_out [m]
        real, intent(in),       dimension(:,:)   :: t       ! temperature in layer between p_out and p [K]
        real, intent(in),       dimension(:,:)   :: qv      ! water vapor in layer between p_out and p [kg/kg]

        p_out = p * exp( -dz / (Rd / gravity * ( t * ( 1 + 0.608 * qv ) )))

        ! note: derived from WRF formulate to compute height from pressure
        ! z(k) = z(k-1) - &
        !     R_d / g * 0.5 * ( t(k) * ( 1 + 0.608 * qv(k) ) +   &
        !                     t(k-1) * ( 1 + 0.608 * qv(k-1) ) ) * &
        !     LOG ( p(k) / p(k-1) )

    end subroutine compute_p_offset

    ! !>----------------------------------------------------------
    ! !! Compute the height of level z_out based on the pressure at z, z_out,
    ! !! and the temperature and humidity in between them
    ! !!
    ! !! Input temperature is real temperature in Kelvin  [K]
    ! !! Input humidity is mixing ratio                   [kg/kg]
    ! !! Input pressures are in Pascals                   [Pa]
    ! !! Heights (input and output) are in meters         [m]
    ! !!
    ! !!----------------------------------------------------------
    ! subroutine compute_z_offset(z_out, z, p0, p1, t, qv)
    !     implicit none
    !     real, intent(inout),    dimension(:,:)   :: z_out   !
    !     real, intent(in),       dimension(:,:)   :: z       !
    !     real, intent(in),       dimension(:,:)   :: p0      !
    !     real, intent(in),       dimension(:,:)   :: p1      !
    !     real, intent(in),       dimension(:,:)   :: t       !
    !     real, intent(in),       dimension(:,:)   :: qv      !
    !
    !     z_out = z - &
    !             Rd / gravity * ( t * ( 1 + 0.608 * qv ) ) * &
    !             LOG ( p1 / p0 )
    !
    !     ! note: WRF formulate to compute height from pressure
    !     ! z(k) = z(k-1) - &
    !     !     R_d / g * 0.5 * ( t(k) * ( 1 + 0.608 * qv(k) ) +   &
    !     !                     t(k-1) * ( 1 + 0.608 * qv(k-1) ) ) * &
    !     !     LOG ( p(k) / p(k-1) )
    !
    ! end subroutine compute_z_offset

    !>----------------------------------------------------------
    !! Convert relative humidity, temperature, and pressure to water vapor mixing ratio
    !!
    !! Input temperature is real temperature in Kelvin  [K]
    !! Input relative humidity is fractional            [0-1]
    !! Input pressure s in Pascals                      [Pa]
    !! Output mixing ratio is in kg / kg                [kg/kg]
    !!
    !!----------------------------------------------------------
    pure elemental function rh_to_mr(input_rh, t, p) result(mr)
        implicit none
        real, intent(in) :: input_rh, t, p
        real :: mr
        real :: es, e, rh

        rh = min(1.0, max(0.0, input_rh))

        ! saturated vapor pressure
        ! convert temperature to saturated vapor pressure (in Pa)
        es = 611.2 * exp(17.67 * (t - 273.15) / (t - 29.65))

        ! convert relative humidity to vapor pressure
        e = rh * es

        ! finally convert vapor pressure to mixing ratio
        mr = 0.62197 * e / (p - e)
    end function


    !>----------------------------------------------------------
    !! Convert temperature, specific humidity, and pressure to relative humidity
    !!
    !! Input temperature is real temperature in Kelvin  [K]
    !! Input specific humidity is in kg / kg            [kg/kg]
    !! Input pressure s in Pascals                      [Pa]
    !! Output relative humidity is fractional           [0-1]
    !!
    !!----------------------------------------------------------
    pure elemental function relative_humidity(t,qv,p)
        implicit none
        real               :: relative_humidity
        real,   intent(in) :: t
        real,   intent(in) :: qv, p
        real               :: mr, e, es

        ! convert specific humidity to mixing ratio
        mr = qv / (1-qv)
        ! convert mixing ratio to vapor pressure
        e = mr * p / (0.62197+mr)
        ! convert temperature to saturated vapor pressure
        es = 611.2 * exp(17.67 * (t - 273.15) / (t - 29.65))
        ! finally return relative humidity
        relative_humidity = e / es

        ! because it is an approximation things could go awry and rh outside or reasonable bounds could break something else.
        ! alternatively air could be supersaturated (esp. on boundary cells) but cloud fraction calculations will break.
        relative_humidity = min(1.0, max(0.0, relative_humidity))

    end function relative_humidity



    !>----------------------------------------------------------
    !! Calculate direction [0-2*pi) from u and v wind speeds
    !!
    !!----------------------------------------------------------
    pure function calc_direction(u,v) result(direction)
        implicit none
        real, intent(in) :: u,v
        real :: direction

        if (v<0) then
            direction = atan(u/v) + pi
        elseif (v==0) then
            if (u>0) then
                direction=pi/2.0
            else
                direction=pi*1.5
            endif
        else
            if (u>=0) then
                direction = atan(u/v)
            else
                direction = atan(u/v) + (2*pi)
            endif
        endif

    end function calc_direction

    !>----------------------------------------------------------
    !! Calculate the strength of the u wind field given a direction [0-2*pi] and magnitude
    !!
    !!----------------------------------------------------------
    pure elemental function calc_speed(u, v) result(speed)
        implicit none
        real, intent(in) :: u,v
        real :: speed

        speed = sqrt(u**2 + v**2)
    end function calc_speed

    !>----------------------------------------------------------
    !! Calculate the strength of the u wind field given a direction [0-2*pi] and magnitude
    !!
    !!----------------------------------------------------------
    pure elemental function calc_u(direction, magnitude) result(u)
        implicit none
        real, intent(in) :: direction, magnitude
        real :: u

        u = sin(direction) * magnitude
    end function calc_u

    !>----------------------------------------------------------
    !! Calculate the strength of the v wind field given a direction [0-2*pi] and magnitude
    !!
    !!----------------------------------------------------------
    pure elemental function calc_v(direction, magnitude) result(v)
        implicit none
        real, intent(in) :: direction, magnitude
        real :: v

        v = cos(direction) * magnitude
    end function calc_v

    !>----------------------------------------------------------
    !! Calculate the saturated adiabatic lapse rate from a T/Moisture input
    !!
    !! return the moist / saturated adiabatic lapse rate for a given
    !! Temperature and mixing ratio (really MR could be calculated as f(T))
    !! from http://glossary.ametsoc.org/wiki/Saturation-adiabatic_lapse_rate
    !!
    !!----------------------------------------------------------
    pure elemental function calc_sat_lapse_rate(T,mr) result(sat_lapse)
        implicit none
        real, intent(in) :: T,mr  ! inputs T in K and mr in kg/kg
        real :: L
        real :: sat_lapse

        L=LH_vaporization ! short cut for imported parameter
        sat_lapse = gravity*((1 + (L*mr) / (Rd*T))          &
                    / (cp + (L*L*mr*(Rd/Rw)) / (Rd*T*T) ))
    end function calc_sat_lapse_rate

    !>----------------------------------------------------------
    !! Calculate the moist brunt vaisala frequency (Nm^2)
    !! formula from Durran and Klemp, 1982 after Lalas and Einaudi 1974
    !!
    !!----------------------------------------------------------
    pure elemental function calc_moist_stability(t_top, t_bot, z_top, z_bot, qv_top, qv_bot, qc) result(BV_freq)
        implicit none
        real, intent(in) :: t_top, t_bot, z_top, z_bot, qv_top, qv_bot, qc
        real :: t,qv, dz, sat_lapse
        real :: BV_freq

        t  = ( t_top +  t_bot)/2
        qv = (qv_top + qv_bot)/2
        dz = ( z_top - z_bot)
        sat_lapse = calc_sat_lapse_rate(t,qv)

        BV_freq = (gravity/t) * ((t_top-t_bot)/dz + sat_lapse) * &
                  (1 + (LH_vaporization*qv)/(Rd*t)) - (gravity/(1+qv+qc) * (qv_top-qv_bot)/dz)
    end function calc_moist_stability

    !>----------------------------------------------------------
    !! Calculate the dry brunt vaisala frequency (Nd^2)
    !!
    !!----------------------------------------------------------
    pure elemental function calc_dry_stability(th_top, th_bot, z_top, z_bot) result(BV_freq)
        implicit none
        real, intent(in) :: th_top, th_bot, z_top, z_bot
        real :: BV_freq

        BV_freq = gravity * (log(th_top)-log(th_bot)) / (z_top - z_bot)
    end function calc_dry_stability

    !>----------------------------------------------------------
    !! Calculate either moist or dry brunt vaisala frequency
    !!
    !!----------------------------------------------------------
    pure function calc_stability(th_top, th_bot, pii_top, pii_bot, z_top, z_bot, qv_top, qv_bot, qc) result(BV_freq)
        implicit none
        real, intent(in) :: th_top, th_bot, pii_top, pii_bot, z_top, z_bot, qv_top, qv_bot, qc
        real :: BV_freq

        if (qc<1e-7) then
            if (variable_N) then
                BV_freq = calc_dry_stability(th_top, th_bot, z_top, z_bot)
            else
                BV_freq = N_squared
            endif
        else
            if (variable_N) then
                BV_freq = calc_moist_stability(th_top*pii_top, th_bot*pii_bot, z_top, z_bot, qv_top, qv_bot, qc)
            else
                BV_freq = N_squared/10.0 ! might be better as max(1e-7,N_squared-(1e-4))
            endif
        endif

    end function calc_stability

    !>----------------------------------------------------------
    !! Calculate the non-dimensional Froude number for flow over a barrier
    !!
    !! Used to identify topographically blocked flow (Fr < 0.75-1.25)
    !!
    !!----------------------------------------------------------
    pure function calc_froude(brunt_vaisalla_frequency, barrier_height, wind_speed) result(froude)
        implicit none
        real, intent(in) :: brunt_vaisalla_frequency    ! [ 1 / s ]
        real, intent(in) :: barrier_height              ! [ m ]
        real, intent(in) :: wind_speed                  ! [ m / s ]
        real :: froude                                  ! []
        real :: denom

        denom = (barrier_height * brunt_vaisalla_frequency)

        if (denom==0) then
            froude = 100 ! anything over ~5 is effectively infinite anyway
        else
            froude = wind_speed / denom
        endif

    end function calc_froude

    !>----------------------------------------------------------
    !! Compute the fraction of blocking winds to apply
    !!
    !!----------------------------------------------------------
    function blocking_fraction(froude) result(fraction)
        implicit none
        real, intent(in) :: froude
        real :: fraction

        fraction = (max_froude-froude) * froude_gain
        fraction = min(max( fraction, 0.), 1.)

    end function blocking_fraction


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
    elemental function sat_mr(temperature,pressure)
    ! Calculate the saturated mixing ratio at a temperature (K), pressure (Pa)
        implicit none
        real,intent(in) :: temperature,pressure
        real :: e_s,a,b
        real :: sat_mr

        ! from http://www.dtic.mil/dtic/tr/fulltext/u2/778316.pdf
        !   Lowe, P.R. and J.M. Ficke., 1974: THE COMPUTATION OF SATURATION VAPOR PRESSURE
        !       Environmental Prediction Research Facility, Technical Paper No. 4-74
        ! which references:
        !   Murray, F. W., 1967: On the computation of saturation vapor pressure.
        !       Journal of Applied Meteorology, Vol. 6, pp. 203-204.
        ! Also notes a 6th order polynomial and look up table as viable options.
        if (temperature < 273.15) then
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

    !> -------------------------------
    !!
    !! Convert p [Pa] at shifting it to a given elevatiom [m]
    !!
    !! -------------------------------
    elemental function pressure_at_elevation(sealevel_pressure, elevation) result(pressure)
        implicit none
        real, intent(in) :: sealevel_pressure, elevation
        real :: pressure

        pressure = sealevel_pressure * (1 - 2.25577E-5 * elevation)**5.25588

    end function

    !>------------------------------------------------------------
    !!  Adjust the pressure field for the vertical shift between the low and high-res domains
    !!
    !!  Ideally this should include temperature... but it isn't entirely clear
    !!  what it would mean to do that, what temperature do you use? Current time-step even though you are adjusting future time-step P?
    !!  Alternatively, could adjust input pressure to SLP with future T then adjust back to elevation with current T?
    !!  Currently if T is supplied, it uses the mean of the high and low-res T to split the difference.
    !!  Equations from : http://www.wmo.int/pages/prog/www/IMOP/meetings/SI/ET-Stand-1/Doc-10_Pressure-red.pdf
    !!  excerpt from CIMO Guide, Part I, Chapter 3 (Edition 2008, Updated in 2010) equation 3.2
    !!  http://www.meteormetrics.com/correctiontosealevel.htm
    !!
    !! @param pressure  The pressure field to be adjusted
    !! @param z_lo      The 3D vertical coordinate of the input pressures
    !! @param z_hi      The 3D vertical coordinate of the computed/adjusted pressures
    !! @param lowresT   OPTIONAL 3D temperature field of the input pressures
    !! @param lowresT   OPTIONAL 3D temperature field of the computed/adjusted pressures
    !! @retval pressure The pressure field after adjustment
    !!
    !!------------------------------------------------------------
    subroutine update_pressure(pressure,z_lo,z_hi, lowresT, hiresT)
        implicit none
        real,dimension(:,:,:), intent(inout) :: pressure
        real,dimension(:,:,:), intent(in) :: z_lo,z_hi
        real,dimension(:,:,:), intent(in), optional :: lowresT, hiresT

        ! local variables 1D arrays operate on a complete x row at a time
        real,dimension(:),allocatable::slp !sea level pressure [Pa]
        ! vapor pressure, change in height, change in temperature with height and mean temperature
        real,dimension(:),allocatable:: dz, tmean !, e, dTdz
        integer :: nx, ny, nz, i, j, nz_lo
        nx = size(pressure,1)
        nz = size(pressure,2)
        nz_lo = size(z_lo, 2)
        ny = size(pressure,3)

        if (present(lowresT)) then
            ! OpenMP parallelization directives
            !$omp parallel shared(pressure, z_lo,z_hi, lowresT, hiresT) &
            !$omp private(i,j, dz, tmean) firstprivate(nx,ny,nz)  !! private(e, dTdz)

            ! create the temporary variables needed internally (must be inside the parallel region)
            allocate(dz(nx))
            allocate(tmean(nx))

            !$omp do
            do j=1,ny
                ! is an additional loop over z more cache friendly?
                do i=1,nz
                    ! vapor pressure
!                     e = qv(:,:,j) * pressure(:,:,j) / (0.62197+qv(:,:,j))

                    ! change in elevation (note reverse direction from "expected" because the formula is an SLP reduction)
                    dz   = (z_lo(:,min(i, nz_lo),j) - z_hi(:,i,j))

                    ! lapse rate (not sure if this should be positive or negative)
                    ! dTdz = (loresT(:,:,j) - hiresT(:,:,j)) / dz
                    ! mean temperature between levels
                    if (present(hiresT)) then
                        tmean= (hiresT(:,i,j) + lowresT(:,min(i, nz_lo),j)) / 2
                    else
                        tmean= lowresT(:,min(i, nz_lo),j)
                    endif

                    ! Actual pressure adjustment
                    ! slp= ps*np.exp(((g/R)*Hp) / (ts - a*Hp/2.0 + e*Ch))
                    pressure(:,i,j) = pressure(:,i,j) * exp( ((gravity/Rd) * dz) / tmean )   !&
                    !                     (tmean + (e * 0.12) ) ) ! alternative

                    ! alternative formulation M=0.029, R=8.314?
                    ! p= p0*(t0/(t0+dtdz*z))**((g*M)/(R*dtdz))
                    ! do i=1,nz
                    !     pressure(:,i,j) = pressure(:,i,j)*(t0/(tmean(:,i)+dTdz(:,i)*z))**((g*M)/(R*dtdz))
                    ! enddo
                enddo
            enddo
            !$omp end do

            deallocate(dz, tmean)
            !$omp end parallel
        else

            ! OpenMP parallelization directives
            !$omp parallel shared(pressure, z_lo,z_hi) &
            !$omp private(slp,i,j) firstprivate(nx,ny,nz)

            ! allocate thread local data
            allocate(slp(nx))
            !$omp do
            do j=1,ny
                do i=1,nz
                    ! slp = pressure(:,i,j) / (1 - 2.25577E-5 * z_lo(:,i,j))**5.25588
                    pressure(:,i,j) = pressure(:,i,j) * (1 - 2.25577e-5 * (z_hi(:,i,j)-z_lo(:,min(i, nz_lo),j)))**5.25588
                enddo
            enddo
            !$omp end do
            deallocate(slp)
            !$omp end parallel
        endif
    end subroutine update_pressure


    !> -------------------------------
    !!
    !! Compute exner function to convert potential_temperature to temperature
    !!
    !! -------------------------------
    elemental function exner_function(pressure) result(exner)
        implicit none
        real, intent(in) :: pressure
        real :: exner

        associate(po=>100000) !, Rd=>287.058, cp=>1003.5)
            exner = (pressure / po) ** (Rd/cp)

        end associate
    end function


    !> -------------------------------
    !!
    !! Initialize module level variables with configuration options
    !!
    !! -------------------------------
    subroutine init_atm_utilities(options)
        implicit none
        type(options_t) :: options

        N_squared   = options%lt_options%N_squared
        variable_N  = options%lt_options%variable_N

        max_froude  = options%block_options%block_fr_max
        min_froude  = options%block_options%block_fr_min
        froude_gain = 1 / max(max_froude-min_froude, 0.001)

    end subroutine init_atm_utilities

!+---+-----------------------------------------------------------------+
!..Cloud fraction scheme by G. Thompson (NCAR-RAL), not intended for
!.. combining with any cumulus or shallow cumulus parameterization
!.. scheme cloud fractions.  This is intended as a stand-alone for
!.. cloud fraction and is relatively good at getting widespread stratus
!.. and stratoCu without caring whether any deep/shallow Cu param schemes
!.. is making sub-grid-spacing clouds/precip.  Under the hood, this
!.. scheme follows Mocko and Cotton (1995) in applicaiton of the
!.. Sundqvist et al (1989) scheme but using a grid-scale dependent
!.. RH threshold, one each for land v. ocean points based on
!.. experiences with HWRF testing.
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+

    SUBROUTINE cal_cldfra3(CLDFRA, qv, qc, qi, qs, dz,                &
        &                 p, t, XLAND, gridkm,                             &
        &                 modify_qvapor, max_relh,                         &
        &                 kts,kte)
   !
        USE module_mp_thompson   , ONLY : rsif, rslf
        IMPLICIT NONE
   !
        INTEGER, INTENT(IN):: kts, kte
        LOGICAL, INTENT(IN):: modify_qvapor
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: qv, qc, qi, cldfra
        REAL, DIMENSION(kts:kte), INTENT(IN):: p, t, dz, qs
        REAL, INTENT(IN):: gridkm, XLAND, max_relh

   !..Local vars.
        REAL:: RH_00L, RH_00O, RH_00
        REAL:: entrmnt=0.5
        INTEGER:: k
        REAL:: TC, qvsi, qvsw, RHUM, delz
        REAL, DIMENSION(kts:kte):: qvs, rh, rhoa

   !+---+

   !..Initialize cloud fraction, compute RH, and rho-air.

        DO k = kts,kte
            CLDFRA(K) = 0.0
            qvsw = rslf(P(k), t(k))
            qvsi = rsif(P(k), t(k))

            tc = t(k) - 273.15
            if (tc .ge. -12.0) then
                qvs(k) = qvsw
            elseif (tc .lt. -35.0) then
                qvs(k) = qvsi
            else
                qvs(k) = qvsw - (qvsw-qvsi)*(-12.0-tc)/(-12.0+35.)
            endif

            rh(k) = MAX(0.01, qv(k)/qvs(k))
            rhoa(k) = p(k)/(287.0*t(k))
        ENDDO


   !..First cut scale-aware. Higher resolution should require closer to
   !.. saturated grid box for higher cloud fraction.  Simple functions
   !.. chosen based on Mocko and Cotton (1995) starting point and desire
   !.. to get near 100% RH as grid spacing moves toward 1.0km, but higher
   !.. RH over ocean required as compared to over land.

        DO k = kts,kte

            delz = MAX(100., dz(k))
            RH_00L = 0.65 + SQRT(1./(25.0+gridkm*gridkm*delz*0.01))
            RH_00O = 0.81 + SQRT(1./(50.0+gridkm*gridkm*delz*0.01))
            RHUM = rh(k)

            if (qc(k).gt.1.E-7 .or. qi(k).ge.1.E-7                         &
        &                    .or. (qs(k).gt.1.E-6 .and. t(k).lt.273.)) then
               CLDFRA(K) = 1.0
               qvs(k) = qv(k)
            else

                IF ((XLAND-1.5).GT.0.) THEN                                  !--- Ocean
                    RH_00 = RH_00O
                ELSE                                                         !--- Land
                    RH_00 = RH_00L
                ENDIF

                tc = t(k) - 273.15
                if (tc .lt. -12.0) RH_00 = RH_00L

                if (tc .ge. 20.0) then
                    CLDFRA(K) = 0.0
                elseif (tc .ge. -12.0) then
                    RHUM = MIN(rh(k), 1.0)
                    CLDFRA(K) = MAX(0., 1.0-SQRT((1.005-RHUM)/(1.005-RH_00)))
                else
                    if (max_relh.gt.1.12 .or. (.NOT.(modify_qvapor)) ) then
   !..For HRRR model, the following look OK.
                        RHUM = MIN(rh(k), 1.45)
                        RH_00 = RH_00 + (1.45-RH_00)*(-12.0-tc)/(-12.0+100.)
                        if (RH_00 .ge. 1.5) then
                            WRITE (*,*) ' FATAL: RH_00 too large (1.5): ', RH_00, RH_00L, tc
                        endif
                        CLDFRA(K) = MAX(0., 1.0-SQRT((1.5-RHUM)/(1.5-RH_00)))
                    else
   !..but for the GFS model, RH is way lower.
                        RHUM = MIN(rh(k), 1.05)
                        RH_00 = RH_00 + (1.05-RH_00)*(-12.0-tc)/(-12.0+100.)
                        if (RH_00 .ge. 1.05) then
                            WRITE (*,*) ' FATAL: RH_00 too large (1.05): ', RH_00, RH_00L, tc
                        endif
                        CLDFRA(K) = MAX(0., 1.0-SQRT((1.05-RHUM)/(1.05-RH_00)))
                    endif
                endif
                if (CLDFRA(K).gt.0.) CLDFRA(K) = MAX(0.01, MIN(CLDFRA(K),0.9))
            endif
        ENDDO

        call find_cloudLayers(qvs, cldfra, T, P, Dz, entrmnt,             &
        &                      qc, qi, qs, kts,kte)

   !..Do a final total column adjustment since we may have added more than 1mm
   !.. LWP/IWP for multiple cloud decks.

        call adjust_cloudFinal(cldfra, qc, qi, rhoa, dz, kts,kte)
        if (modify_qvapor) then
            DO k = kts,kte
                if (cldfra(k).gt.0.20 .and. cldfra(k).lt.1.0) then
                  qv(k) = qvs(k)
                endif
            ENDDO
        endif

    END SUBROUTINE cal_cldfra3

!+---+-----------------------------------------------------------------+
!..From cloud fraction array, find clouds of multi-level depth and compute
!.. a reasonable value of LWP or IWP that might be contained in that depth,
!.. unless existing LWC/IWC is already there.

    SUBROUTINE find_cloudLayers(qvs1d, cfr1d, T1d, P1d, Dz1d, entrmnt,&
            &                            qc1d, qi1d, qs1d, kts,kte)
       !
        IMPLICIT NONE
       !
        INTEGER, INTENT(IN):: kts, kte
        REAL, INTENT(IN):: entrmnt
        REAL, DIMENSION(kts:kte), INTENT(IN):: qs1d,qvs1d,T1d,P1d,Dz1d
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: cfr1d, qc1d, qi1d

       !..Local vars.
        REAL, DIMENSION(kts:kte):: theta
        REAL:: theta1, theta2, delz
        INTEGER:: k, k2, k_tropo, k_m12C, k_cldb, k_cldt, kbot
        LOGICAL:: in_cloud

       !+---+

        k_m12C = 0
        DO k = kte, kts, -1
            theta(k) = T1d(k)*((100000.0/P1d(k))**(287.05/1004.))
            if (T1d(k)-273.16 .gt. -12.0 .and. P1d(k).gt.10100.0) k_m12C = MAX(k_m12C, k)
        ENDDO
        if (k_m12C .le. kts) k_m12C = kts

        if (k_m12C.gt.kte-3) then
            WRITE (*,*) 'DEBUG-GT: WARNING, no possible way neg12C can occur this high up: ', k_m12C
            do k = kte, kts, -1
                WRITE (*,*) 'DEBUG-GT,  k,  P, T : ', k,P1d(k)*0.01,T1d(k)-273.15
            enddo
            write(*,*) ('FATAL ERROR, problem in temperature profile.')
        endif

       !..Find tropopause height, best surrogate, because we would not really
       !.. wish to put fake clouds into the stratosphere.  The 10/1500 ratio
       !.. d(Theta)/d(Z) approximates a vertical line on typical SkewT chart
       !.. near typical (mid-latitude) tropopause height.  Since messy data
       !.. could give us a false signal of such a transition, do the check over
       !.. three K-level change, not just a level-to-level check.  This method
       !.. has potential failure in arctic-like conditions with extremely low
       !.. tropopause height, as would any other diagnostic, so ensure resulting
       !.. k_tropo level is above 700hPa.

        DO k = kte-3, kts, -1
            theta1 = theta(k)
            theta2 = theta(k+2)
            delz = dz1d(k) + dz1d(k+1) + dz1d(k+2)
            if ( ((((theta2-theta1)/delz) .lt. 10./1500. ) .AND.       &
            &                 (P1d(k).gt.8500.)) .or. (P1d(k).gt.70000.) ) then
                goto 86
            endif
        ENDDO
    86  continue
        k_tropo = MAX(kts+2, MIN(k+2, kte-1))

        !if (k_tropo.gt.kte-2) then
        !    WRITE (*,*) 'DEBUG-GT: CAUTION, tropopause appears to be very high up: ', k_tropo
        !    do k = kte, kts, -1
        !        WRITE (*,*) 'DEBUG-GT,   P, T : ', k,P1d(k)*0.01,T1d(k)-273.16
        !    enddo
        !endif

       !..Eliminate possible fractional clouds above supposed tropopause.
        DO k = k_tropo+1, kte
            if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.1.0) then
                cfr1d(k) = 0.
            endif
        ENDDO

       !..We would like to prevent fractional clouds below LCL in idealized
       !.. situation with deep well-mixed convective PBL, that otherwise is
       !.. likely to get clouds in more realistic capping inversion layer.

        kbot = kts+2
        DO k = kbot, k_m12C
            if ( (theta(k)-theta(k-1)) .gt. 0.025E-3*Dz1d(k)) EXIT
        ENDDO
        kbot = MAX(kts+1, k-2)
        DO k = kts, kbot
            if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.1.0) cfr1d(k) = 0.
        ENDDO

       !..Starting below tropo height, if cloud fraction greater than 1 percent,
       !.. compute an approximate total layer depth of cloud, determine a total
       !.. liquid water/ice path (LWP/IWP), then reduce that amount with tuning
       !.. parameter to represent entrainment factor, then divide up LWP/IWP
       !.. into delta-Z weighted amounts for individual levels per cloud layer.

        k_cldb = k_tropo
        in_cloud = .false.
        k = k_tropo
        DO WHILE (.not. in_cloud .AND. k.gt.k_m12C+1)
            k_cldt = 0
            if (cfr1d(k).ge.0.01) then
                in_cloud = .true.
                k_cldt = MAX(k_cldt, k)
            endif
            if (in_cloud) then
                DO k2 = k_cldt-1, k_m12C, -1
                    if (cfr1d(k2).lt.0.01 .or. k2.eq.k_m12C) then
                        k_cldb = k2+1
                        goto 87
                    endif
                ENDDO
        87      continue
                in_cloud = .false.
            endif
            if ((k_cldt - k_cldb + 1) .ge. 2) then
                call adjust_cloudIce(cfr1d, qi1d, qs1d, qvs1d, T1d, Dz1d,   &
            &                           entrmnt, k_cldb,k_cldt,kts,kte)
                k = k_cldb
            elseif ((k_cldt - k_cldb + 1) .eq. 1) then
                if (cfr1d(k_cldb).gt.0.and.cfr1d(k_cldb).lt.1.)             &
            &               qi1d(k_cldb)=0.05*qvs1d(k_cldb)
                k = k_cldb
            endif
                k = k - 1
        ENDDO


        k_cldb = k_m12C + 3
        in_cloud = .false.
        k = k_m12C + 2
        DO WHILE (.not. in_cloud .AND. k.gt.kbot)
            k_cldt = 0
            if (cfr1d(k).ge.0.01) then
                in_cloud = .true.
                k_cldt = MAX(k_cldt, k)
            endif
            if (in_cloud) then
                DO k2 = k_cldt-1, kbot, -1
                    if (cfr1d(k2).lt.0.01 .or. k2.eq.kbot) then
                        k_cldb = k2+1
                        goto 88
                    endif
                ENDDO
        88      continue
                in_cloud = .false.
            endif
            if ((k_cldt - k_cldb + 1) .ge. 2) then
                call adjust_cloudH2O(cfr1d, qc1d, qvs1d, T1d, Dz1d,         &
            &                           entrmnt, k_cldb,k_cldt,kts,kte)
                k = k_cldb
            elseif ((k_cldt - k_cldb + 1) .eq. 1) then
                if (cfr1d(k_cldb).gt.0.and.cfr1d(k_cldb).lt.1.)             &
            &                qc1d(k_cldb)=0.05*qvs1d(k_cldb)
                k = k_cldb
            endif
            k = k - 1
        ENDDO

    END SUBROUTINE find_cloudLayers

!+---+-----------------------------------------------------------------+

    SUBROUTINE adjust_cloudIce(cfr,qi,qs,qvs,T,dz,entr, k1,k2,kts,kte)
                !
        IMPLICIT NONE
                !
        INTEGER, INTENT(IN):: k1,k2, kts,kte
        REAL, INTENT(IN):: entr
        REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, qs, qvs, T, dz
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: qi
        REAL:: iwc, max_iwc, tdz, this_iwc, this_dz
        INTEGER:: k

        tdz = 0.
        do k = k1, k2
            tdz = tdz + dz(k)
        enddo
        max_iwc = ABS(qvs(k2)-qvs(k1))

        do k = k1, k2
            max_iwc = MAX(1.E-6, max_iwc - (qi(k)+qs(k)))
        enddo
        max_iwc = MIN(1.E-3, max_iwc)

        this_dz = 0.0
        do k = k1, k2
            if (k.eq.k1) then
                this_dz = this_dz + 0.5*dz(k)
            else
                this_dz = this_dz + dz(k)
            endif
            this_iwc = max_iwc*this_dz/tdz
            iwc = MAX(1.E-6, this_iwc*(1.-entr))
            if (cfr(k).gt.0.0.and.cfr(k).lt.1.0.and.T(k).ge.203.16) then
                qi(k) = qi(k) + cfr(k)*cfr(k)*iwc
            endif
        enddo

    END SUBROUTINE adjust_cloudIce

    !+---+-----------------------------------------------------------------

    SUBROUTINE adjust_cloudH2O(cfr, qc, qvs,T,dz,entr, k1,k2,kts,kte)
                !
        IMPLICIT NONE
                !
        INTEGER, INTENT(IN):: k1,k2, kts,kte
        REAL, INTENT(IN):: entr
        REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, qvs, T, dz
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: qc
        REAL:: lwc, max_lwc, tdz, this_lwc, this_dz
        INTEGER:: k

        tdz = 0.
        do k = k1, k2
            tdz = tdz + dz(k)
        enddo
        max_lwc = ABS(qvs(k2)-qvs(k1))

        do k = k1, k2
            max_lwc = MAX(1.E-6, max_lwc - qc(k))
        enddo
        max_lwc = MIN(1.E-3, max_lwc)
        this_dz = 0.0
        do k = k1, k2
            if (k.eq.k1) then
                this_dz = this_dz + 0.5*dz(k)
            else
                this_dz = this_dz + dz(k)
            endif
            this_lwc = max_lwc*this_dz/tdz
            lwc = MAX(1.E-6, this_lwc*(1.-entr))
            if (cfr(k).gt.0.0.and.cfr(k).lt.1.0.and.T(k).ge.253.16) then
                qc(k) = qc(k) + cfr(k)*cfr(k)*lwc
            endif
        enddo

    END SUBROUTINE adjust_cloudH2O

    !+---+-----------------------------------------------------------------+

    !..Do not alter any grid-explicitly resolved hydrometeors, rather only
    !.. the supposed amounts due to the cloud fraction scheme.

    SUBROUTINE adjust_cloudFinal(cfr, qc, qi, Rho,dz, kts,kte)

        IMPLICIT NONE
                !
        INTEGER, INTENT(IN):: kts,kte
        REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, Rho, dz
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: qc, qi
        REAL:: lwp, iwp, xfac
        INTEGER:: k

        lwp = 0.
        iwp = 0.
        do k = kts, kte
            if (cfr(k).gt.0.0) then
                lwp = lwp + qc(k)*Rho(k)*dz(k)
                iwp = iwp + qi(k)*Rho(k)*dz(k)
            endif
        enddo

        if (lwp .gt. 1.5) then
            xfac = 1.5/lwp
            do k = kts, kte
                if (cfr(k).gt.0.0 .and. cfr(k).lt.1.0) then
                    qc(k) = qc(k)*xfac
                endif
            enddo
        endif

        if (iwp .gt. 1.5) then
            xfac = 1.5/iwp
                do k = kts, kte
                    if (cfr(k).gt.0.0 .and. cfr(k).lt.1.0) then
                        qi(k) = qi(k)*xfac
                    endif
                enddo
        endif

    END SUBROUTINE adjust_cloudFinal

end module mod_atm_utilities
