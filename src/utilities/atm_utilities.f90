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
    pure function calc_speed(u, v) result(speed)
        implicit none
        real, intent(in) :: u,v
        real :: speed

        speed = sqrt(u**2 + v**2)
    end function calc_speed

    !>----------------------------------------------------------
    !! Calculate the strength of the u wind field given a direction [0-2*pi] and magnitude
    !!
    !!----------------------------------------------------------
    pure function calc_u(direction, magnitude) result(u)
        implicit none
        real, intent(in) :: direction, magnitude
        real :: u

        u = sin(direction) * magnitude
    end function calc_u

    !>----------------------------------------------------------
    !! Calculate the strength of the v wind field given a direction [0-2*pi] and magnitude
    !!
    !!----------------------------------------------------------
    pure function calc_v(direction, magnitude) result(v)
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
    pure function calc_sat_lapse_rate(T,mr) result(sat_lapse)
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
    pure function calc_moist_stability(t_top, t_bot, z_top, z_bot, qv_top, qv_bot, qc) result(BV_freq)
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
    pure function calc_dry_stability(th_top, th_bot, z_top, z_bot) result(BV_freq)
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

    !> -------------------------------
    !!
    !! Compute exner function to convert potential_temperature to temperature
    !!
    !! -------------------------------
    elemental function exner_function(pressure) result(exner)
        implicit none
        real, intent(in) :: pressure
        real :: exner

        associate(po=>100000, Rd=>287.058, cp=>1003.5)
            exner = (pressure / po) ** (Rd/cp)

        end associate
    end function


    subroutine init_atm_utilities(options)
        implicit none
        type(options_t) :: options

        N_squared   = options%lt_options%N_squared
        variable_N  = options%lt_options%variable_N

        max_froude  = options%block_options%block_fr_max
        min_froude  = options%block_options%block_fr_min
        froude_gain = 1 / max(max_froude-min_froude, 0.001)

    end subroutine init_atm_utilities

end module mod_atm_utilities