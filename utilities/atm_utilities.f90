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

    use icar_constants, only : pi, gravity, Rd, Rw, cp, LH_vaporization
    use data_structures

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


    subroutine init_atm_utilities(options)
        implicit none
        type(options_type) :: options

        N_squared   = options%lt_options%N_squared
        variable_N  = options%lt_options%variable_N

        max_froude  = options%block_options%block_fr_max
        min_froude  = options%block_options%block_fr_min
        froude_gain = 1 / max(max_froude-min_froude, 0.001)

    end subroutine init_atm_utilities

end module mod_atm_utilities
