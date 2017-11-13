!>------------------------------------------------------------
!!  date / time delta module
!!
!!  Defines a time_delta_t object
!!
!!  Permits subtraction of two time objects resulting in a time_delta object
!!  time_delta can be accessed in days, hours, minutes, or seconds
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module time_delta_object
    implicit none
    private

    type, public :: time_delta_t
        private
        double precision :: delta = 0
    contains
        procedure, public :: days
        procedure, public :: hours
        procedure, public :: minutes
        procedure, public :: seconds

        procedure, public :: set_time_delta
    end type time_delta_t

contains

    ! below this is the definition of the time object methods
    subroutine set_time_delta(self, days, hours, minutes, seconds)
        implicit none
        class(time_delta_t) :: self
        double precision, intent(in), optional :: days, hours, minutes, seconds
        double precision :: dt

        dt=0
        if (present(days)) dt = dt + days
        if (present(hours)) dt = dt + hours / 24.0
        if (present(minutes)) dt = dt + minutes / 1440.0
        if (present(seconds)) dt = dt + seconds / 86400.0

        self%delta = dt

    end subroutine set_time_delta


    function days(self)
        implicit none
        class(time_delta_t), intent(in) :: self
        double precision :: days

        days = self%delta

    end function days

    function hours(self)
        implicit none
        class(time_delta_t), intent(in) :: self
        double precision :: hours

        hours = self%delta * 24

    end function hours

    function minutes(self)
        implicit none
        class(time_delta_t), intent(in) :: self
        double precision :: minutes

        minutes = self%delta * 1440

    end function minutes

    function seconds(self)
        implicit none
        class(time_delta_t), intent(in) :: self
        double precision :: seconds

        seconds = self%delta * 86400

    end function seconds

end module time_delta_object
