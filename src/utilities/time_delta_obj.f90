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
    use icar_constants, only : MAXSTRINGLENGTH
    use iso_fortran_env, only: real64

    implicit none
    private

    type, public :: time_delta_t
        private
        real(real64) :: delta = 0
    contains
        procedure, public :: days
        procedure, public :: hours
        procedure, public :: minutes
        procedure, public :: seconds
        procedure, public :: as_string
        generic,   public :: set         => set_time_delta_d
        generic,   public :: set         => set_time_delta_f
        generic,   public :: set         => set_time_delta_i

        procedure, private :: set_time_delta_d
        procedure, private :: set_time_delta_f
        procedure, private :: set_time_delta_i
    end type time_delta_t

contains
    ! below this is the definition of the time delta object methods

    subroutine set_time_delta_d(self, seconds, days, hours, minutes)
        implicit none
        class(time_delta_t) :: self
        real(real64), intent(in) :: seconds
        real(real64), intent(in), optional :: days, hours, minutes
        real(real64) :: dt

        dt=0
        if (present(days)) dt = dt + days
        if (present(hours)) dt = dt + hours / 24.0
        if (present(minutes)) dt = dt + minutes / 1440.0

        dt = dt + seconds / 86400.0

        self%delta = dt

    end subroutine set_time_delta_d

    subroutine set_time_delta_f(self, seconds, days, hours, minutes)
        implicit none
        class(time_delta_t) :: self
        real, intent(in) :: seconds
        real, intent(in), optional :: days, hours, minutes
        real(real64) :: dt

        dt=0
        if (present(days)) dt = dt + days
        if (present(hours)) dt = dt + hours / 24.0D0
        if (present(minutes)) dt = dt + minutes / 1440.0D0

        dt = dt + seconds / 86400.0D0

        self%delta = dt

    end subroutine set_time_delta_f

    subroutine set_time_delta_i(self, seconds, days, hours, minutes)
        implicit none
        class(time_delta_t) :: self
        integer, intent(in) :: seconds
        integer, intent(in), optional :: days, hours, minutes
        real(real64) :: dt

        dt=0
        if (present(days)) dt = dt + days
        if (present(hours)) dt = dt + hours / 24.0D0
        if (present(minutes)) dt = dt + minutes / 1440.0D0

        dt = dt + seconds / 86400.0D0

        self%delta = dt

    end subroutine set_time_delta_i


    function days(self)
        implicit none
        class(time_delta_t), intent(in) :: self
        real(real64) :: days

        days = self%delta

    end function days

    function hours(self)
        implicit none
        class(time_delta_t), intent(in) :: self
        real(real64) :: hours

        hours = self%delta * 24

    end function hours

    function minutes(self)
        implicit none
        class(time_delta_t), intent(in) :: self
        real(real64) :: minutes

        minutes = self%delta * 1440

    end function minutes

    function seconds(self)
        implicit none
        class(time_delta_t), intent(in) :: self
        real(real64) :: seconds

        seconds = self%delta * 86400

    end function seconds

    function as_string(self) result(pretty_string)
        implicit none
        class(time_delta_t), intent(in) :: self
        character(len=MAXSTRINGLENGTH)  :: pretty_string

        if (abs(self%seconds()) < 1) then
            write(pretty_string,"(F8.4,A)") self%seconds(), " seconds"
        elseif (abs(self%seconds()) <= 60) then
            write(pretty_string,"(F6.2,A)") self%seconds(), " seconds"
        elseif (abs(self%minutes()) <= 60) then
            write(pretty_string,"(F6.2,A)") self%minutes(), " minutes"
        elseif (abs(self%hours()) <= 24) then
            write(pretty_string,"(F6.2,A)") self%hours(), " hours"
        elseif (abs(self%days()) <= 365) then
            write(pretty_string,"(F10.2,A)") self%days(), " days"
        else
            write(pretty_string,"(F10.2,A,F6.1,A)") self%days(), " days (roughly ",self%days()/365.25, " years)"
        endif

    end function as_string

end module time_delta_object
