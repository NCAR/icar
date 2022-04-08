!>------------------------------------------------------------
!!  date / time module
!!
!!  Defines a Time_type object
!!
!!  Contains various utilities for working with dates and times.
!!  Utilities are capable of handling multiple model calendars.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module time_object
    use time_delta_object, only : time_delta_t
    use iso_fortran_env, only: real128

    implicit none

    private

    integer, parameter :: MAXSTRINGLENGTH = 1024
    integer, parameter, public :: GREGORIAN=0, NOLEAP=1, THREESIXTY=2, NOCALENDAR=-1
    integer, parameter, public :: NON_VALID_YEAR = -9999

    !>------------------------------------------------------------
    !!  date / time Object
    !!
    !!  Time_type object stores date and time in a convenient object.
    !!  datetime is stored internally as both a modified julian day (or days since a known date)
    !!  and as a Y/M/D h:m:s integers
    !!
    !!------------------------------------------------------------
    type, public :: Time_type
        integer :: year_zero = 1800  ! reference date
        integer :: month_zero= 1
        integer :: day_zero  = 1
        integer :: hour_zero = 0
        integer :: calendar
        integer, dimension(13) :: month_start
        integer :: year, month, day, hour, minute, second

        real(real128) :: current_date_time = 0

      contains
        procedure, public  :: date        => calendar_date
        procedure, public  :: mjd         => get_mjd
        procedure, public  :: seconds     => get_seconds
        procedure, public  :: day_of_year => calc_day_of_year
        procedure, public  :: year_fraction=>calc_year_fraction
        procedure, public  :: date_to_mjd => date_to_mjd
        procedure, public  :: as_string   => as_string
        procedure, public  :: equals      => equals_with_precision
        procedure, public  :: units       => units
        procedure, public  :: broadcast   => bcast

        generic,   public  :: init        => time_init_c
        generic,   public  :: init        => time_init_i
        generic,   public  :: set         => set_from_string
        generic,   public  :: set         => set_from_date
        generic,   public  :: set         => set_from_mjd
        procedure, private :: time_init_c
        procedure, private :: time_init_i
        procedure, private :: set_from_string
        procedure, private :: set_from_date
        procedure, private :: set_from_mjd
        procedure, private :: set_calendar

        ! operator overloading to permit comparison of time objects as (t1 < t2), (t1 == t2) etc.
        procedure, private :: equal
        generic,   public  :: operator(==)   => equal
        procedure, private :: not_equal
        generic,   public  :: operator(/=)   => not_equal
        procedure, private :: greater_than
        generic,   public  :: operator(>)    => greater_than
        procedure, private :: greater_or_eq
        generic,   public  :: operator(>=)   => greater_or_eq
        procedure, private :: less_than
        generic,   public  :: operator(<)    => less_than
        procedure, private :: less_or_eq
        generic,   public  :: operator(<=)   => less_or_eq

        procedure, private :: addition
        generic,   public  :: operator(+)    => addition
        procedure, private :: difference
        generic,   public  :: operator(-)    => difference

    end type Time_type

interface

    !>------------------------------------------------------------
    !!  Initialize the time object
    !!
    !!  Set the object calendar and base year
    !!
    !!------------------------------------------------------------
    module subroutine time_init_c(this, calendar_name, year_zero, month_zero, day_zero, hour_zero)
        implicit none
        class(Time_type) :: this
        character(len=*), intent(in) :: calendar_name
        integer, intent(in), optional :: year_zero, month_zero, day_zero, hour_zero

    end subroutine time_init_c

    !>------------------------------------------------------------
    !!  Initialize the time object
    !!
    !!  Set the object calendar and base year
    !!
    !!------------------------------------------------------------
    module subroutine time_init_i(this, calendar, year_zero, month_zero, day_zero, hour_zero)
        implicit none
        class(Time_type) :: this
        integer, intent(in) :: calendar
        integer, intent(in), optional :: year_zero, month_zero, day_zero, hour_zero

    end subroutine time_init_i

    !>------------------------------------------------------------
    !!  Set the calendar from a given name
    !!
    !!  Looks at a calendar_name string to identify matches to known calendars
    !!  Known calendars are : 'gregorian', 'standard', '365-day', 'noleap', '360-day'
    !!  But only the first 5 characters are required.
    !!
    !!  Sets the object calendar attribute.
    !!
    !!------------------------------------------------------------
    module subroutine set_calendar(this, calendar_name)
        implicit none
        class(Time_type) :: this
        character(len=*), intent(in) :: calendar_name

    end subroutine set_calendar


    !>------------------------------------------------------------
    !!  Return the current date number (seconds since reference time)
    !!
    !!  For a gregorian calendar, if no year is specified, this will be
    !!  a modified Julian day (*86400).  For all other calendars it will be
    !!  seconds since startyear/startmonth/startday
    !!
    !!------------------------------------------------------------
    module function get_seconds(this) result(seconds)
        implicit none
        class(Time_type) :: this
        real(real128) :: seconds

    end function get_seconds


    !>------------------------------------------------------------
    !!  Return the current date number (days since initial year)
    !!
    !!  For a gregorian calendar, if no year is specified, this will be
    !!  a modified Julian day.  For all other calendars it will be days since 0-1-1
    !!
    !!------------------------------------------------------------
    module function get_mjd(this) result(mjd)
        implicit none
        class(Time_type) :: this
        real(real128) :: mjd

    end function get_mjd


    !>------------------------------------------------------------
    !!  Convert a Year, Month, Day, hour, minute, second into a single number
    !!
    !!  This number will be a Modified Julian Day for the default gregorian calendar
    !!  or it will be a number of days since whatever the initial year was
    !!
    !!------------------------------------------------------------
    module function date_to_mjd(this, year, month, day, hour, minute, second)
        implicit none
        class(Time_type), intent(in) :: this
        integer, intent(in) :: year, month, day, hour, minute, second
        real(real128) :: date_to_mjd

    end function date_to_mjd

    !>------------------------------------------------------------
    !!  Calculate the year, month, day, hour, minute second for the current date_time object
    !!
    !!  Y, M, D, h, m, s are modified on return
    !!  arguably, seconds should be a real number, not an integer...
    !!
    !!------------------------------------------------------------
    module subroutine calendar_date(this, year, month, day, hour, minute, second)
        implicit none
        class(Time_type), intent(in) :: this
        integer, intent(out) :: year, month, day, hour, minute, second

    end subroutine calendar_date


    !>------------------------------------------------------------
    !!  Return the day of the year corresponding to the current date_time
    !!
    !!  Calculate the day of the year from a "modified julian day" or other days since date
    !!
    !!------------------------------------------------------------
    module function calc_day_of_year(this, lon)
        implicit none
        real                        :: calc_day_of_year
        class(Time_type), intent(in):: this
        real,             intent(in), optional :: lon

    end function calc_day_of_year


    !>------------------------------------------------------------
    !!  Return the fractional way through the year corresponding to the current date_time
    !!
    !!  "Jan  1, 00:00" is 0
    !!  "Dec 31, 24:00" is 1.0
    !!
    !!------------------------------------------------------------
    module function calc_year_fraction(this, lon)
        implicit none
        real                        :: calc_year_fraction
        class(Time_type)            :: this
        real, intent(in), optional  :: lon

    end function calc_year_fraction



    !>------------------------------------------------------------
    !!  Set the current date based on an input string
    !!
    !!  Convert an input date string in the form YYYY/MM/DD or YYYY/MM/DD hh:mm:ss
    !!
    !!------------------------------------------------------------
    module subroutine set_from_string(this, date)
        implicit none
        class(Time_type), intent(inout) :: this
        character (len=*), intent(in) :: date

    end subroutine set_from_string

    !>------------------------------------------------------------
    !!  Set the current date based on an input set of integer date components
    !!
    !!  Set the datetime object to a given input date (Year, Month, Day, Hour, Minute, second)
    !!
    !!------------------------------------------------------------
    module subroutine set_from_date(this, year, month, day, hour, minute, second)
        implicit none
        class(Time_type), intent(inout) :: this
        integer, intent(in) :: year, month, day
        integer, intent(in), optional :: hour, minute, second

    end subroutine set_from_date

    !>------------------------------------------------------------
    !!  Set the current date based on an input day number (e.g. Modified Julian day)
    !!
    !!  Set the datetime object to a given input datetime number
    !!
    !!------------------------------------------------------------
    module subroutine set_from_mjd(this, days)
        implicit none
        class(Time_type), intent(inout) :: this
        real(real128), intent(in) :: days

    end subroutine set_from_mjd

    !>------------------------------------------------------------
    !!  Create a formated "units" string for this object
    !!
    !!  For example "days since 1858-11-17 00:00:00"
    !!
    !!------------------------------------------------------------
    module function units(this)
        implicit none
        class(Time_type), intent(in)   :: this
        character(len=MAXSTRINGLENGTH) :: units

    end function units

    module subroutine bcast(this, source, first_image, last_image)
        implicit none
        class(Time_type), intent(inout) :: this
        integer,          intent(in)    :: source
        integer,          intent(in)    :: first_image
        integer,          intent(in)    :: last_image

    end subroutine


    !>------------------------------------------------------------
    !!  Convert the date object into a string in the 0-filled format : "YYYY/MM/DD hh:mm:ss"
    !!
    !!------------------------------------------------------------
    module function as_string(this, input_format) result(pretty_string)
        implicit none
        class(Time_type), intent(in) :: this
        character(len=*), intent(in), optional :: input_format
        character(len=MAXSTRINGLENGTH) :: pretty_string

    end function as_string

    !>------------------------------------------------------------
    !!  Test that time 1 is greater than time 2
    !!
    !!  Used to implement the [>] operator
    !!  If the two input times have the same refernce calender, it just returns
    !!  > on the internal time representation between the two
    !!  Otherwise it retrieves the year, month, day, ... from both and checks for
    !!  > on all elements sequentially (year, then month, then day...)
    !!
    !!------------------------------------------------------------
    module function greater_than(t1, t2)
        implicit none
        class(Time_type), intent(in) :: t1, t2
        logical :: greater_than

    end function greater_than

    !>------------------------------------------------------------
    !!  Test that time 1 is greater than or equal to time 2
    !!
    !!  Used to implement the [>=] operator
    !!  If the two input times have the same refernce calender, it just returns
    !!  >= on the internal time representation between the two
    !!  Otherwise it retrieves the year, month, day, ... from both and checks for
    !!  >= on all elements sequentially (year, then month, then day...)
    !!
    !!------------------------------------------------------------
    module function greater_or_eq(t1, t2)
        implicit none
        class(Time_type), intent(in) :: t1, t2
        logical :: greater_or_eq

    end function greater_or_eq


    !>------------------------------------------------------------
    !!  Test that time 1 is equal to time 2
    !!
    !!  Used to implement the [==] operator
    !!  If the two input times have the same refernce calender, it just returns
    !!  == on the internal time representation between the two
    !!  Otherwise it retrieves the year, month, day, ... from both and checks for
    !!  equality on all elements
    !!
    !!------------------------------------------------------------
    module function equal(t1, t2)
        implicit none
        class(Time_type), intent(in) :: t1, t2
        logical :: equal

    end function equal


    !>------------------------------------------------------------
    !!  Test that time 1 is not equal to time 2
    !!
    !!  Used to implement the [/=] operator
    !!  Simply returns not ==
    !!
    !!------------------------------------------------------------
    module function not_equal(t1, t2)
        implicit none
        class(Time_type), intent(in) :: t1, t2
        logical :: not_equal

    end function not_equal


    !>------------------------------------------------------------
    !!  Test that time 1 is less than or equal to time 2
    !!
    !!  Used to implement the [<=] operator
    !!  Simply returns not >
    !!
    !!------------------------------------------------------------
    module function less_or_eq(t1, t2)
        implicit none
        class(Time_type), intent(in) :: t1, t2
        logical :: less_or_eq

    end function less_or_eq

    !>------------------------------------------------------------
    !!  Test that time 1 is less than time 2
    !!
    !!  Used to implement the [<] operator
    !!  Simply returns not >=
    !!
    !!------------------------------------------------------------
    module function less_than(t1, t2)
        implicit none
        class(Time_type), intent(in) :: t1, t2
        logical :: less_than

    end function less_than

    !>------------------------------------------------------------
    !!  Compare that two time objects are equal to within a specified precision
    !!
    !!  If precision is not specified, just returns the == operator result
    !!  If precision is specified, returns true if abs(t1-t2) <= precision
    !!
    !!------------------------------------------------------------
    module function equals_with_precision(t1,t2, precision) result(equals)
        implicit none
        class(Time_type),    intent(in) :: t1, t2
        type(time_delta_t),  intent(in), optional :: precision
        logical :: equals

    end function equals_with_precision


    !>------------------------------------------------------------
    !!  Subtract two times and return a time_delta object
    !!
    !!------------------------------------------------------------
    module function difference(t1, t2) result(dt)
        implicit none
        class(Time_type), intent(in) :: t1, t2
        type(time_delta_t) :: dt

    end function difference

    !>------------------------------------------------------------
    !!  Add a given time delta to a time object
    !!
    !!  returns a new time object
    !!
    !!------------------------------------------------------------
    module function addition(t1, dt) result(t2)
        implicit none
        class(Time_type),   intent(in) :: t1
        type(time_delta_t), intent(in) :: dt
        type(Time_type) :: t2

    end function addition

end interface
end module time_object
