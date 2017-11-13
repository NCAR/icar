module time_io

    use data_structures
    use icar_constants
    use time_object, only : Time_type
    use string,      only : get_integer
    use io_routines, only : io_read, io_read_attribute

    implicit none

contains

    function time_gain_from_units(units) result(gain)
        implicit none
        character(len=*), intent(in) :: units
        double precision :: gain

        if ((units(1:4)=="days").or.(units(1:4)=="Days")) then
            gain = 1.0D0
        else if ((units(1:4)=="hour").or.(units(1:4)=="Hour")) then
            gain = 24.0D0
        else if ((units(1:3)=="min").or.(units(1:3)=="Min")) then
            gain = 1440.0D0
        else if ((units(1:3)=="sec").or.(units(1:3)=="Sec")) then
            gain = 86400.0D0
        else
            write(*,*) trim(units)
            stop "Error: unknown units"
        endif

    end function time_gain_from_units

    function year_from_units(units) result(year)
        implicit none
        character(len=*), intent(in) :: units
        integer :: year

        integer :: since_loc, year_loc

        since_loc = index(units,"since")

        year_loc = index(units(since_loc:)," ")
        year_loc = year_loc+since_loc

        year = get_integer(units(year_loc:year_loc+3))

    end function year_from_units

    function month_from_units(units) result(month)
        implicit none
        character(len=*), intent(in) :: units
        integer :: month

        integer :: since_loc, month_loc

        since_loc = index(units,"since")

        month_loc = index(units(since_loc:)," ") + 4
        month_loc = month_loc + since_loc

        month = get_integer(units(month_loc:month_loc+1))

    end function month_from_units

    function day_from_units(units) result(day)
        implicit none
        character(len=*), intent(in) :: units
        integer :: day

        integer :: since_loc, day_loc

        since_loc = index(units,"since")

        day_loc = index(units(since_loc:)," ") + 7
        day_loc = day_loc + since_loc

        day = get_integer(units(day_loc:day_loc+1))

    end function day_from_units

    subroutine read_times(filename, varname, times, timezone_offset)
        implicit none
        character(len=*),   intent(in) :: filename, varname
        type(Time_type),    intent(inout), dimension(:) :: times
        double precision, optional :: timezone_offset

        double precision, allocatable, dimension(:) :: temp_times
        integer :: time_idx, error
        integer :: start_year, start_month, start_day
        character(len=MAXSTRINGLENGTH) :: calendar, units
        double precision :: calendar_gain

        ! first read the time variable (presumebly a 1D double precision array)
        call io_read(filename, varname, temp_times)

        ! attempt to read the calendar attribute from the time variable
        call io_read_attribute(filename,"calendar", calendar, var_name=varname, error=error)
        ! if time attribute it not present, set calendar to one specified in the config file
        if (error/=0) then
            print*, "WARNING: assuming standard/gregorian calendar for file "//trim(filename)
            calendar = "standard"
        endif

        ! attempt to read the units for this time variable
        call io_read_attribute(filename, "units", units, var_name=varname, error=error)

        ! if units attribute was present, then read information from it.
        if (error==0) then
            start_year    = year_from_units(units)
            start_month   = month_from_units(units)
            start_day     = day_from_units(units)
            ! based off of the string "Days since" (or "seconds" or...)
            calendar_gain = time_gain_from_units(units)
        else
            stop "Time variable does not have units attribute"
        endif

        ! puts the units to seconds since ...
        ! in case it is in units of e.g. "hours since" or "days since"
        temp_times = temp_times * calendar_gain
        if (present(timezone_offset)) then
            temp_times = temp_times + timezone_offset / 24.0
        endif

        do time_idx = 1, size(temp_times,1)

            call times(time_idx)%init(calendar, start_year, start_month, start_day)
            call times(time_idx)%set(days=temp_times(time_idx))

        end do

        deallocate(temp_times)

    end subroutine read_times

end module time_io
