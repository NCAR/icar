module time_io

    use data_structures
    use icar_constants
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use string,             only : get_integer
    use io_routines,        only : io_read, io_read_attribute
    use iso_fortran_env, only: real64, real128

    implicit none

contains

    function find_timestep_in_file(filename, time_var, time, time_at_step, precision, error) result(step)
        implicit none
        character(len=*),  intent(in) :: filename
        character(len=*),  intent(in) :: time_var
        type(Time_type),   intent(in) :: time
        type(Time_type),   intent(inout), optional :: time_at_step
        type(time_delta_t),intent(in),    optional :: precision
        integer,           intent(inout), optional :: error
        integer :: step

        type(Time_type), allocatable :: times_in_file(:)
        type(time_delta_t) :: max_dt
        integer :: i,n
        logical :: found

        call max_dt%set(seconds=1.0)
        if (present(precision)) max_dt = precision
        if (present(error)) error=0

        ! read the times for all timesteps in the specified file
        call read_times(filename, time_var, times_in_file)

        step = -1
        found= .False.
        n    = size(times_in_file)
        ! loop through times looking for a time that matches the input time to within
        ! a specified maximum delta t
        do i = 1, n
            if (.not.found) then
                if (times_in_file(i)%equals(time, precision=max_dt)) then
                    step = i
                    found=.True.
                elseif (times_in_file(i) > time) then
                    step = i-1
                    found=.True.
                endif
            endif
        enddo

        if (step < 1) then
            if (present(error)) then
                error = 1
            else
                write(*,*) "ERROR: Unable to find requested date in file."
                write(*,*) "Filename: ",trim(filename)
                write(*,*) "  time  : ",trim(time%as_string())
                write(*,*) "First time in file : ", trim(times_in_file(1)%as_string())
                write(*,*) " Last time in file : ", trim(times_in_file(n)%as_string())
                stop "Unable to find date in file"
            endif
        endif

        if (present(time_at_step)) time_at_step = times_in_file(step)

        deallocate(times_in_file)

    end function find_timestep_in_file


    function time_gain_from_units(units) result(gain)
        implicit none
        character(len=*), intent(in) :: units
        real(real128) :: gain

        if ((units(1:4)=="days").or.(units(1:4)=="Days")) then
            gain = 1.0Q0
        else if ((units(1:4)=="hour").or.(units(1:4)=="Hour")) then
            gain = 24.0Q0
        else if ((units(1:3)=="min").or.(units(1:3)=="Min")) then
            gain = 1440.0Q0
        else if ((units(1:3)=="sec").or.(units(1:3)=="Sec")) then
            gain = 86400.0Q0
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

        month_loc = index(units(since_loc:)," ") + 5
        month_loc = month_loc + since_loc

        month = get_integer(units(month_loc:month_loc+1))

    end function month_from_units

    function day_from_units(units) result(day)
        implicit none
        character(len=*), intent(in) :: units
        integer :: day

        integer :: since_loc, day_loc

        since_loc = index(units,"since")

        day_loc = index(units(since_loc:)," ") + 8
        day_loc = day_loc + since_loc

        day = get_integer(units(day_loc:day_loc+1))

    end function day_from_units

    function hour_from_units(units, error) result(hour)
        implicit none
        character(len=*), intent(in) :: units
        integer, intent(out), optional :: error
        integer :: hour

        integer :: since_loc, hour_loc
        if (present(error)) error = 0

        since_loc = index(units,"since")

        hour_loc = index(units(since_loc:)," ") + 11
        hour_loc = hour_loc + since_loc

        ! default return value if hours can't be read from the units attribute (e.g. they aren't present)
        hour = 0

        if( hour_loc+1 <= len(units) ) then
           if (trim(units(hour_loc:hour_loc+1)) /= "") then
               hour = get_integer(units(hour_loc:hour_loc+1))
           endif
        else
           if (present(error)) error = 1
        endif

    end function hour_from_units


    subroutine read_times(filename, varname, times, timezone_offset, curstep)
        implicit none
        character(len=*),   intent(in) :: filename, varname
        type(Time_type),    intent(inout), allocatable, dimension(:) :: times
        real(real128),      intent(in), optional :: timezone_offset
        integer,            intent(in), optional :: curstep

        real(real64),  allocatable, dimension(:) :: temp_times_64
        real(real128), allocatable, dimension(:) :: temp_times_128
        integer :: time_idx, error
        integer :: start_year, start_month, start_day, start_hour
        character(len=MAXSTRINGLENGTH) :: calendar, units
        real(real128) :: calendar_gain

        ! first read the time variable (presumebly a 1D real(real64) array)
        if (present(curstep)) then
            call io_read(trim(filename), trim(varname), temp_times_64, curstep=curstep)
        else
            call io_read(trim(filename), trim(varname), temp_times_64)
        endif

        ! attempt to read the calendar attribute from the time variable
        call io_read_attribute(trim(filename),"calendar", calendar, var_name=trim(varname), error=error)
        ! if time attribute it not present, set calendar to one specified in the config file
        if (error/=0) then
            if (this_image()==1) write(*,*) "WARNING: assuming standard/gregorian calendar for file "//trim(filename)
            calendar = "standard"
        endif

        ! attempt to read the units for this time variable
        call io_read_attribute(trim(filename), "units", units, var_name=trim(varname), error=error)

        ! if units attribute was present, then read information from it.
        if (error==0) then
            start_year    = year_from_units(units)
            start_month   = month_from_units(units)
            start_day     = day_from_units(units)
            start_hour    = hour_from_units(units)
            ! based off of the string "Days since" (or "seconds" or...)
            calendar_gain = time_gain_from_units(units)
        else
            stop "Time variable does not have units attribute"
        endif

        ! converts the input units to "days since ..."
        ! in case it is in units of e.g. "hours since" or "seconds since"
        allocate(temp_times_128(size(temp_times_64)))
        temp_times_128 = temp_times_64 / calendar_gain

        if (present(timezone_offset)) then
            temp_times_128 = temp_times_128 + timezone_offset / 24.0
        endif

        if (allocated(times)) deallocate(times)
        allocate(times(size(temp_times_128)))

        do time_idx = 1, size(temp_times_128,1)
            call times(time_idx)%init(calendar, start_year, start_month, start_day, start_hour)
            call times(time_idx)%set(days=temp_times_128(time_idx))
        end do

        deallocate(temp_times_64, temp_times_128)

    end subroutine read_times

end module time_io
