!>------------------------------------------------------------
!!  Module to test the calendar routines
!!
!!  Loops through 2100? years checking that conversions two and from an
!!  internal date-time (e.g. Modified Julian Day) and YMD hms are consistent
!!  test multiple calendars
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module calendar_test_module
    use time
    integer, parameter :: STRING_LENGTH = 255
    real, parameter :: MAX_ERROR = 1e-5 ! allow less than 1 second error (over a 2100 yr period)
contains
    logical function calendar_test(calendar_name,error)
        character(len=STRING_LENGTH), intent(in) :: calendar_name
        character(len=STRING_LENGTH), intent(out) :: error

        double precision :: mjd_input, mjd_output
        double precision :: min_mjd, max_mjd, mjd_step
        integer :: year, month, day, hour, minute, second

        error=""
        calendar_test=.True.

        call time_init(calendar_name)
        min_mjd=365.0*1.d0
        max_mjd=365.0*2100.d0
        mjd_step=0.1

        MJDLOOP: do mjd_input = min_mjd, max_mjd, mjd_step
            ! test that input and output Modified Julian Days stay the same
            call calendar_date(mjd_input, year, month, day, hour, minute, second)
            mjd_output=date_to_mjd(year, month, day, hour, minute, second)

            ! test that month and day values are at least realistic (rare to fail)
            if ((day<1).or.(day>31)) then
                calendar_test=.False.
                write(error,'(6I6, " Error:",f10.6)') year, month, day, hour, minute, second, (mjd_output-mjd_input)
                print*, "Testing range:"
                print*, "   min_mjd, max_mjd, mjd_step"
                print*, min_mjd, max_mjd, mjd_step
                print*, "   Corresponding First Date"
                call calendar_date(min_mjd, year, month, day, hour, minute, second)
                print*, year, month, day, hour, minute, second
                print*, "   Corresponding Last Date"
                call calendar_date(max_mjd, year, month, day, hour, minute, second)
                print*, year, month, day, hour, minute, second
                print*, " "
                EXIT MJDLOOP
            endif
            if ((month<1).or.(month>12)) then
                calendar_test=.False.
                write(error,'(6I6, " Error:",f10.6)') year, month, day, hour, minute, second, (mjd_output-mjd_input)
                print*, "Testing range:"
                print*, "   min_mjd, max_mjd, mjd_step"
                print*, min_mjd, max_mjd, mjd_step
                print*, "   Corresponding First Date"
                call calendar_date(min_mjd, year, month, day, hour, minute, second)
                print*, year, month, day, hour, minute, second
                print*, "   Corresponding Last Date"
                call calendar_date(max_mjd, year, month, day, hour, minute, second)
                print*, year, month, day, hour, minute, second
                print*, " "
                EXIT MJDLOOP
            endif

            ! this is the main test it is likely to fail
            if (abs(mjd_output-mjd_input)>MAX_ERROR) then
                calendar_test=.False.
                write(error,'(6I6, " Error:",f10.6)') year, month, day, hour, minute, second, (mjd_output-mjd_input)
                print*, "Testing range:"
                print*, "   min_mjd, max_mjd, mjd_step"
                print*, min_mjd, max_mjd, mjd_step
                print*, "   Corresponding First Date"
                call calendar_date(min_mjd, year, month, day, hour, minute, second)
                print*, year, month, day, hour, minute, second
                print*, "   Corresponding Last Date"
                call calendar_date(max_mjd, year, month, day, hour, minute, second)
                print*, year, month, day, hour, minute, second
                print*, " "
                EXIT MJDLOOP
            endif
        end do MJDLOOP

    end function calendar_test

    subroutine detailed_tests(calendar_name)
        character(len=STRING_LENGTH), intent(in) :: calendar_name

        double precision :: mjd_input, mjd_output
        double precision :: min_mjd, max_mjd, mjd_step
        integer :: year, month, day, hour, minute, second

        call time_init(calendar_name)
        min_mjd=365.0*1.d0
        max_mjd=365.0*2.d0
        mjd_step=1
        print*, "Detailed testoutput:"
        print*, "Testing range:"
        print*, "   min_mjd, max_mjd, mjd_step"
        print*, min_mjd, max_mjd, mjd_step
        print*, "   Corresponding First Date"
        call calendar_date(min_mjd, year, month, day, hour, minute, second)
        print*, year, month, day, hour, minute, second
        print*, "   Corresponding Last Date"
        call calendar_date(max_mjd, year, month, day, hour, minute, second)
        print*, year, month, day, hour, minute, second

        MJDLOOP: do mjd_input = min_mjd, max_mjd, mjd_step
            ! test that input and output Modified Julian Days stay the same
            call calendar_date(mjd_input, year, month, day, hour, minute, second)
            mjd_output=date_to_mjd(year, month, day, hour, minute, second)
            print*, mjd_input, mjd_output, mjd_output - mjd_input
            print*, year, month, day, hour, minute, second
        end do MJDLOOP


    end subroutine

end module calendar_test_module

program test_calendar
    use calendar_test_module
    integer, parameter :: NCALENDARS=5
    character(len=STRING_LENGTH),dimension(NCALENDARS) :: calendars_to_test
    character(len=STRING_LENGTH) :: calendar_error
    character(len=STRING_LENGTH) :: options

    integer :: current_calendar
    integer :: error
    logical :: file_exists

    calendars_to_test=[character(len=STRING_LENGTH) :: "gregorian","standard","365-day","noleap","360-day"]
    options=""

    if (command_argument_count()>0) then
        call get_command_argument(1,options, status=error)
    endif

    if (trim(options)=="help") then
        print*, "calendar_test [detailed] [help]"
        print*, "   detailed: print an entire year of dates for visual inspection"
        print*, "   help: print this message"
        stop
    endif

    do current_calendar=1,NCALENDARS
        print*, "Testing : ", trim(calendars_to_test(current_calendar))
        if (calendar_test(calendars_to_test(current_calendar),error=calendar_error)) then
            print*, " PASSED "
            if (trim(options)=="detailed") then
                call detailed_tests(calendars_to_test(current_calendar))
            endif
        else
            print*, " FAILED ", trim(calendars_to_test(current_calendar))
            print*, trim(calendar_error)
            call detailed_tests(calendars_to_test(current_calendar))
        endif
    end do
end program test_calendar
