!>------------------------------------------------------------
!!  date / time module
!!
!!  Contains various utilities for working with dates and times.
!!  Utilities are capable of handling multiple model calendars.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module time
    implicit none
    ! define calendars
    integer, parameter :: GREGORIAN=0, NOLEAP=1, THREESIXTY=2
    integer, parameter :: YEAR_ZERO=1800  ! starting year for noleap and 360day calendars
    integer :: calendar
    integer, dimension(13) :: month_start
    private
    public :: time_init, date_to_mjd, calendar_date, calc_day_of_year, parse_date, calc_year_fraction
    public :: calendar
    public :: GREGORIAN, NOLEAP, THREESIXTY, YEAR_ZERO
contains

    subroutine time_init(calendar_name)
        implicit none
        character(len=*), intent(in) :: calendar_name
        integer :: i

        ! zero based month_starts (will have 1 added below)
        month_start=[0,31,59,90,120,151,181,212,243,273,304,334,365]

        if (trim(calendar_name)=="gregorian") then
            calendar=GREGORIAN
        else if (trim(calendar_name)=="standard") then
            calendar=GREGORIAN
        else if (trim(calendar_name)=="365-day") then
            calendar=NOLEAP
        else if (trim(calendar_name)=="noleap") then
            calendar=NOLEAP
        else if (trim(calendar_name)=="360-day") then
            calendar=THREESIXTY
        else
            write(*,*) "Unknown Calendar: ", trim(calendar_name)
            write(*,*) "Acceptable names = "
            write(*,*) "  gregorian, standard, 365-day, noleap, 360-day"
            write(*,*) " "
            stop
        endif

        if (calendar==THREESIXTY) then
            do i=0,12
                month_start(i+1)=i*30
            end do
        endif
        do i=0,12
            month_start(i+1)=month_start(i+1)+1
        end do


    end subroutine time_init

    !   algorithms from Wikipedia: http://en.wikipedia.org/wiki/Julian_day
    !   originally from Richards, E. G. (2013). Calendars. In S. E. Urban & P. K. Seidelmann, eds.
    !                   Explanatory Supplement to the Astronomical Almanac, 3rd ed. (pp. 585–624).
    !                   Mill Valley, Calif.: University Science Books. ISBN 978-1-89138-985-6
    !                   p617-9
    function date_to_mjd(year, month, day, hour, minute, second)
        implicit none
        integer, intent(in) :: year, month, day, hour, minute, second
        double precision :: date_to_mjd
        double precision :: d,m,y
        integer :: a,b

        if (calendar==GREGORIAN) then
            a = (14-month)/12
            y = year+4800-a
            m = month+12*a-3
            ! Gregorian calendar
            b = day + floor((153*m+2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045
            ! Julian calendar
            ! b = day + floor(153*m+2/5) + 365*y + floor(y/4) - 32083
            date_to_mjd = b + (((second/60d+0)+minute)/60d+0 + hour-12)/24.0 - 2400000.5
        else if (calendar==NOLEAP) then
            date_to_mjd = (year-YEAR_ZERO)*365 + month_start(month)-1 + day-1 + (hour + (minute+second/60d+0)/60d+0)/24d+0
        else if (calendar==THREESIXTY) then
            date_to_mjd = (year-YEAR_ZERO)*360 + month_start(month)-1 + day-1 + (hour + (minute+second/60d+0)/60d+0)/24d+0
        end if

    end function date_to_mjd

    ! compute the year, month, day, hour, minute, second corresponding
    ! to the input modified julian day (mjd)
    ! note mjd for NOLEAP and 360day calendars is not a true MJD
    ! arguably, seconds should be a real number, not an integer...
    subroutine calendar_date(inputmjd, year, month, day, hour, minute, second)
        implicit none
        double precision, intent(in) :: inputmjd
        integer, intent(out) :: year, month, day, hour, minute, second
        integer :: y=4716,j=1401,m=2,n=12,r=4,p=1461
        integer :: v=3,u=5,s=153,w=2,B=274277,C=-38
        integer ::f,e,g,h, jday
        double precision :: day_fraction,mjd
        mjd = inputmjd+1d-5 ! add less than one second
        if (calendar==GREGORIAN) then
            jday=nint(mjd+2400000.5)
            f=jday+j+(((4*jday+B)/146097)*3)/4+C
            e=r*f+v
            g=mod(e,p)/r
            h=u*g+w
            day=mod(h,s)/u+1
            month=mod(h/s+m,n)+1
            year=e/p-y+(n+m-month)/n
        else if (calendar==NOLEAP) then
            year=floor(mjd/365)
            day_fraction=mjd - year*365+1
            do f=1,12
                if (day_fraction>month_start(f)) then
                    month=f
                endif
            end do
            day = floor(day_fraction - month_start(month))+1
            year=year+YEAR_ZERO
        else if (calendar==THREESIXTY) then
            year=floor(mjd/360)
            day_fraction=mjd - year*360+1
            do f=1,12
                if (day_fraction>month_start(f)) then
                    month=f
                endif
            end do
            day = floor(day_fraction - month_start(month))+1
            year=year+YEAR_ZERO
        end if

        day_fraction=mod(mjd,1.0)
        hour=floor(day_fraction*24)

        day_fraction=day_fraction*24-hour
        minute=floor(day_fraction*60)

        day_fraction=day_fraction*60-minute
        second = nint((day_fraction-(24d0*60*1d-5))*60)

    end subroutine

    ! calculate the day of the year from a "modified julian day"
    ! note mjd for NOLEAP and 360day calendars is not a true MJD
    function calc_day_of_year(mjd)
        implicit none
        real :: calc_day_of_year
        double precision, intent(in) :: mjd

        integer :: year, month, day, hour, minute, second

        if (calendar==GREGORIAN) then
            call calendar_date(mjd,year, month, day, hour, minute, second)
            calc_day_of_year = mjd - date_to_mjd(year, 1,1,0,0,0)
        else if (calendar==NOLEAP) then
            calc_day_of_year = mod(mjd,365.0)
        else if (calendar==THREESIXTY) then
            calc_day_of_year = mod(mjd,360.0)
        endif
    end function calc_day_of_year


    function calc_year_fraction(mjd)
        implicit none
        real :: calc_year_fraction
        double precision, intent(in) :: mjd
        double precision :: year_start

        integer :: year, month, day, hour, minute, second

        if (calendar==GREGORIAN) then
            call calendar_date(mjd,year, month, day, hour, minute, second)
            year_start = date_to_mjd(year, 1,1,0,0,0)
            calc_year_fraction = (mjd - year_start) / (date_to_mjd(year+1, 1,1,0,0,0) - year_start)
        else if (calendar==NOLEAP) then
            calc_year_fraction = calc_day_of_year(mjd) / 365.0
        else if (calendar==THREESIXTY) then
            calc_year_fraction = calc_day_of_year(mjd) / 360.0
        endif
    end function calc_year_fraction

    ! convert an input date string in the form YYYY/MM/DD or YYYY/MM/DD hh:mm:ss
    ! into integer variables
    subroutine parse_date(date, year, month, day, hour, min, sec)
      implicit none
      character (len=*), intent(in) :: date
      integer, intent(out) :: sec, min, hour, day, month, year

      read(date(9:10),*) day
      read(date(6:7),*) month
      read(date(1:4),*) year
      if(len_trim(date) <= 11) then
         sec = 0
         min = 0
         hour = 0
      else
         read(date(18:19), *) sec
         read(date(15:16), *) min
         read(date(12:13), *) hour
      endif
    end subroutine parse_date


end module time
