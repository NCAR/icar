module time
	implicit none
	private
	public :: date_to_mjd, calendar_date, parse_date
contains
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

		a = (14-month)/12
		y = year+4800-a
		m = month+12*a-3
		! Gregorian calendar
		b = day + floor((153*m+2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045
		! Julian calendar
! 		b = day + floor(153*m+2/5) + 365*y + floor(y/4) - 32083
		date_to_mjd = b + (((second/60d+0)+minute)/60d+0 + hour-12)/24.0 - 2400000.5

	end function date_to_mjd

	subroutine calendar_date(mjd, year, month, day, hour, minute, second)
		implicit none
		double precision, intent(in) :: mjd
		integer, intent(out) :: year, month, day, hour, minute, second
		integer :: y=4716,j=1401,m=2,n=12,r=4,p=1461
		integer :: v=3,u=5,s=153,w=2,B=274277,C=-38
		integer ::f,e,g,h, jday
		double precision :: day_fraction
		
		jday=nint(mjd+2400000.5)
		f=jday+j+(((4*jday+B)/146097)*3)/4+C
		e=r*f+v
		g=mod(e,p)/r
		h=u*g+w
		day=mod(h,s)/u+1
		month=mod(h/s+m,n)+1
		year=e/p-y+(n+m-month)/n
		
		day_fraction=mod(mjd,1.0)
		hour=floor(day_fraction*24+1e-5)
		
		day_fraction=day_fraction*24-hour
		minute=floor(day_fraction*60+1e-2)
		
		day_fraction=day_fraction*60-minute
		second=nint(day_fraction*60)
		
	end subroutine


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
