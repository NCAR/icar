!>------------------------------------------------------------
!! Module provides geographic interpolation procedures
!! Various functions used for spatial interpolation from low-res
!! forcing grid to high-res model grid.
!!
!! <pre>
!!  Entry points:
!!      geo_LUT     : creates a geographic look uptable to convert one
!!                    grid to another
!!      geo_interp  : interpolates from one grid to another and
!!                    loops over the third dimension
!!      geo_interp2d: interpolates from one grid to another
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module geo
    use data_structures
    use icar_constants, only : kMAINTAIN_LON, kDATELINE_CENTERED, kPRIME_CENTERED, kGUESS_LON

    implicit none

    private

    public::point_in_poly
    public::point_is_on_line
    public::geo_LUT      ! Create a geographic Look up table
    public::geo_interp   ! apply geoLUT to interpolate in 2d for a 3d grid
    public::geo_interp2d ! apply geoLUT to interpolate in 2d for a 2d grid
    public::standardize_coordinates

contains
    !>------------------------------------------------------------
    !! Calculate the weights to use for bilinear interpolation between surrounding points x, y
    !!
    !! Note that this routine should only be used when your input forcing grid is regular in lat / lon
    !! If it is irregular (e.g. a lambert conformal conic projection) then use the triangulation routine
    !! This routine will end up with points located incorrectly and can create artifacts in edge cases.
    !!
    !! Takes a point (yi,xi) and a set of 4 surrounding points (y[4],x[4]) as input
    !!
    !! @param   real    yi  input y location to interpolate to.
    !! @param   real    y   Array of surrounding y-locations to interpolate from.
    !! @param   real    xi  input x location to interpolate to.
    !! @param   real    x   Array of surrounding x-locations to interpolate from.
    !! @retval  real    weights
    !!
    !!------------------------------------------------------------
    function bilin_weights(yi,y,xi,x)
        implicit none
        real,intent(in)::yi,y(0:3),xi,x(0:3)
        real::x0,x1,x2,x3,y5,y6,f1,f2
        real, dimension(4) ::bilin_weights

        ! handle the special case if x(1) and x(0) are identical
        if ((x(1)-x(0))==0) then
            x0=1
        else
            ! compute the linear interpolation between x0 and x1
            x0=abs((xi-x(0))/(x(1)-x(0)))
        endif
        ! the weights for the other x are just the complement
        x1=1-x0

        ! handle the special case if x(3) and x(2) are identical
        if ((x(3)-x(2))==0) then
            x2=1
        else
            ! compute the linear interpolation between x2 and x3
            x2=abs((xi-x(2))/(x(3)-x(2)))
        endif
        ! the weights for the other x are just the complement
        x3=1-x2

        ! now compute the y interpolation weights
        ! first find the y locations after performing the x interpolation
        y5=y(0)*x1+y(1)*x0
        y6=y(2)*x3+y(3)*x2
        ! now find the interpolation weights between those interpolated y points
        if ((y6-y5)==0) then
            ! the special case again
            f1=1
        else
            ! the standard linear interpolation case
            f1=(yi-y5)/(y6-y5)
        endif
        ! the complement
        f2=1-f1

        ! note that the final weights are a mixture of the x and y weights
        bilin_weights=(/x1*f2,x0*f2,x3*f1,x2*f1/)
    end function bilin_weights

    !>------------------------------------------------------------
    !! Calculate the weights to use for triangular interpolation between surrounding points x, y
    !!
    !! Takes a point (yi,xi) and a set of 4 surrounding points (y[4],x[4]) as input
    !! Of the four surrounding points, the first two should be verticies of the triangle
    !! and the last vertex will be the average of all four points.
    !!
    !! Based on Barycentric coordinates: https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    !!  see also explanation here: https://codeplea.com/triangular-interpolation
    !!
    !! @param   real    yi  input y location to interpolate to.
    !! @param   real    y   Array of surrounding y-locations to interpolate from.
    !! @param   real    xi  input x location to interpolate to.
    !! @param   real    x   Array of surrounding x-locations to interpolate from.
    !! @retval  real    weights
    !!
    !!------------------------------------------------------------
    function tri_weights(yi,y,xi,x)
        implicit none
        real,intent(in)::yi,y(4),xi,x(4)
        real::x0,y0
        real, dimension(4) ::tri_weights

        real :: w1,w2,w3,denom

        x0 = sum(x)/4
        y0 = sum(y)/4

        associate(y1 => y(1), y2 => (y(2)), y3 => y0, &
                  x1 => x(1), x2 => (x(2)), x3 => x0)

        denom = ((y2-y3) * (x1-x3) + (x3-x2) * (y1-y3))
        if (denom==0) then
            write(*,*) "ERROR: Triangular interpolation broken"
            write(*,*) xi, yi
            write(*,*) "Triangle vertices"
            write(*,*) x1, y1
            write(*,*) x2, y2
            write(*,*) x3, y3
            error stop "Denominator in triangulation is broken"
        endif

        ! This is the core of the algorithm weights for Barycentric coordinates.
        ! compute weights for each x,y point
        w1 = ((y2-y3) * (xi-x3) + (x3-x2) * (yi-y3)) &
                           / denom

        w2 = ((y3-y1) * (xi-x3) + (x1-x3) * (yi-y3)) &
                           / denom

        w3 = 1 - w1 - w2

        if (minval([w1,w2,w3]) < -1e-2) then
            write(*,*) "ERROR: Point not located in bounding triangle"
            write(*,*) xi, yi
            write(*,*) "Triangle vertices"
            write(*,*) x1, y1
            write(*,*) x2, y2
            write(*,*) x3, y3
            write(*,*) w1, w2, w3
            error stop "Triangulation is broken"
        endif
        w1=max(0.,w1); w2=max(0.,w2); w3=max(0.,w3);

        if (abs((w1+w2+w3)-1)>1e-2) then
            write(*,*) "Error, w1+w2+w3 != 1"
            write(*,*) w1,w2,w3
            write(*,*) "Point"
            write(*,*) xi, yi
            write(*,*) "Triangle vertices"
            write(*,*) x1, y1
            write(*,*) x2, y2
            write(*,*) x3, y3
        endif
        denom = 1.0 / (w1+w2+w3)
        w1 = w1 * denom
        w2 = w2 * denom
        w3 = w3 * denom

        end associate

        tri_weights = [w1, w2, w3, 0.0]

    end function tri_weights

    !>------------------------------------------------------------
    !! Calculate the weights to use for inverse distance interpolation between surrounding points x, y
    !!
    !! Takes a point (yi,xi) and a set of 4 surrounding points (y[4],x[4]) as input
    !!
    !! @param   real    yi  input y location to interpolate to.
    !! @param   real    y   Array of surrounding y-locations to interpolate from.
    !! @param   real    xi  input x location to interpolate to.
    !! @param   real    x   Array of surrounding x-locations to interpolate from.
    !! @retval  real    weights
    !!
    !!------------------------------------------------------------
    function idw_weights(yi,y,xi,x)
        implicit none
        real, intent(in)::yi, y(:), xi, x(:)

        real, allocatable, dimension(:) :: idw_weights
        real, allocatable, dimension(:) :: distances
        integer :: n
        real    :: total_inverse_distance

        n = size(x)
        allocate(distances(n))
        allocate(idw_weights(n))

        distances = 1.0/sqrt((y-yi)**2 + (x-xi)**2)

        total_inverse_distance = sum(distances)

        idw_weights = distances / total_inverse_distance

    end function idw_weights

    !>------------------------------------------------------------
    !! Find the minimum x step
    !!
    !!------------------------------------------------------------
    integer function minxw(xw,longrid,xpos,ypos,lon)
        implicit none
        integer, intent(in) :: xw, xpos, ypos
        real,    intent(in) :: longrid(:,:), lon
        real    :: curdist, dx
        integer :: ime

        ime = ubound(longrid,1)

        curdist = abs(lon - longrid(xpos,ypos))
        if (xpos < ime) then
            dx = abs(longrid(xpos,ypos) - longrid(xpos+1,ypos))
        else
            dx = abs(longrid(xpos,ypos) - longrid(xpos-1,ypos))
        endif
        minxw = max(xw, ceiling(curdist/dx))
        if (minxw > xw) then
            if (minxw > 4) then
                minxw = ceiling(minxw*0.8)
            endif
        endif
    end function minxw

    !>------------------------------------------------------------
    !!  Find the minimum next step size in the y direction
    !!  yw = minyw(yw,lo%lat,xc,yc,lat)
    !!
    !!------------------------------------------------------------
    integer function minyw(yw,latgrid,xpos,ypos,lat)
        integer, intent(in) :: yw, xpos, ypos
        real,    intent(in) :: latgrid(:,:), lat
        real    :: curdist, dx
        integer :: jme

        jme = ubound(latgrid, 2)

        curdist = abs(lat - latgrid(xpos,ypos))
        if (ypos < jme) then
            dx = abs(latgrid(xpos, ypos) - latgrid(xpos, ypos+1))
        else
            dx = abs(latgrid(xpos,ypos)-latgrid(xpos,ypos-1))
        endif
        minyw = max(yw, ceiling(curdist/dx))
        if (minyw > yw) then
            if (minyw > 4) then
                minyw = ceiling(minyw*0.8)
            endif
        endif
    end function minyw

    !>------------------------------------------------------------
    !! Find a location lat,lon in lo%lat,lon grids
    !!  Assumes that the lat/lon grids are semi-regular (dx/dy are not constant but they are nearly so)
    !!  Calculates dx/dy at the middle of the lat/lon grids, then calculates the location of lat,lon
    !!  input point using those coordinates.
    !!  Next proceeds to calculate a new position based on the dx/dy at that position until the new
    !!  position is within 1 gridcell of the current position
    !!  Once that "within 1" position is found, search all cells in a 3x3 grid for the minimum distance
    !!  return a position datatype that includes the found x/y location
    !!
    !! currently tries up to three methods to find the location
    !!   1) assume relatively even dxdy and compute location (iterate a <20times)
    !!   2) use a log(n) search (divide and conquer) also iterate <20times)
    !!   a) if 1 or 2 are approximately successful:
    !!       search a small region around the "best" point to find the real best point
    !!   3) in the diabolical case where 1 and 2 fail (rare) just search every single location (n^2)
    !!   Could add a 2.5) downhill search algorithm...
    !!
    !!  @param  lo              low resolution interpolable object to search for new position
    !!  @param  lat             latitude of new position to find
    !!  @param  lon             longitude of new position to find
    !!  @param  lastpos         previously found position to begin search from
    !!  @retval find_location   position that most closely matches input lat/lon
    !!
    !!------------------------------------------------------------
    type(position) function find_location(lo,lat,lon,lastpos)
        implicit none
        type(interpolable_type),intent(inout)::lo      ! input interpolable object to search
        real,intent(in)::lat,lon                        ! input position to find
        type(position),intent(in)::lastpos              ! position to start search from
        ! locals
        real::mindist, curdist,dx,dy,xsign,ysign,x,y    ! temporary variables to use during the search
        integer::nx,ny,xc,yc,xw,yw,iterations,xstep,ystep   ! more temporary variables
        integer :: ims,ime,jms,jme

        nx = size(lo%lat, 1)
        ny = size(lo%lat, 2)

        ims = lbound(lo%lat,1)
        ime = ubound(lo%lat,1)
        jms = lbound(lo%lat,2)
        jme = ubound(lo%lat,2)

        ! calculate dx/dy of the grid
        dx = lo%lon(ims+1,jms) - lo%lon(ims,jms)
        dy = lo%lat(ims,jms+1) - lo%lat(ims,jms)

        ! Handle weird edge case for WRF ideal simulations
        if (dx == 0) then
            if (.not.lo%dx_errors_printed) then
                write(*,*) "ERROR : geo_find_location : DX = 0 !  Check your inputfile Longitude data"
                write(*,*) "  Attempting to continue by assuming the grids are the same or can wrap arround"
                lo%dx_errors_printed=.True.
            endif
            find_location%x = -1
            find_location%y = -1
            return
        endif
        if (dy == 0) then
            if (.not.lo%dy_errors_printed) then
                write(*,*) "ERROR : geo_find_location : DY = 0 !  Check your inputfile Latitude data"
                write(*,*) "  Attempting to continue by assuming the grids are the same or can wrap arround"
                lo%dy_errors_printed=.True.
            endif
            find_location%x = -1
            find_location%y = -1
            return
        endif
        ! current/starting position = the middle of the grid this is the default that should be given
        !       xc=nx/2
        !       yc=ny/2
        xc = lastpos%x
        yc = lastpos%y

        x = lo%lon(xc, yc)
        y = lo%lat(xc, yc)

        ! steps to take = the difference the between the current point and the input point / dx
        xstep = (lon-x)/dx
        ystep = (lat-y)/dy

        ! CASE 1 assume a quasi regular dx/dy and calculate new location
        !  while we need to step by more than one grid cell, iterate
        !  if the grid is highly regular, we will only iterate 1-2x, highly irregular might require more iterations
        !  in the diabolical case, this could fail?
        !  in most cases this succeeds with a few iterations.
        iterations = 0
        do while (((abs(xstep)>1).or.(abs(ystep)>1)).and.iterations<20)
            iterations = iterations+1
            ! update the current x/y locations
            ! force it to be <nx-1 so we can calculate dx from (xc+1)-xc
            xc = max(ims,min(ime-1,xc+xstep))
            yc = max(jms,min(jme-1,yc+ystep))
            x = lo%lon(xc,yc)
            y = lo%lat(xc,yc)
            ! calculate a new dx/dy and x/y step
            if (xc<nx) then
                dx = lo%lon(xc+1,yc) - lo%lon(xc,yc)
            else
                dx = lo%lon(xc,yc) - lo%lon(xc-1,yc)
            endif
            if(yc<ny) then
                dy = lo%lat(xc,yc+1) - lo%lat(xc,yc)
            else
                dy = lo%lat(xc,yc) - lo%lat(xc,yc-1)
            endif
            xstep = NINT( (lon-x) / dx )
            ystep = NINT( (lat-y) / dy )
        enddo
        ! because one or both steps could actually be 1...
        ! this is deliberate so we can find the edge of the array if necessary
        xc = max(ims, min(ime, xc+xstep ) )
        yc = max(jms, min(jme, yc+ystep ) )

        ! CASE 2 use a log(n search)
        ! in case we hit some pathologically varying dx case
        ! use a "straightforward" log(n) search
        if (iterations >= 20) then
            nx = size(lo%lat, 1)
            ny = size(lo%lat, 2)

            xc = nx/2
            yc = ny/2
            xw = xc
            yw = yc
            ! use xsign and ysign in case lat/lon don't increase in a positive index direction
            xsign = sign(1.0, lo%lon(ims+1,jms) - lo%lon(ims,jms))
            ysign = sign(1.0, lo%lat(ims,jms+1) - lo%lat(ims,jms))

            ! use a O(log(n)) search instead of O(nxm)
            ! start at the halfway point and find the best direction to take in both directions
            iterations = 0
            do while ( (( xw > 2 ).or.( yw > 2 )).and.( iterations < 20) )
                iterations = iterations+1
                xc = min( max(xc, ims), ime)
                yc = min( max(yc, jms), jme)
                ! figure out which direction to step, then step half of the last step distance
                if (lo%lat(xc,yc) > lat) then
                    if (yw == 2) then
                        yw = 1
                    else
                        yw = yw/2 + 1
                    endif
                    yc = yc - ysign * yw
                endif
                yc = min( max(yc, jms), jme)
                if (lo%lat(xc,yc) < lat) then
                    if (yw == 2) then
                        yw = 1
                    else
                        yw = yw/2 + 1
                    endif
                    yc = yc + ysign * yw
                endif
                yc = min( max(yc, jms), jme)
                ! in case lat is exactly equal to lat(xc,yc)
                if (lo%lat(xc,yc) == lat) then
                    yw = 0
                endif
                if (lo%lon(xc,yc) > lon) then
                    if (xw == 2) then
                        xw = 1
                    else
                        xw = xw/2 + 1
                    endif
                    xc = xc - xsign * xw
                endif
                xc = min( max(xc, ims), ime)
                if (lo%lon(xc,yc) < lon) then
                    if (xw == 2) then
                        xw = 1
                    else
                        xw = xw/2 + 1
                    endif
                    xc = xc + xsign * xw
                endif
                xc = min( max(xc, ims), ime)
                ! in case lon is exactly equal to lon(xc,yc)
                if (lo%lon(xc,yc) == lon) then
                    xw = 0
                endif
                xw = minxw(xw, lo%lon, xc, yc, lon)
                yw = minyw(yw, lo%lat, xc, yc, lat)
            enddo
        endif

        ! once we have a "good" location we need to find the actual minimum (we could be above or below it by 1 or 2)
        if (iterations < 20) then
            mindist = 9999.9
            do xw = max(ims, xc-15), min(ime, xc+15)
                do yw = max(jms, yc-15), min(jme, yc+15)
                    curdist = sqrt( (lo%lat(xw,yw) - lat)**2 + (lo%lon(xw,yw) - lon)**2 )
                    if (curdist < mindist) then
                        mindist = curdist
                        x = xw
                        y = yw
                    endif
                enddo
            enddo
        else
            ! CASE 3 search every possible grid cell!
            ! naive search
!           write(*,*) "using n^2 search..."
!           note which point this is for debugging purposes, we really shouldn't get here (often)
!           write(*,*) lat,lon,lo%lat(xc,yc),lo%lon(xc,yc)
!           write(*,*) xc,yc,nx,ny
            mindist=9999.9
            do xw = ims, ime
                do yw = jms, jme
                    curdist = sqrt((lo%lat(xw,yw) - lat)**2 + (lo%lon(xw,yw) - lon)**2)
                    if (curdist < mindist) then
                        mindist = curdist
                        x = xw
                        y = yw
                    endif
                enddo
            enddo
        endif


        find_location%x = x
        find_location%y = y
    end function find_location

    !>------------------------------------------------------------
    !! Determine if a point is on a line defined by two points
    !!
    !! WARNING: assumes y is between y0 and y1 and x < max(x0,x1)
    !!
    !! x,y are the real point location
    !! x0,y0 defines one end of a line
    !! x1,y1 defines the other end of the line
    !!
    !! precision optionally defines the tolerance permitted to consider around the line
    !!
    !!------------------------------------------------------------
    function point_is_on_line(x,y,x0,y0,x1,y1, precision, error) result(online)
        ! WARNING: assumes y is between y0 and y1 and x < max(x0,x1)
        implicit none
        real, intent(in) :: x,y,x0,y0,x1,y1
        real, intent(in), optional :: precision
        real, intent(out),optional :: error
        logical :: online

        real :: internal_precision
        double precision :: slope, offset

        internal_precision = 1e-4
        if (present(precision)) internal_precision = precision

        online = .False.
        ! can we abort testing whether we are actually on a border?
        if (x >= min(x0,x1)) then
            ! we COULD be on a border
            if (x0 /= x1) then
                ! find the equation of the line and check if we are on it.
                ! this is performed in double precision space to be as precise as possible
                ! so that rounding errors are minimized
                slope = (DBLE(y1) - y0) / (DBLE(x1) - x0)
                offset = y0 - (slope * x0)
                ! note that due to rounding errors we could be "close" but not on it...
                if (abs( ((slope * x) + offset) - y) < internal_precision) then
                    ! we are on a border, return true
                    online = .True.
                    if (present(error)) error = abs( ((slope * x) + offset) - y)
                endif
            else
                ! if x0=x1 and x is between x0,x1 and y is between y0,y1 then
                if (x == x0) then ! this is assumed from the input requirements, also assumes y is between y0 and y1
                    online = .True.
                    if (present(error)) error = 0
                endif
            endif
        endif

    end function

    !>------------------------------------------------------------
    !! Determine if a point intersected a vertex
    !!
    !! y is the real point location
    !! y0 defines one end of a line
    !! y1 defines the other end of the line
    !!
    !! precision optionally defines the tolerance permitted to consider around the line
    !!
    !! sets top_vertex, bottom_vertex to True if the point intersected the top or bottom ends of the line (which might be the same)
    !! Flips the inside/outside boolean if it intersected one or both verticies
    !!
    !!------------------------------------------------------------
    subroutine check_vertex_hits(y, y0, y1, inside, top_vertex, bottom_vertex, precision, error)
        implicit none
        real,    intent(in)    :: y, y0, y1
        logical, intent(inout) :: top_vertex, bottom_vertex, inside
        real,    intent(in)    :: precision
        real,    intent(out), optional :: error

        if (present(error)) error = 0

        ! if the segment just parallels the line, then don't reset top_vertex, bottom_vertex history
        ! and don't flip inside/out
        if (y0==y1) return

        ! check if we hit the y0 vertex with our ray
        if (abs(y-y0) < precision) then
            if (present(error)) error = abs(y-y0)
            ! if so, check if we are above or below the other vertex segment
            if (y<=y1) then
                ! we are below the other vertex
                if (top_vertex) then
                    ! if we already found the other vertex (from testing the previous segment)
                    ! so reverse our inside outside to avoid double counting this vertex hit
                    top_vertex=.False.
                    inside = .not.inside
                else
                    ! otherwise, store the fact that we hit a bottom vertex for the next segment
                    bottom_vertex=.True.
                endif
            else if (y>=y1) then
                if (bottom_vertex) then
                    bottom_vertex=.False.
                    inside=.not.inside
                else
                    top_vertex=.True.
                endif
            endif
        endif

        if (abs(y-y1) < precision) then
            if (present(error)) error = abs(y-y1)

            if (y<=y0) then
                if (top_vertex) then
                    top_vertex=.False.
                    inside = .not.inside
                else
                    bottom_vertex=.True.
                endif
            else if (y>=y0) then
                if (bottom_vertex) then
                    bottom_vertex=.False.
                    inside=.not.inside
                else
                    top_vertex=.True.
                endif
            endif
        endif
    end subroutine check_vertex_hits

    !>------------------------------------------------------------
    !! Determine if a point is inside a given polygon or not
    !!
    !!   Polygon is a list of (x,y) pairs. This fuction returns True or False.
    !!   The algorithm is called the Ray Casting Method.
    !!
    !! loosely based off: http://www.ariel.com.au/a/python-point-int-poly.html
    !! Note additional discussion here: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
    !!
    !!------------------------------------------------------------
    function point_in_poly(x, y, poly, precision, error) result(inside)
        implicit none
        real, intent(in) :: x, y
        real, intent(in) :: poly(:,:)
        real, intent(in), optional :: precision
        real, intent(out), optional :: error
        logical :: inside
        logical :: top_vertex, bottom_vertex

        real :: internal_precision
        logical :: point_on_line
        integer :: n, i
        real :: x0,y0,x1,y1
        double precision :: slope, x_line
        real :: err, err_temp, line_err_temp

        internal_precision = 1e-4
        if (present(precision)) internal_precision = precision
        n = size(poly,2)
        err = 9999
        if (present(error)) error = 0

        ! by default we assume we are not in the polygon
        point_on_line = .False.
        inside = .False.
        top_vertex = .False.
        bottom_vertex = .False.

        ! check if point is a vertex
        do i=1,n
            ! check if the point is the next vertex
            if ((x==poly(1,i)).and.(y==poly(2,i))) then
                ! if the point is a vertex, we are inside the polygon, so return true
                inside=.True.
                return
            endif
        enddo

        x0 = poly(1,1)
        y0 = poly(2,1)
        ! now step through all edges checking to see if the number of edges that have to be crossed is odd or even
        do i=1,n
            x1 = poly(1,mod(i,n)+1)
            y1 = poly(2,mod(i,n)+1)
            ! the first if statements test to see if we can abort early
            if (y >= min(y0,y1)) then
                if (y <= max(y0,y1)) then
                    if (x <= max(x0,x1)) then

                        if (point_is_on_line(x,y,x0,y0,x1,y1, internal_precision, error=line_err_temp)) then
                            point_on_line = .True.
                            err = min(err, line_err_temp)
                            if (present(error)) error = err
                        endif

                        if (y0 /= y1) then
                            ! note slope is reversed from "normal" because we are looking across x
                            ! performed in double precision to match what is done in point_is_on_line (roughly)
                            slope = (dble(x1)-x0)/(dble(y1)-y0)
                            x_line = (y-dble(y0)) * slope + x0
                        else
                            x_line = -1e10 ! only inside if x0==x1
                        endif

                        ! if we cross one line then we are inside, but if we cross two we are outside, and so on
                        if ((x0 == x1) .or. (dble(x) <= x_line)) then
                            inside = .not. inside
                            call check_vertex_hits(y,y0,y1,inside, top_vertex, bottom_vertex, internal_precision, error=err_temp)
                            err = max(err, err_temp)
                        else
                            top_vertex = .False.
                            bottom_vertex = .False.
                        endif
                    endif
                endif
            endif

            x0 = x1
            y0 = y1
        enddo

        if (inside) err = 0
        if (point_on_line) inside = .True.
        if (present(error)) error = err

        ! inside will be set to the correct value
    end function point_in_poly


    function test_surrounding(lat,lon, lo, positions, nx, ny, n, error) result(surrounded)
        implicit none
        real,                    intent(in) :: lat, lon
        type(interpolable_type),intent(in) :: lo
        type(fourpos),           intent(inout) :: positions
        integer,                 intent(in) :: nx, ny, n
        real,                    intent(out), optional :: error
        logical :: surrounded
        integer :: i
        real :: err

        real :: polygon(2,n)

        err = 0

        ! first enforce that all points passed in are within the bounds of the lat/lon arrays
        do i=1,n
            positions%x(i) = min(max(positions%x(i),1), nx)
            positions%y(i) = min(max(positions%y(i),1), ny)
        enddo

        ! then setup a polygon containing the lat/lon data from the input
        do i=1,n
            polygon(1,i) = lo%lon( positions%x(i), positions%y(i))
            polygon(2,i) = lo%lat( positions%x(i), positions%y(i))
        enddo

        ! finally test that the provided lat/lon point falls within this polygon
        surrounded = point_in_poly(lon, lat, polygon, error=err)

        if (present(error)) error = err

    end function test_surrounding

    function test_triangle(lat,lon, lo, positions, nx, ny, error) result(surrounded)
        implicit none
        real,                    intent(in) :: lat, lon
        type(interpolable_type),intent(in) :: lo
        type(fourpos),           intent(inout) :: positions
        integer,                 intent(in) :: nx, ny
        real,                    intent(out), optional :: error
        logical :: surrounded
        integer :: i
        real :: err

        real :: polygon(2,3)

        err = 0

        ! assumes all points passed in are within the bounds of the lat/lon arrays!
        ! this is enforced in test_surrounding which should have passed by now
        ! then setup a polygon containing the lat/lon data from the input
        do i=1,2
            polygon(1,i) = lo%lon( positions%x(i), positions%y(i))
            polygon(2,i) = lo%lat( positions%x(i), positions%y(i))
        enddo

        ! set up the third point in the triangle as the average of all 4 surrounding points
        polygon(:,3) = 0
        do i=1,4
            polygon(1,3) = polygon(1,3) + lo%lon( positions%x(i), positions%y(i))
            polygon(2,3) = polygon(2,3) + lo%lat( positions%x(i), positions%y(i))
        enddo
        polygon(:,3) = polygon(:,3) / 4

        ! finally test that the provided lat/lon point falls within this polygon
        surrounded = point_in_poly(lon, lat, polygon, error=err)

        if (present(error)) error = err

    end function test_triangle


    !>------------------------------------------------------------
    !!  Given a closest position, return the 4 points surrounding the lat/lon position in lo%lat/lon
    !!
    !!------------------------------------------------------------
    type(fourpos) function find_surrounding(lo,lat,lon,pos,nx,ny)
        implicit none
        type(interpolable_type),intent(in)::lo
        real,intent(in)::lat,lon
        type(position),intent(in)::pos
        integer,intent(in) :: nx,ny
        integer :: i, j
        integer :: xdeltas(4), ydeltas(4), dx, dy
        integer :: lowerx, upperx, lowery, uppery
        real :: best_err, current_err, tri1_err, tri2_err
        logical :: tri_1, tri_2
        type(fourpos) :: search_point

        best_err = 1e10

        xdeltas = [-1,-1,1,1]
        ydeltas = [-1,1,-1,1]

        find_surrounding%x = -999
        find_surrounding%y = -999

        lowerx = lbound(lo%lat,1)
        upperx = ubound(lo%lat,1)
        lowery = lbound(lo%lat,2)
        uppery = ubound(lo%lat,2)

        do i = 1, 4
            dx = xdeltas(i)
            dy = ydeltas(i)

            ! check that the current search box exists in the input data
            if (((pos%x + dx) > upperx) .or. ((pos%x + dx) < lowerx) .or. ((pos%y + dy) > uppery) .or. ((pos%y + dy) < lowery)) then
                cycle
            endif

            search_point%x = [pos%x, pos%x+dx, pos%x+dx, pos%x  ]
            search_point%y = [pos%y, pos%y,   pos%y+dy, pos%y+dy]

            if (test_surrounding(lat, lon, lo, search_point, nx, ny, 4, error=current_err)) then
                if (current_err < best_err) then
                    best_err = current_err
                    ! if it is not in this triangle, it must be in the other.

                    tri_1 = test_triangle(lat, lon, lo, search_point, nx, ny, error=tri1_err)

                    search_point%x = [pos%x, pos%x,   pos%x+dx, pos%x+dx ]
                    search_point%y = [pos%y, pos%y+dy, pos%y,   pos%y+dy]

                    tri_2 = test_triangle(lat, lon, lo, search_point, nx, ny, error=tri2_err)

                    if (tri_1) then
                        search_point%x = [pos%x, pos%x+dx, pos%x+dx, pos%x  ]
                        search_point%y = [pos%y, pos%y,   pos%y+dy, pos%y+dy]

                        if (tri_2) then
                            if (tri2_err < tri1_err) then
                                search_point%x = [pos%x, pos%x,   pos%x+dx, pos%x+dx ]
                                search_point%y = [pos%y, pos%y+dy, pos%y,   pos%y+dy]
                            endif
                        endif
                    elseif (.not.tri_2) then
                        ! in edge cases the x,y point is technically closest, but the lat, lon point won't fall in a triangle that contains it
                        ! so search the other possible triangles in this 4 grid cell box too.
                        search_point%x = [pos%x,    pos%x+dx, pos%x+dx, pos%x]
                        search_point%y = [pos%y+dy, pos%y+dy, pos%y,    pos%y]

                        tri_1 = test_triangle(lat, lon, lo, search_point, nx, ny, error=tri1_err)

                        if (.not.tri_1) then
                            search_point%x = [pos%x+dx, pos%x+dx, pos%x,    pos%x]
                            search_point%y = [pos%y,    pos%y+dy, pos%y+dy, pos%y]

                            tri_2 = test_triangle(lat, lon, lo, search_point, nx, ny, error=tri2_err)
                        endif

                    endif
                    find_surrounding = search_point


                    if (.not.test_triangle(lat, lon, lo, search_point, nx, ny)) then
                        write(*,*) "Warning point in box, but not triangle"
                        write(*,*) search_point%x
                        write(*,*) search_point%y
                        write(*,*) lat, lon
                        do j=1,4
                            write(*,*) lo%lat(search_point%x(j), search_point%y(j)), &
                                       lo%lon(search_point%x(j), search_point%y(j))
                        enddo
                    endif

                endif

            endif
        enddo

        if (find_surrounding%x(1) > -999) return

        write(*,*) "ERROR: Failed to find point", lon, lat, " in a quadrant near:"
        write(*,*) lo%lon(pos%x, pos%y), lo%lat(pos%x, pos%y)
        write(*,*) pos
        write(*,*) nx, ny
        error stop "Failed to find point"
        ! enforce that surround points fall within the bounds of the full domain

    end function find_surrounding

    !>------------------------------------------------------------
    !!  Compute the geographic look up table from LOw resolution grid to HIgh resolution grid
    !!
    !!------------------------------------------------------------
    subroutine geo_LUT(hi, lo)
        implicit none
        type(interpolable_type), intent(in)    :: hi
        type(interpolable_type), intent(inout) :: lo
        type(fourpos) :: xy
        type(position) :: curpos, lastpos
        integer :: nx, ny, i, j, k, lo_nx, lo_ny
        integer :: ims,ime,jms,jme
        real, dimension(4) :: lat, lon

        nx    = size(hi%lat,1)
        ny    = size(hi%lat,2)
        lo_nx = size(lo%lat,1)
        lo_ny = size(lo%lat,2)

        ims = lbound(hi%lat,1)
        ime = ubound(hi%lat,1)
        jms = lbound(hi%lat,2)
        jme = ubound(hi%lat,2)

        allocate(lo%geolut%x(4,ims:ime, jms:jme))
        allocate(lo%geolut%y(4,ims:ime, jms:jme))
        allocate(lo%geolut%w(4,ims:ime, jms:jme))

        curpos%x = 1
        curpos%y = 1

        do j=jms,jme
            curpos%x = 1

            ! use a brute force approach to find a starting
            lastpos  = find_location(lo, hi%lat(ims,j), hi%lon(ims,j), curpos)
            if (lastpos%x < 1) then
                ! something broke, try assuming that the grids are the same
                ! should put in some check here to make sure this doesn't happen in a real case
                ! lastpos%x = 1
                ! lastpos%y = 1
                write(*,*) "Error in Geographic interpolation, check input lat / lon grids", ims,j
                error stop
            endif

            do i = ims, ime
                curpos = find_location(lo, hi%lat(i,j), hi%lon(i,j), lastpos)

                if (curpos%x < 1) then
                    ! something broke, try assuming that the grids are the same possibly wrapping (for ideal)
                    write(*,*) "Error in Geographic interpolation, check input lat / lon grids", i,j
                    curpos%x = mod(i-1, lo_nx) + 1
                    curpos%y = mod(j-1, lo_ny) + 1
                    ! all "surrounding" grid points are "this" grid point
                    lo%geolut%x(:,i,j) = curpos%x
                    lo%geolut%y(:,i,j) = curpos%y
                    lo%geolut%w(:,i,j) = 0.25
                else
                    ! Found a good point, now find the other 3 of the surrounding 4 points
                    xy = find_surrounding(lo, hi%lat(i,j), hi%lon(i,j), curpos, lo_nx, lo_ny)
                    lo%geolut%x(:,i,j) = xy%x
                    lo%geolut%y(:,i,j) = xy%y
                    ! load those latitutes and longitudes into 1D arrays to calculate weights
                    do k=1,4
                        lat(k) = lo%lat(xy%x(k), xy%y(k))
                        lon(k) = lo%lon(xy%x(k), xy%y(k))
                    enddo
                    ! and calculate the weights to apply to each gridcell
                    lo%geolut%w(:,i,j) = tri_weights(hi%lat(i,j), lat, hi%lon(i,j), lon)
                    ! lo%geolut%w(:,i,j) = idw_weights(hi%lat(i,j), lat, hi%lon(i,j), lon)
                    ! lo%geolut%w(:,i,j) = bilin_weights(hi%lat(i,j), lat, hi%lon(i,j), lon)

                endif

            enddo
        enddo

    end subroutine geo_LUT

    !>------------------------------------------------------------
    !!  Interpolate boundaries of fieldout to fieldin using geolut
    !!
    !!------------------------------------------------------------
    subroutine boundary_interpolate(fieldout, fieldin, geolut)
        implicit none
        real,intent(inout)::fieldout(:,:,:)
        real,intent(in)::fieldin(:,:,:)
        type(geo_look_up_table),intent(in)::geolut
        integer::nx,nz,ny
        integer :: ims,ime,jms,jme,kms,kme
        integer:: i,j,k,l,localx,localy
        real::localw
        real :: local_center

        nx=size(fieldout,1)
        nz=size(fieldout,2)
        ny=size(fieldout,3)

        ims = lbound(fieldout,1)
        ime = ubound(fieldout,1)
        jms = lbound(fieldout,2)
        jme = ubound(fieldout,2)
        kms = lbound(fieldout,3)
        kme = ubound(fieldout,3)

        ! use the geographic lookup table generated earlier to
        ! compute a bilinear interpolation from lo to hi
        ! first loop over all x elements on the y ends
        do k=kms,kme,ny-1
            do j=jms,jme
                do i=ims,ime
                    fieldout(i,j,k)=0
                    local_center = 0
                    do l=1,4
                        localx=geolut%x(l,i,k)
                        localy=geolut%y(l,i,k)
                        local_center = local_center + fieldin(localx,j,localy)
                        if (l < 3) then
                            localw=geolut%w(l,i,k)
                            fieldout(i,j,k) = fieldout(i,j,k) + fieldin(localx,j,localy) * localw
                        endif
                    enddo
                    localw = geolut%w(3,i,k)
                    fieldout(i,j,k) = fieldout(i,j,k) + local_center/4 * localw
                    ! fieldout(i,j,k)=0
                    ! do l=1,4
                    !     localx=geolut%x(l,i,k)
                    !     localy=geolut%y(l,i,k)
                    !     localw=geolut%w(l,i,k)
                    !     fieldout(i,j,k)=fieldout(i,j,k)+fieldin(localx,j,localy)*localw
                    ! enddo
                enddo
            enddo
        enddo
        ! then loop over all y elements on the x ends
        do k=kms,kme
            do j=jms,jme
                do i=ims,ime,nx-1
                    fieldout(i,j,k)=0
                    local_center = 0
                    do l=1,4
                        localx=geolut%x(l,i,k)
                        localy=geolut%y(l,i,k)
                        local_center = local_center + fieldin(localx,j,localy)
                        if (l < 3) then
                            localw=geolut%w(l,i,k)
                            fieldout(i,j,k) = fieldout(i,j,k) + fieldin(localx,j,localy) * localw
                        endif
                    enddo
                    localw = geolut%w(3,i,k)
                    fieldout(i,j,k) = fieldout(i,j,k) + local_center/4 * localw
                    ! fieldout(i,j,k)=0
                    ! do l=1,4
                    !     localx=geolut%x(l,i,k)
                    !     localy=geolut%y(l,i,k)
                    !     localw=geolut%w(l,i,k)
                    !     fieldout(i,j,k)=fieldout(i,j,k)+fieldin(localx,j,localy)*localw
                    ! enddo
                enddo
            enddo
        enddo

    end subroutine boundary_interpolate

    !>------------------------------------------------------------
    !!  Interpolate fieldout to fieldin using geolut.
    !!  if boundary_only is true, call boundary_interpolate instead
    !!  loops over y,z,x but geolut is only defined over x,y (for now)
    !!
    !!------------------------------------------------------------
    subroutine geo_interp(fieldout, fieldin, geolut, boundary_only)
        implicit none
        real,                   intent(inout) :: fieldout(:,:,:)
        real,                   intent(in)    :: fieldin(:,:,:)
        type(geo_look_up_table),intent(in)    :: geolut
        logical,                intent(in),   optional :: boundary_only

        integer :: ims,ime, jms,jme, kms,kme
        integer:: i,j,k,l, localx,localy, nz_input
        logical :: boundaries
        real :: localw
        real :: local_center

        ! This can probably all be done on the 1-n assumption, but this keeps the code more flexible in the future
        ims = lbound(fieldout,1)
        ime = ubound(fieldout,1)
        jms = lbound(fieldout,2)
        jme = ubound(fieldout,2)
        kms = lbound(fieldout,3)
        kme = ubound(fieldout,3)

        nz_input = size(fieldin,2)

        boundaries = .False.
        if (present(boundary_only)) boundaries = boundary_only

        ! if we are only processing the boundary, then make x and y increments be the size of the array
        ! so we only hit the edges of the array
        if (boundaries) then
            call boundary_interpolate(fieldout, fieldin, geolut)
        else
        ! use the geographic lookup table generated earlier to
        ! compute a bilinear interpolation from lo to hi
        ! if we are doing the interior too, just iterate over all x and y
            do k=kms,kme
                do j=jms,jme
                    do i=ims,ime
                        fieldout(i,j,k)=0
                        local_center = 0
                        ! This loop is used if doing a triangular bi-linear interpolation on the grid
                        ! Need to see all 4 points to compute the "local_center" value for the triangulation
                        do l=1,4
                            localx=geolut%x(l,i,k)
                            localy=geolut%y(l,i,k)
                            local_center = local_center + fieldin(localx, min(j,nz_input), localy)
                            ! This is the actual interpolation adding the value from the two raw grid points
                            if (l < 3) then

                                localw=geolut%w(l,i,k)
                                fieldout(i,j,k) = fieldout(i,j,k) + fieldin(localx, min(j,nz_input), localy) * localw
                            endif
                        enddo
                        ! This is the actual interpolation adding the value from the center average of all four grid points
                        localw = geolut%w(3,i,k)
                        fieldout(i,j,k) = fieldout(i,j,k) + local_center/4 * localw

                        ! This loop is used if doing a straight bi-linear interpolation on the grid
                        ! do l=1,4
                        !     localx=geolut%x(l,i,k)
                        !     localy=geolut%y(l,i,k)
                        !     localw=geolut%w(l,i,k)
                        !     fieldout(i,j,k)=fieldout(i,j,k)+fieldin(localx,j,localy)*localw
                        ! enddo
                    enddo
                enddo
            enddo
        endif
    end subroutine geo_interp

    !>------------------------------------------------------------
    !!  Interpolate fieldout to fieldin using geolut.
    !!
    !!------------------------------------------------------------
    subroutine geo_interp2d(fieldout, fieldin, geolut)
        implicit none
        real, dimension(:,:),   intent(inout) :: fieldout
        real, dimension(:,:),   intent(in)    :: fieldin
        type(geo_look_up_table),intent(in)    :: geolut

        integer :: i,k,l,ny,nx,localx,localy
        integer :: ims,ime,jms,jme
        real :: localw
        real :: local_center

        nx = size(fieldout,1)
        ny = size(fieldout,2)

        ims = lbound(fieldout,1)
        ime = ubound(fieldout,1)
        jms = lbound(fieldout,2)
        jme = ubound(fieldout,2)

        ! use the geographic lookup table generated earlier to
        ! compute a bilinear interpolation from lo to hi
        do k = jms, jme
            do i = ims, ime
                fieldout(i,k)=0
                local_center = 0
                ! This loop is used if doing a triangular bi-linear interpolation on the grid
                ! Need to see all 4 points to compute the "local_center" value for the triangulation
                do l = 1, 4
                    localx = geolut%x(l,i,k)
                    localy = geolut%y(l,i,k)
                    local_center = local_center + fieldin(localx,localy)
                    ! This is the actual interpolation adding the value from the two raw grid points
                    if (l < 3) then
                        localw = geolut%w(l,i,k)
                        fieldout(i,k) = fieldout(i,k) + fieldin(localx,localy) * localw
                    endif
                enddo
                ! This is the actual interpolation adding the value from the center average of all four grid points
                localw = geolut%w(3,i,k)
                fieldout(i,k) = fieldout(i,k) + local_center/4 * localw

                ! This loop is used if doing a straight bi-linear interpolation on the grid
                ! do l=1,4
                !     localx=geolut%x(l,i,k)
                !     localy=geolut%y(l,i,k)
                !     localw=geolut%w(l,i,k)
                !     fieldout(i,k)=fieldout(i,k)+fieldin(localx,localy)*localw
                ! enddo
            enddo
        enddo

    end subroutine


    !> ----------------------------------------------------------------------------
    !!  Conver longitudes into a common format for geographic interpolation
    !!
    !!  First convert longitudes that may be -180-180 into 0-360 range first.
    !!  If data straddle the 0 degree longitude, convert back to -180-180 so interpolation will be smooth
    !!
    !!  @param[inout]   long_data  2D array of longitudes (either -180 - 180 or 0 - 360)
    !!
    !! ----------------------------------------------------------------------------
    subroutine standardize_coordinates(domain, longitude_system)
        implicit none
        type(interpolable_type), intent(inout) :: domain
        integer,                  intent(in), optional :: longitude_system

        real, dimension(:,:), allocatable :: temporary_geo_data
        integer :: nx, ny, i
        integer :: ims,ime,jms,jme

        ! if the lat, lon data were given as 1D variables, we need to make them 2D for interpolability
        if (size(domain%lat,2)==1) then
            nx = size(domain%lon,1)
            ny = size(domain%lat,1)

            ims = lbound(domain%lon,1)
            ime = ubound(domain%lon,1)
            jms = lbound(domain%lat,1)
            jme = ubound(domain%lat,1)

            allocate(temporary_geo_data(ims:ime, jms:jme))
            do i = jms,jme
                temporary_geo_data(:,i) = domain%lon(:,1)
            end do

            deallocate(domain%lon)
            allocate(domain%lon(ims:ime, jms:jme))
            domain%lon = temporary_geo_data

            do i = ims, ime
                temporary_geo_data(i,:) = domain%lat(:,1)
            end do
            deallocate(domain%lat)
            allocate(domain%lat(ims:ime, jms:jme))
            domain%lat = temporary_geo_data

        endif

        if (present(longitude_system)) then
            if (longitude_system == kMAINTAIN_LON) then
                ! break
            elseif (longitude_system == kDATELINE_CENTERED) then
                where(domain%lon<0) domain%lon = 360+domain%lon
            elseif (longitude_system == kPRIME_CENTERED) then
                where(domain%lon>180) domain%lon = domain%lon - 360
            elseif (longitude_system == kGUESS_LON) then
                ! try to guess which coordinate system to use
                ! first convert all domains into -180 to +180 coordinates
                where(domain%lon>180) domain%lon = domain%lon - 360
                ! also convert from a -180 to 180 coordinate system into a 0-360 coordinate system if closer to the dateline
                if ((minval(domain%lon) < -150) .or. (maxval(domain%lon) > 150)) then
                    where(domain%lon<0) domain%lon = 360+domain%lon
                endif
            else
                write(*,*) "Unknown longitude system of coordinates:", longitude_system
                write(*,*) " Set options%parameters%longitude_system to either:"
                write(*,*) " Prime meridian centered (0 to 360): ", kPRIME_CENTERED
                write(*,*) " Dateline meridian centered (-180 to 180): ", kDATELINE_CENTERED
                error stop
            endif
        endif

    end subroutine standardize_coordinates

end module geo
