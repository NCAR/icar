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
    implicit none

    private

    public::geo_LUT      ! Create a geographic Look up table
    public::geo_interp   ! apply geoLUT to interpolate in 2d for a 3d grid
    public::geo_interp2d ! apply geoLUT to interpolate in 2d for a 2d grid
    public::standardize_coordinates

contains
    !>------------------------------------------------------------
    !! Calculate the weights to use for bilinear interpolation between surrounding points x, y
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
    !! Find the minimum x
    !! Takes a point (yi,xi) and a set of 4 surrounding points (y[4],x[4]) as input
    !!
    !! Example : xw = minxw(xw,lo%lon,xc,yc,lon)
    !!
    !! @param   real    yi  input y location to interpolate to.
    !! @param   real    y   Array of surrounding y-locations to interpolate from.
    !! @param   real    xi  input x location to interpolate to.
    !! @param   real    x   Array of surrounding x-locations to interpolate from.
    !! @retval  real    weights
    !!
    !!------------------------------------------------------------
    integer function minxw(xw,longrid,xpos,ypos,lon)
        implicit none
        integer, intent(in)::xw,xpos,ypos
        real,intent(in)::longrid(:,:),lon
        real::curdist,dx
        integer::nx

        nx=size(longrid,1)
        curdist=abs(lon-longrid(xpos,ypos))
        if (xpos<nx) then
            dx=abs(longrid(xpos,ypos)-longrid(xpos+1,ypos))
        else
            dx=abs(longrid(xpos,ypos)-longrid(xpos-1,ypos))
        endif
        minxw=max(xw,ceiling(curdist/dx))
        if (minxw>xw) then
            if (minxw>4) then
                minxw=ceiling(minxw*0.8)
            endif
        endif
    end function minxw
    !>------------------------------------------------------------
    !!  Find the minimum next step size in the y direction
    !!  yw = minyw(yw,lo%lat,xc,yc,lat)
    !!
    !!------------------------------------------------------------
    integer function minyw(yw,latgrid,xpos,ypos,lat)
        integer, intent(in)::yw,xpos,ypos
        real,intent(in)::latgrid(:,:),lat
        real::curdist,dx
        integer::ny

        ny=size(latgrid,2)
        curdist=abs(lat-latgrid(xpos,ypos))
        if (ypos<ny) then
            dx=abs(latgrid(xpos,ypos)-latgrid(xpos,ypos+1))
        else
            dx=abs(latgrid(xpos,ypos)-latgrid(xpos,ypos-1))
        endif
        minyw=max(yw,ceiling(curdist/dx))
        if (minyw>yw) then
            if (minyw>4) then
                minyw=ceiling(minyw*0.8)
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
        class(interpolable_type),intent(inout)::lo      ! input interpolable object to search
        real,intent(in)::lat,lon                        ! input position to find
        type(position),intent(in)::lastpos              ! position to start search from
        ! locals
        real::mindist, curdist,dx,dy,xsign,ysign,x,y    ! temporary variables to use during the search
        integer::nx,ny,xc,yc,xw,yw,iterations,xstep,ystep   ! more temporary variables

        nx = size(lo%lat, 1)
        ny = size(lo%lat, 2)
        ! calcualte dx/dy at the middle of the grid
        dx = lo%lon(2,1) - lo%lon(1,1)
        dy = lo%lat(1,2) - lo%lat(1,1)

        ! Handle weird edge case for WRF ideal simulations
        if (dx == 0) then
            if (.not.lo%dx_errors_printed) then
                print*, "ERROR : geo_find_location : DX = 0 !  Check your inputfile Longitude data"
                print*, "  Attempting to continue by assuming the grids are the same or can wrap arround"
                lo%dx_errors_printed=.True.
            endif
            find_location%x = -1
            find_location%y = -1
            return
        endif
        if (dy == 0) then
            if (.not.lo%dy_errors_printed) then
                print*, "ERROR : geo_find_location : DY = 0 !  Check your inputfile Latitude data"
                print*, "  Attempting to continue by assuming the grids are the same or can wrap arround"
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

        x = lo%lon(xc,yc)
        y = lo%lat(xc,yc)

        ! steps to take = the difference the between the current point and the input point / dx
        xstep = (lon-x)/dx
        ystep = (lat-y)/dy

        ! CASE 1 assume a quasi regular dx/dy and caluate new location
        !  while we need to step by more than one grid cell, iterate
        !  if the grid is highly regular, we will only iterate 1-2x, highly irregular might require more iterations
        !  in the diabolical case, this could fail?
        !  in most cases this succeeds with a few iterations.
        iterations = 0
        do while (((abs(xstep)>1).or.(abs(ystep)>1)).and.iterations<20)
            iterations = iterations+1
            ! update the current x/y locations
            ! force it to be <nx-1 so we can calculate dx from (xc+1)-xc
            xc = max(1,min(nx-1,xc+xstep))
            yc = max(1,min(ny-1,yc+ystep))
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
        xc = max(1, min(nx, xc+xstep ) )
        yc = max(1, min(ny, yc+ystep ) )

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
            ysign = sign(1.0, lo%lat(1,2) - lo%lat(1,1) )
            xsign = sign(1.0, lo%lon(2,1) - lo%lon(1,1) )

            ! use a O(log(n)) search instead of O(nxm)
            ! start at the halfway point and find the best direction to take in both directions
            iterations = 0
            do while ( (( xw > 2 ).or.( yw > 2 )).and.( iterations < 20) )
                iterations = iterations+1
                xc = min( max(xc, 1), nx)
                yc = min( max(yc, 1), ny)
                ! figure out which direction to step, then step half of the last step distance
                if (lo%lat(xc,yc) > lat) then
                    if (yw == 2) then
                        yw = 1
                    else
                        yw = yw/2 + 1
                    endif
                    yc = yc - ysign * yw
                endif
                yc = min( max(yc, 1), ny)
                if (lo%lat(xc,yc) < lat) then
                    if (yw == 2) then
                        yw = 1
                    else
                        yw = yw/2 + 1
                    endif
                    yc = yc + ysign * yw
                endif
                yc = min( max(yc, 1), ny)
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
                xc = min( max(xc, 1), nx)
                if (lo%lon(xc,yc) < lon) then
                    if (xw == 2) then
                        xw = 1
                    else
                        xw = xw/2 + 1
                    endif
                    xc = xc + xsign * xw
                endif
                xc = min( max(xc, 1), nx)
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
            do xw = max(1, xc-15), min(nx, xc+15)
                do yw = max(1, yc-15), min(ny, yc+15)
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
            do xw=1,nx
                do yw=1,ny
                    curdist=sqrt((lo%lat(xw,yw)-lat)**2 + (lo%lon(xw,yw)-lon)**2)
                    if (curdist<mindist) then
                        mindist=curdist
                        x=xw
                        y=yw
                    endif
                enddo
            enddo
        endif


        find_location%x = x
        find_location%y = y
    end function find_location


    !>------------------------------------------------------------
    !!  Given a closest position, return the 4 points surrounding the lat/lon position in lo%lat/lon
    !!
    !!------------------------------------------------------------
    type(fourpos) function find_surrounding(lo,lat,lon,pos,nx,ny)
        implicit none
        class(interpolable_type),intent(in)::lo
        real,intent(in)::lat,lon
        type(position),intent(in)::pos
        integer,intent(in) :: nx,ny
        integer :: i

        if ((lo%lat(pos%x,pos%y)-lat) > 0) then
            find_surrounding%y = (/pos%y,pos%y,pos%y-1,pos%y-1/)
        else
            find_surrounding%y = (/pos%y,pos%y,pos%y+1,pos%y+1/)
        endif

        if ((lo%lon(pos%x,pos%y)-lon) > 0) then
            find_surrounding%x = (/pos%x,pos%x-1,pos%x,pos%x-1/)
        else
            find_surrounding%x = (/pos%x,pos%x+1,pos%x,pos%x+1/)
        endif

        ! enforce that surround points fall within the bounds of the full domain
        do i=1,4
            find_surrounding%x(i) = min(max(find_surrounding%x(i),1), nx)
            find_surrounding%y(i) = min(max(find_surrounding%y(i),1), ny)
        enddo

    end function find_surrounding

    !>------------------------------------------------------------
    !!  Compute the geographic look up table from LOw resolution grid to HIgh resolution grid
    !!
    !!------------------------------------------------------------
    subroutine geo_LUT(hi, lo)
        implicit none
        class(interpolable_type), intent(in)    :: hi
        class(interpolable_type), intent(inout) :: lo
        type(fourpos) :: xy
        type(position) :: curpos, lastpos
        integer :: nx, ny, i, j, k, lo_nx, lo_ny
        real, dimension(4) :: lat, lon

        nx    = size(hi%lat,1)
        ny    = size(hi%lat,2)
        lo_nx = size(lo%lat,1)
        lo_ny = size(lo%lat,2)

        allocate(lo%geolut%x(4,nx,ny))
        allocate(lo%geolut%y(4,nx,ny))
        allocate(lo%geolut%w(4,nx,ny))

        curpos%x = 1
        curpos%y = 1

        do j=1,ny
            curpos%x = 1
            ! use a brute force approach to find a starting
            lastpos  = find_location(lo, hi%lat(1,j), hi%lon(1,j), curpos)
            if (lastpos%x<1) then
                ! something broke, try assuming that the grids are the same
                ! should put in some check here to make sure this doesn't happen in a real case
                lastpos%x = 1
                lastpos%y = j
                write(*,*) "Error in Geographic interpolation, check input lat / lon grids", 1,j
            endif

            do i=1, nx
                curpos = find_location(lo, hi%lat(i,j), hi%lon(i,j), lastpos)
                if (curpos%x<1) then
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
                    lo%geolut%w(:,i,j) = bilin_weights(hi%lat(i,j), lat, hi%lon(i,j), lon)
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
        integer:: i,j,k,l,localx,localy
        real::localw

        nx=size(fieldout,1)
        nz=size(fieldout,2)
        ny=size(fieldout,3)
        ! use the geographic lookup table generated earlier to
        ! compute a bilinear interpolation from lo to hi
        ! first loop over all x elements on the y ends
        do k=1,ny,ny-1
            do j=1,nz
                do i=1,nx
                    fieldout(i,j,k)=0
                    do l=1,4
                        localx=geolut%x(l,i,k)
                        localy=geolut%y(l,i,k)
                        localw=geolut%w(l,i,k)
                        fieldout(i,j,k)=fieldout(i,j,k)+fieldin(localx,j,localy)*localw
                    enddo
                enddo
            enddo
        enddo
        ! then loop over all y elements on the x ends
        do k=1,ny
            do j=1,nz
                do i=1,nx,nx-1
                    fieldout(i,j,k)=0
                    do l=1,4
                        localx=geolut%x(l,i,k)
                        localy=geolut%y(l,i,k)
                        localw=geolut%w(l,i,k)
                        fieldout(i,j,k)=fieldout(i,j,k)+fieldin(localx,j,localy)*localw
                    enddo
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
    subroutine geo_interp(fieldout,fieldin,geolut,boundary_only)
        implicit none
        real,intent(inout)::fieldout(:,:,:)
        real,intent(in)::fieldin(:,:,:)
        type(geo_look_up_table),intent(in)::geolut
        logical,intent(in)::boundary_only
        integer::nx,nz,ny
        integer:: i,j,k,l,localx,localy
        real::localw

        nx=size(fieldout,1)
        nz=size(fieldout,2)
        ny=size(fieldout,3)

        ! if we are only processing the boundary, then make x and y increments be the size of the array
        ! so we only hit the edges of the array
        if (boundary_only) then
            call boundary_interpolate(fieldout, fieldin, geolut)
        else
        ! use the geographic lookup table generated earlier to
        ! compute a bilinear interpolation from lo to hi
        ! if we are doing the interior too, just iterate over all x and y
            do k=1,ny
                do j=1,nz
                    do i=1,nx
                        fieldout(i,j,k)=0
                        do l=1,4
                            localx=geolut%x(l,i,k)
                            localy=geolut%y(l,i,k)
                            localw=geolut%w(l,i,k)
                            fieldout(i,j,k)=fieldout(i,j,k)+fieldin(localx,j,localy)*localw
                        enddo
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
        real, dimension(:,:),intent(inout) :: fieldout
        real, dimension(:,:),intent(in) :: fieldin
        type(geo_look_up_table),intent(in)::geolut
        integer::i,k,l,ny,nx,localx,localy
        real::localw

        nx=size(fieldout,1)
        ny=size(fieldout,2)
        ! use the geographic lookup table generated earlier to
        ! compute a bilinear interpolation from lo to hi
        do k=1,ny
            do i=1,nx
                fieldout(i,k)=0
                do l=1,4
                    localx=geolut%x(l,i,k)
                    localy=geolut%y(l,i,k)
                    localw=geolut%w(l,i,k)
                    fieldout(i,k)=fieldout(i,k)+fieldin(localx,localy)*localw
                enddo
            enddo
        enddo

    end subroutine  geo_interp2d


    subroutine standardize_coordinates(domain)
        implicit none
        class(interpolable_type), intent(inout) :: domain

        real, dimension(:,:), allocatable :: temporary_geo_data
        integer :: nx, ny, i

        ! if the lat, lon data were given as 1D variables, we need to make them 2D for interpolability
        if (size(domain%lat,2)==1) then
            nx = size(domain%lon,1)
            ny = size(domain%lat,1)

            allocate(temporary_geo_data(nx,ny))
            do i = 1, ny
                temporary_geo_data(:,i) = domain%lon(:,1)
            end do
            deallocate(domain%lon)
            allocate(domain%lon(nx,ny))
            domain%lon = temporary_geo_data

            do i = 1, nx
                temporary_geo_data(i,:) = domain%lat(:,1)
            end do
            deallocate(domain%lat)
            allocate(domain%lat(nx,ny))
            domain%lat = temporary_geo_data

        endif

        ! also convert from a -180 to 180 coordinate system into a 0-360 coordinate system if necessary
        where(domain%lon<0) domain%lon = 360+domain%lon

    end subroutine standardize_coordinates

end module geo
