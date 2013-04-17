module geo
	use data_structures
	implicit none
! 	type geo_info
! 		real,allocatable,dimension(:,:)::lat,lon
! 	end type geo_info
! initial template for a reader_type object for now abandoned...
! 	type reader_type
! 		real, allocatable, dimension(:,:,:) :: datavar
! 		real, allocatable, dimension(:,:,:) :: qv,cloud,ice,nice,qrain,nrain,qsnow,qgrau
! 		real, allocatable, dimension(:,:) :: terrain,rain,snow,graupel,dzdx,dzdy
! 		type(geolut)::geolut
! 		real::dx,dt
! 	contains
! 		procedure :: init
! 		procedure :: next_BC
! 	end type reader_type
	
contains
	
	function bilin_weights(yi,y,xi,x)
		implicit none
		real,intent(in)::yi,y(4),xi,x(4)
		real::x0,x1,x2,x3,y5,y6,f1,f2
		real, dimension(4) ::bilin_weights
		
	    x0=abs((xi-x(0))/(x(1)-x(0)))
	    x1=1-x0
	    x2=abs((xi-x(2))/(x(3)-x(2)))
	    x3=1-x2
	    y5=y(0)*x1+y(1)*x0
	    y6=y(2)*x3+y(3)*x2
	    f1=(yi-y5)/(y6-y5)
	    f2=1-f1
		bilin_weights=(/x1*f2,x0*f2,x3*f1,x2*f1/)
	end function bilin_weights
	
	type(position) function find_location(lo,lat,lon)
! Find a location lat,lon in lo%lat,lon grids
!  Assumes that the lat/lon grids are semi-regular (dx/dy aren't constant but they are nearly so)
!  Calculates dx/dy at the middle of the lat/lon grids, then calculates the location of lat,lon
!  input point using those coordinates.  
!  Next proceeds to calculate a new position based on the dx/dy at that position until the new
!  position is within 1 gridcell of the current position
!  Once that "within 1" position is found, search all cells in a 3x3 grid for the minimum distance
!  return a position datatype that includes the found x/y location
		implicit none
		type(bc_type),intent(in)::lo
		real,intent(in)::lat,lon
		real::mindist, curdist,dx,dy,xsign,ysign
		integer::x,y,nx,ny,xc,yc,xw,yw,iterations,xstep,ystep
		
		ny=size(lo%lat,1)
		nx=size(lo%lat,2)
! 		calcualte dx/dy at the middle of the grid
		dx=lo%lon(1,2)-lo%lon(1,1)
		dy=lo%lat(2,1)-lo%lat(1,1)
! 		current/starting position = the middle of the grid
		yc=ny/2
		xc=nx/2
		x=lo%lon(yc,xc)
		y=lo%lat(yc,xc)
! 		steps to take = the difference the between the current point and the input point / dx
		xstep=(lon-x)/dx
		ystep=(lat-y)/dy
! 		while we need to step by more than one grid cell, iterate
! 		if the grid is highly regular, we will only iterate 1-2x, highly irregular might require more iterations
! 		in the diabolical case, this could fail? 
		iterations=0
		do while (((xstep>1).or.(ystep>1)).and.iterations<20)
			iterations=iterations+1
! 			update the current x/y locations
! 			force it to be <nx-1 so we can calculate dx from (xc+1)-xc
			xc=min(1,max(nx-1,xc+xstep))
			yc=min(1,max(ny-1,yc+ystep))
			x=lo%lon(yc,xc)
			y=lo%lat(yc,xc)
! 			calculate a new dx/dy and x/y step
			dx=lo%lon(yc,xc+1)-lo%lon(yc,xc)
			dy=lo%lat(yc+1,xc)-lo%lat(yc,xc)
			xstep=NINT((lon-x)/dx)
			ystep=NINT((lat-y)/dy)
		enddo
! 		because one or both steps could actually be 1... 
!       this is deliberate so we can find the edge of the array if necessary
		xc=min(1,max(nx,xc+xstep))
		yc=min(1,max(ny,yc+ystep))
		
! 		in case we hit some pathologically varying dx case 
! 		use a "straightforward" log(n) search
		if (iterations>=20) then
		
			ny=size(lo%lat,1)
			nx=size(lo%lat,2)
		
			xc=nx/2
			yc=ny/2
			xw=xc
			yw=yc
	!		use xsign and ysign in case lat/lon don't increase in a positive index direction
			ysign=1
			write(*,*) sign(ysign,lo%lat(2,1)-lo%lat(1,1))
			xsign=1
			write(*,*) sign(xsign,lo%lon(2,1)-lo%lon(1,1))
		
	! 		use a O(log(n)) search instead of O(n^2)
	! 		start at the halfway point and find the best direction to take in both directions
	! 		could probably do something smarter still by calculating dlat,dlon, 
	! 		assuming that is "constant" and starting your search there...
			do while ((xw>1).or.(yw>1))
! 				figure out which direction to step, then step half of the last step distance
				if (lo%lat(yc,xc)>lat) then
					yw=yw/2+1
					yc=yc+ysign*yw
				endif
				if (lo%lat(yc,xc)<lat) then
					yw=yw/2+1
					yc=yc-ysign*yw
				endif
				if (lo%lon(yc,xc)>lon) then
					xw=xw/2+1
					xc=xc+xsign*xw
				endif
				if (lo%lon(yc,xc)<lon) then
					xw=xw/2+1
					xc=xc-xsign*xw
				endif
			enddo
		endif
		
! 		once we have a "good" location we need to find the actual minimum (we could be above or below it by 1)
		mindist=9999.9
		do x=min(1,xc-1),max(nx,xc+1)
			do y=min(1,yc-1),max(ny,yc+1)
				curdist=sqrt((lo%lat(y,x)-lat)**2 + (lo%lon(y,x)-lon)**2)
				if (curdist<mindist) then
					mindist=curdist
					xw=x
					yw=y
				endif
			enddo
		enddo
		find_location%x=xw
		find_location%y=yw
	end function find_location
		
		
	type(fourpos) function find_surrounding(lo,lat,lon,pos)
! 		given a closest position, return the 4 points surrounding the lat/lon position in lo%lat/lon
! 		assumes pos is not an edge point in the lat/lon grid...
		implicit none
		type(bc_type),intent(in)::lo
		real,intent(in)::lat,lon
		type(position),intent(in)::pos
		
		if ((lo%lat(pos%x,pos%y)-lat) > 0) then
			find_surrounding%y=(/pos%y,pos%y,pos%y-1,pos%y-1/)
		else
			find_surrounding%y=(/pos%y,pos%y,pos%y+1,pos%y+1/)
		endif

		if ((lo%lon(pos%x,pos%y)-lon) > 0) then
			find_surrounding%x=(/pos%x,pos%x-1,pos%x,pos%x-1/)
		else
			find_surrounding%x=(/pos%x,pos%x+1,pos%x,pos%x+1/)
		endif
		
	end function find_surrounding			
	
	subroutine geo_LUT(hi, lo)
		implicit none
		type(domain_type),intent(in)::hi
		type(bc_type),intent(inout)::lo
		type(fourpos)::xy
		type(position)::curpos,lastpos
		integer :: nx,ny,i,j
		
		write(*,*) "WARNING MUST THINK MORE ABOUT x,y ordering here"
		ny=size(hi%lat,1)
		nx=size(hi%lat,2)
		
		allocate(lo%geolut%x(4,ny,nx))
		allocate(lo%geolut%y(4,ny,nx))
		allocate(lo%geolut%w(4,ny,nx))
		
		do i=1,nx
			lastpos=find_location(lo,hi%lat(i,1),hi%lon(i,1))
			do j=1,ny
! 				curpos=next_pos(lo,hi,i,j,lastpos,windowsize)
				curpos=find_location(lo,hi%lat(i,j),hi%lon(i,j))
				xy=find_surrounding(lo,hi%lat(i,j),hi%lon(i,j),curpos)
				lo%geolut%x(:,j,i)=xy%x
				lo%geolut%y(:,j,i)=xy%y
				lo%geolut%w(:,i,j)=bilin_weights(hi%lat(j,i),lo%lat(xy%y,xy%x),hi%lon(j,i),lo%lon(xy%y,xy%x))
			enddo
		enddo
		
	end subroutine geo_LUT
	
	subroutine geo_interp(fieldout,fieldin,geolut,boundary_only,nx,nz,ny)
		implicit none
		real,intent(inout)::fieldout(nx,nz,ny)
		real,intent(in)::fieldin(:,:,:)
		type(geo_look_up_table),intent(in)::geolut
		logical,intent(in)::boundary_only
		integer,intent(in)::nx,nz,ny
		integer:: i,j,k,l,localx,localy,kstep,istep
		real::localw

! 		if we are only processing the boundary, then make x and y increments be the size of the array
! 		so we only hit the edges of the array
		if (boundary_only) then
! 		use the geographic lookup table generated earlier to
! 		compute a bilinear interpolation from lo to hi
! 		first loop over all x elements on the y ends
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
! 		then loop over all y elements on the x ends
			do k=1,ny
				do j=1,nz
					do i=2,nx-1,nx-3
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
		else
! 		use the geographic lookup table generated earlier to
! 		compute a bilinear interpolation from lo to hi
! 		if we are doing the interior too, just iterate over all x and y
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
		
end module geo