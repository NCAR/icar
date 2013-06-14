module geo
	use data_structures
	implicit none
	
	private
	public::geo_LUT
	public::geo_interp,geo_interp2d
	
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
		real,intent(in)::yi,y(0:3),xi,x(0:3)
		real::x0,x1,x2,x3,y5,y6,f1,f2
		real, dimension(4) ::bilin_weights
		if ((x(1)-x(0))==0) then
			x0=1
		else
		    x0=abs((xi-x(0))/(x(1)-x(0)))
		endif
	    x1=1-x0
		if ((x(3)-x(2))==0) then
			x2=1
		else
		    x2=abs((xi-x(2))/(x(3)-x(2)))
		endif
	    x3=1-x2
	    y5=y(0)*x1+y(1)*x0
	    y6=y(2)*x3+y(3)*x2
		if ((y6-y5)==0) then
			f1=1
		else
		    f1=(yi-y5)/(y6-y5)
		endif
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
		
		nx=size(lo%lat,1)
		ny=size(lo%lat,2)
! 		calcualte dx/dy at the middle of the grid
		dx=lo%lon(2,1)-lo%lon(1,1)
		dy=lo%lat(1,2)-lo%lat(1,1)
! 		current/starting position = the middle of the grid
		xc=nx/2
		yc=ny/2
		x=lo%lon(xc,yc)
		y=lo%lat(xc,yc)
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
			xc=max(1,min(nx-1,xc+xstep))
			yc=max(1,min(ny-1,yc+ystep))
			x=lo%lon(xc,yc)
			y=lo%lat(xc,yc)
! 			calculate a new dx/dy and x/y step
			dx=lo%lon(xc+1,yc)-lo%lon(xc,yc)
			dy=lo%lat(xc,yc+1)-lo%lat(xc,yc)
			xstep=NINT((lon-x)/dx)
			ystep=NINT((lat-y)/dy)
		enddo
! 		because one or both steps could actually be 1... 
!       this is deliberate so we can find the edge of the array if necessary
		xc=max(1,min(nx,xc+xstep))
		yc=max(1,min(ny,yc+ystep))
		
! 		in case we hit some pathologically varying dx case 
! 		use a "straightforward" log(n) search
		if (iterations>=20) then
			write(*,*) "   Using log(n) search for :",lat,lon
			nx=size(lo%lat,1)
			ny=size(lo%lat,2)
		
			xc=nx/2
			yc=ny/2
			xw=xc
			yw=yc
	!		use xsign and ysign in case lat/lon don't increase in a positive index direction
			ysign=sign(1.0,lo%lat(1,2)-lo%lat(1,1))
			xsign=sign(1.0,lo%lon(2,1)-lo%lon(1,1))
		
	! 		use a O(log(n)) search instead of O(nxm)
	! 		start at the halfway point and find the best direction to take in both directions
			do while ((xw>1).or.(yw>1))
! 				figure out which direction to step, then step half of the last step distance
				if (lo%lat(xc,yc)>lat) then
					if (yw==2) then
						yw=1
					else
						yw=yw/2+1
					endif
					yc=yc-ysign*yw
				endif
				if (lo%lat(xc,yc)<lat) then
					if (yw==2) then
						yw=1
					else
						yw=yw/2+1
					endif
					yc=yc+ysign*yw
				endif
				if (lo%lon(xc,yc)>lon) then
					if (xw==2) then
						xw=1
					else
						xw=xw/2+1
					endif
					xc=xc-xsign*xw
				endif
				if (lo%lon(xc,yc)<lon) then
					if (xw==2) then
						xw=1
					else
						xw=xw/2+1
					endif
					xc=xc+xsign*xw
				endif
			enddo
		endif
		
! 		once we have a "good" location we need to find the actual minimum (we could be above or below it by 1 or 2)
		mindist=9999.9
		do x=max(1,xc-5),min(nx,xc+5)
			do y=max(1,yc-5),min(ny,yc+5)
				curdist=sqrt((lo%lat(x,y)-lat)**2 + (lo%lon(x,y)-lon)**2)
				if (curdist<mindist) then
					mindist=curdist
					xw=x
					yw=y
				endif
			enddo
		enddo

! naive version that may actually be fast enough...
! 		mindist=9999.9
! 		do x=1,nx
! 			do y=1,ny
! 				curdist=sqrt((lo%lat(x,y)-lat)**2 + (lo%lon(x,y)-lon)**2)
! 				if (curdist<mindist) then
! 					mindist=curdist
! 					xw=x
! 					yw=y
! 				endif
! 			enddo
! 		enddo

		
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
		integer :: nx,ny,i,j,k
		real,dimension(4) :: lat,lon
		
		nx=size(hi%lat,1)
		ny=size(hi%lat,2)
		
		allocate(lo%geolut%x(4,nx,ny))
		allocate(lo%geolut%y(4,nx,ny))
		allocate(lo%geolut%w(4,nx,ny))
		
		do j=1,ny
! 			lastpos=find_location(lo,hi%lat(1,j),hi%lon(1,j))
			do i=1,nx
! 				curpos=next_pos(lo,hi,i,j,lastpos,windowsize)
				curpos=find_location(lo,hi%lat(i,j),hi%lon(i,j))
				xy=find_surrounding(lo,hi%lat(i,j),hi%lon(i,j),curpos)
				lo%geolut%x(:,i,j)=xy%x
				lo%geolut%y(:,i,j)=xy%y
				do k=1,4
					lat(k)=lo%lat(xy%x(k),xy%y(k))
					lon(k)=lo%lon(xy%x(k),xy%y(k))
				enddo
				lo%geolut%w(:,i,j)=bilin_weights(hi%lat(i,j),lat,hi%lon(i,j),lon)
			enddo
		enddo
		
	end subroutine geo_LUT
	
	subroutine geo_interp(fieldout,fieldin,geolut,boundary_only)
		implicit none
		real,intent(inout)::fieldout(:,:,:)
		real,intent(in)::fieldin(:,:,:)
		type(geo_look_up_table),intent(in)::geolut
		logical,intent(in)::boundary_only
		integer::nx,nz,ny
		integer:: i,j,k,l,localx,localy,kstep,istep
		real::localw
		
		nx=size(fieldout,1)
		nz=size(fieldout,2)
		ny=size(fieldout,3)
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
	
	subroutine geo_interp2d(fieldout, fieldin, geolut)
		real, dimension(:,:),intent(out) :: fieldout
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

	
end module geo