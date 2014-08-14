module vertical_interpolation
	use data_structures
	implicit none
	
	private
	public vLUT
	public vinterp
contains
	
	function weights(zin,ztop,zbot)
		! Compute weights to interpolate between ztop and zbottom to get to zin
		implicit none
		real, intent(in) :: zin,ztop,zbot
		real, dimension(2) :: weights
		real :: zrange
		
		weights(1)=(zin-zbot)/(ztop-zbot)
		weights(2)=1-weights(1)
	end function weights
	
	function find_match(zin,z,guess)
		! Find the two points that border the input z in a column of z values
		! 	if zin < z(1) find_match(1)=-1
		! 	if zin > z(n) then find_match(1)=9999
		! 	else zin>=z(find_match(1)) and zin<z(find_match(2))
		! guess is an optional value that allows one to start searching z for zin in a good location
		
		implicit none
		real, intent(in) :: zin
		real, intent(in),dimension(:) :: z
		integer, optional, intent(inout)::guess
		real,dimension(2) :: find_match
		integer::n,i,endpt
		
		n=size(z)
		find_match(1)=-1
		
		if (.not.present(guess)) then
			guess=n/2
		endif
		
		if (z(guess)>=zin) then
			! then we should search downward from z(guess)
			endpt=1
			i=guess
			do while ((i>=endpt).and.(find_match(1)==-1))
				if (z(i)<=zin) then
					find_match(1)=i
					find_match(2)=i+1
				endif
				i=i-1
			end do
		else
			! then we should search upward from z(guess)
			endpt=n
			i=guess
			do while ((i<=endpt).and.(find_match(1)==-1))
				if (z(i)>zin) then
					find_match(1)=i-1
					find_match(2)=i
				endif
				i=i+1
			end do
			
			if (find_match(1)==-1) then
				find_match(1)=-2
			endif
		endif
		
	end function find_match
	
	subroutine vLUT(hi,lo)
		! Compute the vertical interpolation look up table from a LOw-resolution forcing grid 
		! to a HIgh-resolution model grid
		! NOTE that the low-resolution grid has already been interpolated horizontally to the hi grid
		! and the low-resolution grid may actually have more vertical resolution than the hi grid
		!
		! The output vlut is stored in the "lo-res" forcing data structure (as with geolut)
		!
		implicit none
		class(interpolable_type), intent(in)    :: hi
		class(interpolable_type), intent(inout) :: lo
		integer::nx,ny,nz,i,j,k,guess,lo_nz
		integer,dimension(2) :: curpos
		real,dimension(2) :: curweights
		
		nx=size(hi%z,1)
		nz=size(hi%z,2)
		ny=size(hi%z,3)

		lo_nz=size(lo%z,2)
		
		allocate(lo%vert_lut%z(2,nx,nz,ny))
		allocate(lo%vert_lut%w(2,nx,nz,ny))
		do j=1,ny
			do i=1,nx
				guess=1
				do k=1,nz
					
					curpos=find_match(hi%z(i,k,j),lo%z(i,:,j),guess=guess)
					if (curpos(1)>0) then
						! matched within the grid
						curweights=weights(hi%z(i,k,j),lo%z(i,curpos(1),j),lo%z(i,curpos(2),j))
					elseif (curpos(1)==-1) then
						! matched below the grid
						curpos(1)=1
						curpos(2)=1
						curweights=0.5
					elseif (curpos(1)==-2) then
						! matched above the grid
						curpos(1)=lo_nz
						curpos(2)=lo_nz
						curweights=0.5
					else
						write(*,*) "find_match Failed to return appropriate position"
						write(*,*) " at grid location:"
						write(*,*) i,k,j
						write(*,*) "z to match = ",hi%z(i,k,j)
						write(*,*) "from Z-column="
						write(*,*) lo%z(i,:,j)
						stop
					endif
					lo%vert_lut%z(:,i,k,j)=curpos
					lo%vert_lut%w(:,i,k,j)=curweights
! 					if ((k==1).and.(curpos(1)>20)) then
! 						print*, "i=",i
! 						print*, "j=",j
! 						print*, "pos=",curpos
! 						print*, "weights=",curweights
! 						print*, "z to match = ",hi%z(i,k,j)
! 						print*, "from Z-column:"
! 						print*, lo%z(i,16:25,j)
! 					endif
					guess=curpos(2)
				enddo !k=1,z
			enddo !i=1,x
		enddo !j=1,y
		
	end subroutine vLUT
	
	subroutine vinterp_boundary(hi,lo,vlut)
		implicit none
		real,dimension(:,:,:), intent(inout) :: hi
		real,dimension(:,:,:), intent(in)    :: lo
		class(vert_look_up_table),intent(in) :: vlut
		integer :: i,j,k,nx,ny,nz

		nx=size(hi,1)
		nz=size(hi,2)
		ny=size(hi,3)
		
		j=1
		do k=1,nz
			do i=1,nx
				hi(i,k,j)=lo(i,vlut%z(1,i,k,j),j)*vlut%w(1,i,k,j) + lo(i,vlut%z(2,i,k,j),j)*vlut%w(2,i,k,j)
			enddo
		enddo
		j=ny
		do k=1,nz
			do i=1,nx
				hi(i,k,j)=lo(i,vlut%z(1,i,k,j),j)*vlut%w(1,i,k,j) + lo(i,vlut%z(2,i,k,j),j)*vlut%w(2,i,k,j)
			enddo
		enddo
		
		i=1
		do j=1,ny
			do k=1,nz
				hi(i,k,j)=lo(i,vlut%z(1,i,k,j),j)*vlut%w(1,i,k,j) + lo(i,vlut%z(2,i,k,j),j)*vlut%w(2,i,k,j)
			enddo
		enddo
		i=nx
		do j=1,ny
			do k=1,nz
				hi(i,k,j)=lo(i,vlut%z(1,i,k,j),j)*vlut%w(1,i,k,j) + lo(i,vlut%z(2,i,k,j),j)*vlut%w(2,i,k,j)
			enddo
		enddo
		
	end subroutine vinterp_boundary
	
	subroutine vinterp(hi,lo,vlut,boundary_only)
		implicit none
		real,dimension(:,:,:), intent(inout) :: hi
		real,dimension(:,:,:), intent(in)    :: lo
		class(vert_look_up_table),intent(in) :: vlut
		logical,optional,intent(in)::boundary_only
		integer :: i,j,k,nx,ny,nz
		
		if (present(boundary_only)) then
			if (boundary_only) then
				call vinterp_boundary(hi,lo,vlut)
				return
			endif
		endif
		
		nx=size(hi,1)
		nz=size(hi,2)
		ny=size(hi,3)
		
		do j=1,ny
			do k=1,nz
				do i=1,nx
					hi(i,k,j)=lo(i,vlut%z(1,i,k,j),j)*vlut%w(1,i,k,j) + lo(i,vlut%z(2,i,k,j),j)*vlut%w(2,i,k,j)
				enddo
			enddo
		enddo
		
	end subroutine vinterp
	
end module vertical_interpolation