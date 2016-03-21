!>------------------------------------------------------------
!!  Module to provide vertical interpolation
!!  includes setting up a vertical Look Up Table (vLUT)
!!  and performing vertical interpolation (vinterp)
!!
!!  Similar in concept to the geo module  
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module vertical_interpolation
    use data_structures
    implicit none
    
    private
    public vLUT
    public vLUT_forcing
    public vinterp
contains
    
    function weights(zin,ztop,zbot)
        ! Compute weights to interpolate between ztop and zbottom to get to zin
        implicit none
        real, intent(in) :: zin,ztop,zbot
        real, dimension(2) :: weights
        real :: zrange
        
        if (ztop==zbot) then
            weights(1)=0.5
            weights(2)=0.5
        else
            weights(1)=(zin-zbot)/(ztop-zbot)
            weights(2)=1-weights(1)
        endif
    end function weights
    
    ! Find the two points that border the input z in a column of z values
    !   if zin < z(1) find_match(1)=-1
    !   if zin > z(n) then find_match(1)=9999
    !   else zin>=z(find_match(1)) and zin<z(find_match(2))
    ! guess is an optional value that allows one to start searching z for zin in a good location
    function find_match(zin,z,guess)
        
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
                    if (i==n) then
                        find_match(2)=i
                    endif
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
    
    ! Compute the vertical interpolation look up table from a LOw-resolution forcing grid 
    ! to a HIgh-resolution model grid
    ! NOTE that the low-resolution grid has already been interpolated horizontally to the hi grid
    ! and the low-resolution grid may actually have more vertical resolution than the hi grid
    !
    ! The output vlut is stored in the "lo-res" forcing data structure (as with geolut)
    !
    subroutine vLUT(hi,lo)
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
                    guess=curpos(2)
                enddo !k=1,z
            enddo !i=1,x
        enddo !j=1,y
        
    end subroutine vLUT

    subroutine vLUT_forcing(hi,lo)
        ! identical to vLUT above, but for forcing data
        ! only change is that the vertical axis is the last axis
        ! instead of the middle axis
        implicit none
        class(interpolable_type), intent(in)    :: hi
        class(interpolable_type), intent(inout) :: lo
        integer::nx,ny,nz,i,j,k,guess,lo_nz
        integer,dimension(2) :: curpos
        real,dimension(2) :: curweights
        
        nx=size(hi%z,1)
        ny=size(hi%z,2)  ! difference from vLUT
        nz=size(hi%z,3)  ! difference from vLUT

        lo_nz=size(lo%z,3)   ! difference from vLUT
        
        if (allocated(lo%vert_lut%z)) then
            deallocate(lo%vert_lut%z, lo%vert_lut%w)
        endif
        allocate(lo%vert_lut%z(2,nx,ny,nz))   ! difference from vLUT
        allocate(lo%vert_lut%w(2,nx,ny,nz))   ! difference from vLUT
        do k=1,ny  ! difference from vLUT
            do i=1,nx
                guess=1
                do j=1,nz  ! difference from vLUT
                    
                    curpos=find_match(hi%z(i,k,j),lo%z(i,k,:),guess=guess)  ! difference from vLUT
                    if (curpos(1)>0) then
                        ! matched within the grid
                        curweights=weights(hi%z(i,k,j),lo%z(i,k,curpos(1)),lo%z(i,k,curpos(2)))  ! difference from vLUT
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
                        write(*,*) lo%z(i,k,:)  ! difference from vLUT
                        stop
                    endif
                    lo%vert_lut%z(:,i,k,j)=curpos
                    lo%vert_lut%w(:,i,k,j)=curweights
                    guess=curpos(2)
                enddo !j=1,z  ! difference from vLUT
            enddo !i=1,x
        enddo !k=1,y  ! difference from vLUT
        
    end subroutine vLUT_forcing


    
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
    
    subroutine vinterp(hi, lo, vlut, boundary_only, axis)
        implicit none
        real,dimension(:,:,:),    intent(inout) :: hi
        real,dimension(:,:,:),    intent(in)    :: lo
        class(vert_look_up_table),intent(in)    :: vlut
        logical, optional,        intent(in)    :: boundary_only
        integer, optional,        intent(in)    :: axis
        integer :: i,j,k, nx,ny,nz, zaxis
        
        if (present(boundary_only)) then
            if (boundary_only) then
                call vinterp_boundary(hi,lo,vlut)
                return
            endif
        endif
        
        if (present(axis)) then
            zaxis=axis
        else
            zaxis=2
        endif
        
        if (zaxis==2) then
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
            
        elseif (zaxis==3) then
            ! Wind arrays often have different x and y dimensions from the mass grid
            ! so use the lesser of the two (not a perfect interpolation for wind, but should capture most of it)
            ! axis=3 means we are doing a time_varying z interpolation
            nx=min(size(hi,1), size(vlut%z,2))
            ny=min(size(hi,2), size(vlut%z,3))
            nz=min(size(hi,3), size(vlut%z,4))
            do j=1,nz
                do k=1,ny
                    do i=1,nx
                        hi(i,k,j)=lo(i,k,vlut%z(1,i,k,j))*vlut%w(1,i,k,j) + lo(i,k,vlut%z(2,i,k,j))*vlut%w(2,i,k,j)
                    enddo
                enddo
            enddo
        else
            write(*,*) "Vertical interpolation over the first axis not supported yet"
            write(*,*) "  if needed, update vinterp.f90"
            stop
        endif
        
    end subroutine vinterp
    
end module vertical_interpolation
