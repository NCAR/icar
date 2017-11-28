submodule(grid_interface) grid_implementation
    use assertions_mod,       only : assert, assertions
    implicit none


contains

    module function get_dims(this) result(dims)
        class(grid_t) :: this
        integer :: dims(3)

        dims(1) = this%ime - this%ims + 1
        dims(2) = this%kme - this%kms + 1
        dims(3) = this%jme - this%jms + 1
    end function

    ! note currently this is supplied in the domain object, however it should be moved to the grid object along with associated my_ny, etc routines
    ! anything that defines the domain grid should be in here...

    !> -------------------------------
    !! Decompose the domain into as even a set of tiles as possible in two dimensions
    !!
    !! Searches through possible numbers of x and y tiles that multiple evenly to
    !! give the total number of images requested.
    !!
    !! For each x/y split compute the number of grid cells in both dimensions in each tile
    !! return the split that provides the closest match between the number of x and y grid cells
    !!
    !! -------------------------------
    ! module subroutine domain_decomposition(this, nx, ny, nimages)
    !     class(domain_t), intent(inout) :: this
    !     integer,         intent(in)    :: nx, ny, nimages
    !     integer :: ysplit, xsplit, xs, ys, i
    !     real :: best, current, x, y
    !
    !     xsplit = 1
    !     ysplit = nimages
    !     xs = xsplit
    !     ys = ysplit
    !
    !     x = (nx/real(xsplit))
    !     y = (ny/real(ysplit))
    !
    !     if (y > x) then
    !         best = abs(1 - ( y / x ))
    !     else
    !         best = abs(1 - ( x / y ))
    !     endif
    !
    !     do i=nimages,1,-1
    !         if (mod(nimages,i)==0) then
    !             ysplit = i
    !             xsplit = nimages / i
    !
    !             x = (nx/float(xsplit))
    !             y = (ny/float(ysplit))
    !
    !             if (y > x) then
    !                 current = abs(1 - ( y / x ))
    !             else
    !                 current = abs(1 - ( x / y ))
    !             endif
    !
    !             if (current < best) then
    !                 best = current
    !                 xs = xsplit
    !                 ys = ysplit
    !             endif
    !         endif
    !     enddo
    !
    !     this%ximages = xs
    !     this%yimages = ys
    !
    !     this%ximg = mod(this_image(),  this%ximages)
    !     this%yimg =     this_image() / this%ximages
    !
    !     if (assertions) call assert((xs*ys) == nimages, "Number of tiles does not sum to number of images")
    !
    ! end subroutine domain_decomposition


end submodule
