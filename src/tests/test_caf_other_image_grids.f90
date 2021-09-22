program other_image_grids

    use grid_interface, only : grid_t
    implicit none

    type(grid_t) :: grid
    integer,     allocatable :: test_data(:,:)[:]
    integer :: these_dims(6)

    integer :: i
    integer :: nx = 1024
    integer :: ny = 1234
    integer :: nz = 13
    logical :: passed[*] = .True.

    allocate(test_data(6, num_images())[*])

    do i = 1, num_images()
        call grid%set_grid_dimensions(nx, ny, nz, for_image=i)
        test_data(:,i) = [grid%ims, grid%ime, &
                          grid%kms, grid%kme, &
                          grid%jms, grid%jme]
    end do

    call grid%set_grid_dimensions(nx, ny, nz)

    these_dims = [grid%ims, grid%ime, &
                  grid%kms, grid%kme, &
                  grid%jms, grid%jme]

    sync all

    do i = 1, num_images()
        if (any( test_data(:,this_image())[i] /= these_dims )) then
            print*, "TEST_FAILED:"
            print*, this_image(), these_dims
            print*, i, test_data(:,this_image())[i]
            passed = .False.
        endif
    enddo

    sync all

    if (this_image() == 1) then
        do i = 1, num_images()
            passed = passed .and. passed[i]
        enddo
        if (passed) print*, "Test Passed"
    endif


end program other_image_grids
