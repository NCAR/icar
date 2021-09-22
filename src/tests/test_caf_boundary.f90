program test_caf_boundary

    use options_interface,  only : options_t
    use boundary_interface, only : boundary_t
    use variable_interface, only : variable_t

    implicit none

    type(options_t) :: options
    type(boundary_t):: boundary

    type(variable_t):: var
    character(len=1024) :: name
    integer :: i

    ! read in the config file (gets boundary forcing filenames and variable names)
    call options%init()

    ! Initialize the boundary structure (reads data and distributes to all images)
    call boundary%init(options)

    ! loop over images and variables printing out the values in one gridcell to make sure they are not junk (and that it doesn't crash with bounds checking on)
    ! do i=1,num_images()
    !     if (i==this_image()) then
    !
    !         associate(list => boundary%variables)
    !         ! iterate through the variables in the variables dictionary
    !         call list%reset_iterator()
    !         do while (list%has_more_elements())
    !             ! get the next variable in the structure
    !             var = list%next(name)
    !             ! print out a data point from either the 2D or 3D array depending on what type of data these are.
    !             if (var%three_d) print*, this_image(), trim(name), var%data_3d(2,2,2)
    !             if (var%two_d)   print*, this_image(), trim(name), var%data_2d(2,2)
    !         end do
    !
    !         end associate
    !
    !     endif
    !     sync all
    ! enddo


    call boundary%update_forcing(options)

    ! do i=1,num_images()
    !     if (i==this_image()) then
    !
    !         associate(list => boundary%variables)
    !         ! iterate through the variables in the variables dictionary
    !         call list%reset_iterator()
    !         do while (list%has_more_elements())
    !             ! get the next variable in the structure
    !             var = list%next(name)
    !             ! print out a data point from either the 2D or 3D array depending on what type of data these are.
    !             if (var%three_d) print*, this_image(), trim(name), var%data_3d(2,2,2)
    !             if (var%two_d)   print*, this_image(), trim(name), var%data_2d(2,2)
    !         end do
    !
    !         end associate
    !
    !     endif
    !     sync all
    ! enddo

    sync all

    if (this_image()==1) print*, "Done"


end program
