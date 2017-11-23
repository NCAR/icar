program test_options

    use options_interface,  only : options_t

    implicit none

    type(options_t) :: options
    integer :: i

    call options%init()

    ! if (this_image()==1) then
    !     print*, this_image(), trim(options%comment)
    !     print*, this_image(), trim(options%parameters%comment)
    ! endif
    !
    sync all

    do i=1,num_images()
        if (this_image()==i) then
            print*, this_image(), trim(options%comment)
            print*, this_image(), trim(options%parameters%comment)
        endif
        sync all
    end do

end program test_options
