program test_caf

    implicit none

    integer, parameter :: MAX_LIST_SIZE=8

    type dataset
        integer :: i
        real :: list(MAX_LIST_SIZE)
        real, allocatable :: list2(:)
    end type dataset

    type(dataset), allocatable :: test[:]
    integer :: i

    allocate(test[*])
    allocate(test%list2(8))

    if (this_image()==1) then
        print*, this_image(), num_images()
        test%i=1024
        ! allocate(test%list2(8))
        test%list = 10
        test%list2 = 7
        print*, size(test%list2), test%list2(1)
    endif

    sync all
    if (this_image()>1) then

        test = test[1]

        ! test%list2(:) = test[1]%list2(:)

        print*, this_image(), test%i, size(test%list2), test%list2(1)
        print*, this_image(), test%i, size(test%list), test%list(1)
    endif

    print*, this_image()

    deallocate(test)

    print*, this_image(), " Deallocated"
end program test_caf
