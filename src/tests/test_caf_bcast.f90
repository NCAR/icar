program caf_bcast_test

    use co_util, only : broadcast

    implicit none

    call test_co_4dd()

    call test_4dr()

    print*, "passed"


contains

    subroutine test_4dr
        implicit none

        real, allocatable :: array_data(:,:,:,:)
        real, allocatable :: test_data(:,:,:,:)

        integer :: nx,ny,nz,nt, i, j, k, l
        integer :: source_image

        nx = 3
        ny = 5
        nz = 6
        nt = 9

        allocate(array_data(nx,ny,nz,nt))
        allocate(test_data(nx,ny,nz,nt))

        array_data = 1

        ! all images calculate the test data for comparison after the broadcast
        do l=1,nt
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        test_data(i,j,k,l) = l/10.0+(k+(j+(i*10))*10)
                    enddo
                enddo
            enddo
        enddo

        ! pick the image to use as a source image and copy the test data into the array to be broadcast
        source_image = max(1,(num_images()/2 + 1))
        if (this_image() == source_image) then
            array_data = test_data
        endif

        ! broadcast array from source image to all images between 1 and num_images (ie all images)
        call broadcast(array_data, source_image, 1, num_images(), .true.)

        do l=1,nt
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        if (array_data(i,j,k,l)/=test_data(i,j,k,l)) then
                            print*, this_image(), " failed on :",i,j,k,l
                            print*, array_data(i,j,k,l), " should be",test_data(i,j,k,l)
                            stop
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine



    subroutine test_co_4dd
        implicit none

        real(kind=8), allocatable :: array_data(:,:,:,:)[:]
        real(kind=8), allocatable :: test_data(:,:,:,:)

        integer :: nx,ny,nz,nt, i, j, k, l
        integer :: source_image

        nx = 3
        ny = 5
        nz = 6
        nt = 9

        allocate(array_data(nx,ny,nz,nt)[*])
        allocate(test_data(nx,ny,nz,nt))

        array_data = 1

        ! all images calculate the test data for comparison after the broadcast
        do l=1,nt
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        test_data(i,j,k,l) = l/10.0+(k+(j+(i*10))*10)
                    enddo
                enddo
            enddo
        enddo

        ! pick the image to use as a source image and copy the test data into the array to be broadcast
        source_image = max(1,(num_images()/2 + 1))
        if (this_image() == source_image) then
            array_data = test_data
        endif

        ! broadcast array from source image to all images between 1 and num_images (ie all images)
        call broadcast(array_data, source_image, 1, num_images())

        do l=1,nt
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        if (array_data(i,j,k,l)/=test_data(i,j,k,l)) then
                            print*, this_image(), " failed on :",i,j,k,l
                            print*, array_data(i,j,k,l), " should be",test_data(i,j,k,l)
                            stop
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine

end program
