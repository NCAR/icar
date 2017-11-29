module co_util
    ! Note that gfortran 7 does not support coarray "put"s (in which data are sent to a remote image)
    ! as such the current code uses "get"s in which data are fetched from a remote image.
    ! puts have potential to run faster than gets because they can operate more asynchronously
    ! It does not matter for some or all of these routines as they are already doing a sync before continuing.

    implicit none

contains
    recursive subroutine co_bcast(coarray, source_image, first_image, last_image)
        implicit none
        real(kind=8), intent(inout) :: coarray(:,:,:,:)[*]
        integer, intent(in) :: source_image, first_image, last_image
        integer :: dest_image

        if (first_image==last_image) return

        if (source_image/=first_image) then
            dest_image=first_image

            if (this_image()==source_image) then
                ! This is a "put" and needs to come before the sync
                ! coarray(:,:,:,:)[dest_image] = coarray
                sync images(dest_image)
            elseif (this_image()== dest_image) then
                ! This is a "get" and needs to come after the sync
                sync images(source_image)
                coarray = coarray(:,:,:,:)[source_image]
            endif

            if (this_image()<source_image) then
                call co_bcast(coarray, dest_image, dest_image, source_image-1)
            else
                call co_bcast(coarray, source_image, source_image, last_image)
            endif
        else
            dest_image = ((last_image-first_image)+1)/2 + first_image

            if (this_image()==source_image) then
                ! This is a "put" and needs to come before the sync
                ! coarray(:,:,:,:)[dest_image] = coarray
                sync images(dest_image)
            elseif (this_image()== dest_image) then
                ! This is a "get" and needs to come after the sync
                sync images(source_image)
                coarray = coarray(:,:,:,:)[source_image]
            endif

            if (this_image()<dest_image) then
                call co_bcast(coarray, source_image, first_image, dest_image-1)
            else
                call co_bcast(coarray, dest_image, dest_image, last_image)
            endif
        endif

    end subroutine co_bcast

    subroutine print_in_image_order(input)
        implicit none
        real, intent(in) :: input(:,:)

        integer :: i, j
        integer :: nx, xstep
        integer :: ny, ystep

        nx = size(input,1)
        ny = size(input,2)

        ! print 8 rows with 8 columns each from an arbitrary sized array
        ystep = ny / 8
        xstep = nx / 8

        do i=1,num_images()
            if (this_image()==i) then
                write(*,*) this_image()

                ! loop through the 8 rows printing every 8th column
                do j=lbound(input,2),ubound(input,2),ystep
                    write(*,*) input(::xstep,j)
                enddo
            endif
            ! after a given image has printed it's array contents, sync all to make all images wait
            sync all
        end do

    end subroutine print_in_image_order

end module co_util
