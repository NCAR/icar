module co_util
    ! Note that gfortran 7 does not support coarray "put"s (in which data are sent to a remote image)
    ! as such the current code uses "get"s in which data are fetched from a remote image.
    ! puts have potential to run faster than gets because they can operate more asynchronously
    ! It does not matter for some or all of these routines as they are already doing a sync before continuing.

    use icar_constants, only : kMAX_STRING_LENGTH

    implicit none

    !>------------------------------------------------------------
    !! Generic interface to various broadcast routines
    !!------------------------------------------------------------
    interface broadcast
        module procedure co_bcast_4dd, co_bcast_3dd, co_bcast_2dd, co_bcast_1dd, & ! double co-arrays
                         co_bcast_4dr, co_bcast_3dr, co_bcast_2dr, co_bcast_1dr, & ! real co-arrays
                            bcast_4dr,    bcast_3dr,    bcast_2dr,    bcast_1dr, & ! real arrays
                         co_bcast_i,      bcast_i,   co_bcast_1di,    bcast_1di, & ! int scalars and 1d arrays and co-arrays
                         co_bcast_r,      bcast_r, & ! 1dr is taken care of above with multi-dim reals
                         co_bcast_c,      bcast_c,   co_bcast_1dc,    bcast_1dc, & ! character scalars and 1d arrays and co-arrays
                         co_bcast_l,      bcast_l                                  ! logical scalars and co-array (scalars)
    end interface

contains
    ! Note, this first routine is included completely, all variants of this are almost identical routines based on include files
    ! the only parts left in are the parts that differ between routines
    recursive subroutine co_bcast_4dr(coarray, source_image, first_image, last_image)
        implicit none
        real, intent(inout) :: coarray(:,:,:,:)[*]
        integer, intent(in) :: source_image, first_image, last_image
        integer :: dest_image

        !call co_broadcast(coarray,source_image)
        !return

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
                call broadcast(coarray, dest_image, dest_image, source_image-1)
            else
                call broadcast(coarray, source_image, source_image, last_image)
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
                call broadcast(coarray, source_image, first_image, dest_image-1)
            else
                call broadcast(coarray, dest_image, dest_image, last_image)
            endif
        endif

    end subroutine

    recursive subroutine co_bcast_3dr(coarray, source_image, first_image, last_image)
        implicit none
        real, intent(inout) :: coarray(:,:,:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:,:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:,:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_2dr(coarray, source_image, first_image, last_image)
        implicit none
        real, intent(inout) :: coarray(:,:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_1dr(coarray, source_image, first_image, last_image)
        implicit none
        real, intent(inout) :: coarray(:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine



    recursive subroutine co_bcast_4dd(coarray, source_image, first_image, last_image)
        implicit none
        real(kind=8), intent(inout) :: coarray(:,:,:,:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:,:,:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:,:,:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_3dd(coarray, source_image, first_image, last_image)
        implicit none
        real(kind=8), intent(inout) :: coarray(:,:,:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:,:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:,:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_2dd(coarray, source_image, first_image, last_image)
        implicit none
        real(kind=8), intent(inout) :: coarray(:,:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:,:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_1dd(coarray, source_image, first_image, last_image)
        implicit none
        real(kind=8), intent(inout) :: coarray(:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_1di(coarray, source_image, first_image, last_image)
        implicit none
        integer, intent(inout) :: coarray(:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine


    recursive subroutine co_bcast_1dc(coarray, source_image, first_image, last_image)
        implicit none
        character(len=*), intent(inout) :: coarray(:)[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray(:)[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray(:)[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine


    recursive subroutine co_bcast_c(coarray, source_image, first_image, last_image)
        implicit none
        character(len=*), intent(inout) :: coarray[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_l(coarray, source_image, first_image, last_image)
        implicit none
        logical, intent(inout) :: coarray[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_r(coarray, source_image, first_image, last_image)
        implicit none
        real, intent(inout) :: coarray[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine

    recursive subroutine co_bcast_i(coarray, source_image, first_image, last_image)
        implicit none
        integer, intent(inout) :: coarray[*]
        INCLUDE 'broadcast_core/bcast_part1.inc'
                coarray = coarray[source_image]
        INCLUDE 'broadcast_core/bcast_part2.inc'
                coarray = coarray[source_image]
        INCLUDE 'broadcast_core/bcast_part3.inc'
    end subroutine


    subroutine bcast_4dr(data_array, source_image, first_image, last_image, create_co_array)
        implicit none
        real,         intent(inout) :: data_array(:,:,:,:)
        integer,      intent(in)    :: source_image, first_image, last_image
        logical,      intent(in)    :: create_co_array

        real, allocatable   :: coarray(:,:,:,:)[:]
        integer :: n1,n2,n3,n4

        !call co_broadcast(data_array, source_image)
        !return

        n1 = size(data_array, 1)
        n2 = size(data_array, 2)
        n3 = size(data_array, 3)
        n4 = size(data_array, 4)

        allocate(coarray(n1,n2,n3,n4)[*], source=data_array)

        call broadcast(coarray, source_image, first_image, last_image)

        data_array=coarray

        deallocate(coarray)

    end subroutine


    subroutine bcast_3dr(data_array, source_image, first_image, last_image, create_co_array)
        implicit none
        real,         intent(inout) :: data_array(:,:,:)
        integer,      intent(in)    :: source_image, first_image, last_image
        logical,      intent(in)    :: create_co_array

        real, allocatable   :: coarray(:,:,:)[:]
        integer :: n1,n2,n3

        !call co_broadcast(data_array, source_image)
        !return

        n1 = size(data_array, 1)
        n2 = size(data_array, 2)
        n3 = size(data_array, 3)

        allocate(coarray(n1,n2,n3)[*], source=data_array)

        call broadcast(coarray, source_image, first_image, last_image)

        data_array=coarray

        deallocate(coarray)

    end subroutine


    subroutine bcast_2dr(data_array, source_image, first_image, last_image, create_co_array)
        implicit none
        real,         intent(inout) :: data_array(:,:)
        integer,      intent(in)    :: source_image, first_image, last_image
        logical,      intent(in)    :: create_co_array

        real, allocatable   :: coarray(:,:)[:]
        integer :: n1,n2

        !call co_broadcast(data_array, source_image)
        !return

        n1 = size(data_array, 1)
        n2 = size(data_array, 2)

        allocate(coarray(n1,n2)[*], source=data_array)

        call broadcast(coarray, source_image, first_image, last_image)

        data_array=coarray

        deallocate(coarray)

    end subroutine

    subroutine bcast_1dr(data_array, source_image, first_image, last_image, create_co_array)
        implicit none
        real,         intent(inout) :: data_array(:)
        integer,      intent(in)    :: source_image, first_image, last_image
        logical,      intent(in)    :: create_co_array

        real, allocatable   :: coarray(:)[:]
        integer :: n1

        !call co_broadcast(data_array, source_image)
        !return

        n1 = size(data_array, 1)

        allocate(coarray(n1)[*], source=data_array)

        call broadcast(coarray, source_image, first_image, last_image)

        data_array=coarray

        deallocate(coarray)

    end subroutine

    subroutine bcast_1di(data_array, source_image, first_image, last_image, create_co_array)
        implicit none
        integer,      intent(inout) :: data_array(:)
        integer,      intent(in)    :: source_image, first_image, last_image
        logical,      intent(in)    :: create_co_array

        integer, allocatable   :: coarray(:)[:]
        integer :: n1

        !call co_broadcast(data_array, source_image)
        !return

        n1 = size(data_array, 1)

        allocate(coarray(n1)[*], source=data_array)

        call broadcast(coarray, source_image, first_image, last_image)

        data_array=coarray

        deallocate(coarray)

    end subroutine

    subroutine bcast_1dc(data_array, source_image, first_image, last_image, create_co_array)
        implicit none
        character(len=*),   intent(inout) :: data_array(:)
        integer,            intent(in)    :: source_image, first_image, last_image
        logical,            intent(in)    :: create_co_array

        character(len=kMAX_STRING_LENGTH), allocatable   :: coarray(:)[:]
        integer :: n1

        !call co_broadcast(data_array, source_image)
        !return

        n1 = size(data_array, 1)

        allocate(coarray(n1)[*], source=data_array)

        call broadcast(coarray, source_image, first_image, last_image)

        data_array=coarray

        deallocate(coarray)

    end subroutine

    subroutine bcast_r(scalar, source_image, first_image, last_image, create_co_array)
        implicit none
        real,         intent(inout) :: scalar
        integer,      intent(in)    :: source_image, first_image, last_image
        logical,      intent(in)    :: create_co_array

        real, allocatable :: coscalar[:]

        allocate(coscalar[*])
        coscalar = scalar

        !call co_broadcast(coscalar, source_image)
        !return

        call broadcast(coscalar, source_image, first_image, last_image)

        scalar = coscalar
        deallocate(coscalar)

    end subroutine

    subroutine bcast_i(scalar, source_image, first_image, last_image, create_co_array)
        implicit none
        integer,      intent(inout) :: scalar
        integer,      intent(in)    :: source_image, first_image, last_image
        logical,      intent(in)    :: create_co_array

        integer, allocatable :: coscalar[:]

        allocate(coscalar[*])
        coscalar = scalar

        !call co_broadcast(coscalar, source_image)
        !return

        call broadcast(coscalar, source_image, first_image, last_image)

        scalar = coscalar
        deallocate(coscalar)
    end subroutine


    subroutine bcast_l(scalar, source_image, first_image, last_image, create_co_array)
        implicit none
        logical, intent(inout) :: scalar
        integer, intent(in)    :: source_image, first_image, last_image
        logical, intent(in)    :: create_co_array

        logical, allocatable :: coscalar[:]

        allocate(coscalar[*])
        coscalar = scalar

        !call co_broadcast(coscalar, source_image)
        !return

        call broadcast(coscalar, source_image, first_image, last_image)

        scalar = coscalar
        deallocate(coscalar)
    end subroutine


    subroutine bcast_c(scalar, source_image, first_image, last_image, create_co_array)
        implicit none
        character(len=*), intent(inout) :: scalar
        integer,          intent(in)    :: source_image, first_image, last_image
        logical,          intent(in)    :: create_co_array

        character(len=kMAX_STRING_LENGTH), allocatable :: coscalar[:]

        allocate(coscalar[*])
        coscalar = scalar

        !call co_broadcast(coscalar, source_image)
        !return

        call broadcast(coscalar, source_image, first_image, last_image)

        scalar = coscalar
        deallocate(coscalar)
    end subroutine

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
