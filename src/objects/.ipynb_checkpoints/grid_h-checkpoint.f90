module grid_interface

    use icar_constants, only : kDEFAULT_HALO_SIZE, kMAX_DIM_LENGTH

    implicit none

    private
    public :: grid_t

    type grid_t
        integer :: yimg,    ximg
        integer :: yimages, ximages
        integer :: ims, ime
        integer :: jms, jme
        integer :: kms, kme
        integer :: ns_halo_nx, ew_halo_ny, halo_nz, halo_size
        integer :: nx_global, ny_global
        integer :: nx, ny, nz

        integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                    its,ite, jts,jte, kts,kte    ! for the data tile to process   (t)

        logical :: is2d, is3d
        character(len=kMAX_DIM_LENGTH), allocatable :: dimensions(:)

    contains
        procedure :: get_dims
        procedure :: domain_decomposition
        procedure :: set_grid_dimensions

    end type

interface
    module function get_dims(this) result(dims)
        implicit none
        class(grid_t), intent(in) :: this
        integer, allocatable :: dims(:)
    end function

    module subroutine domain_decomposition(this, nx, ny, nimages, ratio, for_image)
        implicit none
        class(grid_t),  intent(inout) :: this
        integer,        intent(in)    :: nx, ny, nimages
        real,           intent(in), optional :: ratio
        integer,        intent(in), optional :: for_image
    end subroutine

    module subroutine set_grid_dimensions(this, nx, ny, nz, nx_extra, ny_extra, halo_width, for_image)
        implicit none
        class(grid_t),   intent(inout) :: this
        integer,         intent(in)    :: nx, ny, nz
        integer,         intent(in), optional :: nx_extra, ny_extra, halo_width, for_image
    end subroutine
end interface
end module
