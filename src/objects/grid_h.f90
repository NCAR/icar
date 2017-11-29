module grid_interface
    implicit none

    private
    public :: grid_t

    type grid_t
        integer :: yimg,    ximg
        integer :: yimages, ximages
        integer :: ims, ime
        integer :: jms, jme
        integer :: kms, kme
        integer :: ns_halo_nx, ew_halo_ny, halo_nz
        integer :: nx_global, ny_global
        integer :: nx, ny, nz

        integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                    its,ite, jts,jte, kts,kte    ! for the data tile to process   (t)

    contains
        procedure :: get_dims
        procedure :: domain_decomposition
        procedure :: get_grid_dimensions

    end type

interface
    module function get_dims(this) result(dims)
        implicit none
        class(grid_t) :: this
        integer :: dims(3)
    end function

    module subroutine domain_decomposition(this, nx, ny, nimages, ratio)
        implicit none
        class(grid_t),  intent(inout) :: this
        integer,        intent(in)    :: nx, ny, nimages
        real,           intent(in), optional :: ratio
    end subroutine

    module subroutine get_grid_dimensions(this, nx, ny, nz, nx_extra, ny_extra)
        implicit none
        class(grid_t),   intent(inout) :: this
        integer,         intent(in)    :: nx, ny, nz
        integer,         intent(in), optional :: nx_extra, ny_extra
    end subroutine
end interface
end module
