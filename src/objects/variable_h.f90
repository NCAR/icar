module variable_interface
    use icar_constants,      only : kMAX_DIM_LENGTH
    use grid_interface,      only : grid_t
    use meta_data_interface, only : meta_data_t

    implicit none

    ! defines a variable type that can store data and attributes
    ! have to think about how to handle multiple variable types (int, 2d, etc)
    ! could add multiple "local" variables or create multiple variable types...
    type, extends(meta_data_t) :: variable_t
        real, pointer :: data_3d(:,:,:) => null()
        real, pointer :: data_2d(:,:)   => null()

        logical                         :: unlimited_dim = .False.
        logical                         :: three_d = .False.
        logical                         :: two_d = .False.

        integer :: n_dimensions
        integer,                        allocatable :: dim_len(:)
        character(len=kMAX_DIM_LENGTH), allocatable :: dimensions(:)

        ! note these are used for netcdf output
        integer, allocatable    :: dim_ids(:)
        integer                 :: var_id = -1

    contains
        procedure, public  :: broadcast
        procedure, public  :: init_grid
        procedure, public  :: init_dims
        generic,   public  :: initialize => init_grid
        generic,   public  :: initialize => init_dims

    ! inherited from meta_data
    !     procedure, public : add_attribute

    end type

    interface

        module subroutine broadcast(this, source, start, end)
            implicit none
            class(variable_t),  intent(inout) :: this
            integer,            intent(in)    :: source, start, end
        end subroutine


        module subroutine init_grid(this, grid)
            implicit none
            class(variable_t),  intent(inout) :: this
            type(grid_t),       intent(in)    :: grid
        end subroutine

        module subroutine init_dims(this, dims)
            implicit none
            class(variable_t),  intent(inout) :: this
            integer,            intent(in)    :: dims(:)
        end subroutine

    end interface

end module
