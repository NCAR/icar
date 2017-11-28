module variable_interface
    use icar_constants
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

        ! Note this dummy coarray is only here because types that inherit from variable_t have coarrays
        ! so the parent type is required to...
        ! integer, allocatable :: dummy[:]

    ! contains
    !     procedure, public : add_attribute
    !     procedure, public : set_attribute
    !     procedure, public : get_attribute
    !
    !     procedure, public : set_dimension
    !     procedure, public : get_dimension

    end type variable_t

    ! interface
    !     module subroutine add_attribute(this, attr_name, attr_value)
    !         implicit none
    !         type(variable_t), intent(inout) :: this
    !         character(len=kMAX_ATTR_LENGTH), intent(in) :: attr_name
    !         character(len=kMAX_ATTR_LENGTH), intent(in) :: attr_value
    !     end subroutine
    !
    !
    !     module subroutine set_dimension(this, dim, name, length)
    !         implicit none
    !         type(variable_t), intent(inout) :: this
    !         character(len=kMAX_ATTR_LENGTH), intent(in) :: attr_name
    !         character(len=kMAX_ATTR_LENGTH), intent(in) :: attr_value
    !     end subroutine
    !
    ! end interface
    !  complains that new_variable procedure doesn't exist
    ! interface variable_t
    !     module procedure new_variable
    ! end interface
    !
    ! interface
    !     module function new_variable(name)
    !         implicit none
    !         character(len=kMAX_NAME_LENGTH), intent(in) :: name
    !         type(variable_t) :: new_variable
    !     end function
    ! end interface
    !
    ! complains that it can't have a coarray in the output
    ! interface variable_t
    !     module function new_variable(name)
    !         implicit none
    !         character(len=kMAX_NAME_LENGTH), intent(in) :: name
    !         type(variable_t) :: new_variable
    !     end function
    ! end interface

end module
