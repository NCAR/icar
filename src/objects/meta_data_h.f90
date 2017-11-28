module meta_data_interface
    use icar_constants

    implicit none

    private
    public :: meta_data_t

    type meta_data_t
    !   private
        character(len=kMAX_NAME_LENGTH) :: name
        integer :: n_attrs = 0
        character(len=kMAX_ATTR_LENGTH), allocatable :: attribute_names(:)
        character(len=kMAX_ATTR_LENGTH), allocatable :: attribute_values(:)

    contains

        procedure, public  :: add_attribute
        ! procedure, public  :: get_attribute
    end type

    interface

        module subroutine add_attribute(this, input_name, input_value)
            implicit none
            class(meta_data_t), intent(inout) :: this
            character(len=*),   intent(in)    :: input_name
            character(len=*),   intent(in)    :: input_value
        end subroutine

        ! module subroutine get_attribute(this, name, value)
        !     implicit none
        !     class(output_t),   intent(inout)  :: this
        !     character(len=kMAX_ATTR_LENGTH), intent(in) :: name
        !     character(len=kMAX_ATTR_LENGTH), intent(value) :: value
        ! end subroutine

    end interface
end module
