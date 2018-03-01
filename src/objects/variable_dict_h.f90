!>------------------------------------------------
!! Defines the interface for a variable dictionary
!!
!! Dictionary can be accessed by string keys
!! Dictionary stores variable_t types
!!
!!------------------------------------------------

module variable_dict_interface
    ! variable type to store... this could be made unlimited but that complicates use
    use variable_interface,     only : variable_t
    use icar_constants

    !>------------------------------------------------
    !! Defines the object that is actually stored in the dictionary array
    !!
    !! This object stores the key-value pair and can also store an associated variable
    !!
    !!------------------------------------------------
    type var_dict_element
        character(len=kMAX_NAME_LENGTH) :: name         ! serves as the dictionary key
        type(variable_t)                :: var          ! store the primary dictionary variable for the key name
        type(variable_t)                :: domain_var   ! store the associated dictionary variable
    end type

    !>------------------------------------------------
    !! Defines the dictionary
    !!
    !! Primarily and array of dictionary elements
    !!
    !! Methods are defined to add to the dictionary and retrieve from it
    !! Note that there is no way to remove objects from the dictionary at present
    !! Also that the routines to work with an associated variable (domain_var) are not implemented yet
    !!
    !!------------------------------------------------
    type var_dict_t
        ! this is the primary data structure that holds the dictionary data
        type(var_dict_element), allocatable :: var_list(:)

        ! used to iterate through variables in dictionary
        integer :: current_variable

        ! keep track of how many variables we have stored
        integer :: n_vars = 0
        ! keep track of how many variables we are able to store (== size(var_list))
        integer :: max_vars = 0

        ! keep track of wheather or not this dictionary has been initialized yet or not
        logical :: initialized = .False.

    contains
        procedure :: reset_iterator     ! reset internal counters to make it possible to iterate over elements
        procedure :: has_more_elements  ! test if there are more elements to iterate over still
        procedure :: next               ! continue iterating through array elements
        procedure :: get_var            ! get a variable for a given key
        procedure :: add_var            ! store a variable with a given key
        procedure :: get_domain_var     ! get the associated domain variable for a key
        procedure :: set_domain_var     ! set the associated domain variable for a key

        procedure :: init               ! initialize the dictionary (e.g. allocate the var_list)
    end type

interface

    !>------------------------------------------------
    !! Additional documentation provided in the associated object implementation
    !!
    !!------------------------------------------------

    module subroutine init(this)
        implicit none
        class(var_dict_t),   intent(inout)  :: this
    end subroutine


    !>-------------------------
    !! Module subroutines for managing the dictionary as an iterator
    !!
    !!-------------------------
    module subroutine reset_iterator(this)
        implicit none
        class(var_dict_t),   intent(inout)  :: this
    end subroutine

    module function has_more_elements(this) result(boolean)
        implicit none
        class(var_dict_t),   intent(in) :: this
        logical :: boolean
    end function

    module function next(this, name, err) result(var_data)
        implicit none
        class(var_dict_t),   intent(inout)  :: this
        character(len=*),    intent(out),   optional :: name
        integer,             intent(out),   optional :: err
        type(variable_t)                    :: var_data
    end function


    !>-------------------------
    !! Primary subroutines to add and retrieve elements
    !!
    !!-------------------------
    module function get_var(this, varname, err) result(var_data)
        implicit none
        class(var_dict_t),   intent(in) :: this
        character(len=*),    intent(in) :: varname
        integer,             intent(out),   optional :: err
        type(variable_t)                :: var_data
    end function

    module subroutine add_var(this, varname, var_data, domain_var, save_state, err)
        implicit none
        class(var_dict_t),   intent(inout)  :: this
        character(len=*),    intent(in)     :: varname
        type(variable_t),    intent(in)     :: var_data
        type(variable_t),    intent(in), optional :: domain_var
        logical,             intent(in), optional :: save_state
        integer,             intent(out),optional :: err
    end subroutine

    module function get_domain_var(this, varname) result(var_data)
        implicit none
        class(var_dict_t),   intent(in) :: this
        character(len=*),    intent(in) :: varname
        type(variable_t),    pointer    :: var_data
    end function

    module subroutine set_domain_var(this, varname, var_data)
        implicit none
        class(var_dict_t),   intent(inout)  :: this
        character(len=*),    intent(in)     :: varname
        type(variable_t),    intent(in)     :: var_data
    end subroutine

end interface

end module
