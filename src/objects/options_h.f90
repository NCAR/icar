module options_interface

    use icar_constants,             only : kMAX_STRING_LENGTH, kMAX_STORAGE_VARS
    use options_types,              only : parameter_options_type, physics_type, mp_options_type, lt_options_type, &
                                           block_options_type, adv_options_type, lsm_options_type, bias_options_type

    implicit none

    private
    public :: options_t

    type :: options_t
        character(len=kMAX_STRING_LENGTH) :: comment

        integer :: vars_to_advect(   kMAX_STORAGE_VARS ) = 0
        integer :: vars_to_allocate( kMAX_STORAGE_VARS ) = 0
        integer :: vars_for_restart( kMAX_STORAGE_VARS ) = 0

        type(parameter_options_type),   allocatable :: parameters[:]

        ! defines which physics package to be used.
        type(physics_type),             allocatable :: physics[:]

        ! physics parameterization options
        type(mp_options_type),          allocatable :: mp_options[:]

        type(lt_options_type),          allocatable :: lt_options[:]

        type(block_options_type),       allocatable :: block_options[:]

        type(adv_options_type),         allocatable :: adv_options[:]

        type(lsm_options_type),         allocatable :: lsm_options[:]

        type(bias_options_type),        allocatable :: bias_options[:]

    contains

        procedure, public  :: init
        procedure, public  :: alloc_vars
        procedure, public  :: restart_vars
        procedure, public  :: advect_vars
    end type


interface

    module subroutine init(this)
        implicit none
        class(options_t),   intent(inout)  :: this
    end subroutine

    module subroutine alloc_vars(this, input_vars, var_idx, error)
        implicit none
        class(options_t), intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error
    end subroutine

    module subroutine restart_vars(this, input_vars, var_idx, error)
        implicit none
        class(options_t), intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error
    end subroutine

    module subroutine advect_vars(this, input_vars, var_idx, error)
        implicit none
        class(options_t), intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error
    end subroutine

end interface

end module
