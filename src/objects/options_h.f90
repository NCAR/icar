module options_interface

    use icar_constants,             only : kMAX_STRING_LENGTH
    use options_types,              only : parameter_options_type, physics_type, mp_options_type, lt_options_type, &
                                           block_options_type, adv_options_type, lsm_options_type, bias_options_type

    implicit none

    private
    public :: options_t

    type :: options_t
        character(len=kMAX_STRING_LENGTH) :: comment

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
    end type


interface

    module subroutine init(this)
        implicit none
        class(options_t),   intent(inout)  :: this

    end subroutine

end interface

end module
