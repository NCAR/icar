module boundary_interface
    use icar_constants
    use options_interface,        only : options_t
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t
    use meta_data_interface,      only : meta_data_t
    use time_object,              only : Time_type
    use time_delta_object,        only : time_delta_t
    use data_structures,          only : interpolable_type
    use icar_constants

    implicit none

    private
    public :: boundary_t

    ! ------------------------------------------------
    ! boundary conditions type, must be linearizable so we can remove low res linear wind field
    ! ------------------------------------------------
    type :: boundary_t
        type(meta_data_t)    :: info

        ! list of input files
        character (len=kMAX_FILE_LENGTH), allocatable :: file_list(:)
        !   manage file pointer and position in file for boundary conditions
        integer :: curfile
        integer :: curstep

        type(Time_type)                   :: current_time   ! the date/time of the forcing data in memory
        type(time_delta_t)                :: forcing_dt     ! the time step in between two forcing steps
        character(len=kMAX_STRING_LENGTH) :: time_var       ! the name of the input time variable [optional]

        type(var_dict_t)                  :: variables      ! a dictionary with all forcing data

        ! boundary data coordinate system
        real, dimension(:,:),   allocatable :: lat, lon
        real, dimension(:,:,:), allocatable :: z

        type(interpolable_type) :: geo
        type(interpolable_type) :: geo_u
        type(interpolable_type) :: geo_v
        type(interpolable_type) :: original_geo

    contains

        procedure :: init
        procedure :: init_external
        procedure :: update_forcing

        procedure :: distribute_update
        procedure :: distribute_initial_conditions

        ! procedure :: find_start_time
        procedure :: init_local
        procedure :: init_local2

    end type

    interface

    ! Set default component values
    module subroutine init(this, options)
        implicit none
        class(boundary_t), intent(inout) :: this
        type(options_t),   intent(inout) :: options
    end subroutine

    module subroutine init_external(this, options)
        implicit none
        class(boundary_t), intent(inout) :: this
        type(options_t),   intent(inout) :: options
    end subroutine

    module subroutine init_local2(this, options, file_list, var_list, dim_list, start_time, &
                                 lat_ext, lon_ext, zvar_ext, time_ext) !, p_var, ps_var)
        implicit none
        class(boundary_t),               intent(inout)          :: this
        type(options_t),                 intent(inout)          :: options
        character(len=kMAX_NAME_LENGTH), intent(in)             :: file_list(:)
        character(len=kMAX_NAME_LENGTH), intent(in)             :: var_list (:)
        integer,                         intent(in)             :: dim_list (:)
        type(Time_type),                 intent(in), optional             :: start_time
        character(len=kMAX_NAME_LENGTH), intent(in)             :: lat_ext
        character(len=kMAX_NAME_LENGTH), intent(in)             :: lon_ext
        character(len=kMAX_NAME_LENGTH), intent(in), optional   :: zvar_ext
        character(len=kMAX_NAME_LENGTH), intent(in), optional   :: time_ext
        ! character(len=kMAX_NAME_LENGTH), intent(in)     :: p_var
        ! character(len=kMAX_NAME_LENGTH), intent(in)     :: ps_var
    end subroutine

    module subroutine init_local(this, options, file_list, var_list, dim_list, start_time, &
                                 lat_var, lon_var, z_var, time_var, p_var, ps_var)
        implicit none
        class(boundary_t),               intent(inout)  :: this
        type(options_t),                 intent(inout)  :: options
        character(len=kMAX_NAME_LENGTH), intent(in)     :: file_list(:)
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list (:)
        integer,                         intent(in)     :: dim_list (:)
        type(Time_type),                 intent(in)     :: start_time
        character(len=kMAX_NAME_LENGTH), intent(in)     :: lat_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: lon_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: z_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: time_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: p_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: ps_var
    end subroutine

    module subroutine update_forcing(this, options)
        implicit none
        class(boundary_t), intent(inout) :: this
        type(options_t),   intent(inout) :: options
    end subroutine

    module subroutine distribute_update(this)
        implicit none
        class(boundary_t), intent(inout) :: this
    end subroutine

    module subroutine distribute_initial_conditions(this)
        implicit none
        class(boundary_t), intent(inout) :: this
    end subroutine

  end interface

end module
