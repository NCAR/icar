submodule(variable_interface) variable_implementation
  implicit none


contains

    !> -------------------------------
    !! Initialize a variable object from a given grid
    !!
    !! Allocates 2d/3d data structure as appropriate
    !!
    !! -------------------------------
    module subroutine init_grid(this, grid)
        implicit none
        class(variable_t),  intent(inout) :: this
        type(grid_t),       intent(in)    :: grid
        integer :: err

        this%dimensions = grid%dimensions
        this%dim_len    = grid%get_dims()

        this%two_d   = grid%is2d
        this%three_d = grid%is3d

        if (grid%is2d) then
            this%n_dimensions = 2
            if (associated(this%data_2d)) deallocate(this%data_2d)
            allocate(this%data_2d(grid%ims:grid%ime,    &
                                  grid%jms:grid%jme), stat=err)
            if (err /= 0) stop "variable:grid:2d: Allocation request denied"
        endif

        if (grid%is3d) then
            this%n_dimensions = 3
            if (associated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(grid%ims:grid%ime,    &
                                  grid%kms:grid%kme,    &
                                  grid%jms:grid%jme), stat=err)
            if (err /= 0) stop "variable:grid:3d: Allocation request denied"
        endif

    end subroutine

    !> -------------------------------
    !! Initialize a variable object from a given array of dimension sizes
    !!
    !! Allocates 2d/3d data structure as appropriate
    !!
    !! -------------------------------
    module subroutine init_dims(this, dims)
        implicit none
        class(variable_t),  intent(inout) :: this
        integer,            intent(in)    :: dims(:)

        integer :: err

        this%dim_len    = dims

        this%two_d   = size(dims) == 2
        this%three_d = size(dims) == 3

        if (this%two_d) then
            this%n_dimensions = 2
            this%dimensions = ['x','y']
            if (associated(this%data_2d)) deallocate(this%data_2d)
            allocate(this%data_2d(dims(1), dims(2)), stat=err)
            if (err /= 0) stop "variable:dims:2d: Allocation request denied"
        endif

        if (this%three_d) then
            this%n_dimensions = 3
            this%dimensions = ['x','y','z']
            if (associated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(dims(1), dims(2), dims(3)), stat=err)
            if (err /= 0) stop "variable:dims:3d: Allocation request denied"
        endif

    end subroutine


    module subroutine broadcast(this, source, start, end)
        implicit none
        class(variable_t),  intent(inout) :: this
        integer,            intent(in)    :: source, start, end



        ! all the components of the variable that may need to be broadcast
        ! real, pointer :: data_3d(:,:,:) => null()
        ! real, pointer :: data_2d(:,:)   => null()
        !
        ! logical                         :: unlimited_dim = .False.
        ! logical                         :: three_d = .False.
        ! logical                         :: two_d = .False.
        !
        ! integer :: n_dimensions
        ! integer,                        allocatable :: dim_len(:)
        ! character(len=kMAX_DIM_LENGTH), allocatable :: dimensions(:)
        !
        ! ! note these are used for netcdf output
        ! integer, allocatable    :: dim_ids(:)
        ! integer                 :: var_id = -1

        ! also keep in mind, any inherited meta_data elements must be broadcast too
        ! type attribute_t
        !     character(len=kMAX_NAME_LENGTH) :: name
        !     character(len=kMAX_ATTR_LENGTH) :: value
        ! end type
        !
        ! character(len=kMAX_NAME_LENGTH) :: name
        ! integer :: n_attrs = 0
        !
        ! type(attribute_t), allocatable :: attributes(:)

    end subroutine


end submodule
