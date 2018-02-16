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

        this%dimensions = grid%dimensions
        this%dim_len    = grid%get_dims()

        this%two_d   = grid%is2d
        this%three_d = grid%is3d

        if (grid%is2d) then
            this%n_dimensions = 2
            if (associated(this%data_2d)) deallocate(this%data_2d)
            allocate(this%data_2d(grid%ims:grid%ime,    &
                                  grid%jms:grid%jme))
        endif

        if (grid%is3d) then
            this%n_dimensions = 3
            if (associated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(grid%ims:grid%ime,    &
                                  grid%kms:grid%kme,    &
                                  grid%jms:grid%jme))
        endif

    end subroutine

    module subroutine init_dims(this, dims)
        implicit none
        class(variable_t),  intent(inout) :: this
        integer,            intent(in)    :: dims(:)

        this%dim_len    = dims

        this%two_d   = size(dims) == 2
        this%three_d = size(dims) == 3

        if (this%two_d) then
            this%n_dimensions = 2
            this%dimensions = ['x','y']
            if (associated(this%data_2d)) deallocate(this%data_2d)
            allocate(this%data_2d(dims(1), dims(2)))
        endif

        if (this%three_d) then
            this%n_dimensions = 3
            this%dimensions = ['x','y','z']
            if (associated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(dims(1), dims(2), dims(3)))
        endif

    end subroutine


end submodule
