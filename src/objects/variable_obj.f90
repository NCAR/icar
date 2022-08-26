submodule(variable_interface) variable_implementation
    use icar_constants, only : kREAL, kDOUBLE
    use co_util,        only : broadcast
    implicit none


contains

    !> -------------------------------
    !! Initialize a variable object from a given grid
    !!
    !! Allocates 2d/3d data structure as appropriate
    !!
    !! -------------------------------
    module subroutine init_grid(this, grid, forcing_var, force_boundaries, dtype)
        implicit none
        class(variable_t),  intent(inout) :: this
        type(grid_t),       intent(in)    :: grid
        character(len=*),   intent(in), optional :: forcing_var
        logical,            intent(in), optional :: force_boundaries
        integer,            intent(in), optional :: dtype

        integer :: err

        this%dtype = kREAL
        if (present(dtype)) this%dtype = dtype

        this%dimensions = grid%dimensions
        this%dim_len    = grid%get_dims()

        this%two_d   = grid%is2d
        this%three_d = grid%is3d

        this%forcing_var = ""
        if (present(forcing_var)) this%forcing_var = forcing_var

        this%force_boundaries = .True.
        if (present(force_boundaries)) this%force_boundaries = force_boundaries

        if (grid%is2d) then
            this%n_dimensions = 2
            if (associated(this%data_2d)) deallocate(this%data_2d)
            if (this%dtype == kREAL) then
                allocate(this%data_2d(grid%ims:grid%ime,    &
                                      grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:2d: Allocation request failed"

                this%data_2d = 0
            elseif (this%dtype == kDOUBLE) then
                allocate(this%data_2dd(grid%ims:grid%ime,    &
                                       grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:2d: Allocation request failed"

                this%data_2dd = 0
            endif

            if (trim(this%forcing_var) /= "") then
                allocate(this%dqdt_2d(grid%ims:grid%ime,    &
                                      grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:dqdt_2d: Allocation request failed"

                this%dqdt_2d = 0
            endif

        endif

        if (grid%is3d) then
            this%n_dimensions = 3
            if (associated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(grid%ims:grid%ime,    &
                                  grid%kms:grid%kme,    &
                                  grid%jms:grid%jme), stat=err)
            if (err /= 0) stop "variable:grid:3d: Allocation request failed"

            this%data_3d = 0

            if (trim(this%forcing_var) /= "") then
                allocate(this%dqdt_3d(grid%ims:grid%ime,    &
                                      grid%kms:grid%kme,    &
                                      grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:dqdt_3d: Allocation request failed"

                this%dqdt_3d = 0
            endif

        endif

    end subroutine

    !> -------------------------------
    !! Initialize a variable object from a given array of dimension sizes
    !!
    !! Allocates 2d/3d data structure as appropriate
    !!
    !! -------------------------------
    module subroutine init_dims(this, dims, forcing_var, force_boundaries)
        implicit none
        class(variable_t),  intent(inout) :: this
        integer,            intent(in)    :: dims(:)
        character(len=*),   intent(in), optional :: forcing_var
        logical,            intent(in), optional :: force_boundaries

        integer :: err

        this%dim_len    = dims

        this%two_d   = size(dims) == 2
        this%three_d = size(dims) == 3

        this%forcing_var = ""
        if (present(forcing_var)) this%forcing_var = forcing_var

        this%force_boundaries = .True.
        if (present(force_boundaries)) this%force_boundaries = force_boundaries

        if (this%two_d) then
            this%n_dimensions = 2
            this%dimensions = ['x','y']
            if (associated(this%data_2d)) deallocate(this%data_2d)
            allocate(this%data_2d(dims(1), dims(2)), stat=err)
            if (err /= 0) stop "variable:dims:2d: Allocation request denied"
            this%data_2d = 0

            if (trim(this%forcing_var) /= "") then
                if (associated(this%dqdt_2d)) deallocate(this%dqdt_2d)
                allocate(this%dqdt_2d(dims(1), dims(2)), stat=err)
                if (err /= 0) stop "variable:dims:dqdt_2d: Allocation request denied"

                this%dqdt_2d = 0
            endif
        endif

        if (this%three_d) then
            this%n_dimensions = 3
            this%dimensions = ['x','y','z']
            if (associated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(dims(1), dims(2), dims(3)), stat=err)
            if (err /= 0) stop "variable:dims:3d: Allocation request denied"

            this%data_3d = 0

            if (trim(this%forcing_var) /= "") then
                if (associated(this%dqdt_3d)) deallocate(this%dqdt_3d)
                allocate(this%dqdt_3d(dims(1), dims(2), dims(3)), stat=err)
                if (err /= 0) stop "variable:dims:dqdt_3d: Allocation request denied"

                this%dqdt_3d = 0
            endif
        endif

    end subroutine


    module subroutine bcast_var(this, source, start_img, end_img)
        implicit none
        class(variable_t),  intent(inout) :: this
        integer,            intent(in)    :: source
        integer,            intent(in),   optional :: start_img, end_img

        integer :: first, last, attr_array_size
        character(len=kMAX_STRING_LENGTH) :: error_message
        integer :: i

        if (present(start_img)) first = start_img
        if (present(end_img))   last  = end_img


        ! all the components of the variable that may need to be broadcast
        ! First we simply broadcast all the "easy" scalar values
        call broadcast(this%unlimited_dim, source, first, last, create_co_array=.True.)
        call broadcast(this%n_dimensions,  source, first, last, create_co_array=.True.)
        call broadcast(this%three_d,       source, first, last, create_co_array=.True.)
        call broadcast(this%two_d,         source, first, last, create_co_array=.True.)
        call broadcast(this%forcing_var,   source, first, last, create_co_array=.True.)

        ! these attributes are inherited from meta_data parent class and must be broadcast too
        call broadcast(this%name,          source, first, last, create_co_array=.True.)
        call broadcast(this%n_attrs,       source, first, last, create_co_array=.True.)

        ! we have to figure out how big the attribute array is as n_attrs is the number stored in it, not the memory allocated for it
        if (this_image()==source) then
            if (allocated(this%attributes)) then
                attr_array_size = size(this%attributes)
            else
                attr_array_size = 1
                allocate( this%attributes( attr_array_size ) )
            endif
        endif

        call broadcast(attr_array_size,    source, first, last, create_co_array=.True.)

        if (.not.allocated(this%attributes)) allocate(this%attributes(attr_array_size))
        do i=1, this%n_attrs
            call broadcast(this%attributes(i)%name,  source, first, last, create_co_array=.True.)
            call broadcast(this%attributes(i)%value, source, first, last, create_co_array=.True.)
        enddo


        ! Handle anything that potentially has 2 or 3 dimensions separately
        ! First handle the if 3D case
        if (this%three_d) then
            if (size(this%dim_len) /= 3) deallocate(this%dim_len)
            if (.not.allocated(this%dim_len))       allocate(this%dim_len(3))
            call broadcast(this%dim_len, source, first, last, create_co_array=.True.)

            if (.not.allocated(this%dimensions))    allocate(this%dimensions(3))
            call broadcast(this%dimensions, source, first, last, create_co_array=.True.)

            if (.not.associated(this%data_3d))      allocate(this%data_3d(this%dim_len(1), this%dim_len(2), this%dim_len(3)))
            if (any(shape(this%data_3d) /= this%dim_len)) then
                write(error_message, '(A,A,A,I6,A,3I5,A,3I5)') &
                        "ERROR: variable ", trim(this%name), " has the wrong shape on image: ", &
                        this_image(), " has:",shape(this%data_3d)," should have:",this%dim_len
                write(*,*) trim(error_message)
                stop
            endif
            call broadcast(this%data_3d, source, first, last, create_co_array=.True.)

        ! Then handle the if 2D case
        else if (this%two_d) then
            if (allocated(this%dim_len) .and. (size(this%dim_len) /= 2)) deallocate(this%dim_len)
            if (.not.allocated(this%dim_len))       allocate(this%dim_len(2))
            call broadcast(this%dim_len, source, first, last, create_co_array=.True.)

            if (.not.allocated(this%dimensions))    allocate(this%dimensions(2))
            call broadcast(this%dimensions, source, first, last, create_co_array=.True.)

            if (.not.associated(this%data_2d))      allocate(this%data_2d(this%dim_len(1), this%dim_len(2)))
            if (any(shape(this%data_2d) /= this%dim_len)) then
                write(error_message, '(A,A,A,I6,A,2I5,A,2I5)') &
                        "ERROR: variable ", trim(this%name), " has the wrong shape on image: ", &
                        this_image(), " has:",shape(this%data_2d)," should have:",this%dim_len
                write(*,*) trim(error_message)
                stop
            endif
            call broadcast(this%data_2d, source, first, last, create_co_array=.True.)
        endif

    end subroutine


end submodule
