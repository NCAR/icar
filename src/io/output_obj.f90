submodule(output_interface) output_implementation
  implicit none


contains

    module subroutine set_domain(this, domain)
        class(output_t),  intent(inout)  :: this
        type(domain_t),   intent(in)     :: domain
        integer :: i

        if (.not.this%is_initialized) call this%init()

        do i=1,domain%info%n_attrs
            call this%add_attribute(domain%info%attribute_names(i), domain%info%attribute_values(i))
        enddo

    end subroutine


    module subroutine add_to_output(this, variable)
        class(output_t),   intent(inout)  :: this
        type(variable_t),  intent(in)     :: variable

        if (.not.this%is_initialized) call this%init()

        if (this%n_variables == size(this%variables)) call this%increase_var_capacity()

        this%n_variables = this%n_variables + 1
        this%variables(this%n_variables) = variable

    end subroutine


    module subroutine save_file(this, filename)
        class(output_t), intent(inout)  :: this
        character(len=*), intent(in) :: filename
        integer :: err

        if (.not.this%is_initialized) call this%init()

        ! open file
        this%filename = filename
        err = nf90_open(filename, NF90_WRITE, this%ncfile_id)
        if (err /= NF90_NOERR) then
            call check( nf90_create(filename, NF90_CLOBBER, this%ncfile_id), "Opening:"//trim(filename))
            this%creating=.True.
        endif

        ! define variables or find variable IDs (and dimensions)
        call setup_variables(this)

        ! add global attributes such as the image number
        call add_global_attributes(this)

        if (this%creating) then
            ! End define mode. This tells netCDF we are done defining metadata.
            call check( nf90_enddef(this%ncfile_id), "end define mode" )
            this%creating=.false.
        endif

        ! store output
        call save_data(this)

        ! close file
        call check(nf90_close(this%ncfile_id), "Closing file "//trim(filename))
    end subroutine

    subroutine add_global_attributes(this)
        implicit none
        class(output_t), intent(inout)  :: this
        integer :: i

        if (this%n_attrs > 0) then
            do i=1,this%n_attrs
                call check( nf90_put_att(   this%ncfile_id,             &
                                            NF90_GLOBAL,                &
                                            trim(this%attribute_names(i)),    &
                                            trim(this%attribute_values(i))),  &
                                            "global attr:"//trim(this%attribute_names(i)))
            enddo
        endif

        call check(nf90_put_att(this%ncfile_id, NF90_GLOBAL, "image", this_image()))

    end subroutine add_global_attributes

    subroutine setup_variables(this)
        implicit none
        class(output_t), intent(inout) :: this
        integer :: i

        ! iterate through variables creating or setting up variable IDs if they exist, also dimensions
        do i=1,this%n_variables
            ! create all dimensions or find dimension IDs if they exist already
            call setup_dims_for_var(this, this%variables(i))

            call setup_variable(this, this%variables(i))
        end do

    end subroutine setup_variables


    subroutine save_data(this)
        implicit none
        class(output_t), intent(in) :: this
        integer :: i

        do i=1,this%n_variables
            associate(var => this%variables(i))
                if (var%three_d) then
                    call check( nf90_put_var(this%ncfile_id, var%var_id,  var%data_3d),   &
                                "saving:"//trim(var%name) )
                elseif (var%two_d) then
                    call check( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d),   &
                                "saving:"//trim(var%name) )
                endif
            end associate
        end do

    end subroutine save_data

    subroutine setup_dims_for_var(this, var)
        implicit none
        class(output_t),    intent(inout) :: this
        class(variable_t),  intent(inout) :: var
        integer :: i, err

        if (allocated(var%dim_ids)) deallocate(var%dim_ids)

        if (var%three_d) then
            allocate(var%dim_ids(3))
        elseif (var%two_d) then
            allocate(var%dim_ids(2))
        endif

        do i = 1, size(var%dim_ids)

            ! Try to find the dimension ID if it exists already.
            err = nf90_inq_dimid(this%ncfile_id, trim(var%dimensions(i)), var%dim_ids(i))

            ! probably the dimension doesn't exist in the file, so we will create it.
            if (err/=NF90_NOERR) then

                ! 4d = time dimension, only write if it has one
                if (i == 4) then
                    if (var%unlimited_dim) then
                        call check( nf90_def_dim(this%ncfile_id, trim(var%dimensions(i)), NF90_UNLIMITED, &
                                    var%dim_ids(i) ), "def_dim"//var%dimensions(i) )
                    endif
                else
                    call check( nf90_def_dim(this%ncfile_id, var%dimensions(i), var%dim_len(i),       &
                                var%dim_ids(i) ), "def_dim"//var%dimensions(i) )
                endif
            endif
        end do

    end subroutine setup_dims_for_var

    subroutine setup_variable(this, var)
        implicit none
        class(output_t),   intent(inout) :: this
        class(variable_t), intent(inout) :: var
        integer :: i, err

        err = nf90_inq_varid(this%ncfile_id, var%name, var%var_id)

        ! if the variable was not found in the netcdf file then we will define it.
        if (err /= NF90_NOERR) then
            if (var%unlimited_dim) then
                call check( nf90_def_var(this%ncfile_id, var%name, NF90_REAL, var%dim_ids, var%var_id), &
                            "Defining variable:"//trim(var%name) )
            else
                call check( nf90_def_var(this%ncfile_id, var%name, NF90_REAL, var%dim_ids(1:3), var%var_id), &
                            "Defining variable:"//trim(var%name) )
            endif

            ! setup attributes
            do i=1,size(var%attribute_names)
                call check( nf90_put_att(this%ncfile_id,                &
                                         var%var_id,                    &
                                         trim(var%attribute_names(i)),        &
                                         trim(var%attribute_values(i))),      &
                            "saving attribute"//trim(var%attribute_names(i)))
            enddo
        endif

    end subroutine setup_variable

    module subroutine init(this)
        implicit none
        class(output_t),   intent(inout)  :: this

        allocate(this%variables(kINITIAL_VAR_SIZE))
        this%n_variables = 0
        this%n_dims      = 0
        this%is_initialized = .True.


    end subroutine

    module subroutine increase_var_capacity(this)
        implicit none
        class(output_t),   intent(inout)  :: this
        class(variable_t), allocatable :: new_variables(:)

        ! assert allocated(this%variables)
        allocate(new_variables, source=this%variables)
        ! new_variables = this%variables

        deallocate(this%variables)

        allocate(this%variables(size(new_variables)*2))
        this%variables(:size(new_variables)) = new_variables

        deallocate(new_variables)

    end subroutine


    !>------------------------------------------------------------
    !! Simple error handling for common netcdf file errors
    !!
    !! If status does not equal nf90_noerr, then print an error message and STOP
    !! the entire program.
    !!
    !! @param   status  integer return code from nc_* routines
    !! @param   extra   OPTIONAL string with extra context to print in case of an error
    !!
    !!------------------------------------------------------------
    subroutine check(status,extra)
        implicit none
        integer, intent ( in) :: status
        character(len=*), optional, intent(in) :: extra

        ! check for errors
        if(status /= nf90_noerr) then
            ! print a useful message
            print *, trim(nf90_strerror(status))
            if(present(extra)) then
                ! print any optionally provided context
                write(*,*) trim(extra)
            endif
            ! STOP the program execution
            stop "Stopped"
        end if
    end subroutine check


end submodule
