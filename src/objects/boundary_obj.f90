!>------------------------------------------------------------
!!  Handles reading boundary conditions from the forcing file(s)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
submodule(boundary_interface) boundary_implementation

    use array_utilities,        only : interpolate_in_z
    use io_routines,            only : io_getdims, io_read, io_maxDims, io_variable_is_present
    use time_io,                only : read_times, find_timestep_in_file
    use co_util,                only : broadcast
    use string,                 only : str
    use mod_atm_utilities,      only : rh_to_mr, compute_3d_p, compute_3d_z, exner_function
    use geo,                    only : standardize_coordinates
    use vertical_interpolation, only : vLUT, vinterp

    implicit none
contains

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file into image 1
    !!
    !! Distributes initial conditions to all other images
    !!
    !!------------------------------------------------------------
    module subroutine init(this, options)
        class(boundary_t), intent(inout) :: this
        type(options_t),   intent(inout) :: options

        character(len=kMAX_NAME_LENGTH), allocatable :: vars_to_read(:)
        integer,                         allocatable :: var_dimensions(:)

        this%file_list = options%parameters%boundary_files

        ! the parameters option type can't contain allocatable arrays because it is a coarray
        ! so we need to allocate the vars_to_read and var_dimensions outside of the options type
        call setup_variable_lists(options%parameters%vars_to_read, options%parameters%dim_list, vars_to_read, var_dimensions)

        ! Read through forcing variable names stored in "options"
        ! needs to read each one to find the grid information for it
        ! then create grid and initialize a variable...
        ! also need to explicitly save lat and lon data
        ! if (this_image() == 1) then
            call this%init_local(options,                           &
                                 this%file_list,                    &
                                 vars_to_read, var_dimensions,      &
                                 options%parameters%start_time,     &
                                 options%parameters%latvar,         &
                                 options%parameters%lonvar,         &
                                 options%parameters%zvar,           &
                                 options%parameters%time_var,       &
                                 options%parameters%pvar,           &
                                 options%parameters%psvar           &
                                 )

        ! endif
        ! call this%distribute_initial_conditions()


        call setup_boundary_geo(this, options%parameters%longitude_system)

    end subroutine

    !>-----------------------------------------------------------
    !! initialize structure for additional external files (i.e. SWE, snowheight etc)
    !!
    !>-----------------------------------------------------------
    module subroutine init_external(this, options)
        class(boundary_t), intent(inout) :: this  ! the additional starting conditions
        type(options_t),   intent(inout) :: options

        character(len=kMAX_NAME_LENGTH), allocatable :: vars_to_read(:)
        integer,                         allocatable :: var_dimensions(:)
        character(len=kMAX_NAME_LENGTH)  :: lat_ext
        character(len=kMAX_NAME_LENGTH)  :: lon_ext
        character(len=kMAX_NAME_LENGTH)  :: zvar_ext
        ! character(len=kMAX_NAME_LENGTH), allocatable :: ext_list(:)
        ! integer, allocatable :: dim_list(:)

        allocate(this%file_list(1))
        ! allocate(ext_list(1))
        ! allocate(dim_list(1))

        this%file_list(1) = options%parameters%external_files

        call setup_variable_lists(options%parameters%ext_var_list, options%parameters%ext_dim_list, vars_to_read, var_dimensions)  ! this simply copies o%p%ext_vars to vars_to_read and same for dims (+some checks)


        if ( trim(options%parameters%z_ext) /= "") zvar_ext = options%parameters%z_ext

        call this%init_local2(options,                          &
                             this%file_list,                    &
                             vars_to_read, var_dimensions,      &
                             start_time = options%parameters%start_time,     &
                             lat_ext = options%parameters%lat_ext,        &  !lat_ext,
                             lon_ext = options%parameters%lon_ext,        &  !lon_ext      &
                             !zvar_ext = options%parameters%zvar_ext,      &  !,- include if statement for this??
                             time_ext = options%parameters%time_ext &
                             )

        call setup_boundary_geo(this, options%parameters%longitude_system)


    end subroutine


    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the external file
    !!
    !!   For now this is only available for static data (no time component)
    !!
    !!------------------------------------------------------------
    module subroutine init_local2(this, options, file_list, var_list, dim_list, start_time, &
                                 lat_ext, lon_ext, zvar_ext, time_ext)!, swe_ext, tsnow_ext)
        class(boundary_t),               intent(inout)          :: this
        type(options_t),                 intent(inout)          :: options
        character(len=kMAX_NAME_LENGTH), intent(in)             :: file_list(:)
        character(len=kMAX_NAME_LENGTH), intent(in)             :: var_list (:)
        integer,                         intent(in)             :: dim_list (:)
        type(Time_type),                 intent(in), optional   :: start_time
        character(len=kMAX_NAME_LENGTH), intent(in)             :: lat_ext
        character(len=kMAX_NAME_LENGTH), intent(in)             :: lon_ext
        character(len=kMAX_NAME_LENGTH), intent(in), optional   :: zvar_ext
        character(len=kMAX_NAME_LENGTH), intent(in), optional   :: time_ext

        real, allocatable :: temp_z(:,:,:)

        integer :: i, nx, ny, nz

        if (present(time_ext) .and. (trim(time_ext)/="") ) then
            ! figure out while file and timestep contains the requested start_time
            ! if (this_image()==1) print * ," looking for time_ext"
            call set_curfile_curstep(this, start_time, file_list, time_ext)
            call read_bc_time(this%current_time, file_list(this%curfile), time_ext, this%curstep)
        else
            this%curstep=1 !?
            this%curfile=1 !?
        endif


        !  read in latitude and longitude coordinate data
        call io_read(file_list(this%curfile), lat_ext, this%lat, this%curstep)
        call io_read(file_list(this%curfile), lon_ext, this%lon, this%curstep)


        ! if ( (size(this%lat,2)==1) .and. (size(this%lon,2)==1) ) then  !(size(this%lat,1)/=size(this%lon,1)) (shape(this%lat) /= shape(this%lon)) .and.
        !     this%lat = reshape(this%lat, (/size(this%lat,1), size(this%lon,1)/), pad=this%lat)
        !     this%lon = reshape(this%lon, (/size(this%lat,1), size(this%lon,1)/), pad=this%lon)
        ! endif

        ! read in the height coordinate of the input data
        ! if (present(zvar_ext) .and. (trim(zvar_ext)/="") ) then
        !     if (this_image()==1)  print *, " zvar_ext ", trim(zvar_ext)
        !     if (.not. options%parameters%compute_z ) then ! .and. present(zvar_ext)
        !         ! call io_read(file_list(this%curfile), z_var,   temp_z,   this%curstep)
        !         call io_read(file_list(1), zvar_ext,   temp_z,   1)
        !         nx = size(temp_z,1)
        !         ny = size(temp_z,2)
        !         nz = size(temp_z,3)
        !         if (allocated(this%z)) deallocate(this%z)
        !         allocate(this%z(nx,nz,ny))
        !         this%z = reshape(temp_z, shape=[nx,nz,ny], order=[1,3,2])
        !     else
        !         call io_read(file_list(this%curfile), file_list(1),   temp_z,   this%curstep)
        !         nx = size(temp_z,1)
        !         ny = size(temp_z,2)
        !         ! nz = size(temp_z,3)
        !         if (allocated(this%z)) deallocate(this%z)
        !         allocate(this%z(nx,nz,ny))
        !     endif
        ! else
            ! print *, "  ext var = 2D"
            if ( (size(this%lat,2)==1) .and. (size(this%lon,2)==1) ) then
                if (this_image()==1) print*, "  external conditions provided on regular 1D grid"
                nx = size(this%lon,1)
                ny = size(this%lat,1)
                nz = 0
            else
                nx = size(this%lat,1)
                ny = size(this%lat,2)
                nz = 0
            endif

        ! endif

        ! ! call assert(size(var_list) == size(dim_list), "list of variable dimensions must match list of variables")

        do i=1, size(var_list)

            call add_var_to_dict(this%variables, file_list(this%curfile), var_list(i), dim_list(i), this%curstep, [nx, nz, ny])
            ! if (this_image()==1)  print*, i," var_list(i): ",trim(var_list(i)), " dim_list(i): ", dim_list(i), " this%curstep ", this%curstep

        end do

    end subroutine


    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file
    !!
    !!------------------------------------------------------------
    module subroutine init_local(this, options, file_list, var_list, dim_list, start_time, &
                                 lat_var, lon_var, z_var, time_var, p_var, ps_var)
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

        type(variable_t)  :: test_variable
        real, allocatable :: temp_z(:,:,:)


        integer :: i, nx, ny, nz

        ! figure out while file and timestep contains the requested start_time
        call set_curfile_curstep(this, start_time, file_list, time_var)
        call read_bc_time(this%current_time, file_list(this%curfile), time_var, this%curstep)

        !  read in latitude and longitude coordinate data
        call io_read(file_list(this%curfile), lat_var, this%lat, this%curstep)
        call io_read(file_list(this%curfile), lon_var, this%lon, this%curstep)

        ! read in the height coordinate of the input data
        if (.not. options%parameters%compute_z) then
            ! call io_read(file_list(this%curfile), z_var,   temp_z,   this%curstep)
            call io_read(file_list(1), z_var,   temp_z,   1)
            nx = size(temp_z,1)
            ny = size(temp_z,2)
            nz = size(temp_z,3)

            if (allocated(this%z)) deallocate(this%z)
            allocate(this%z(nx,nz,ny))

            this%z = reshape(temp_z, shape=[nx,nz,ny], order=[1,3,2])
        else
            call io_read(file_list(this%curfile), p_var,   temp_z,   this%curstep)
            nx = size(temp_z,1)
            ny = size(temp_z,2)
            nz = size(temp_z,3)

            if (allocated(this%z)) deallocate(this%z)
            allocate(this%z(nx,nz,ny))

        endif


        ! call assert(size(var_list) == size(dim_list), "list of variable dimensions must match list of variables")

        do i=1, size(var_list)

            call add_var_to_dict(this%variables, file_list(this%curfile), var_list(i), dim_list(i), this%curstep, [nx, nz, ny])

        end do

        call update_computed_vars(this, options, update=options%parameters%time_varying_z)

    end subroutine

    !>------------------------------------------------------------
    !! Setup the main geo structure in the boundary structure
    !!
    !!------------------------------------------------------------
    subroutine setup_boundary_geo(this, longitude_system)
        implicit none
        type(boundary_t), intent(inout) :: this
        integer,          intent(in)    :: longitude_system


        if (allocated(this%geo%lat)) deallocate(this%geo%lat)
        allocate( this%geo%lat, source=this%lat)

        if (allocated(this%geo%lon)) deallocate(this%geo%lon)
        allocate( this%geo%lon, source=this%lon)

        call standardize_coordinates(this%geo, longitude_system)

        this%geo_u = this%geo
        this%geo_v = this%geo


        if ( allocated(this%z) )  then
            ! geo%z will be interpolated from this%z to the high-res grids for vinterp in domain... not a great separation
            ! here we save the original z dataset so that it can be used to interpolate varying z through time.
            if (allocated(this%original_geo%z)) deallocate(this%original_geo%z)
            allocate( this%original_geo%z, source=this%z)
        endif

    end subroutine



    !>------------------------------------------------------------
    !! Reads and adds a variable into the variable dictionary
    !!
    !! Given a filename, varname, number of dimensions (2,3) and timestep
    !! Read the timestep of varname from filename and stores the result in a variable structure
    !!
    !! Variable is then added to a master variable dictionary
    !!
    !!------------------------------------------------------------
    subroutine add_var_to_dict(var_dict, file_name, var_name, ndims, timestep, dims)
        implicit none
        type(var_dict_t), intent(inout) :: var_dict
        character(len=*), intent(in)    :: file_name
        character(len=*), intent(in)    :: var_name
        integer,          intent(in)    :: ndims
        integer,          intent(in)    :: timestep
        integer,          intent(in)    :: dims(3)

        real, allocatable :: temp_2d_data(:,:)
        real, allocatable :: temp_3d_data(:,:,:)
        type(variable_t)  :: new_variable
        integer           :: nx,ny,nz


        if (ndims==2) then
            call io_read(file_name, var_name, temp_2d_data, timestep)
            ! print*, "    file_name, var_name ", file_name, trim(var_name)

            call new_variable%initialize( shape( temp_2d_data ) )
            new_variable%data_2d = temp_2d_data

            call var_dict%add_var(var_name, new_variable)
            ! print*, "    added boundary variable ", trim(var_name), " to dict with shape ", shape( temp_2d_data )
            ! do not deallocate data arrays because they are pointed to inside the var_dict now
            ! deallocate(new_variable%data_2d)

        elseif (ndims==3) then
            call io_read(file_name, var_name, temp_3d_data, timestep)

            nx = size(temp_3d_data, 1)
            ny = size(temp_3d_data, 2)
            nz = size(temp_3d_data, 3)


            call new_variable%initialize( [nx,nz,ny] )
            new_variable%data_3d = reshape(temp_3d_data, shape=[nx,nz,ny], order=[1,3,2])

            call var_dict%add_var(var_name, new_variable)

            ! do not deallocate data arrays because they are pointed to inside the var_dict now
            ! deallocate(new_variable%data_3d)

        ! these variables are computed (e.g. pressure from height or height from pressure)
        elseif (ndims==-3) then
            call new_variable%initialize( dims )
            new_variable%computed = .True.

            call var_dict%add_var(var_name, new_variable)
        endif

    end subroutine

    !>------------------------------------------------------------
    !! Reads a new set of forcing data for the next time step
    !!
    !!------------------------------------------------------------
    module subroutine update_forcing(this, options)
        class(boundary_t), intent(inout) :: this
        type(options_t),   intent(inout) :: options

        real, allocatable :: data3d(:,:,:), data2d(:,:)
        type(variable_t)  :: var
        type(variable_t)  :: pvar, zvar, tvar
        character(len=kMAX_NAME_LENGTH) :: name
        integer :: nx, ny, nz, err


        ! if (this_image()==1) then
            call update_forcing_step(this, this%file_list, options%parameters%time_var)

            call read_bc_time(this%current_time, this%file_list(this%curfile), options%parameters%time_var, this%curstep)

            associate(list => this%variables)

            ! loop through the list of variables that need to be read in
            call list%reset_iterator()
            do while (list%has_more_elements())
                ! get the next variable in the structure
                var = list%next(name)

                ! note that pressure (and maybe z eventually?) can be computed from the other so they may not be read
                if (var%computed) then
                    cycle
                elseif (var%three_d) then
                    ! because the data arrays are pointers, this should update the data stored in this%variables
                    call io_read(this%file_list(this%curfile), name, data3d, this%curstep)

                    nx = size(data3d, 1)
                    ny = size(data3d, 2)
                    nz = size(data3d, 3)

                    ! need to vinterp this dataset to the original vertical levels (if necessary)

                    var%data_3d(:,:,:) = reshape(data3d, shape=[nx,nz,ny], order=[1,3,2])

                else if (var%two_d) then
                    call io_read(this%file_list(this%curfile), name, data2d, this%curstep)
                    var%data_2d(:,:) = data2d(:,:)
                endif

            end do

            ! after reading all variables that can be read, not compute any remaining variables (e.g. z from p+ps)
            call update_computed_vars(this, options, update=.True.)

            ! if the vertical levels of the forcing data change over time, they need to be interpolated to the original levels here.
            if (options%parameters%time_varying_z) then
                call interpolate_original_levels(this, options)
            endif

            end associate
        ! endif

        ! call this%distribute_update()

    end subroutine

    subroutine interpolate_original_levels(this, options)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options

        type(variable_t)        :: input_z, var
        type(interpolable_type) :: input_geo
        real, allocatable :: temp_3d(:,:,:)
        character(len=kMAX_NAME_LENGTH) :: name


        associate(list => this%variables)

        input_z = list%get_var(options%parameters%zvar)

        if (options%parameters%z_is_geopotential) then
            input_z%data_3d = input_z%data_3d / gravity
        endif

        ! if (options%parameters%z_is_on_interface) then
        !     call interpolate_in_z(input_z%data_3d)
        ! endif

        allocate(input_geo%z, source=input_z%data_3d)

        ! create a vertical interpolation look up table for the current time step
        call vLUT(this%original_geo, input_geo)


        ! loop through the list of variables that were read in and might need to be interpolated in 3D
        call list%reset_iterator()
        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next(name)

            if (var%three_d) then
                ! need to vinterp this dataset to the original vertical levels (if necessary)

                temp_3d = var%data_3d
                call vinterp(var%data_3d, temp_3d, input_geo%vert_lut)

            endif

        end do
        end associate

    end subroutine


    subroutine update_computed_vars(this, options, update)
        implicit none
        class(boundary_t),   intent(inout)   :: this
        type(options_t),     intent(in)      :: options
        logical,             intent(in),    optional :: update

        integer           :: err
        type(variable_t)  :: var, pvar, zvar, tvar, qvar
        logical :: update_internal

        integer :: nx,ny,nz
        real, allocatable :: temp_z(:,:,:)
        character(len=kMAX_NAME_LENGTH) :: name

        update_internal = .False.
        if (present(update)) update_internal = update

        associate(list => this%variables)

        if (options%parameters%qv_is_relative_humidity) then
            call compute_mixing_ratio_from_rh(list, options)
        endif

        if (options%parameters%qv_is_spec_humidity) then
            call compute_mixing_ratio_from_sh(list, options)
        endif

        ! because z is not updated over time, we don't want to reapply this every time, only in the initialization
        if (.not. update_internal) then
            if (options%parameters%z_is_geopotential) then
                this%z = this%z / gravity
            endif

            if (options%parameters%z_is_on_interface) then
                call interpolate_in_z(this%z)
            endif
        endif

        if (options%parameters%t_offset /= 0) then
            tvar = list%get_var(options%parameters%tvar)
            tvar%data_3d = tvar%data_3d + options%parameters%t_offset
        endif

        ! loop through the list of variables that need to be read in
        call list%reset_iterator()
        do while (list%has_more_elements())

            ! get the next variable in the structure
            var = list%next(name)

            if (var%computed) then

                if (trim(name) == trim(options%parameters%zvar)) then
                    call compute_z_update(this, list, options)
                endif

                if (trim(name) == trim(options%parameters%pvar)) then
                    call compute_p_update(this, list, options, var)
                endif

            endif
        end do

        if (.not.options%parameters%t_is_potential) then

            tvar = list%get_var(options%parameters%tvar)
            pvar = list%get_var(options%parameters%pvar)

            tvar%data_3d = tvar%data_3d / exner_function(pvar%data_3d)
        endif

        end associate

    end subroutine update_computed_vars


    subroutine compute_mixing_ratio_from_rh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: pvar, tvar, qvar

        tvar = list%get_var(options%parameters%tvar)
        pvar = list%get_var(options%parameters%pvar)
        qvar = list%get_var(options%parameters%qvvar)

        if (maxval(qvar%data_3d) > 2) then
            qvar%data_3d = qvar%data_3d/100.0
        endif

        if (options%parameters%t_is_potential) then
            qvar%data_3d = rh_to_mr(qvar%data_3d, tvar%data_3d * exner_function(pvar%data_3d), pvar%data_3d)
        else
            qvar%data_3d = rh_to_mr(qvar%data_3d, tvar%data_3d, pvar%data_3d)
        endif


    end subroutine compute_mixing_ratio_from_rh


    subroutine compute_mixing_ratio_from_sh(list, options)
        implicit none
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: qvar

        qvar = list%get_var(options%parameters%qvvar)

        qvar%data_3d = qvar%data_3d / (1 + qvar%data_3d)

    end subroutine compute_mixing_ratio_from_sh


    subroutine compute_z_update(this, list, options)
        implicit none
        class(boundary_t),  intent(inout)   :: this
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options

        integer           :: err
        type(variable_t)  :: var, pvar, zvar, tvar, qvar
        real, allocatable, target :: real_t(:,:,:)

        real, pointer :: t(:,:,:)

        qvar = list%get_var(options%parameters%qvvar)
        tvar = list%get_var(options%parameters%tvar)
        zvar = list%get_var(options%parameters%hgtvar)
        var = list%get_var(options%parameters%pvar, err)

        pvar = list%get_var(options%parameters%pslvar, err)

        if (options%parameters%t_is_potential) then
            ! stop "Need real air temperature to compute height"
            allocate(real_t, mold=tvar%data_3d)
            real_t = exner_function(var%data_3d) * tvar%data_3d
            t => real_t
        else
            t => tvar%data_3d
        endif

        if (err == 0) then
            call compute_3d_z(var%data_3d, pvar%data_2d, this%z, t, qvar%data_3d)

        else
            pvar = list%get_var(options%parameters%psvar, err)
            if (err == 0) then
                call compute_3d_z(var%data_3d, pvar%data_2d, this%z, t, qvar%data_3d, zvar%data_2d)
            else
                write(*,*) "ERROR reading surface pressure or sea level pressure, variables not found"
                error stop
            endif
        endif
        zvar = list%get_var(options%parameters%zvar)
        zvar%data_3d = this%z

    end subroutine compute_z_update

    subroutine compute_p_update(this, list, options, pressure_var)
        implicit none
        class(boundary_t),  intent(inout)   :: this
        type(var_dict_t),   intent(inout)   :: list
        type(options_t),    intent(in)      :: options
        type(variable_t),   intent(inout)   :: pressure_var

        integer           :: err
        type(variable_t)  :: pvar, zvar, tvar, qvar

        if (options%parameters%t_is_potential) stop "Need real air temperature to compute pressure"

        qvar = list%get_var(options%parameters%qvvar)
        tvar = list%get_var(options%parameters%tvar)
        zvar = list%get_var(options%parameters%hgtvar)

        pvar = list%get_var(options%parameters%pslvar, err)

        if (err == 0) then
            call compute_3d_p(pressure_var%data_3d, pvar%data_2d, this%z, tvar%data_3d, qvar%data_3d, zvar%data_2d)

        else
            pvar = list%get_var(options%parameters%psvar, err)

            if (err == 0) then
                call compute_3d_p(pressure_var%data_3d, pvar%data_2d, this%z, tvar%data_3d, qvar%data_3d)
            else
                write(*,*) "ERROR reading surface pressure or sea level pressure, variables not found"
                error stop
            endif
        endif


    end subroutine compute_p_update


    !>------------------------------------------------------------
    !! Sends the udpated forcing data to all other images
    !!
    !!------------------------------------------------------------
    module subroutine distribute_update(this)
        class(boundary_t), intent(inout) :: this

        type(variable_t)  :: var

        associate(list => this%variables)

        call list%reset_iterator()
        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next()

            if (var%three_d) then
                call broadcast(var%data_3d, 1, 1, num_images(), create_co_array = .True.)

            else if (var%two_d) then
                call broadcast(var%data_2d, 1, 1, num_images(), create_co_array = .True.)
            endif

        enddo

        end associate

        call this%current_time%broadcast(1, 1, num_images())

    end subroutine


    !>------------------------------------------------------------
    !! broadcast a 2d array that may need to be allocated on remote images
    !!
    !! (this needs to be moved to co_util)
    !!------------------------------------------------------------
    subroutine distribute_2d_array(arr, source, first, last)
        implicit none
        real, allocatable, intent(inout) :: arr(:,:)
        integer,           intent(in)    :: source, first, last

        integer :: dims(2)

        ! because this array is probably not allocated, and may not be the correct shape if it is, we have to
        ! broadcast the shape of the array first.
        if (this_image()==source) dims = shape(arr)
        call broadcast(dims, source, first, last, create_co_array=.True.)

        if (allocated(arr)) then
            if (any(dims /= shape(arr))) deallocate(arr)
        endif
        if (.not.allocated(arr)) allocate( arr( dims(1), dims(2) ) )
        call broadcast(arr, source, first, last, create_co_array=.True.)

    end subroutine

    !>------------------------------------------------------------
    !! broadcast a 3d array that may need to be allocated on remote images
    !!
    !! (this needs to be moved to co_util)
    !!------------------------------------------------------------
    subroutine distribute_3d_array(arr, source, first, last)
        implicit none
        real, allocatable, intent(inout) :: arr(:,:,:)
        integer,           intent(in)    :: source, first, last

        integer :: dims(3)

        ! because this array is probably not allocated, and may not be the correct shape if it is, we have to
        ! broadcast the shape of the array first.
        if (this_image()==source) dims = shape(arr)
        call broadcast(dims, source, first, last, create_co_array=.True.)

        if (allocated(arr)) then
            if (any(dims /= shape(arr))) deallocate(arr)
        endif
        if (.not.allocated(arr)) allocate( arr( dims(1), dims(2), dims(3) ) )
        call broadcast(arr, source, first, last, create_co_array=.True.)

    end subroutine

    !>------------------------------------------------------------
    !! Distribute the initial conditions from image 1 into all other images
    !!
    !!------------------------------------------------------------
    module subroutine distribute_initial_conditions(this)
      class(boundary_t), intent(inout) :: this

      integer                         :: number_of_variables
      type(variable_t)                :: temp_variable
      integer                         :: i
      character(len=kMAX_NAME_LENGTH) :: name

      ! needs to distribute lat, lon, time and all vars in var_dict
      ! broadcast(variable, from:image, to: first_image, last_image)
      call distribute_2d_array(this%lat, 1, 1, num_images())
      call distribute_2d_array(this%lon, 1, 1, num_images())
      call distribute_3d_array(this%z,   1, 1, num_images())

      if (this_image() == 1) then
          call this%variables%reset_iterator()
          number_of_variables = this%variables%n_vars
      endif

      call broadcast(number_of_variables, 1, 1, num_images(), create_co_array=.True.)

      do i=1, number_of_variables
          ! if we are the sending image,
          if (this_image() == 1) then
              temp_variable = this%variables%next(name)
          endif

          ! broadcasting is handled within the variable_t object (though it calls broadcast on its components)
          ! we can't make a coarray of a variable_t because it has a pointer to data in it.
          call temp_variable%broadcast(1, 1, num_images())
          call broadcast(name, 1, 1, num_images(), create_co_array=.True.)

          if (this_image() /= 1) then
              call this%variables%add_var(name, temp_variable, save_state=.True.)
          endif
      enddo

      call this%current_time%broadcast(1, 1, num_images())

    end subroutine


    !>------------------------------------------------------------
    !! Setup the vars_to_read and var_dimensions arrays given a master set of variables
    !!
    !! Count the number of variables specified, then allocate and store those variables in a list just their size.
    !! The master list will have all variables, but not all will be set
    !!------------------------------------------------------------
    subroutine setup_variable_lists(master_var_list, master_dim_list, vars_to_read, var_dimensions)
        implicit none
        character(len=kMAX_NAME_LENGTH), intent(in)                 :: master_var_list(:)
        integer,                         intent(in)                 :: master_dim_list(:)
        character(len=kMAX_NAME_LENGTH), intent(inout), allocatable :: vars_to_read(:)
        integer,                         intent(inout), allocatable :: var_dimensions(:)

        integer :: n_valid_vars
        integer :: i, curvar, err

        n_valid_vars = 0
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                n_valid_vars = n_valid_vars + 1
            endif
        enddo

        allocate(vars_to_read(  n_valid_vars), stat=err)
        if (err /= 0) stop "vars_to_read: Allocation request denied"

        allocate(var_dimensions(  n_valid_vars), stat=err)
        if (err /= 0) stop "var_dimensions: Allocation request denied"

        curvar = 1
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                vars_to_read(curvar) = master_var_list(i)
                var_dimensions(curvar) = master_dim_list(i)
                ! if (this_image()==1) print *, "in variable list: ", vars_to_read(curvar)
                curvar = curvar + 1
            endif
        enddo
    end subroutine

    !>------------------------------------------------------------
    !! Find the time step in the input forcing to start the model on
    !!
    !! The model start date (start_time) may not be the same as the first forcing
    !! date (initial_time).  Convert the difference between the two into forcing
    !! steps by dividing by the time delta between forcing steps (in_dt) after
    !! converting in_dt from seconds to days.
    !!
    !! @param  options  model options structure
    !! @retval step     integer number of steps into the forcing sequence
    !!
    !!------------------------------------------------------------


    !>------------------------------------------------------------
    !! Figure out how many time steps are in a file based on a specified variable
    !!
    !! By default assumes that the variable has three dimensions.  If not, var_space_dims must be set
    !!------------------------------------------------------------
    function get_n_timesteps(filename, varname, var_space_dims) result(steps_in_file)
        implicit none
        character(len=*), intent(in) :: filename, varname
        integer,          intent(in), optional :: var_space_dims
        integer :: steps_in_file

        integer :: dims(io_maxDims)
        integer :: space_dims

        space_dims=3
        if (present(var_space_dims)) space_dims = var_space_dims

        call io_getdims(filename, varname, dims)

        if (dims(1) == space_dims) then
            steps_in_file = 1
        else
            steps_in_file = dims(dims(1)+1)
        endif

    end function


    !>------------------------------------------------------------
    !! Set the boundary data structure to the correct time step / file in the list of files
    !!
    !! Reads the time_var from each file successively until it finds a timestep that matches time
    !!------------------------------------------------------------
    subroutine set_curfile_curstep(bc, time, file_list, time_var)
        implicit none
        type(boundary_t),   intent(inout) :: bc
        type(Time_type),    intent(in) :: time
        character(len=*),   intent(in) :: file_list(:)
        character(len=*),   intent(in) :: time_var


        integer :: error

        ! these are module variables that should be correctly set when the subroutine returns
        bc%curfile = 1
        bc%curstep = 1
        error = 1
        bc%curfile = 0
        do while ( (error/=0) .and. (bc%curfile < size(file_list)) )
            bc%curfile = bc%curfile + 1
            bc%curstep = find_timestep_in_file(file_list(bc%curfile), time_var, time, error=error)
        enddo

        if (error==1) then
            stop "Ran out of files to process while searching for matching time variable!"
        endif

    end subroutine


    !>------------------------------------------------------------
    !! Update the curstep and curfile (increments curstep and curfile if necessary)
    !!
    !!------------------------------------------------------------
    subroutine update_forcing_step(bc, file_list, time_var)
        implicit none
        type(boundary_t),   intent(inout) :: bc
        character(len=*),   intent(in) :: file_list(:)
        character(len=*),   intent(in) :: time_var


        integer :: steps_in_file

        bc%curstep = bc%curstep + 1 ! this may be all we have to do most of the time

        ! check that we haven't stepped passed the end of the current file
        steps_in_file = get_n_timesteps(file_list(bc%curfile), time_var, 0)

        if (steps_in_file < bc%curstep) then
            ! if we have, use the next file
            bc%curfile = bc%curfile + 1
            ! and the first timestep in the next file
            bc%curstep = 1

            ! if we have run out of input files, stop with an error message
            if (bc%curfile > size(file_list)) then
                stop "Ran out of files to process while searching for matching time variable!"
            endif

        endif
    end subroutine


    !>------------------------------------------------------------
    !!  Read in the time step from a boundary conditions file if available
    !!
    !!  if not time_var is specified, nothing happens
    !!
    !! Should update this to save the times so they don't all have to be read again every timestep
    !!
    !! @param model_time    Double Scalar to hold time data
    !! @param filename      Name of the NetCDF file to read.
    !! @param varname       Name of the time variable to read from <filename>.
    !! @param curstep       The time step in <filename> to read.
    !!
    !!------------------------------------------------------------
    subroutine read_bc_time(model_time, filename, time_var, curstep)
        implicit none
        type(Time_type),    intent(inout) :: model_time
        character(len=*),   intent(in)    :: filename, time_var
        integer,            intent(in)    :: curstep

        type(Time_type), dimension(:), allocatable :: times

        call read_times(filename, time_var, times, curstep=curstep)
        model_time = times(1)
        deallocate(times)

    end subroutine read_bc_time


end submodule
