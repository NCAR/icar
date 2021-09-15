submodule(output_interface) output_implementation
  use output_metadata,          only : get_metadata
  implicit none

contains

    module subroutine set_domain(this, domain)
        class(output_t),  intent(inout)  :: this
        type(domain_t),   intent(in)     :: domain
        integer :: i

        if (.not.this%is_initialized) call this%init()

        do i=1,domain%info%n_attrs
            call this%add_attribute(domain%info%attributes(i)%name, domain%info%attributes(i)%value)
        enddo

    end subroutine


    module subroutine add_to_output(this, variable)
        class(output_t),   intent(inout) :: this
        type(variable_t),  intent(in)     :: variable

        if (.not.this%is_initialized) call this%init()

        if (associated(variable%data_2d).or.associated(variable%data_3d)) then

            if (this%n_variables == size(this%variables)) call this%increase_var_capacity()

            this%n_variables = this%n_variables + 1

            this%variables(this%n_variables) = variable
        endif

    end subroutine


    module subroutine save_file(this, filename, current_step, time)
        class(output_t),  intent(inout) :: this
        character(len=*), intent(in)    :: filename
        integer,          intent(in)    :: current_step
        type(Time_type),  intent(in)    :: time
        integer :: err

        if (.not.this%is_initialized) call this%init()

        ! open file
        this%filename = filename
        err = nf90_open(filename, NF90_WRITE, this%ncfile_id)
        if (err /= NF90_NOERR) then
            call check( nf90_create(filename, NF90_CLOBBER, this%ncfile_id), "Opening:"//trim(filename))
            this%creating=.True.
        else
            ! in case we need to add a new variable when setting up variables
            call check(nf90_redef(this%ncfile_id), "Setting redefine mode for: "//trim(filename))
        endif

        ! define variables or find variable IDs (and dimensions)
        call setup_variables(this, time)

        if (this%creating) then
            ! add global attributes such as the image number, domain dimension, creation time
            call add_global_attributes(this)

        endif
        ! End define mode. This tells netCDF we are done defining metadata.
        call check( nf90_enddef(this%ncfile_id), "end define mode" )

        ! store output
        call save_data(this, current_step, time)

        this%creating = .false.
        ! close file
        call check(nf90_close(this%ncfile_id), "Closing file "//trim(filename))
    end subroutine

    module subroutine add_variables(this, var_list, domain)
        class(output_t), intent(inout)  :: this
        integer,         intent(in)     :: var_list(:)
        type(domain_t),  intent(in)     :: domain

        if (0<var_list( kVARS%u) )                          call this%add_to_output( get_metadata( kVARS%u                            , domain%u%data_3d))
        if (0<var_list( kVARS%v) )                          call this%add_to_output( get_metadata( kVARS%v                            , domain%v%data_3d))
        if (0<var_list( kVARS%w) )                          call this%add_to_output( get_metadata( kVARS%w                            , domain%w%data_3d))
        if (0<var_list( kVARS%w) )                          call this%add_to_output( get_metadata( kVARS%w_real                       , domain%w_real%data_3d))
        if (0<var_list( kVARS%nsquared) )                   call this%add_to_output( get_metadata( kVARS%nsquared                     , domain%nsquared%data_3d))
        if (0<var_list( kVARS%water_vapor) )                call this%add_to_output( get_metadata( kVARS%water_vapor                  , domain%water_vapor%data_3d))
        if (0<var_list( kVARS%potential_temperature) )      call this%add_to_output( get_metadata( kVARS%potential_temperature        , domain%potential_temperature%data_3d))
        if (0<var_list( kVARS%cloud_water) )                call this%add_to_output( get_metadata( kVARS%cloud_water                  , domain%cloud_water_mass%data_3d))
        if (0<var_list( kVARS%cloud_number_concentration))  call this%add_to_output( get_metadata( kVARS%cloud_number_concentration   , domain%cloud_number%data_3d))
        if (0<var_list( kVARS%cloud_ice) )                  call this%add_to_output( get_metadata( kVARS%cloud_ice                    , domain%cloud_ice_mass%data_3d))
        if (0<var_list( kVARS%ice_number_concentration))    call this%add_to_output( get_metadata( kVARS%ice_number_concentration     , domain%cloud_ice_number%data_3d))
        if (0<var_list( kVARS%rain_in_air) )                call this%add_to_output( get_metadata( kVARS%rain_in_air                  , domain%rain_mass%data_3d))
        if (0<var_list( kVARS%rain_number_concentration))   call this%add_to_output( get_metadata( kVARS%rain_number_concentration    , domain%rain_number%data_3d))
        if (0<var_list( kVARS%snow_in_air) )                call this%add_to_output( get_metadata( kVARS%snow_in_air                  , domain%snow_mass%data_3d))
        if (0<var_list( kVARS%snow_number_concentration) )  call this%add_to_output( get_metadata( kVARS%snow_number_concentration    , domain%snow_number%data_3d))
        if (0<var_list( kVARS%graupel_in_air) )             call this%add_to_output( get_metadata( kVARS%graupel_in_air               , domain%graupel_mass%data_3d))
        if (0<var_list( kVARS%graupel_number_concentration))call this%add_to_output( get_metadata( kVARS%graupel_number_concentration , domain%graupel_number%data_3d))
        if (0<var_list( kVARS%precipitation) )              call this%add_to_output( get_metadata( kVARS%precipitation                , domain%accumulated_precipitation%data_2d))
        if (0<var_list( kVARS%convective_precipitation) )   call this%add_to_output( get_metadata( kVARS%convective_precipitation     , domain%accumulated_convective_pcp%data_2d))
        if (0<var_list( kVARS%snowfall) )                   call this%add_to_output( get_metadata( kVARS%snowfall                     , domain%accumulated_snowfall%data_2d))
        if (0<var_list( kVARS%graupel) )                    call this%add_to_output( get_metadata( kVARS%graupel                      , domain%graupel%data_2d))
        if (0<var_list( kVARS%pressure) )                   call this%add_to_output( get_metadata( kVARS%pressure                     , domain%pressure%data_3d))
        if (0<var_list( kVARS%temperature) )                call this%add_to_output( get_metadata( kVARS%temperature                  , domain%temperature%data_3d))
        if (0<var_list( kVARS%exner) )                      call this%add_to_output( get_metadata( kVARS%exner                        , domain%exner%data_3d))
        if (0<var_list( kVARS%z) )                          call this%add_to_output( get_metadata( kVARS%z                            , domain%z%data_3d))
        if (0<var_list( kVARS%dz_interface) )               call this%add_to_output( get_metadata( kVARS%dz_interface                 , domain%dz_interface%data_3d))
        if (0<var_list( kVARS%z_interface) )                call this%add_to_output( get_metadata( kVARS%z_interface                  , domain%z_interface%data_3d))
        if (0<var_list( kVARS%dz) )                         call this%add_to_output( get_metadata( kVARS%dz                           , domain%dz_mass%data_3d))
        if (0<var_list( kVARS%density) )                    call this%add_to_output( get_metadata( kVARS%density                      , domain%density%data_3d))
        if (0<var_list( kVARS%pressure_interface) )         call this%add_to_output( get_metadata( kVARS%pressure_interface           , domain%pressure_interface%data_3d))
        if (0<var_list( kVARS%cloud_fraction) )             call this%add_to_output( get_metadata( kVARS%cloud_fraction               , domain%cloud_fraction%data_2d))
        if (0<var_list( kVARS%shortwave) )                  call this%add_to_output( get_metadata( kVARS%shortwave                    , domain%shortwave%data_2d))
        if (0<var_list( kVARS%longwave) )                   call this%add_to_output( get_metadata( kVARS%longwave                     , domain%longwave%data_2d))
        if (0<var_list( kVARS%vegetation_fraction) )        call this%add_to_output( get_metadata( kVARS%vegetation_fraction          , domain%vegetation_fraction%data_3d))
        if (0<var_list( kVARS%lai) )                        call this%add_to_output( get_metadata( kVARS%lai                          , domain%lai%data_2d))
        if (0<var_list( kVARS%canopy_water) )               call this%add_to_output( get_metadata( kVARS%canopy_water                 , domain%canopy_water%data_2d))
        if (0<var_list( kVARS%snow_water_equivalent) )      call this%add_to_output( get_metadata( kVARS%snow_water_equivalent        , domain%snow_water_equivalent%data_2d))
        if (0<var_list( kVARS%snow_height) )      call this%add_to_output( get_metadata( kVARS%snow_height        , domain%snow_height%data_2d))
        if (0<var_list( kVARS%skin_temperature) )           call this%add_to_output( get_metadata( kVARS%skin_temperature             , domain%skin_temperature%data_2d))
        if (0<var_list( kVARS%soil_water_content) )         call this%add_to_output( get_metadata( kVARS%soil_water_content           , domain%soil_water_content%data_3d))
        if (0<var_list( kVARS%soil_temperature) )           call this%add_to_output( get_metadata( kVARS%soil_temperature             , domain%soil_temperature%data_3d))
        if (0<var_list( kVARS%latitude) )                   call this%add_to_output( get_metadata( kVARS%latitude                     , domain%latitude%data_2d))
        if (0<var_list( kVARS%longitude) )                  call this%add_to_output( get_metadata( kVARS%longitude                    , domain%longitude%data_2d))
        if (0<var_list( kVARS%u_latitude) )                 call this%add_to_output( get_metadata( kVARS%u_latitude                   , domain%u_latitude%data_2d))
        if (0<var_list( kVARS%u_longitude) )                call this%add_to_output( get_metadata( kVARS%u_longitude                  , domain%u_longitude%data_2d))
        if (0<var_list( kVARS%v_latitude) )                 call this%add_to_output( get_metadata( kVARS%v_latitude                   , domain%v_latitude%data_2d))
        if (0<var_list( kVARS%v_longitude) )                call this%add_to_output( get_metadata( kVARS%v_longitude                  , domain%v_longitude%data_2d))
        if (0<var_list( kVARS%terrain) )                    call this%add_to_output( get_metadata( kVARS%terrain                      , domain%terrain%data_2d))
        if (0<var_list( kVARS%sensible_heat) )              call this%add_to_output( get_metadata( kVARS%sensible_heat                , domain%sensible_heat%data_2d))
        if (0<var_list( kVARS%latent_heat) )                call this%add_to_output( get_metadata( kVARS%latent_heat                  , domain%latent_heat%data_2d))
        if (0<var_list( kVARS%u_10m) )                      call this%add_to_output( get_metadata( kVARS%u_10m                        , domain%u_10m%data_2d))
        if (0<var_list( kVARS%v_10m) )                      call this%add_to_output( get_metadata( kVARS%v_10m                        , domain%v_10m%data_2d))
        if (0<var_list( kVARS%temperature_2m) )             call this%add_to_output( get_metadata( kVARS%temperature_2m               , domain%temperature_2m%data_2d))
        if (0<var_list( kVARS%humidity_2m) )                call this%add_to_output( get_metadata( kVARS%humidity_2m                  , domain%humidity_2m%data_2d))
        if (0<var_list( kVARS%surface_pressure) )           call this%add_to_output( get_metadata( kVARS%surface_pressure             , domain%surface_pressure%data_2d))
        if (0<var_list( kVARS%longwave_up) )                call this%add_to_output( get_metadata( kVARS%longwave_up                  , domain%longwave_up%data_2d))
        if (0<var_list( kVARS%ground_heat_flux) )           call this%add_to_output( get_metadata( kVARS%ground_heat_flux             , domain%ground_heat_flux%data_2d))
        if (0<var_list( kVARS%soil_deep_temperature) )      call this%add_to_output( get_metadata( kVARS%soil_deep_temperature        , domain%soil_deep_temperature%data_2d))
        if (0<var_list( kVARS%soil_totalmoisture) )         call this%add_to_output( get_metadata( kVARS%soil_totalmoisture           , domain%soil_totalmoisture%data_2d))
        if (0<var_list( kVARS%roughness_z0) )               call this%add_to_output( get_metadata( kVARS%roughness_z0                 , domain%roughness_z0%data_2d))


    end subroutine

    subroutine add_global_attributes(this)
        implicit none
        class(output_t), intent(inout)  :: this
        integer :: i

        character(len=19)       :: todays_date_time
        integer,dimension(8)    :: date_time
        character(len=49)       :: date_format
        character(len=5)        :: UTCoffset
        character(len=64)       :: err
        integer                 :: ncid

        ncid = this%ncfile_id

        err="Creating global attributes"
        call check( nf90_put_att(ncid,NF90_GLOBAL,"Conventions","CF-1.6"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"title","Intermediate Complexity Atmospheric Research (ICAR) model output"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"institution","National Center for Atmospheric Research"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"history","Created:"//todays_date_time//UTCoffset), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"references", &
                    "Gutmann et al. 2016: The Intermediate Complexity Atmospheric Model (ICAR). J.Hydrometeor. doi:10.1175/JHM-D-15-0155.1, 2016."), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"contact","Ethan Gutmann : gutmann@ucar.edu"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"git",VERSION), trim(err))

        if (this%n_attrs > 0) then
            do i=1,this%n_attrs
                call check( nf90_put_att(   this%ncfile_id,             &
                                            NF90_GLOBAL,                &
                                            trim(this%attributes(i)%name),    &
                                            trim(this%attributes(i)%value)),  &
                                            "global attr:"//trim(this%attributes(i)%name))
            enddo
        endif

        call date_and_time(values=date_time,zone=UTCoffset)
        date_format='(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)'
        write(todays_date_time,date_format) date_time(1:3),date_time(5:7)

        call check(nf90_put_att(this%ncfile_id, NF90_GLOBAL,"history","Created:"//todays_date_time//UTCoffset), "global attr")
        call check(nf90_put_att(this%ncfile_id, NF90_GLOBAL, "image", this_image()))

    end subroutine add_global_attributes

    subroutine setup_variables(this, time)
        implicit none
        class(output_t), intent(inout) :: this
        type(Time_type), intent(in)    :: time
        integer :: i

        ! iterate through variables creating or setting up variable IDs if they exist, also dimensions
        do i=1,this%n_variables
            ! create all dimensions or find dimension IDs if they exist already
            call setup_dims_for_var(this, this%variables(i))

            call setup_variable(this, this%variables(i))
        end do

        call setup_time_variable(this, time)

    end subroutine setup_variables

    subroutine setup_time_variable(this, time)
        implicit none
        class(output_t), intent(inout) :: this
        type(Time_type), intent(in)    :: time
        integer :: err
        character(len=kMAX_NAME_LENGTH) :: calendar

        associate(var => this%time)
        var%name = "time"
        var%dimensions = [ "time" ]
        var%n_dimensions = 1

        select case (time%calendar)
            case(GREGORIAN)
                calendar = "proleptic_gregorian"
            case(NOLEAP)
                calendar = "noleap"
            case(THREESIXTY)
                calendar = "360-day"
            case default
                calendar = "standard"
        end select


        err = nf90_inq_varid(this%ncfile_id, var%name, var%var_id)

        ! if the variable was not found in the netcdf file then we will define it.
        if (err /= NF90_NOERR) then

            if (allocated(var%dim_ids)) deallocate(var%dim_ids)
            allocate(var%dim_ids(1))

            ! Try to find the dimension ID if it exists already.
            err = nf90_inq_dimid(this%ncfile_id, trim(var%dimensions(1)), var%dim_ids(1))

            ! if the dimension doesn't exist in the file, create it.
            if (err/=NF90_NOERR) then
                call check( nf90_def_dim(this%ncfile_id, trim(var%dimensions(1)), NF90_UNLIMITED, &
                            var%dim_ids(1) ), "def_dim"//var%dimensions(1) )
            endif

            call check( nf90_def_var(this%ncfile_id, var%name, NF90_DOUBLE, var%dim_ids(1), var%var_id), "Defining time" )
            call check( nf90_put_att(this%ncfile_id, var%var_id,"standard_name","time"))
            call check( nf90_put_att(this%ncfile_id, var%var_id,"calendar",trim(calendar)))
            call check( nf90_put_att(this%ncfile_id, var%var_id,"units",time%units()))
            call check( nf90_put_att(this%ncfile_id, var%var_id,"UTCoffset","0"))

        endif
        end associate

    end subroutine setup_time_variable


    subroutine save_data(this, current_step, time)
        implicit none
        class(output_t), intent(in) :: this
        integer,         intent(in) :: current_step
        type(Time_type), intent(in) :: time
        integer :: i
        integer :: dim_3d(3)

        integer :: start_three_D_t(4) = [1,1,1,1]
        integer :: start_two_D_t(3)  = [1,1,1]
        start_three_D_t(4) = current_step
        start_two_D_t(3)   = current_step

        do i=1,this%n_variables
            associate(var => this%variables(i))
                if (var%three_d) then
                    dim_3d = var%dim_len
                    if (var%unlimited_dim) then
                        call check( nf90_put_var(this%ncfile_id, var%var_id,  reshape(var%data_3d, shape=dim_3d, order=[1,3,2]), start_three_D_t),   &
                                    "saving:"//trim(var%name) )
                    elseif (this%creating) then
                        call check( nf90_put_var(this%ncfile_id, var%var_id,  reshape(var%data_3d, shape=dim_3d, order=[1,3,2]) ),   &
                                    "saving:"//trim(var%name) )
                    endif

                elseif (var%two_d) then
                    if (var%unlimited_dim) then
                        call check( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d, start_two_D_t),   &
                                "saving:"//trim(var%name) )
                    elseif (this%creating) then
                        call check( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d),   &
                                "saving:"//trim(var%name) )
                    endif
                endif
            end associate
        end do

        call check( nf90_put_var(this%ncfile_id, this%time%var_id, time%mjd(), [current_step] ),   &
                    "saving:"//trim(this%time%name) )


    end subroutine save_data

    subroutine setup_dims_for_var(this, var)
        implicit none
        class(output_t),    intent(inout) :: this
        type(variable_t),   intent(inout) :: var
        integer :: i, err

        if (allocated(var%dim_ids)) deallocate(var%dim_ids)

        allocate(var%dim_ids(var%n_dimensions))

        do i = 1, size(var%dim_ids)

            ! Try to find the dimension ID if it exists already.
            err = nf90_inq_dimid(this%ncfile_id, trim(var%dimensions(i)), var%dim_ids(i))

            ! probably the dimension doesn't exist in the file, so we will create it.
            if (err/=NF90_NOERR) then

                ! assume that the last dimension should be the unlimited dimension (generally a good idea...)
                if (var%unlimited_dim .and. (i==size(var%dim_ids))) then
                    call check( nf90_def_dim(this%ncfile_id, trim(var%dimensions(i)), NF90_UNLIMITED, &
                                var%dim_ids(i) ), "def_dim"//var%dimensions(i) )
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
        type(variable_t),  intent(inout) :: var
        integer :: i, err

        err = nf90_inq_varid(this%ncfile_id, var%name, var%var_id)

        ! if the variable was not found in the netcdf file then we will define it.
        if (err /= NF90_NOERR) then
            call check( nf90_def_var(this%ncfile_id, var%name, NF90_REAL, var%dim_ids, var%var_id), &
                        "Defining variable:"//trim(var%name) )

            ! setup attributes
            do i=1,size(var%attributes)
                call check( nf90_put_att(this%ncfile_id,                &
                                         var%var_id,                    &
                                         trim(var%attributes(i)%name),        &
                                         trim(var%attributes(i)%value)),      &
                            "saving attribute"//trim(var%attributes(i)%name))
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
        type(variable_t),  allocatable :: new_variables(:)

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
