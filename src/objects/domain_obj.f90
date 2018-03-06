submodule(domain_interface) domain_implementation
    use assertions_mod,       only : assert, assertions
    use grid_interface,       only : grid_t
    use options_interface,    only : options_t
    use mod_atm_utilities,    only : exner_function, update_pressure
    use icar_constants,       only : kVARS
    use string,               only : str
    use microphysics,         only : mp_simple_var_request
    use co_util,              only : broadcast
    use io_routines,          only : io_read, io_write
    use geo,                  only : geo_lut, geo_interp, geo_interp2d, standardize_coordinates
    use array_utilities,      only : array_offset_x, array_offset_y, smooth_array
    use vertical_interpolation,only : vinterp, vLUT

    implicit none

    interface setup
        module procedure setup_var, setup_exch
    end interface

    ! primary public routines : init, get_initial_conditions, halo_send, halo_retrieve, or halo_exchange
contains


    !> -------------------------------
    !! Initialize the size of the domain
    !!
    !! -------------------------------
    module subroutine init(this, options)
        class(domain_t), intent(inout) :: this
        class(options_t),intent(inout) :: options

        this%dx = options%parameters%dx

        call this%var_request(options)

        call read_domain_shape(this, options)

        call create_variables(this, options)

        call initialize_variables(this, options)

        call setup_meta_data(this)

    end subroutine


    !> -------------------------------
    !! Set up the initial conditions for the domain
    !!
    !! This includes setting up all of the geographic interpolation now that we have the forcing grid
    !! and interpolating the first time step of forcing data on to the high res domain grids
    !!
    !! -------------------------------
    module subroutine get_initial_conditions(this, forcing, options)
      implicit none
      class(domain_t),  intent(inout) :: this
      class(boundary_t),intent(inout) :: forcing
      class(options_t), intent(in)    :: options

      ! create geographic lookup table for domain
      call setup_geo_interpolation(this, forcing)

      ! for all variables with a forcing_var /= "", get forcing, interpolate to local domain
      call interpolate_datasets(this, forcing)

      this%model_time = forcing%current_time

    end subroutine


    !> -------------------------------
    !! Send the halos from all exchangable objects to their neighbors
    !!
    !! -------------------------------
    module subroutine halo_send(this)
      class(domain_t), intent(inout) :: this
      if (associated(this%water_vapor%data_3d))           call this%water_vapor%send()
      if (associated(this%potential_temperature%data_3d)) call this%potential_temperature%send()
      if (associated(this%cloud_water_mass%data_3d))      call this%cloud_water_mass%send()
      if (associated(this%cloud_ice_mass%data_3d))        call this%cloud_ice_mass%send()
      if (associated(this%cloud_ice_number%data_3d))      call this%cloud_ice_number%send()
      if (associated(this%rain_mass%data_3d))             call this%rain_mass%send()
      if (associated(this%rain_number%data_3d))           call this%rain_number%send()
      if (associated(this%snow_mass%data_3d))             call this%snow_mass%send()
      if (associated(this%graupel_mass%data_3d))          call this%graupel_mass%send()
    end subroutine

    !> -------------------------------
    !! Get the halos from all exchangable objects from their neighbors
    !!
    !! -------------------------------
    module subroutine halo_retrieve(this)
      class(domain_t), intent(inout) :: this
      if (associated(this%water_vapor%data_3d))           call this%water_vapor%retrieve() ! the retrieve call will sync all
      if (associated(this%potential_temperature%data_3d)) call this%potential_temperature%retrieve(no_sync=.True.)
      if (associated(this%cloud_water_mass%data_3d))      call this%cloud_water_mass%retrieve(no_sync=.True.)
      if (associated(this%cloud_ice_mass%data_3d))        call this%cloud_ice_mass%retrieve(no_sync=.True.)
      if (associated(this%cloud_ice_number%data_3d))      call this%cloud_ice_number%retrieve(no_sync=.True.)
      if (associated(this%rain_mass%data_3d))             call this%rain_mass%retrieve(no_sync=.True.)
      if (associated(this%rain_number%data_3d))           call this%rain_number%retrieve(no_sync=.True.)
      if (associated(this%snow_mass%data_3d))             call this%snow_mass%retrieve(no_sync=.True.)
      if (associated(this%graupel_mass%data_3d))          call this%graupel_mass%retrieve(no_sync=.True.)
    end subroutine

    !> -------------------------------
    !! Send and get the halos from all exchangable objects to/from their neighbors
    !!
    !! -------------------------------
    module subroutine halo_exchange(this)
      class(domain_t), intent(inout) :: this
      call this%halo_send()

      call this%halo_retrieve()
    end subroutine


    !> -------------------------------
    !! Allocate and or initialize all domain variables if they have been requested
    !!
    !! -------------------------------
    subroutine create_variables(this, opt)
        class(domain_t), intent(inout)  :: this
        class(options_t),intent(in)     :: opt
        integer :: i,j

        integer :: ims, ime, jms, jme

        ims = this%grid%ims
        ime = this%grid%ime
        jms = this%grid%jms
        jme = this%grid%jme

        if (this_image()==1) print *,"  Initializing variables"

        if (0<opt%vars_to_allocate( kVARS%u) )                          call setup(this%u,                        this%u_grid,   forcing_var=opt%parameters%uvar,       list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%v) )                          call setup(this%v,                        this%v_grid,   forcing_var=opt%parameters%vvar,       list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%w) )                          call setup(this%w,                        this%grid )
        if (0<opt%vars_to_allocate( kVARS%water_vapor) )                call setup(this%water_vapor,              this%grid,     forcing_var=opt%parameters%qvvar,      list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%potential_temperature) )      call setup(this%potential_temperature,    this%grid,     forcing_var=opt%parameters%tvar,       list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%cloud_water) )                call setup(this%cloud_water_mass,         this%grid,     forcing_var=opt%parameters%qcvar,      list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%cloud_number_concentration))  call setup(this%cloud_number,             this%grid )
        if (0<opt%vars_to_allocate( kVARS%cloud_ice) )                  call setup(this%cloud_ice_mass,           this%grid,     forcing_var=opt%parameters%qivar,      list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%ice_number_concentration))    call setup(this%cloud_ice_number,         this%grid )
        if (0<opt%vars_to_allocate( kVARS%rain_in_air) )                call setup(this%rain_mass,                this%grid,     forcing_var=opt%parameters%qrvar,      list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%rain_number_concentration))   call setup(this%rain_number,              this%grid )
        if (0<opt%vars_to_allocate( kVARS%snow_in_air) )                call setup(this%snow_mass,                this%grid,     forcing_var=opt%parameters%qsvar,      list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%snow_number_concentration) )  call setup(this%snow_number,              this%grid )
        if (0<opt%vars_to_allocate( kVARS%graupel_in_air) )             call setup(this%graupel_mass,             this%grid,     forcing_var=opt%parameters%qgvar,      list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%graupel_number_concentration))call setup(this%graupel_number,           this%grid )
        if (0<opt%vars_to_allocate( kVARS%precipitation) )              call setup(this%accumulated_precipitation,this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snowfall) )                   call setup(this%accumulated_snowfall,     this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%pressure) )                   call setup(this%pressure,                 this%grid,     forcing_var=opt%parameters%pvar,       list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%temperature) )                call setup(this%temperature,              this%grid )
        if (0<opt%vars_to_allocate( kVARS%exner) )                      call setup(this%exner,                    this%grid )
        if (0<opt%vars_to_allocate( kVARS%z) )                          call setup(this%z,                        this%grid )
        if (0<opt%vars_to_allocate( kVARS%dz_interface) )               call setup(this%dz_interface,             this%grid )
        if (0<opt%vars_to_allocate( kVARS%z_interface) )                call setup(this%z_interface,              this%grid )
        if (0<opt%vars_to_allocate( kVARS%dz) )                         call setup(this%dz_mass,                  this%grid )
        if (0<opt%vars_to_allocate( kVARS%density) )                    call setup(this%density,                  this%grid )
        if (0<opt%vars_to_allocate( kVARS%pressure_interface) )         call setup(this%pressure_interface,       this%grid )
        if (0<opt%vars_to_allocate( kVARS%graupel) )                    call setup(this%graupel,                  this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%shortwave) )                  call setup(this%shortwave,                this%grid2d,   forcing_var=opt%parameters%swdown_var,  list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%longwave) )                   call setup(this%longwave,                 this%grid2d,   forcing_var=opt%parameters%lwdown_var,  list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%vegetation_fraction) )        call setup(this%vegetation_fraction,      this%grid_monthly )
        if (0<opt%vars_to_allocate( kVARS%lai) )                        call setup(this%lai,                      this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%canopy_water) )               call setup(this%canopy_water,             this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snow_water_equivalent) )      call setup(this%snow_water_equivalent,    this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%skin_temperature) )           call setup(this%skin_temperature,         this%grid2d,   forcing_var=opt%parameters%sst_var,     list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%soil_water_content) )         call setup(this%soil_water_content,       this%grid_soil,forcing_var=opt%parameters%soil_t_var,  list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%soil_temperature) )           call setup(this%soil_temperature,         this%grid_soil,forcing_var=opt%parameters%soil_vwc_var,list=this%variables_to_force)
        if (0<opt%vars_to_allocate( kVARS%latitude) )                   call setup(this%latitude,                 this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%longitude) )                  call setup(this%longitude,                this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%u_latitude) )                 call setup(this%u_latitude,               this%u_grid2d )
        if (0<opt%vars_to_allocate( kVARS%u_longitude) )                call setup(this%u_longitude,              this%u_grid2d )
        if (0<opt%vars_to_allocate( kVARS%v_latitude) )                 call setup(this%v_latitude,               this%v_grid2d )
        if (0<opt%vars_to_allocate( kVARS%v_longitude) )                call setup(this%v_longitude,              this%v_grid2d )
        if (0<opt%vars_to_allocate( kVARS%terrain) )                    call setup(this%terrain,                  this%grid2d )

        ! integer variable_t types aren't available yet...
        if (0<opt%vars_to_allocate( kVARS%precipitation) ) allocate(this%precipitation_bucket     (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%snowfall) )      allocate(this%snowfall_bucket          (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%land_cover) )    allocate(this%land_cover_type          (ims:ime, jms:jme),          source=0 )

        ! call setup_forcing_variable

    end subroutine

    !> -------------------------------
    !! Setup a regular variable.
    !!
    !! Initializes the variable
    !! including the forcing_variable if it was set
    !! and adds that variable to the list of variables that has forcing data if the list is supplied
    !! and the forcing_var is both present and not blank ("")
    !!
    !! -------------------------------
    subroutine setup_var(var, grid, forcing_var, list)
        implicit none
        type(variable_t),   intent(inout) :: var
        type(grid_t),       intent(in)    :: grid
        character(len=*),   intent(in),   optional :: forcing_var
        type(var_dict_t),   intent(inout),optional :: list

        if (present(forcing_var)) then
            call var%initialize(grid, forcing_var=forcing_var)

            if (present(list)) then
                if (Len(Trim(forcing_var)) /= 0) then
                    call list%add_var(forcing_var, var)
                endif
            endif
        else

            call var%initialize(grid)
        endif

    end subroutine

    !> -------------------------------
    !! Setup an exchangeable variable.
    !!
    !! Initializes the variable
    !! including the forcing_variable if it was set
    !! and adds that variable to the list of variables that has forcing data if the list is supplied
    !! and the forcing_var is both present and not blank ("")
    !!
    !! -------------------------------
    subroutine setup_exch(var, grid, forcing_var, list)
        implicit none
        type(exchangeable_t),   intent(inout) :: var
        type(grid_t),           intent(in)    :: grid
        character(len=*),       intent(in),   optional :: forcing_var
        type(var_dict_t),       intent(inout),optional :: list

        if (present(forcing_var)) then
            call var%initialize(grid, forcing_var=forcing_var)
            if (present(list)) then
                if (Len(Trim(forcing_var)) /= 0) then
                    call list%add_var(forcing_var, var%meta_data)
                endif
            endif
        else

            call var%initialize(grid)
        endif

    end subroutine

    !> ---------------------------------
    !! Load the data in varname from filename into data_array
    !!
    !! The first / master image reads the file from the disk
    !! Other images get the data broadcast from the master image
    !!
    !! ---------------------------------
    subroutine load_data(filename, varname, data_array, grid)
        implicit none
        character(len=*),  intent(in)   :: filename, varname
        real, allocatable, intent(inout):: data_array(:,:)
        type(grid_t),      intent(in)   :: grid

        if (this_image()==1) then
            call io_read(filename, varname, data_array)
        else
            if (allocated(data_array)) deallocate(data_array)
            allocate(data_array(grid%nx_global, grid%ny_global))
        endif

        call broadcast(data_array, 1, 1, num_images(), .true.)

    end subroutine


    !> ---------------------------------
    !! Read the core model variables from disk
    !!
    !! Reads Terrain, lat, lon and u/v lat/lon on the high-res domain grid
    !! Passing data between images and disk is handled by load_data
    !!
    !! ---------------------------------
    subroutine read_core_variables(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        class(options_t),intent(in)     :: options
        real, allocatable :: temporary_data(:,:)

        ! Read the terrain data
        call load_data(options%parameters%init_conditions_file,   &
                       options%parameters%hgt_hi,                 &
                       temporary_data, this%grid)
        this%terrain%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)

        ! Read the latitude data
        call load_data(options%parameters%init_conditions_file,   &
                       options%parameters%lat_hi,                 &
                       temporary_data, this%grid)
        this%latitude%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)

        ! Read the longitude data
        call load_data(options%parameters%init_conditions_file,   &
                       options%parameters%lon_hi,                 &
                       temporary_data, this%grid)
        this%longitude%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)

        ! Read the u-grid longitude data if specified, other wise interpolate from mass grid
        if (options%parameters%ulon_hi /= "") then
            call load_data(options%parameters%init_conditions_file,   &
                           options%parameters%ulon_hi,                &
                           temporary_data, this%u_grid)
            this%u_longitude%data_2d = temporary_data(this%u_grid%ims:this%u_grid%ime, this%u_grid%jms:this%u_grid%jme)
        else
            associate(ulon => this%u_longitude%data_2d, lon=>this%longitude%data_2d, nx=>this%grid%nx)
                ulon(1,:)    = (1.5 * lon(1,:)  - 0.5 * lon(2,:))    ! extrapolate past the end
                ulon(2:nx,:) = (lon(1:nx-1,:) + lon(2:nx,:)) / 2     ! interpolate between points
                ulon(nx+1,:) = (1.5 * lon(nx,:) - 0.5 * lon(nx-1,:)) ! extrapolate past the end
            end associate
        endif

        ! Read the u-grid latitude data if specified, other wise interpolate from mass grid
        if (options%parameters%ulat_hi /= "") then
            call load_data(options%parameters%init_conditions_file,   &
                           options%parameters%ulat_hi,                &
                           temporary_data, this%u_grid)
            this%u_latitude%data_2d = temporary_data(this%u_grid%ims:this%u_grid%ime, this%u_grid%jms:this%u_grid%jme)

        else
            associate(ulat => this%u_latitude%data_2d, lat=>this%latitude%data_2d, nx=>this%grid%nx)
                ulat(1,:)    = (1.5 * lat(1,:)  - 0.5 * lat(2,:))    ! extrapolate past the end
                ulat(2:nx,:) = (lat(1:nx-1,:) + lat(2:nx,:)) / 2     ! interpolate between points
                ulat(nx+1,:) = (1.5 * lat(nx,:) - 0.5 * lat(nx-1,:)) ! extrapolate past the end
            end associate

        endif

        ! Read the v-grid longitude data if specified, other wise interpolate from mass grid
        if (options%parameters%vlon_hi /= "") then
            call load_data(options%parameters%init_conditions_file,   &
                           options%parameters%vlon_hi,                &
                           temporary_data, this%v_grid)
            this%v_longitude%data_2d = temporary_data(this%v_grid%ims:this%v_grid%ime, this%v_grid%jms:this%v_grid%jme)

        else
            associate(vlon => this%v_longitude%data_2d, lon=>this%longitude%data_2d, ny=>this%grid%ny)
                vlon(:,1)    = (1.5 * lon(:,1)  - 0.5 * lon(:,2))  ! extrapolate past the end
                vlon(:,2:ny) = (lon(:,1:ny-1) + lon(:,2:ny)) / 2   ! interpolate between points
                vlon(:,ny+1) = (1.5 * lon(:,ny) - 0.5 * lon(:,ny-1)) ! extrapolate past the end
            end associate
        endif

        ! Read the v-grid latitude data if specified, other wise interpolate from mass grid
        if (options%parameters%vlat_hi /= "") then
            call load_data(options%parameters%init_conditions_file,   &
                           options%parameters%vlat_hi,                &
                           temporary_data, this%v_grid)
            this%v_latitude%data_2d = temporary_data(this%v_grid%ims:this%v_grid%ime, this%v_grid%jms:this%v_grid%jme)

        else
            associate(vlat => this%v_latitude%data_2d, lat=>this%latitude%data_2d, ny=>this%grid%ny)
                vlat(:,1)    = (1.5 * lat(:,1)  - 0.5 * lat(:,2))   ! extrapolate past the end
                vlat(:,2:ny) = (lat(:,1:ny-1) + lat(:,2:ny)) / 2    ! interpolate between points
                vlat(:,ny+1) = (1.5 * lat(:,ny) - 0.5 * lat(:,ny-1))  ! extrapolate past the end
            end associate
        endif

        if (this_image()==1) write(*,*) "  Finished reading core domain variables"

    end subroutine

    !> -------------------------------
    !! Setup a single Geographic structure given a latitude, longitude, and z array
    !!
    !! -------------------------------
    subroutine setup_geo(geo, latitude, longitude, z)
        implicit none
        class(interpolable_type), intent(inout) :: geo
        real,                     intent(in)    :: latitude(:,:)
        real,                     intent(in)    :: longitude(:,:)
        real,                     intent(in)    :: z(:,:,:)

        if (allocated(geo%lat)) deallocate(geo%lat)
        allocate( geo%lat, source=latitude)

        if (allocated(geo%lon)) deallocate(geo%lon)
        allocate( geo%lon, source=longitude)

        if (allocated(geo%z)) deallocate(geo%z)
        allocate( geo%z, source=z)

        ! This makes 2D variables out of lat/lon if they come in as 1D variables
        ! This also puts the longitudes onto a 0-360 if they are -180-180 (important for Alaska)
        ! Though if working in Europe the -180-180 grid is better ideally the optimal value should be checked.
        ! and good luck if you want to work over the poles...
        call standardize_coordinates(geo)

    end subroutine

    !> -------------------------------
    !! Setup the Geographic structures for all relevant grids in the domain
    !!
    !! -------------------------------
    subroutine setup_domain_geo(this)
        implicit none
        type(domain_t), intent(inout) :: this

        real, allocatable :: temp_z(:,:,:)

        call setup_geo(this%geo,   this%latitude%data_2d,   this%longitude%data_2d,   this%z%data_3d)
        call array_offset_x(this%z%data_3d, temp_z)
        call setup_geo(this%geo_u, this%u_latitude%data_2d, this%u_longitude%data_2d, temp_z)
        call array_offset_y(this%z%data_3d, temp_z)
        call setup_geo(this%geo_v, this%v_latitude%data_2d, this%v_longitude%data_2d, temp_z)

        if (allocated(temp_z)) deallocate(temp_z)
    end subroutine setup_domain_geo

    !> -------------------------------
    !! Initialize various domain variables, mostly z, dz, etc.
    !!
    !! -------------------------------
    subroutine initialize_variables(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        class(options_t),intent(in)     :: options

        integer :: i

        call read_core_variables(this, options)

        associate(z                     => this%z%data_3d,                      &
                  z_interface           => this%z_interface%data_3d,            &
                  dz                    => options%parameters%dz_levels,        &
                  dz_mass               => this%dz_mass%data_3d,                &
                  dz_interface          => this%dz_interface%data_3d,           &
                  exner                 => this%exner%data_3d,                  &
                  pressure              => this%pressure%data_3d,               &
                  temperature           => this%temperature%data_3d,            &
                  potential_temperature => this%potential_temperature%data_3d,  &
                  terrain               => this%terrain%data_2d)

            i = this%grid%kms
            dz_mass(:,i,:)      = dz(i) / 2
            z(:,i,:)            = terrain + dz_mass(:,i,:)
            z_interface(:,i,:)  = terrain

            do i = this%grid%kms+1, this%grid%kme
                dz_mass(:,i,:)     = (dz(i) + dz(i-1)) / 2
                dz_interface(:,i,:)= dz(i)
                z(:,i,:)           = z(:,i-1,:)           + dz_mass(:,i,:)
                z_interface(:,i,:) = z_interface(:,i-1,:) + dz_interface(:,i,:)
            enddo

        end associate

        call setup_domain_geo(this)

    end subroutine initialize_variables

    !> -------------------------------
    !! Populare the metadata structure in the domain for later output
    !!
    !! -------------------------------
    subroutine setup_meta_data(this)
        implicit none
        class(domain_t), intent(inout) :: this

        call this%info%add_attribute("ids",str(this%grid%ids))
        call this%info%add_attribute("ide",str(this%grid%ide))
        call this%info%add_attribute("jds",str(this%grid%jds))
        call this%info%add_attribute("jde",str(this%grid%jde))
        call this%info%add_attribute("kds",str(this%grid%kds))
        call this%info%add_attribute("kde",str(this%grid%kde))

        call this%info%add_attribute("ims",str(this%grid%ims))
        call this%info%add_attribute("ime",str(this%grid%ime))
        call this%info%add_attribute("jms",str(this%grid%jms))
        call this%info%add_attribute("jme",str(this%grid%jme))
        call this%info%add_attribute("kms",str(this%grid%kms))
        call this%info%add_attribute("kme",str(this%grid%kme))

        call this%info%add_attribute("its",str(this%grid%its))
        call this%info%add_attribute("ite",str(this%grid%ite))
        call this%info%add_attribute("jts",str(this%grid%jts))
        call this%info%add_attribute("jte",str(this%grid%jte))
        call this%info%add_attribute("kts",str(this%grid%kts))
        call this%info%add_attribute("kte",str(this%grid%kte))


    end subroutine setup_meta_data


    !> -------------------------------
    !! Add variables needed by all domains to the list of requested variables
    !!
    !! -------------------------------
    module subroutine var_request(this, options)
        class(domain_t), intent(inout) :: this
        class(options_t),intent(inout) :: options

        ! List the variables that are required to be allocated for any domain
        call options%alloc_vars(                                                    &
                     [kVARS%z,                      kVARS%z_interface,              &
                      kVARS%dz,                     kVARS%dz_interface,             &
                      kVARS%terrain,                kVARS%pressure,                 &
                      kVARS%exner,                  kVARS%potential_temperature,    &
                      kVARS%latitude,               kVARS%longitude,                &
                      kVARS%u_latitude,             kVARS%u_longitude,              &
                      kVARS%v_latitude,             kVARS%v_longitude               ])

        ! List the variables that are required for any restart
        call options%restart_vars(                                                  &
                     [kVARS%z,                                                      &
                      kVARS%terrain,                                                &
                      kVARS%latitude,               kVARS%longitude,                &
                      kVARS%u_latitude,             kVARS%u_longitude,              &
                      kVARS%v_latitude,             kVARS%v_longitude               ])

    end subroutine var_request

    !> -------------------------------
    !! Read in the shape of the domain required and setup the grid objects
    !!
    !! -------------------------------
    subroutine read_domain_shape(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        class(options_t),intent(in)     :: options

        real, allocatable :: temporary_data(:,:)
        integer :: nx_global, ny_global, nz_global

        ! This doesn't need to read in this variable, it could just request the dimensions
        ! but this is not a performance sensitive part of the code (for now)
        call io_read(options%parameters%init_conditions_file,   &
                     options%parameters%hgt_hi,                 &
                     temporary_data)

        nx_global = size(temporary_data,1)
        ny_global = size(temporary_data,2)
        nz_global = options%parameters%nz

        call this%grid%set_grid_dimensions(         nx_global, ny_global, nz_global)
        call this%u_grid%set_grid_dimensions(       nx_global, ny_global, nz_global, nx_extra = 1)
        call this%v_grid%set_grid_dimensions(       nx_global, ny_global, nz_global, ny_extra = 1)

        call this%grid2d%set_grid_dimensions(       nx_global, ny_global, 0)
        call this%u_grid2d%set_grid_dimensions(     nx_global, ny_global, 0, nx_extra = 1)
        call this%v_grid2d%set_grid_dimensions(     nx_global, ny_global, 0, ny_extra = 1)

        call this%grid_soil%set_grid_dimensions(    nx_global, ny_global, 4)
        call this%grid_monthly%set_grid_dimensions( nx_global, ny_global, 12)

        deallocate(temporary_data)

    end subroutine

    !> -------------------------------
    !! Check that a set of variables is within realistic bounds (i.e. >0)
    !!
    !! Need to add more variables to the list
    !!
    !! -------------------------------
    module subroutine enforce_limits(this)
      class(domain_t), intent(inout) :: this
      if (associated(this%water_vapor%data_3d)           ) where(this%water_vapor%data_3d < 0)             this%water_vapor%data_3d = 0
      if (associated(this%potential_temperature%data_3d) ) where(this%potential_temperature%data_3d < 0)   this%potential_temperature%data_3d = 0
      if (associated(this%cloud_water_mass%data_3d)      ) where(this%cloud_water_mass%data_3d < 0)        this%cloud_water_mass%data_3d = 0
      if (associated(this%cloud_ice_mass%data_3d)        ) where(this%cloud_ice_mass%data_3d < 0)          this%cloud_ice_mass%data_3d = 0
      if (associated(this%cloud_ice_number%data_3d)      ) where(this%cloud_ice_number%data_3d < 0)        this%cloud_ice_number%data_3d = 0
      if (associated(this%rain_mass%data_3d)             ) where(this%rain_mass%data_3d < 0)               this%rain_mass%data_3d = 0
      if (associated(this%rain_number%data_3d)           ) where(this%rain_number%data_3d < 0)             this%rain_number%data_3d = 0
      if (associated(this%snow_mass%data_3d)             ) where(this%snow_mass%data_3d < 0)               this%snow_mass%data_3d = 0
      if (associated(this%graupel_mass%data_3d)          ) where(this%graupel_mass%data_3d < 0)            this%graupel_mass%data_3d = 0

    end subroutine


    !> -------------------------------
    !! Setup the Geographic look up tables for interpolating a given forcing data set to each of the grids
    !!
    !! -------------------------------
    subroutine setup_geo_interpolation(this, forcing)
        implicit none
        class(domain_t),  intent(inout) :: this
        class(boundary_t),intent(inout) :: forcing

        integer :: nx, ny, nz

        ! this%geo and forcing%geo have to be of class interpolable
        ! which means they must contain lat, lon, z, geolut, and vLUT components
        call geo_LUT(this%geo,   forcing%geo)
        call geo_LUT(this%geo_u, forcing%geo_u)
        call geo_LUT(this%geo_v, forcing%geo_v)

        nx = size(this%geo%z, 1)
        nz = size(forcing%z,  2)
        ny = size(this%geo%z, 3)

        allocate(forcing%geo%z(nx, nz, ny))
        call geo_interp(forcing%geo%z, forcing%z, forcing%geo%geolut)
        call vLUT(this%geo,   forcing%geo)

        allocate(forcing%geo_u%z(nx+1, nz, ny))
        call geo_interp(forcing%geo_u%z, forcing%z, forcing%geo_u%geolut)
        call vLUT(this%geo_u, forcing%geo_u)

        allocate(forcing%geo_v%z(nx, nz, ny+1))
        call geo_interp(forcing%geo_v%z, forcing%z, forcing%geo_v%geolut)
        call vLUT(this%geo_v, forcing%geo_v)

    end subroutine

    !> -------------------------------
    !! Loop through all variables for which forcing data have been supplied and interpolate the forcing data to the domain
    !!
    !! -------------------------------
    subroutine interpolate_datasets(this, forcing)
        implicit none
        class(domain_t),  intent(inout) :: this
        class(boundary_t),intent(in)    :: forcing

        ! temporary to hold the variable to be interpolated
        type(variable_t) :: var_to_interpolate
        integer :: nz

        ! make sure the dictionary is reset to point to the first variable
        call this%variables_to_force%reset_iterator()

        ! No iterate through the dictionary as long as there are more elements present
        do while (this%variables_to_force%has_more_elements())
            ! get the next variable
            var_to_interpolate = this%variables_to_force%next()

            ! interpolate
            call interpolate_variable(var_to_interpolate, forcing, &
                                      vert_interp= (trim(var_to_interpolate%forcing_var) /= trim(this%pressure%forcing_var)) )

        enddo

        nz = size(this%geo%z, 2)
        call update_pressure(this%pressure%data_3d, forcing%geo%z(:,:nz,:), this%geo%z)

    end subroutine

    !> -------------------------------
    !! Interpolate one variable by requesting the forcing data from the boundary data structure then
    !! calling the appropriate interpolation routine (2D vs 3D) with the appropriate grid (mass, u, v)
    !!
    !! -------------------------------
    subroutine interpolate_variable(var, forcing, vert_interp)
        implicit none
        class(variable_t), intent(inout) :: var
        class(boundary_t), intent(in)    :: forcing
        logical,           intent(in),   optional :: vert_interp

        type(variable_t) :: input_data
        ! note that 3D variables have a different number of vertical levels, so they have to first be interpolated
        ! to the high res horizontal grid, then vertically interpolated to the actual icar domain
        real, allocatable :: temp_3d(:,:,:)
        logical :: interpolate_vertically
        integer :: nz

        interpolate_vertically=.True.
        if (present(vert_interp)) interpolate_vertically = vert_interp

        input_data = forcing%variables%get_var(var%forcing_var)

        if (var%two_d) then
            call geo_interp2d(var%data_2d, input_data%data_2d, forcing%geo%geolut)

        else if (var%three_d) then
            ! Sequence of if statements to test if this variable needs to be interpolated onto the staggared grids

            ! allocate a temporary variable to hold the horizontally interpolated data before vertical interpolation
            allocate(temp_3d(size(var%data_3d,1), size(input_data%data_3d,2), size(var%data_3d,3) ))

            ! Interpolate to the Mass grid
            if ((size(var%data_3d,1) == size(forcing%geo%geolut%x,2)).and.(size(var%data_3d,3) == size(forcing%geo%geolut%x,3))) then
                call geo_interp(temp_3d, input_data%data_3d, forcing%geo%geolut)

                ! note that pressure (and possibly other variables?) should not be interpolated, it will be adjusted later
                ! really, it whould be interpolated, and the bottom layers (below the forcing model) should be adjusted separately...
                if (interpolate_vertically) then
                    call vinterp(var%data_3d, temp_3d, forcing%geo%vert_lut)
                else
                    nz = size(var%data_3d,2)
                    var%data_3d = temp_3d(:,:nz,:)
                endif

            ! Interpolate to the u staggered grid
            else if ((size(var%data_3d,1) == size(forcing%geo_u%geolut%x,2)).and.(size(var%data_3d,3) == size(forcing%geo_u%geolut%x,3))) then
                call geo_interp(temp_3d, input_data%data_3d, forcing%geo_u%geolut)
                call vinterp(var%data_3d, temp_3d, forcing%geo_u%vert_lut)
                call smooth_array(var%data_3d, windowsize=10, ydim=3)

            ! Interpolate to the v staggered grid
            else if ((size(var%data_3d,1) == size(forcing%geo_v%geolut%x,2)).and.(size(var%data_3d,3) == size(forcing%geo_v%geolut%x,3))) then
                call geo_interp(temp_3d, input_data%data_3d, forcing%geo_v%geolut)
                call vinterp(var%data_3d, temp_3d, forcing%geo_v%vert_lut)
                call smooth_array(var%data_3d, windowsize=10, ydim=3)
            endif

        else
            write(*,*) "ERROR: unknown variable dimensions in domain variable for forcing:", trim(var%forcing_var)
        endif

    end subroutine

    !> -------------------------------
    !! Used to interpolate an exchangeable type, just gets the meta_data structure from it and uses interpolate_variable
    !!
    !! This is not used presently since the meta_data structures are added directly to the variables_to_force dictionary
    !!
    !! -------------------------------
    subroutine interpolate_exchangeable(var, forcing)
        implicit none
        class(exchangeable_t), intent(inout) :: var
        class(boundary_t),     intent(in)    :: forcing

        ! exchangeables all have a meta_data variable_t component with a pointer to the 3D local data
        call interpolate_variable(var%meta_data, forcing)

    end subroutine


end submodule
