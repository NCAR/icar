submodule(domain_interface) domain_implementation
  use assertions_mod,       only : assert, assertions
  use iso_fortran_env,      only : error_unit
  use grid_interface,       only : grid_t
  use options_interface,    only : options_t
  use mod_atm_utilities,    only : exner_function, sat_mr, pressure_at_elevation
  use icar_constants,       only : kVARS
  use string,               only : str
  use microphysics,         only : mp_simple_var_request
  use co_util,              only : broadcast
  use io_routines,          only : io_read, io_write

  implicit none

contains


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

        if (this_image()==1) print *,"Initializing variables"
        if (0<opt%vars_to_allocate( kVARS%u) )                          call this%u%                        initialize( this%u_grid )
        if (0<opt%vars_to_allocate( kVARS%v) )                          call this%v%                        initialize( this%v_grid )
        if (0<opt%vars_to_allocate( kVARS%w) )                          call this%w%                        initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%water_vapor) )                call this%water_vapor%              initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%potential_temperature) )      call this%potential_temperature%    initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%cloud_water) )                call this%cloud_water_mass%         initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%cloud_number_concentration))  call this%cloud_number%             initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%cloud_ice) )                  call this%cloud_ice_mass%           initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%ice_number_concentration))    call this%cloud_ice_number%         initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%rain_in_air) )                call this%rain_mass%                initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%rain_number_concentration))   call this%rain_number%              initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%snow_in_air) )                call this%snow_mass%                initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%snow_number_concentration) )  call this%snow_number%              initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%graupel_in_air) )             call this%graupel_mass%             initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%graupel_number_concentration))call this%graupel_number%           initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%precipitation) )              call this%accumulated_precipitation%initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snowfall) )                   call this%accumulated_snowfall%     initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%pressure) )                   call this%pressure%                 initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%temperature) )                call this%temperature%              initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%exner) )                      call this%exner%                    initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%z) )                          call this%z%                        initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%dz_interface) )               call this%dz_interface%             initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%z_interface) )                call this%z_interface%              initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%dz) )                         call this%dz_mass%                  initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%density) )                    call this%density%                  initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%pressure_interface) )         call this%pressure_interface%       initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%graupel) )                    call this%graupel%                  initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%shortwave) )                  call this%shortwave%                initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%longwave) )                   call this%longwave%                 initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%vegetation_fraction) )        call this%vegetation_fraction%      initialize( this%grid_monthly )
        if (0<opt%vars_to_allocate( kVARS%lai) )                        call this%lai%                      initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%canopy_water) )               call this%canopy_water%             initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snow_water_equivalent) )      call this%snow_water_equivalent%    initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%skin_temperature) )           call this%skin_temperature%         initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%soil_water_content) )         call this%soil_water_content%       initialize( this%grid_soil )
        if (0<opt%vars_to_allocate( kVARS%soil_temperature) )           call this%soil_temperature%         initialize( this%grid_soil )
        if (0<opt%vars_to_allocate( kVARS%latitude) )                   call this%latitude%                 initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%longitude) )                  call this%longitude%                initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%u_latitude) )                 call this%u_latitude%               initialize( this%u_grid2d )
        if (0<opt%vars_to_allocate( kVARS%u_longitude) )                call this%u_longitude%              initialize( this%u_grid2d )
        if (0<opt%vars_to_allocate( kVARS%v_latitude) )                 call this%v_latitude%               initialize( this%v_grid2d )
        if (0<opt%vars_to_allocate( kVARS%v_longitude) )                call this%v_longitude%              initialize( this%v_grid2d )
        if (0<opt%vars_to_allocate( kVARS%terrain) )                    call this%terrain%                  initialize( this%grid2d )

        ! integer variable_t types aren't available yet...
        if (0<opt%vars_to_allocate( kVARS%precipitation) ) allocate(this%precipitation_bucket     (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%snowfall) )      allocate(this%snowfall_bucket          (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%land_cover) )    allocate(this%land_cover_type          (ims:ime, jms:jme),          source=0 )

        sync all
        if (this_image()==1) print *,"Variable Initialization Complete"

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

        if (this_image()==1) print*, "Finished reading core domain variables"

    end subroutine read_core_variables

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
                      kVARS%terrain,                                                &
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

    subroutine read_domain_shape(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        class(options_t),intent(in)     :: options

        real, allocatable :: temporary_data(:,:)
        integer :: nx_global, ny_global, nz_global

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
    !! Initialize the size of the domain using constant dimensions
    !!
    !! -------------------------------
    module subroutine init(this, options)
      class(domain_t), intent(inout) :: this
      class(options_t),intent(inout) :: options

      call mp_simple_var_request(options)

      call this%var_request(options)

      call read_domain_shape(this, options)

      call create_variables(this, options)

      call initialize_variables(this, options)

      call setup_meta_data(this)

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
    !! Send the halos from all exchangable objects to their neighbors
    !!
    !! -------------------------------
    module subroutine halo_send(this)
      class(domain_t), intent(inout) :: this
      call this%water_vapor%send()
      call this%potential_temperature%send()
      call this%cloud_water_mass%send()
      call this%cloud_ice_mass%send()
      call this%cloud_ice_number%send()
      call this%rain_mass%send()
      call this%rain_number%send()
      call this%snow_mass%send()
      call this%graupel_mass%send()
    end subroutine

    !> -------------------------------
    !! Get the halos from all exchangable objects from their neighbors
    !!
    !! -------------------------------
    module subroutine halo_retrieve(this)
      class(domain_t), intent(inout) :: this
      call this%water_vapor%retrieve() ! the retrieve call will sync all
      call this%potential_temperature%retrieve(no_sync=.True.)
      call this%cloud_water_mass%retrieve(no_sync=.True.)
      call this%cloud_ice_mass%retrieve(no_sync=.True.)
      call this%cloud_ice_number%retrieve(no_sync=.True.)
      call this%rain_mass%retrieve(no_sync=.True.)
      call this%rain_number%retrieve(no_sync=.True.)
      call this%snow_mass%retrieve(no_sync=.True.)
      call this%graupel_mass%retrieve(no_sync=.True.)
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

end submodule
