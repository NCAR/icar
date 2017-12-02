submodule(domain_interface) domain_implementation
  use assertions_mod,       only : assert, assertions
  use iso_fortran_env,      only : error_unit
  use grid_interface,       only : grid_t
  use options_interface,    only : options_t
  use mod_atm_utilities,    only : exner_function, sat_mr, pressure_at_elevation
  use icar_constants,       only : kVARS
  use string,               only : str
  use microphysics, only : mp_simple_var_request

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
        if (0<opt%vars_to_allocate( kVARS%u) )                          call this%u%                    initialize( this%u_grid )
        if (0<opt%vars_to_allocate( kVARS%v) )                          call this%v%                    initialize( this%v_grid )
        if (0<opt%vars_to_allocate( kVARS%w) )                          call this%w%                    initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%water_vapor) )                call this%water_vapor%          initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%potential_temperature) )      call this%potential_temperature%initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%cloud_water) )                call this%cloud_water_mass%     initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%cloud_number_concentration))  call this%cloud_number%         initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%cloud_ice) )                  call this%cloud_ice_mass%       initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%ice_number_concentration))    call this%cloud_ice_number%     initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%rain_in_air) )                call this%rain_mass%            initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%rain_number_concentration))   call this%rain_number%          initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%snow_in_air) )                call this%snow_mass%            initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%snow_number_concentration) )  call this%snow_number%          initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%graupel_in_air) )             call this%graupel_mass%         initialize( this%grid )
        if (0<opt%vars_to_allocate( kVARS%graupel_number_concentration))call this%graupel_number%       initialize( this%grid )


        if (0<opt%vars_to_allocate( kVARS%precipitation) ) call this%accumulated_precipitation%     initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%precipitation) ) allocate(this%precipitation_bucket     (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%snowfall) )      call this%accumulated_snowfall%          initialize( this%grid2d )
        if (0<opt%vars_to_allocate( kVARS%snowfall) )      allocate(this%snowfall_bucket          (ims:ime, jms:jme),          source=0)
        if (0<opt%vars_to_allocate( kVARS%pressure) )      call this%pressure%                      initialize(this%grid)
        if (0<opt%vars_to_allocate( kVARS%temperature) )   call this%temperature%                   initialize(this%grid)
        if (0<opt%vars_to_allocate( kVARS%exner) )         call this%exner%                         initialize(this%grid)
        if (0<opt%vars_to_allocate( kVARS%z) )             call this%z%                             initialize(this%grid)
        if (0<opt%vars_to_allocate( kVARS%dz_interface) )  call this%dz_interface%                  initialize(this%grid)
        if (0<opt%vars_to_allocate( kVARS%z_interface) )   call this%z_interface%                   initialize(this%grid)
        if (0<opt%vars_to_allocate( kVARS%dz) )            call this%dz_mass%                       initialize(this%grid)
        if (0<opt%vars_to_allocate( kVARS%density) )       call this%density%                       initialize(this%grid)

        if (0<opt%vars_to_allocate( kVARS%pressure_interface) )       call this%pressure_interface%    initialize(this%grid)
        if (0<opt%vars_to_allocate( kVARS%graupel) )                  call this%graupel%               initialize(this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%shortwave) )                call this%shortwave%             initialize(this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%longwave) )                 call this%longwave%              initialize(this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%land_cover) )               allocate(this%land_cover_type          (ims:ime, jms:jme)         , source=0 )
        if (0<opt%vars_to_allocate( kVARS%vegetation_fraction) )      call this%vegetation_fraction%   initialize(this%grid_monthly)
        if (0<opt%vars_to_allocate( kVARS%lai) )                      call this%lai%                   initialize(this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%canopy_water) )             call this%canopy_water%          initialize(this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%snow_water_equivalent) )    call this%snow_water_equivalent% initialize(this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%skin_temperature) )         call this%skin_temperature%      initialize(this%grid2d)
        if (0<opt%vars_to_allocate( kVARS%soil_water_content) )       call this%soil_water_content%    initialize(this%grid_soil)
        if (0<opt%vars_to_allocate( kVARS%soil_temperature) )         call this%soil_temperature%      initialize(this%grid_soil)

    end subroutine

    !> -------------------------------
    !! Initialize various domain variables, mostly z, dz, etc.
    !!
    !! -------------------------------
    subroutine initialize_variables(this)
        implicit none
        class(domain_t) :: this

        integer :: i

        associate(z                     => this%z%data_3d,                      &
                  z_interface           => this%z_interface%data_3d,            &
                  dz_mass               => this%dz_mass%data_3d,                &
                  dz_interface          => this%dz_interface%data_3d,           &
                  exner                 => this%exner%data_3d,                  &
                  pressure              => this%pressure%data_3d,               &
                  temperature           => this%temperature%data_3d,            &
                  potential_temperature => this%potential_temperature%data_3d,  &
                  terrain               => this%terrain%data_2d)

            i = this%grid%kms
            z(:,i,:) = terrain + dz_mass(:,i,:)
            z_interface(:,i,:) = terrain

            do i = this%grid%kms+1, this%grid%kme
                z(:,i,:)           = z(:,i-1,:)           + dz_mass(:,i,:)
                z_interface(:,i,:) = z_interface(:,i-1,:) + dz_interface(:,i,:)
            enddo
            ! note interface elements need to be specified at an additional level
            i = this%grid%kme + 1
            z_interface(:,i,:) = z_interface(:,i-1,:) + dz_interface(:,i,:)

            exner             = exner_function(pressure)
            temperature       = exner * potential_temperature
        end associate

    end subroutine initialize_variables

    !> -------------------------------
    !! Populare the metadata structure in the domain for later output
    !!
    !! -------------------------------
    subroutine setup_meta_data(this)
        implicit none
        class(domain_t), intent(inout) :: this

        call this%info%add_attribute("ids",str(this%ids))
        call this%info%add_attribute("ide",str(this%ide))
        call this%info%add_attribute("jds",str(this%jds))
        call this%info%add_attribute("jde",str(this%jde))
        call this%info%add_attribute("kds",str(this%kds))
        call this%info%add_attribute("kde",str(this%kde))

        call this%info%add_attribute("ims",str(this%ims))
        call this%info%add_attribute("ime",str(this%ime))
        call this%info%add_attribute("jms",str(this%jms))
        call this%info%add_attribute("jme",str(this%jme))
        call this%info%add_attribute("kms",str(this%kms))
        call this%info%add_attribute("kme",str(this%kme))

        call this%info%add_attribute("its",str(this%its))
        call this%info%add_attribute("ite",str(this%ite))
        call this%info%add_attribute("jts",str(this%jts))
        call this%info%add_attribute("jte",str(this%jte))
        call this%info%add_attribute("kts",str(this%kts))
        call this%info%add_attribute("kte",str(this%kte))

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
                     [kVARS%pressure,               kVARS%pressure_interface,       &
                      kVARS%potential_temperature,  kVARS%temperature,              &
                      kVARS%exner,                  kVARS%density,                  &
                      kVARS%z,                      kVARS%dz_interface,             &
                      kVARS%terrain,                                                &
                      kVARS%latitude,               kVARS%longitude,                &
                      kVARS%u_latitude,             kVARS%u_longitude,              &
                      kVARS%v_latitude,             kVARS%v_longitude               ])

        ! List the variables that are required for any restart
        call options%restart_vars(                                                  &
                     [kVARS%pressure,                                               &
                      kVARS%potential_temperature,                                  &
                      kVARS%z,                                                      &
                      kVARS%terrain,                                                &
                      kVARS%latitude,               kVARS%longitude,                &
                      kVARS%u_latitude,             kVARS%u_longitude,              &
                      kVARS%v_latitude,             kVARS%v_longitude               ])

    end subroutine var_request

    !> -------------------------------
    !! Initialize the size of the domain using constant dimensions
    !!
    !! -------------------------------
    module subroutine init(this, options)
      class(domain_t), intent(inout) :: this
      class(options_t),intent(inout) :: options

      integer :: nx_global, ny_global, nz_global
      ! read input file...
      nx_global=100
      ny_global=200
      nz_global=30

      call this%grid%set_grid_dimensions(nx_global, ny_global, nz_global)
      call this%grid2d%set_grid_dimensions(nx_global, ny_global, 0)
      call this%u_grid%set_grid_dimensions(nx_global, ny_global, nz_global, nx_extra = 1)
      call this%v_grid%set_grid_dimensions(nx_global, ny_global, nz_global, ny_extra = 1)

      call this%grid_soil%set_grid_dimensions(nx_global, ny_global, 4)
      call this%grid_monthly%set_grid_dimensions(nx_global, ny_global, 12)

      call mp_simple_var_request(options)

      call this%var_request(options)

      call create_variables(this, options)

      ! call initialize_variables(this, options)

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
      call this%water_vapor%send()
      call this%potential_temperature%send()
      call this%cloud_water_mass%send()
      call this%cloud_ice_mass%send()
      call this%cloud_ice_number%send()
      call this%rain_mass%send()
      call this%rain_number%send()
      call this%snow_mass%send()
      call this%graupel_mass%send()

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

end submodule
