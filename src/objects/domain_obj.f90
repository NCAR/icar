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


        associate(ims  => this%grid%ims,   &
                  ime  => this%grid%ime,   &
                  jms  => this%grid%jms,   &
                  jme  => this%grid%jme,   &
                  kms  => this%grid%kms,   &
                  kme  => this%grid%kme,   &
                  nsoil=> 4)
            ! used to send data from IO images to all other images
            ! allocate(this%transfer_2d(ims:ime, jms:jme)[*])
            ! allocate(this%transfer_3d(ims:ime, kms:kme, jms:jme)[*])

            if (0<opt%vars_to_allocate( kVARS%precipitation) ) allocate(this%accumulated_precipitation(ims:ime, jms:jme),          source=0.0)
            if (0<opt%vars_to_allocate( kVARS%precipitation) ) allocate(this%precipitation_bucket     (ims:ime, jms:jme),          source=0)
            if (0<opt%vars_to_allocate( kVARS%snowfall) )      allocate(this%accumulated_snowfall     (ims:ime, jms:jme),          source=0.0)
            if (0<opt%vars_to_allocate( kVARS%snowfall) )      allocate(this%snowfall_bucket          (ims:ime, jms:jme),          source=0)
            if (0<opt%vars_to_allocate( kVARS%pressure) )      allocate(this%pressure                 (ims:ime, kms:kme, jms:jme), source=0.0)
            if (0<opt%vars_to_allocate( kVARS%temperature) )   allocate(this%temperature              (ims:ime, kms:kme, jms:jme), source=0.0)
            if (0<opt%vars_to_allocate( kVARS%exner) )         allocate(this%exner                    (ims:ime, kms:kme, jms:jme), source=0.0)
            if (0<opt%vars_to_allocate( kVARS%z) )             allocate(this%z                        (ims:ime, kms:kme, jms:jme), source=0.0)
            if (0<opt%vars_to_allocate( kVARS%dz_interface) )  allocate(this%dz_interface             (ims:ime, kms:kme+1,jms:jme),source=0.0)
            if (0<opt%vars_to_allocate( kVARS%z_interface) )   allocate(this%z_interface              (ims:ime, kms:kme+1,jms:jme),source=0.0)
            if (0<opt%vars_to_allocate( kVARS%dz) )            allocate(this%dz_mass                  (ims:ime, kms:kme, jms:jme), source=0.0)
            if (0<opt%vars_to_allocate( kVARS%density) )       allocate(this%density                  (ims:ime, kms:kme, jms:jme), source=0.0 )

            if (0<opt%vars_to_allocate( kVARS%pressure_interface) )       allocate(this%pressure_interface       (ims:ime, kms:kme, jms:jme), source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%graupel) )                  allocate(this%graupel                  (ims:ime, jms:jme)         , source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%shortwave) )                allocate(this%shortwave                (ims:ime, jms:jme)         , source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%longwave) )                 allocate(this%longwave                 (ims:ime, jms:jme)         , source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%land_cover) )               allocate(this%land_cover_type          (ims:ime, jms:jme)         , source=0 )
            if (0<opt%vars_to_allocate( kVARS%vegetation_fraction) )      allocate(this%vegetation_fraction      (ims:ime, 12, jms:jme)     , source=0.0 ) ! veg fraction could be monthly
            if (0<opt%vars_to_allocate( kVARS%lai) )                      allocate(this%lai                      (ims:ime, jms:jme)         , source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%canopy_water) )             allocate(this%canopy_water             (ims:ime, jms:jme)         , source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%snow_water_equivalent) )    allocate(this%snow_water_equivalent    (ims:ime, jms:jme)         , source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%skin_temperature) )         allocate(this%skin_temperature         (ims:ime, jms:jme)         , source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%soil_water_content) )       allocate(this%soil_water_content       (ims:ime, 1:nsoil, jms:jme), source=0.0 )
            if (0<opt%vars_to_allocate( kVARS%soil_temperature) )         allocate(this%soil_temperature         (ims:ime, 1:nsoil, jms:jme), source=0.0 )

        end associate

    end subroutine

    !> -------------------------------
    !! Initialize various domain variables, mostly z, dz, etc.
    !!
    !! -------------------------------
    subroutine initialize_variables(this)
        implicit none
        class(domain_t) :: this

        integer :: i

        i = this%grid%kms
        this%z(:,i,:) = this%terrain + this%dz_mass(:,i,:)
        this%z_interface(:,i,:) = this%terrain

        do i = this%grid%kms+1, this%grid%kme
            this%z(:,i,:)           = this%z(:,i-1,:)           + this%dz_mass(:,i,:)
            this%z_interface(:,i,:) = this%z_interface(:,i-1,:) + this%dz_interface(:,i,:)
        enddo
        ! note interface elements need to be specified at an additional level
        i = this%grid%kme + 1
        this%z_interface(:,i,:) = this%z_interface(:,i-1,:) + this%dz_interface(:,i,:)

        this%exner             = exner_function(this%pressure)
        this%temperature       = this%exner * this%potential_temperature%local

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
      if (associated(this%water_vapor%local)           ) where(this%water_vapor%local < 0)             this%water_vapor%local = 0
      if (associated(this%potential_temperature%local) ) where(this%potential_temperature%local < 0)   this%potential_temperature%local = 0
      if (associated(this%cloud_water_mass%local)      ) where(this%cloud_water_mass%local < 0)        this%cloud_water_mass%local = 0
      if (associated(this%cloud_ice_mass%local)        ) where(this%cloud_ice_mass%local < 0)          this%cloud_ice_mass%local = 0
      if (associated(this%cloud_ice_number%local)      ) where(this%cloud_ice_number%local < 0)        this%cloud_ice_number%local = 0
      if (associated(this%rain_mass%local)             ) where(this%rain_mass%local < 0)               this%rain_mass%local = 0
      if (associated(this%rain_number%local)           ) where(this%rain_number%local < 0)             this%rain_number%local = 0
      if (associated(this%snow_mass%local)             ) where(this%snow_mass%local < 0)               this%snow_mass%local = 0
      if (associated(this%graupel_mass%local)          ) where(this%graupel_mass%local < 0)            this%graupel_mass%local = 0

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
