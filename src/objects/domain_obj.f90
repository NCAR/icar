submodule(domain_interface) domain_implementation
  use assertions_mod,       only : assert, assertions
  use iso_fortran_env,      only : error_unit
  use grid_interface,       only : grid_t
  use options_interface,    only : options_t
  use mod_atm_utilities,    only : exner_function, sat_mr, pressure_at_elevation
  use string,               only : str

  implicit none

contains


    subroutine initialize_variables(this)
        class(domain_t), intent(inout) :: this
        integer :: i,j

        if (this_image()==1) print *,"Initializing variables"
        call this%u%initialize(                     this%u_grid)
        call this%v%initialize(                     this%v_grid)
        call this%w%initialize(                     this%grid)
        call this%water_vapor%initialize(           this%grid)
        call this%potential_temperature%initialize( this%grid)
        call this%cloud_water_mass%initialize(      this%grid)
        call this%cloud_ice_mass%initialize(        this%grid)
        call this%cloud_ice_number%initialize(      this%grid)
        call this%rain_mass%initialize(             this%grid)
        call this%rain_number%initialize(           this%grid)
        call this%snow_mass%initialize(             this%grid)
        call this%graupel_mass%initialize(          this%grid)


    !   allocate(this%transfer_array_2d(this%nx, this%ny_global)[*])
    !   allocate(this%transfer_array_3d(this%nx, this%nz, this%ny_global)[*])

        associate(ims=>this%grid%ims,   &
                  ime=>this%grid%ime,   &
                  jms=>this%grid%jms,   &
                  jme=>this%grid%jme,   &
                  kms=>this%grid%kms,   &
                  kme=>this%grid%kme    &
            )
          allocate(this%accumulated_precipitation(ims:ime, jms:jme))
          allocate(this%accumulated_snowfall     (ims:ime, jms:jme))
          allocate(this%pressure                 (ims:ime, kms:kme, jms:jme))
          allocate(this%temperature              (ims:ime, kms:kme, jms:jme))
          allocate(this%exner                    (ims:ime, kms:kme, jms:jme))
          allocate(this%z                        (ims:ime, kms:kme, jms:jme))
          allocate(this%dz_interface             (ims:ime, kms:kme, jms:jme))
          allocate(this%z_interface              (ims:ime, kms:kme, jms:jme))
          allocate(this%dz_mass                  (ims:ime, kms:kme, jms:jme))



          do i=kms+1,kme
              this%z(:,i,:)           = this%z(:,i-1,:)           + this%dz_mass(:,i,:)
              this%z_interface(:,i,:) = this%z_interface(:,i-1,:) + this%dz_interface(:,i,:)
              this%pressure(:,i,:)    = pressure_at_elevation(102000.0, this%z(:,i,:))
          enddo

          this%exner             = exner_function(this%pressure)
          this%temperature       = this%exner * this%potential_temperature%local
          this%water_vapor%local = sat_mr(this%temperature,this%pressure)
      end associate

    end subroutine

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
    !! Initialize the size of the domain using constant dimensions
    !!
    !! -------------------------------
    module subroutine init(this, options)
      class(domain_t), intent(inout) :: this
      class(options_t),intent(in)    :: options

      integer :: nx_global, ny_global, nz_global
      ! read input file...

      call this%grid%get_grid_dimensions(nx_global, ny_global, nz_global)
      call this%u_grid%get_grid_dimensions(nx_global, ny_global, nz_global, nx_extra = 1)
      call this%v_grid%get_grid_dimensions(nx_global, ny_global, nz_global, ny_extra = 1)

      call initialize_variables(this)

      call setup_meta_data(this)

    end subroutine


    module subroutine enforce_limits(this)
      class(domain_t), intent(inout) :: this
      where(this%water_vapor%local < 0)             this%water_vapor%local = 0
      where(this%potential_temperature%local < 0)   this%potential_temperature%local = 0
      where(this%cloud_water_mass%local < 0)        this%cloud_water_mass%local = 0
      where(this%cloud_ice_mass%local < 0)          this%cloud_ice_mass%local = 0
      where(this%cloud_ice_number%local < 0)        this%cloud_ice_number%local = 0
      where(this%rain_mass%local < 0)               this%rain_mass%local = 0
      where(this%rain_number%local < 0)             this%rain_number%local = 0
      where(this%snow_mass%local < 0)               this%snow_mass%local = 0
      where(this%graupel_mass%local < 0)            this%graupel_mass%local = 0

    end subroutine

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

    module subroutine halo_retrieve(this)
      class(domain_t), intent(inout) :: this
      call this%water_vapor%retrieve()
      call this%potential_temperature%retrieve(no_sync=.True.)
      call this%cloud_water_mass%retrieve(no_sync=.True.)
      call this%cloud_ice_mass%retrieve(no_sync=.True.)
      call this%cloud_ice_number%retrieve(no_sync=.True.)
      call this%rain_mass%retrieve(no_sync=.True.)
      call this%rain_number%retrieve(no_sync=.True.)
      call this%snow_mass%retrieve(no_sync=.True.)
      call this%graupel_mass%retrieve(no_sync=.True.)
    end subroutine

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

      call this%water_vapor%retrieve()
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
