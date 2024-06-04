submodule(output_interface) output_implementation
    ! use icar_constants,     only : kREAL, kDOUBLE ! these are included in output_interface already
    use output_metadata,    only : get_metadata
    use time_io,            only : get_output_time
    implicit none

contains

    module subroutine init(this, domain, options, file_date_format)
        implicit none
        class(output_t),   intent(inout)  :: this
        type(domain_t),    intent(in)     :: domain
        type(options_t),   intent(in)     :: options
        character(len=49), intent(in)     :: file_date_format
        call this%set_domain(domain)
        call this%add_variables(options%vars_for_restart, domain)
        call this%set_restart_variable(options%parameters%restart_time%as_string(file_date_format))
    end subroutine init

    module subroutine set_domain(this, domain)
        class(output_t),  intent(inout)  :: this
        type(domain_t),   intent(in)     :: domain
        integer :: i

        if (.not.this%is_initialized) call this%init_variables()

        do i=1,domain%info%n_attrs
            call this%add_attribute(domain%info%attributes(i)%name, domain%info%attributes(i)%value)
        enddo

    end subroutine

    module subroutine set_restart_variable(this, restart_time)
        class(output_t),  intent(inout)  :: this
        character(len=25), intent(in) :: restart_time

        if (index(restart_time, "0-00-00_00-00-00") /= 0) then
           this%restarted_from = "Not Restarted"
        else
           this%restarted_from = restart_time
        end if
    end subroutine set_restart_variable


    module subroutine add_to_output(this, variable)
        class(output_t),   intent(inout) :: this
        type(variable_t),  intent(in)     :: variable

        if (.not.this%is_initialized) call this%init_variables()

        if (associated(variable%data_2d).or.associated(variable%data_2dd).or.associated(variable%data_3d)) then

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

        if (.not.this%is_initialized) call this%init_variables()

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

        ! add global attributes such as the image number, domain dimension, creation time
        call add_global_attributes(this)

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
        if (0<var_list( kVARS%precipitation) )              call this%add_to_output( get_metadata( kVARS%precipitation                , domain%accumulated_precipitation%data_2dd))
        if (0<var_list( kVARS%convective_precipitation) )   call this%add_to_output( get_metadata( kVARS%convective_precipitation     , domain%accumulated_convective_pcp%data_2d))
        if (0<var_list( kVARS%snowfall) )                   call this%add_to_output( get_metadata( kVARS%snowfall                     , domain%accumulated_snowfall%data_2dd))
        if (0<var_list( kVARS%graupel) )                    call this%add_to_output( get_metadata( kVARS%graupel                      , domain%graupel%data_2dd))
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
        if (0<var_list( kVARS%shortwave_direct) )           call this%add_to_output( get_metadata( kVARS%shortwave_direct             , domain%shortwave_direct%data_2d))
        if (0<var_list( kVARS%shortwave_diffuse) )          call this%add_to_output( get_metadata( kVARS%shortwave_diffuse            , domain%shortwave_diffuse%data_2d))
        if (0<var_list( kVARS%longwave) )                   call this%add_to_output( get_metadata( kVARS%longwave                     , domain%longwave%data_2d))
        if (0<var_list( kVARS%vegetation_fraction) )        call this%add_to_output( get_metadata( kVARS%vegetation_fraction          , domain%vegetation_fraction%data_3d))
        if (0<var_list( kVARS%vegetation_fraction_max) )    call this%add_to_output( get_metadata( kVARS%vegetation_fraction_max      , domain%vegetation_fraction_max%data_2d))
        if (0<var_list( kVARS%vegetation_fraction_out) )    call this%add_to_output( get_metadata( kVARS%vegetation_fraction_out      , domain%vegetation_fraction_out%data_2d))
        if (0<var_list( kVARS%lai) )                        call this%add_to_output( get_metadata( kVARS%lai                          , domain%lai%data_2d))
        if (0<var_list( kVARS%sai) )                        call this%add_to_output( get_metadata( kVARS%sai                          , domain%sai%data_2d))
        if (0<var_list( kVARS%crop_type) )                  call this%add_to_output( get_metadata( kVARS%crop_type                    , domain%crop_type%data_3d))
        if (0<var_list( kVARS%date_planting) )              call this%add_to_output( get_metadata( kVARS%date_planting                , domain%date_planting%data_2d))
        if (0<var_list( kVARS%date_harvest) )               call this%add_to_output( get_metadata( kVARS%date_harvest                 , domain%date_harvest%data_2d))
        if (0<var_list( kVARS%growing_season_gdd) )         call this%add_to_output( get_metadata( kVARS%growing_season_gdd           , domain%growing_season_gdd%data_2d))
        if (0<var_list( kVARS%irr_frac_total) )             call this%add_to_output( get_metadata( kVARS%irr_frac_total               , domain%irr_frac_total%data_2d))
        if (0<var_list( kVARS%irr_frac_sprinkler) )         call this%add_to_output( get_metadata( kVARS%irr_frac_sprinkler           , domain%irr_frac_sprinkler%data_2d))
        if (0<var_list( kVARS%irr_frac_micro) )             call this%add_to_output( get_metadata( kVARS%irr_frac_micro               , domain%irr_frac_micro%data_2d))
        if (0<var_list( kVARS%irr_frac_flood) )             call this%add_to_output( get_metadata( kVARS%irr_frac_flood               , domain%irr_frac_flood%data_2d))
        if (0<var_list( kVARS%irr_alloc_sprinkler) )        call this%add_to_output( get_metadata( kVARS%irr_alloc_sprinkler          , domain%irr_alloc_sprinkler%data_2d))
        if (0<var_list( kVARS%irr_alloc_micro) )            call this%add_to_output( get_metadata( kVARS%irr_alloc_micro              , domain%irr_alloc_micro%data_2d))
        if (0<var_list( kVARS%irr_alloc_flood) )            call this%add_to_output( get_metadata( kVARS%irr_alloc_flood              , domain%irr_alloc_flood%data_2d))
        if (0<var_list( kVARS%irr_evap_loss_sprinkler) )    call this%add_to_output( get_metadata( kVARS%irr_evap_loss_sprinkler      , domain%irr_evap_loss_sprinkler%data_2d))
        if (0<var_list( kVARS%irr_amt_sprinkler) )          call this%add_to_output( get_metadata( kVARS%irr_amt_sprinkler            , domain%irr_amt_sprinkler%data_2d))
        if (0<var_list( kVARS%irr_amt_micro) )              call this%add_to_output( get_metadata( kVARS%irr_amt_micro                , domain%irr_amt_micro%data_2d))
        if (0<var_list( kVARS%irr_amt_flood) )              call this%add_to_output( get_metadata( kVARS%irr_amt_flood                , domain%irr_amt_flood%data_2d))
        if (0<var_list( kVARS%evap_heat_sprinkler) )        call this%add_to_output( get_metadata( kVARS%evap_heat_sprinkler          , domain%evap_heat_sprinkler%data_2d))
        if (0<var_list( kVARS%mass_ag_grain) )              call this%add_to_output( get_metadata( kVARS%mass_ag_grain                , domain%mass_ag_grain%data_2d))
        if (0<var_list( kVARS%growing_degree_days) )        call this%add_to_output( get_metadata( kVARS%growing_degree_days          , domain%growing_degree_days%data_2d))
        if (0<var_list( kVARS%net_ecosystem_exchange) )     call this%add_to_output( get_metadata( kVARS%net_ecosystem_exchange       , domain%net_ecosystem_exchange%data_2d))
        if (0<var_list( kVARS%gross_primary_prod) )         call this%add_to_output( get_metadata( kVARS%gross_primary_prod           , domain%gross_primary_prod%data_2d))
        if (0<var_list( kVARS%net_primary_prod) )           call this%add_to_output( get_metadata( kVARS%net_primary_prod             , domain%net_primary_prod%data_2d))
        if (0<var_list( kVARS%apar) )                       call this%add_to_output( get_metadata( kVARS%apar                         , domain%apar%data_2d))
        if (0<var_list( kVARS%photosynthesis_total) )       call this%add_to_output( get_metadata( kVARS%photosynthesis_total         , domain%photosynthesis_total%data_2d))
        if (0<var_list( kVARS%stomatal_resist_total) )      call this%add_to_output( get_metadata( kVARS%stomatal_resist_total        , domain%stomatal_resist_total%data_2d))
        if (0<var_list( kVARS%stomatal_resist_sun) )        call this%add_to_output( get_metadata( kVARS%stomatal_resist_sun          , domain%stomatal_resist_sun%data_2d))
        if (0<var_list( kVARS%stomatal_resist_shade) )      call this%add_to_output( get_metadata( kVARS%stomatal_resist_shade        , domain%stomatal_resist_shade%data_2d))
        if (0<var_list( kVARS%gecros_state) )               call this%add_to_output( get_metadata( kVARS%gecros_state                 , domain%gecros_state%data_3d))
        if (0<var_list( kVARS%canopy_water) )               call this%add_to_output( get_metadata( kVARS%canopy_water                 , domain%canopy_water%data_2d))
        if (0<var_list( kVARS%canopy_water_ice) )           call this%add_to_output( get_metadata( kVARS%canopy_water_ice             , domain%canopy_water_ice%data_2d))
        if (0<var_list( kVARS%canopy_water_liquid) )        call this%add_to_output( get_metadata( kVARS%canopy_water_liquid          , domain%canopy_water_liquid%data_2d))
        if (0<var_list( kVARS%canopy_vapor_pressure) )      call this%add_to_output( get_metadata( kVARS%canopy_vapor_pressure        , domain%canopy_vapor_pressure%data_2d))
        if (0<var_list( kVARS%canopy_temperature) )         call this%add_to_output( get_metadata( kVARS%canopy_temperature           , domain%canopy_temperature%data_2d))
        if (0<var_list( kVARS%canopy_fwet) )                call this%add_to_output( get_metadata( kVARS%canopy_fwet                  , domain%canopy_fwet%data_2d))
        if (0<var_list( kVARS%veg_leaf_temperature) )       call this%add_to_output( get_metadata( kVARS%veg_leaf_temperature         , domain%veg_leaf_temperature%data_2d))
        if (0<var_list( kVARS%ground_surf_temperature) )    call this%add_to_output( get_metadata( kVARS%ground_surf_temperature      , domain%ground_surf_temperature%data_2d))
        if (0<var_list( kVARS%frac_within_gap) )            call this%add_to_output( get_metadata( kVARS%frac_within_gap              , domain%frac_within_gap%data_2d))
        if (0<var_list( kVARS%frac_between_gap) )           call this%add_to_output( get_metadata( kVARS%frac_between_gap             , domain%frac_between_gap%data_2d))
        if (0<var_list( kVARS%ground_temperature_bare) )    call this%add_to_output( get_metadata( kVARS%ground_temperature_bare      , domain%ground_temperature_bare%data_2d))
        if (0<var_list( kVARS%ground_temperature_canopy) )  call this%add_to_output( get_metadata( kVARS%ground_temperature_canopy    , domain%ground_temperature_canopy%data_2d))
        if (0<var_list( kVARS%snowfall_ground) )            call this%add_to_output( get_metadata( kVARS%snowfall_ground              , domain%snowfall_ground%data_2d))
        if (0<var_list( kVARS%rainfall_ground) )            call this%add_to_output( get_metadata( kVARS%rainfall_ground              , domain%rainfall_ground%data_2d))
        if (0<var_list( kVARS%snow_water_equivalent) )      call this%add_to_output( get_metadata( kVARS%snow_water_equivalent        , domain%snow_water_equivalent%data_2d))
        if (0<var_list( kVARS%snow_water_eq_prev) )         call this%add_to_output( get_metadata( kVARS%snow_water_eq_prev           , domain%snow_water_eq_prev%data_2d))
        if (0<var_list( kVARS%snow_albedo_prev) )           call this%add_to_output( get_metadata( kVARS%snow_albedo_prev             , domain%snow_albedo_prev%data_2d))
        if (0<var_list( kVARS%snow_temperature) )           call this%add_to_output( get_metadata( kVARS%snow_temperature             , domain%snow_temperature%data_3d))
        if (0<var_list( kVARS%snow_layer_depth) )           call this%add_to_output( get_metadata( kVARS%snow_layer_depth             , domain%snow_layer_depth%data_3d))
        if (0<var_list( kVARS%snow_layer_ice) )             call this%add_to_output( get_metadata( kVARS%snow_layer_ice               , domain%snow_layer_ice%data_3d))
        if (0<var_list( kVARS%snow_layer_liquid_water) )    call this%add_to_output( get_metadata( kVARS%snow_layer_liquid_water      , domain%snow_layer_liquid_water%data_3d))
        if (0<var_list( kVARS%snow_age_factor) )            call this%add_to_output( get_metadata( kVARS%snow_age_factor              , domain%snow_age_factor%data_2d))
        if (0<var_list( kVARS%snow_height) )                call this%add_to_output( get_metadata( kVARS%snow_height                  , domain%snow_height%data_2d))
        if (0<var_list( kVARS%skin_temperature) )           call this%add_to_output( get_metadata( kVARS%skin_temperature             , domain%skin_temperature%data_2d))
        if (0<var_list( kVARS%soil_water_content) )         call this%add_to_output( get_metadata( kVARS%soil_water_content           , domain%soil_water_content%data_3d))
        if (0<var_list( kVARS%eq_soil_moisture) )           call this%add_to_output( get_metadata( kVARS%eq_soil_moisture             , domain%eq_soil_moisture%data_3d))
        if (0<var_list( kVARS%smc_watertable_deep) )        call this%add_to_output( get_metadata( kVARS%smc_watertable_deep          , domain%smc_watertable_deep%data_2d))
        if (0<var_list( kVARS%recharge) )                   call this%add_to_output( get_metadata( kVARS%recharge                     , domain%recharge%data_2d))
        if (0<var_list( kVARS%recharge_deep) )              call this%add_to_output( get_metadata( kVARS%recharge_deep                , domain%recharge_deep%data_2d))
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
        if (0<var_list( kVARS%coeff_momentum_drag) )        call this%add_to_output( get_metadata( kVARS%coeff_momentum_drag          , domain%coeff_momentum_drag%data_2d))
        if (0<var_list( kVARS%coeff_heat_exchange) )        call this%add_to_output( get_metadata( kVARS%coeff_heat_exchange          , domain%coeff_heat_exchange%data_2d))
        if (0<var_list( kVARS%surface_rad_temperature) )    call this%add_to_output( get_metadata( kVARS%surface_rad_temperature      , domain%surface_rad_temperature%data_2d))
        if (0<var_list( kVARS%temperature_2m) )             call this%add_to_output( get_metadata( kVARS%temperature_2m               , domain%temperature_2m%data_2d))
        if (0<var_list( kVARS%humidity_2m) )                call this%add_to_output( get_metadata( kVARS%humidity_2m                  , domain%humidity_2m%data_2d))
        if (0<var_list( kVARS%temperature_2m_veg) )         call this%add_to_output( get_metadata( kVARS%temperature_2m_veg           , domain%temperature_2m_veg%data_2d))
        if (0<var_list( kVARS%temperature_2m_bare) )        call this%add_to_output( get_metadata( kVARS%temperature_2m_bare          , domain%temperature_2m_bare%data_2d))
        if (0<var_list( kVARS%mixing_ratio_2m_veg) )        call this%add_to_output( get_metadata( kVARS%mixing_ratio_2m_veg          , domain%mixing_ratio_2m_veg%data_2d))
        if (0<var_list( kVARS%mixing_ratio_2m_bare) )       call this%add_to_output( get_metadata( kVARS%mixing_ratio_2m_bare         , domain%mixing_ratio_2m_bare%data_2d))
        if (0<var_list( kVARS%surface_pressure) )           call this%add_to_output( get_metadata( kVARS%surface_pressure             , domain%surface_pressure%data_2d))
        if (0<var_list( kVARS%rad_absorbed_total) )         call this%add_to_output( get_metadata( kVARS%rad_absorbed_total           , domain%rad_absorbed_total%data_2d))
        if (0<var_list( kVARS%rad_absorbed_veg) )           call this%add_to_output( get_metadata( kVARS%rad_absorbed_veg             , domain%rad_absorbed_veg%data_2d))
        if (0<var_list( kVARS%rad_absorbed_bare) )          call this%add_to_output( get_metadata( kVARS%rad_absorbed_bare            , domain%rad_absorbed_bare%data_2d))
        if (0<var_list( kVARS%rad_net_longwave) )           call this%add_to_output( get_metadata( kVARS%rad_net_longwave             , domain%rad_net_longwave%data_2d))
        if (0<var_list( kVARS%longwave_up) )                call this%add_to_output( get_metadata( kVARS%longwave_up                  , domain%longwave_up%data_2d))
        if (0<var_list( kVARS%ground_heat_flux) )           call this%add_to_output( get_metadata( kVARS%ground_heat_flux             , domain%ground_heat_flux%data_2d))
        if (0<var_list( kVARS%soil_deep_temperature) )      call this%add_to_output( get_metadata( kVARS%soil_deep_temperature        , domain%soil_deep_temperature%data_2d))
        if (0<var_list( kVARS%evap_canopy) )                call this%add_to_output( get_metadata( kVARS%evap_canopy                  , domain%evap_canopy%data_2d))
        if (0<var_list( kVARS%evap_soil_surface) )          call this%add_to_output( get_metadata( kVARS%evap_soil_surface            , domain%evap_soil_surface%data_2d))
        if (0<var_list( kVARS%transpiration_rate) )         call this%add_to_output( get_metadata( kVARS%transpiration_rate           , domain%transpiration_rate%data_2d))
        if (0<var_list( kVARS%ch_veg) )                     call this%add_to_output( get_metadata( kVARS%ch_veg                       , domain%ch_veg%data_2d))
        if (0<var_list( kVARS%ch_veg_2m) )                  call this%add_to_output( get_metadata( kVARS%ch_veg_2m                    , domain%ch_veg_2m%data_2d))
        if (0<var_list( kVARS%ch_bare) )                    call this%add_to_output( get_metadata( kVARS%ch_bare                      , domain%ch_bare%data_2d))
        if (0<var_list( kVARS%ch_bare_2m) )                 call this%add_to_output( get_metadata( kVARS%ch_bare_2m                   , domain%ch_bare_2m%data_2d))
        if (0<var_list( kVARS%ch_under_canopy) )            call this%add_to_output( get_metadata( kVARS%ch_under_canopy              , domain%ch_under_canopy%data_2d))
        if (0<var_list( kVARS%ch_leaf) )                    call this%add_to_output( get_metadata( kVARS%ch_leaf                      , domain%ch_leaf%data_2d))
        if (0<var_list( kVARS%sensible_heat_veg) )          call this%add_to_output( get_metadata( kVARS%sensible_heat_veg            , domain%sensible_heat_veg%data_2d))
        if (0<var_list( kVARS%sensible_heat_bare) )         call this%add_to_output( get_metadata( kVARS%sensible_heat_bare           , domain%sensible_heat_bare%data_2d))
        if (0<var_list( kVARS%sensible_heat_canopy) )       call this%add_to_output( get_metadata( kVARS%sensible_heat_canopy         , domain%sensible_heat_canopy%data_2d))
        if (0<var_list( kVARS%evap_heat_veg) )              call this%add_to_output( get_metadata( kVARS%evap_heat_veg                , domain%evap_heat_veg%data_2d))
        if (0<var_list( kVARS%evap_heat_bare) )             call this%add_to_output( get_metadata( kVARS%evap_heat_bare               , domain%evap_heat_bare%data_2d))
        if (0<var_list( kVARS%evap_heat_canopy) )           call this%add_to_output( get_metadata( kVARS%evap_heat_canopy             , domain%evap_heat_canopy%data_2d))
        if (0<var_list( kVARS%transpiration_heat) )         call this%add_to_output( get_metadata( kVARS%transpiration_heat           , domain%transpiration_heat%data_2d))
        if (0<var_list( kVARS%ground_heat_veg) )            call this%add_to_output( get_metadata( kVARS%ground_heat_veg              , domain%ground_heat_veg%data_2d))
        if (0<var_list( kVARS%ground_heat_bare) )           call this%add_to_output( get_metadata( kVARS%ground_heat_bare             , domain%ground_heat_bare%data_2d))
        if (0<var_list( kVARS%net_longwave_veg) )           call this%add_to_output( get_metadata( kVARS%net_longwave_veg             , domain%net_longwave_veg%data_2d))
        if (0<var_list( kVARS%net_longwave_bare) )          call this%add_to_output( get_metadata( kVARS%net_longwave_bare            , domain%net_longwave_bare%data_2d))
        if (0<var_list( kVARS%net_longwave_canopy) )        call this%add_to_output( get_metadata( kVARS%net_longwave_canopy          , domain%net_longwave_canopy%data_2d))
        if (0<var_list( kVARS%runoff_surface) )             call this%add_to_output( get_metadata( kVARS%runoff_surface               , domain%runoff_surface%data_2d))
        if (0<var_list( kVARS%runoff_subsurface) )          call this%add_to_output( get_metadata( kVARS%runoff_subsurface            , domain%runoff_subsurface%data_2d))
        if (0<var_list( kVARS%soil_totalmoisture) )         call this%add_to_output( get_metadata( kVARS%soil_totalmoisture           , domain%soil_totalmoisture%data_2d))
        if (0<var_list( kVARS%water_table_depth) )          call this%add_to_output( get_metadata( kVARS%water_table_depth            , domain%water_table_depth%data_2d))
        if (0<var_list( kVARS%water_aquifer) )              call this%add_to_output( get_metadata( kVARS%water_aquifer                , domain%water_aquifer%data_2d))
        if (0<var_list( kVARS%storage_gw) )                 call this%add_to_output( get_metadata( kVARS%storage_gw                   , domain%storage_gw%data_2d))
        if (0<var_list( kVARS%storage_lake) )               call this%add_to_output( get_metadata( kVARS%storage_lake                 , domain%storage_lake%data_2d))
        if (0<var_list( kVARS%roughness_z0) )               call this%add_to_output( get_metadata( kVARS%roughness_z0                 , domain%roughness_z0%data_2d))
        if (0<var_list( kVARS%mass_leaf) )                  call this%add_to_output( get_metadata( kVARS%mass_leaf                    , domain%mass_leaf%data_2d))
        if (0<var_list( kVARS%mass_root) )                  call this%add_to_output( get_metadata( kVARS%mass_root                    , domain%mass_root%data_2d))
        if (0<var_list( kVARS%mass_stem) )                  call this%add_to_output( get_metadata( kVARS%mass_stem                    , domain%mass_stem%data_2d))
        if (0<var_list( kVARS%mass_wood) )                  call this%add_to_output( get_metadata( kVARS%mass_wood                    , domain%mass_wood%data_2d))
        if (0<var_list( kVARS%soil_carbon_fast) )           call this%add_to_output( get_metadata( kVARS%soil_carbon_fast             , domain%soil_carbon_fast%data_2d))
        if (0<var_list( kVARS%soil_carbon_stable) )         call this%add_to_output( get_metadata( kVARS%soil_carbon_stable           , domain%soil_carbon_stable%data_2d))
        if (0<var_list( kVARS%soil_texture_1) )             call this%add_to_output( get_metadata( kVARS%soil_texture_1               , domain%soil_texture_1%data_2d))
        if (0<var_list( kVARS%soil_texture_2) )             call this%add_to_output( get_metadata( kVARS%soil_texture_2               , domain%soil_texture_2%data_2d))
        if (0<var_list( kVARS%soil_texture_3) )             call this%add_to_output( get_metadata( kVARS%soil_texture_3               , domain%soil_texture_3%data_2d))
        if (0<var_list( kVARS%soil_texture_4) )             call this%add_to_output( get_metadata( kVARS%soil_texture_4               , domain%soil_texture_4%data_2d))
        if (0<var_list( kVARS%soil_sand_and_clay) )         call this%add_to_output( get_metadata( kVARS%soil_sand_and_clay           , domain%soil_sand_and_clay%data_3d))
        if (0<var_list( kVARS%re_cloud) )                   call this%add_to_output( get_metadata( kVARS%re_cloud                     , domain%re_cloud%data_3d))
        if (0<var_list( kVARS%re_ice) )                     call this%add_to_output( get_metadata( kVARS%re_ice                       , domain%re_ice%data_3d))
        if (0<var_list( kVARS%re_snow) )                    call this%add_to_output( get_metadata( kVARS%re_snow                      , domain%re_snow%data_3d))
        if (0<var_list( kVARS%out_longwave_rad) )           call this%add_to_output( get_metadata( kVARS%out_longwave_rad             , domain%out_longwave_rad%data_2d))
        if (0<var_list( kVARS%longwave_cloud_forcing) )     call this%add_to_output( get_metadata( kVARS%longwave_cloud_forcing       , domain%longwave_cloud_forcing%data_2d))
        if (0<var_list( kVARS%shortwave_cloud_forcing) )    call this%add_to_output( get_metadata( kVARS%shortwave_cloud_forcing      , domain%shortwave_cloud_forcing%data_2d))
        if (0<var_list( kVARS%cosine_zenith_angle) )        call this%add_to_output( get_metadata( kVARS%cosine_zenith_angle          , domain%cosine_zenith_angle%data_2d))
        if (0<var_list( kVARS%land_emissivity) )            call this%add_to_output( get_metadata( kVARS%land_emissivity              , domain%land_emissivity%data_2d))
        if (0<var_list( kVARS%temperature_interface) )      call this%add_to_output( get_metadata( kVARS%temperature_interface        , domain%temperature_interface%data_3d))
        if (0<var_list( kVARS%tend_swrad) )                 call this%add_to_output( get_metadata( kVARS%tend_swrad                   , domain%tend_swrad%data_3d))
        if (0<var_list( kVARS%t_lake3d) )                   call this%add_to_output( get_metadata( kVARS%t_lake3d                     , domain%t_lake3d%data_3d))
        if (0<var_list( kVARS%lake_icefrac3d) )             call this%add_to_output( get_metadata( kVARS%lake_icefrac3d               , domain%lake_icefrac3d%data_3d))
        if (0<var_list( kVARS%z_lake3d) )                   call this%add_to_output( get_metadata( kVARS%z_lake3d                     , domain%z_lake3d%data_3d))
        if (0<var_list( kVARS%dz_lake3d) )                  call this%add_to_output( get_metadata( kVARS%dz_lake3d                    , domain%dz_lake3d%data_3d))
        if (0<var_list( kVARS%snl2d) )                      call this%add_to_output( get_metadata( kVARS%snl2d                        , domain%snl2d%data_2d))
        if (0<var_list( kVARS%t_grnd2d) )                   call this%add_to_output( get_metadata( kVARS%t_grnd2d                     , domain%t_grnd2d%data_2d))
        if (0<var_list( kVARS%t_soisno3d) )                 call this%add_to_output( get_metadata( kVARS%t_soisno3d                   , domain%t_soisno3d%data_3d))
        if (0<var_list( kVARS%h2osoi_ice3d) )               call this%add_to_output( get_metadata( kVARS%h2osoi_ice3d                 , domain%h2osoi_ice3d%data_3d))
        if (0<var_list( kVARS%h2osoi_liq3d) )               call this%add_to_output( get_metadata( kVARS%h2osoi_liq3d                 , domain%h2osoi_liq3d%data_3d))
        if (0<var_list( kVARS%h2osoi_vol3d) )               call this%add_to_output( get_metadata( kVARS%h2osoi_vol3d                 , domain%h2osoi_vol3d%data_3d))
        if (0<var_list( kVARS%z3d) )                        call this%add_to_output( get_metadata( kVARS%z3d                          , domain%z3d%data_3d))
        if (0<var_list( kVARS%dz3d) )                       call this%add_to_output( get_metadata( kVARS%dz3d                         , domain%dz3d%data_3d))
        if (0<var_list( kVARS%watsat3d) )                   call this%add_to_output( get_metadata( kVARS%watsat3d                     , domain%watsat3d%data_3d))
        if (0<var_list( kVARS%csol3d) )                     call this%add_to_output( get_metadata( kVARS%csol3d                       , domain%csol3d%data_3d))
        if (0<var_list( kVARS%tkmg3d) )                     call this%add_to_output( get_metadata( kVARS%tkmg3d                       , domain%tkmg3d%data_3d))
        if (0<var_list( kVARS%lakemask) )                   call this%add_to_output( get_metadata( kVARS%lakemask                     , domain%lakemask%data_3d))
        if (0<var_list( kVARS%tksatu3d) )                   call this%add_to_output( get_metadata( kVARS%tksatu3d                     , domain%tksatu3d%data_3d))
        if (0<var_list( kVARS%tkdry3d) )                    call this%add_to_output( get_metadata( kVARS%tkdry3d                      , domain%tkdry3d%data_3d))
        if (0<var_list( kVARS%zi3d) )                       call this%add_to_output( get_metadata( kVARS%zi3d                         , domain%zi3d%data_3d))
        if (0<var_list( kVARS%savedtke12d) )                call this%add_to_output( get_metadata( kVARS%savedtke12d                  , domain%savedtke12d%data_2d))
        if (0<var_list( kVARS%lakedepth2d) )                call this%add_to_output( get_metadata( kVARS%lakedepth2d                  , domain%lakedepth2d%data_2d))
        if (0<var_list( kVARS%ivt) )                        call this%add_to_output( get_metadata( kVARS%ivt                          , domain%ivt%data_2d))
        if (0<var_list( kVARS%iwv) )                        call this%add_to_output( get_metadata( kVARS%iwv                          , domain%iwv%data_2d))
        if (0<var_list( kVARS%iwl) )                        call this%add_to_output( get_metadata( kVARS%iwl                          , domain%iwl%data_2d))
        if (0<var_list( kVARS%iwi) )                        call this%add_to_output( get_metadata( kVARS%iwi                          , domain%iwi%data_2d))

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

        ! if not creating the file, only update restarted_from attribute. It
        ! could be opening a previous file whose restart needs to updated
        if (this%creating .eqv. .false.) then
            call check(nf90_put_att(this%ncfile_id, NF90_GLOBAL, "restarted_from", this%restarted_from))
            return
        end if

        ncid = this%ncfile_id
        err="Creating global attributes"
        call check( nf90_put_att(ncid,NF90_GLOBAL,"Conventions","CF-1.6"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"title","Intermediate Complexity Atmospheric Research (ICAR) model output"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"institution","National Center for Atmospheric Research"), trim(err))
        ! initialize todays_date_time variable before use as attribute
        call date_and_time(values=date_time,zone=UTCoffset)
        date_format='(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)'
        write(todays_date_time,date_format) date_time(1:3),date_time(5:7)
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

        call check(nf90_put_att(this%ncfile_id, NF90_GLOBAL, "image", this_image()))
        call check(nf90_put_att(this%ncfile_id, NF90_GLOBAL, "restarted_from", this%restarted_from))

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
                calendar = "360_day"
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
            this%time_units = time%units()
        else
            err = nf90_get_att(this%ncfile_id, var%var_id, "units", this%time_units)
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

        type(Time_type) :: output_time
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
                        if (var%dtype == kREAL) then
                            call check( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d, start_two_D_t),   &
                                    "saving:"//trim(var%name) )
                        elseif (var%dtype == kDOUBLE) then
                            call check( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2dd, start_two_D_t),   &
                                    "saving:"//trim(var%name) )
                        endif
                    elseif (this%creating) then
                        if (var%dtype == kREAL) then
                            call check( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2d),   &
                                    "saving:"//trim(var%name) )
                        elseif (var%dtype == kDOUBLE) then
                            call check( nf90_put_var(this%ncfile_id, var%var_id,  var%data_2dd),   &
                                    "saving:"//trim(var%name) )
                        endif
                    endif
                endif
            end associate
        end do

        output_time = get_output_time(time, units=this%time_units, round_seconds=.True.)

        call check( nf90_put_var(this%ncfile_id, this%time%var_id, dble(output_time%mjd()), [current_step] ),   &
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
            if (var%dtype == kREAL) then
                call check( nf90_def_var(this%ncfile_id, var%name, NF90_REAL, var%dim_ids, var%var_id), &
                            "Defining variable:"//trim(var%name) )
            elseif (var%dtype == kDOUBLE) then
                call check( nf90_def_var(this%ncfile_id, var%name, NF90_DOUBLE, var%dim_ids, var%var_id), &
                            "Defining variable:"//trim(var%name) )
            endif


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

    module subroutine init_variables(this)
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
