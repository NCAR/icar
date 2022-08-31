module output_metadata

    use icar_constants
    use variable_interface,     only : variable_t
    use meta_data_interface,    only : attribute_t
    implicit none

    type(variable_t), allocatable :: var_meta(:)

    !>------------------------------------------------------------
    !! Generic interface to the netcdf read routines
    !!------------------------------------------------------------
    interface get_metadata
        module procedure get_metadata_2d, get_metadata_3d, get_metadata_nod
    end interface


contains

    !>------------------------------------------------------------
    !! Get generic metadata for a no-data object
    !!------------------------------------------------------------
    function get_metadata_nod(var_idx) result(meta_data)
        implicit none
        integer, intent(in) :: var_idx
        type(variable_t) :: meta_data

        if (var_idx > kMAX_STORAGE_VARS) then
            stop "Invalid variable metadata requested"
        endif

        ! initialize the module level var_meta data
        if (.not.allocated(var_meta)) call init_var_meta()

        ! get the var_idx index into var_meta to return
        meta_data = var_meta(var_idx)

        ! set the dimensionality to false
        meta_data%two_d     = .False.
        meta_data%three_d   = .False.

    end function get_metadata_nod

    !>------------------------------------------------------------
    !! Get generic metadata for a two-dimensional variable
    !!
    !! Sets the internal data pointer to point to the input data provided
    !!------------------------------------------------------------
    function get_metadata_2d(var_idx, input_data) result(meta_data)
        implicit none
        integer, intent(in)          :: var_idx
        real,    intent(in), pointer :: input_data(:,:)

        type(variable_t) :: meta_data       ! function result
        integer          :: local_shape(2)  ! store the shape of the input data array

        if (var_idx>kMAX_STORAGE_VARS) then
            stop "Invalid variable metadata requested"
        endif

        if (.not.allocated(var_meta)) call init_var_meta()

        meta_data = var_meta(var_idx)

        if (associated(input_data)) then
            meta_data%data_2d   => input_data
            meta_data%two_d     = .True.
            meta_data%three_d   = .False.
            local_shape(1) = size(input_data, 1)
            local_shape(2) = size(input_data, 2)
            ! for some reason if shape(input_data) is passed as source, then the dim_len bounds are (0:1) instead of 1:2
            allocate(meta_data%dim_len, source=local_shape)
        endif

    end function get_metadata_2d

    !>------------------------------------------------------------
    !! Get generic metadata for a three-dimensional variable
    !!
    !! Sets the internal data pointer to point to the input data provided
    !!------------------------------------------------------------
    function get_metadata_3d(var_idx, input_data) result(meta_data)
        implicit none
        integer, intent(in)          :: var_idx
        real,    intent(in), pointer :: input_data(:,:,:)

        type(variable_t) :: meta_data       ! function result
        integer          :: local_shape(3)  ! store the shape of the input data array

        if (var_idx>kMAX_STORAGE_VARS) then
            stop "Invalid variable metadata requested"
        endif

        ! initialize the module level constant data structure
        if (.not.allocated(var_meta)) call init_var_meta()

        meta_data = var_meta(var_idx)

        if (associated(input_data)) then
            meta_data%data_3d   => input_data
            meta_data%two_d     = .False.
            meta_data%three_d   = .True.
            local_shape(1) = size(input_data, 1)
            local_shape(2) = size(input_data, 3)
            local_shape(3) = size(input_data, 2)
            ! for some reason if shape(input_data) is passed as source, then the dim_len bounds are (0:1) instead of 1:2
            allocate(meta_data%dim_len, source=local_shape)
        endif

    end function get_metadata_3d


    !>------------------------------------------------------------
    !! Get metadata variable name associated with a given index
    !!
    !! Sets the internal data pointer to point to the input data provided
    !!------------------------------------------------------------
    function get_varname(var_idx) result(name)
        implicit none
        integer, intent(in)             :: var_idx
        character(len=kMAX_NAME_LENGTH) :: name       ! function result

        type(variable_t) :: meta_data

        if (var_idx>kMAX_STORAGE_VARS) then
            stop "Invalid variable metadata requested"
        endif

        ! initialize the module level constant data structure
        if (.not.allocated(var_meta)) call init_var_meta()

        meta_data = var_meta(var_idx)

        name = meta_data%name

    end function get_varname

    !>------------------------------------------------------------
    !! Initialize the module level master data structure
    !!
    !!------------------------------------------------------------
    subroutine init_var_meta()
        implicit none
        integer :: i

        character(len=16) :: three_d_u_t_dimensions(4)          = [character(len=16) :: "lon_u","lat_y","level","time"]
        character(len=16) :: three_d_v_t_dimensions(4)          = [character(len=16) :: "lon_x","lat_v","level","time"]
        character(len=16) :: three_d_dimensions(3)              = [character(len=16) :: "lon_x","lat_y","level"]
        character(len=16) :: three_d_t_dimensions(4)            = [character(len=16) :: "lon_x","lat_y","level","time"]
        character(len=16) :: three_d_interface_dimensions(3)    = [character(len=16) :: "lon_x","lat_y","level_i"]
        character(len=16) :: three_d_t_interface_dimensions(4)  = [character(len=16) :: "lon_x","lat_y","level_i","time"]
        character(len=16) :: two_d_dimensions(2)                = [character(len=16) :: "lon_x","lat_y"]
        character(len=16) :: two_d_t_dimensions(3)              = [character(len=16) :: "lon_x","lat_y","time"]
        character(len=16) :: two_d_u_dimensions(2)              = [character(len=16) :: "lon_u","lat_y"]
        character(len=16) :: two_d_v_dimensions(2)              = [character(len=16) :: "lon_x","lat_v"]
        character(len=16) :: three_d_t_soil_dimensions(4)       = [character(len=16) :: "lon_x","lat_y","nsoil","time"]
        character(len=16) :: three_d_t_snow_dimensions(4)       = [character(len=16) :: "lon_x","lat_y","nsnow","time"]
        character(len=16) :: three_d_t_snowsoil_dimensions(4)   = [character(len=16) :: "lon_x","lat_y","nsnowsoil","time"]
        character(len=16) :: three_d_soilcomp_dimensions(3)     = [character(len=16) :: "lon_x","lat_y","nsoil_composition"]
        character(len=16) :: three_d_crop_dimensions(3)         = [character(len=16) :: "lon_x","lat_y","crop"]
        character(len=16) :: three_d_t_gecros_dimensions(4)     = [character(len=16) :: "lon_x","lat_y","gecros","time"]
        character(len=16) :: two_d_month_dimensions(3)          = [character(len=16) :: "lon_x","lat_y","month"]
        character(len=16) :: three_d_t_lake_dimensions(4)           = [character(len=16) :: "lon_x","lat_y","nlevlake","time"]
        character(len=16) :: three_d_t_lake_soisno_dimensions(4)    = [character(len=16) :: "lon_x","lat_y","nlevsoisno","time"] !grid_lake_soisno
        character(len=16) :: three_d_t_lake_soisno_1_dimensions(4)  = [character(len=16) :: "lon_x","lat_y","nlevsoisno_1","time"] 
        character(len=16) :: three_d_t_lake_soi_dimensions(4)       = [character(len=16) :: "lon_x","lat_y","nlevsoi_lake","time"] !grid_lake_soi
        

        if (allocated(var_meta)) deallocate(var_meta)

        if (kVARS%last_var/=kMAX_STORAGE_VARS) then
            stop "ERROR: variable indicies not correctly initialized"
        endif
        ! allocate the actual data structure to be used
        allocate(var_meta(kMAX_STORAGE_VARS))


        !>------------------------------------------------------------
        !!  U  East West Winds
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%u))
            var%name        = "u"
            var%dimensions  = three_d_u_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "grid_eastward_wind"),              &
                               attribute_t("long_name",     "Grid relative eastward wind"),     &
                               attribute_t("units",         "m s-1"),                           &
                               attribute_t("coordinates",   "u_lat u_lon")]
        end associate
        !>------------------------------------------------------------
        !!  V  North South Winds
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%v))
            var%name        = "v"
            var%dimensions  = three_d_v_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "grid_northward_wind"),             &
                               attribute_t("long_name",     "Grid relative northward wind"),    &
                               attribute_t("units",         "m s-1"),                           &
                               attribute_t("coordinates",   "v_lat v_lon")]
        end associate
        !>------------------------------------------------------------
        !!  W  Vertical Winds
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%w))
            var%name        = "w_grid"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "grid_upward_air_velocity"),    &
                               attribute_t("long_name",     "Vertical wind"),                   &
                               attribute_t("description",   "Vertical wind relative to the grid"),&
                               attribute_t("units",         "m s-1"),                           &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        associate(var=>var_meta(kVARS%w_real))
            var%name        = "w"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "upward_air_velocity"),             &
                               attribute_t("long_name",     "Vertical wind"),                   &
                               attribute_t("description",   "Vertical wind including u/v"),     &
                               attribute_t("units",         "m s-1"),                           &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Brunt Vaisala frequency (squared)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%nsquared))
            var%name        = "nsquared"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "square_of_brunt_vaisala_frequency_in_air"),&
                               attribute_t("long_name",     "Burnt Vaisala frequency squared"), &
                               attribute_t("units",         "s-2"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Air Pressure
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%pressure))
            var%name        = "pressure"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_pressure"),                    &
                               attribute_t("long_name",     "Pressure"),                        &
                               attribute_t("units",         "Pa"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Air Pressure on interfaces between mass levels
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%pressure_interface))
            var%name        = "pressure_i"
            var%dimensions  = three_d_t_interface_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_pressure"),                    &
                               attribute_t("long_name",     "Pressure"),                        &
                               attribute_t("units",         "Pa"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Air Pressure on interfaces between mass levels
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%surface_pressure))
            var%name        = "psfc"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_air_pressure"),            &
                               attribute_t("long_name",     "Surface Pressure"),                &
                               attribute_t("units",         "Pa"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Potential Air Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%potential_temperature))
            var%name        = "potential_temperature"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_potential_temperature"),       &
                               attribute_t("long_name",     "Potential Temperature"),           &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Real Air Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%temperature))
            var%name        = "temperature"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_temperature"),                 &
                               attribute_t("long_name",     "Temperature"),                     &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Water Vapor Mixing Ratio
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%water_vapor))
            var%name        = "qv"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "mass_fraction_of_water_vapor_in_air"), &
                               attribute_t("long_name",     "Water Vapor Mixing Ratio"),            &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Cloud water (liquid) mixing ratio
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%cloud_water))
            var%name        = "qc"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "cloud_liquid_water_mixing_ratio"),     &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Cloud water (liquid) number concentration
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%cloud_number_concentration))
            var%name        = "nc"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "number_concentration_of_cloud_droplets_in_air"), &
                               attribute_t("units",         "cm-3"),                                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Cloud ice mixing ratio
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%cloud_ice))
            var%name        = "qi"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "cloud_ice_mixing_ratio"),              &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Cloud ice number concentration
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ice_number_concentration))
            var%name        = "ni"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "number_concentration_of_ice_crystals_in_air"), &
                               attribute_t("units",         "cm-3"),                                            &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Rain water mixing ratio
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%rain_in_air))
            var%name        = "qr"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "mass_fraction_of_rain_in_air"),        &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Rain water number concentration
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%rain_number_concentration))
            var%name        = "nr"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "number_concentration_of_rain_particles_in_air"), &
                               attribute_t("units",         "cm-3"),                                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow in air mixing ratio
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_in_air))
            var%name        = "qs"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "mass_fraction_of_snow_in_air"),        &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow in air number concentration
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_number_concentration))
            var%name        = "ns"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "number_concentration_of_snow_particles_in_air"), &
                               attribute_t("units",         "cm-3"),                                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Graupel mixing ratio
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%graupel_in_air))
            var%name        = "qg"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "mass_fraction_of_graupel_in_air"),     &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Graupel number concentration
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%graupel_number_concentration))
            var%name        = "ng"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "number_concentration_of_graupel_particles_in_air"), &
                               attribute_t("units",         "cm-3"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Precipitation rate at the surface (requires tracking past precipitation amounts)
        !!------------------------------------------------------------
        ! associate(var=>var_meta(kVARS%precip_rate))
        !     var%name        = "precip_rate"
        !     var%dimensions  = two_d_t_dimensions
        !     var%unlimited_dim=.True.
        !     var%attributes  = [attribute_t("standard_name",   "precipitation_flux"),      &
        !                        attribute_t("units",           "kg m-2 s-1"),              &
        !                        attribute_t("coordinates",     "lat lon")]
        ! end associate
        !>------------------------------------------------------------
        !!  Accumulated precipitation at the surface
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%precipitation))
            var%name        = "precipitation"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "precipitation_amount"),                &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Accumulated convective precipitation at the surface
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%convective_precipitation))
            var%name        = "cu_precipitation"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "convective_precipitation_amount"),     &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Accumulated snowfall at the surface
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snowfall))
            var%name        = "snowfall"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "snowfall_amount"),                     &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Accumulated Graupel at the surface
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%graupel))
            var%name        = "graupel"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "graupel_amount"),                      &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Exner function
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%exner))
            var%name        = "exner"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "exner_function_result"),           &
                               attribute_t("units",         "K K-1"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Air Density
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%density))
            var%name        = "density"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_density"),                         &
                               attribute_t("units",         "kg m-3"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Vertical coordinate
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%z))
            var%name        = "z"
            var%dimensions  = three_d_dimensions
            var%attributes  = [attribute_t("standard_name", "height_above_reference_ellipsoid"),    &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Vertical coordinate on the interface between mass levels
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%z_interface))
            var%name        = "z_i"
            var%dimensions  = three_d_interface_dimensions
            var%attributes  = [attribute_t("standard_name", "height_above_reference_ellipsoid"),    &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Vertical layer thickness
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%dz))
            var%name        = "dz"
            var%dimensions  = three_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "layer_thickness"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Thickness of layers between interfaces
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%dz_interface))
            var%name        = "dz_i"
            var%dimensions  = three_d_interface_dimensions
            var%attributes  = [attribute_t("non_standard_name", "layer_thickness"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Cloud cover fraction
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%cloud_fraction))
            var%name        = "clt"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "cloud_area_fraction"),                 &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Effective cloud droplet radius
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%re_cloud))
            var%name        = "re_cloud"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "effective_radius_of_cloud_liquid_water_particles"), &
                               attribute_t("units",         "m"),                                                &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Effective cloud ice radius
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%re_ice))
            var%name        = "re_ice"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "effective_radius_of_stratiform_cloud_ice_particles"), &
                               attribute_t("units",         "m"),                                                  &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Effective snow radius
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%re_snow))
            var%name        = "re_snow"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "effective_radius_of_stratiform_snow_particles"), &
                               attribute_t("units",         "m"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Outgoing longwave radiation
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%out_longwave_rad))
            var%name        = "rlut"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "toa_outgoing_longwave_flux"), &
                               attribute_t("units",         "W m-2"),                      &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Longwave cloud forcing
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%longwave_cloud_forcing))
            var%name        = "lwcf"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "longwave_cloud_forcing"), &
                               attribute_t("units",         "W m-2"),                      &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Shortwave cloud forcing
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%shortwave_cloud_forcing))
            var%name        = "swcf"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "shortwave_cloud_forcing"), &
                               attribute_t("units",         "W m-2"),                       &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Cosine solar zenith angle
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%cosine_zenith_angle))
            var%name        = "cosz"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "cosine_zenith_angle"), &
                               attribute_t("units",         " "),                       &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Tendency from short wave radiation
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%tend_swrad))
            var%name        = "tend_swrad"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "sw_rad_tend"), &
                               attribute_t("units",         " "),               &
                               attribute_t("coordinates",   "lat lon")]
        end associate


        !>------------------------------------------------------------
        !!  Surface emissivity
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%land_emissivity))
            var%name        = "emiss"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_longwave_emissivity"), &
                               attribute_t("units",         " "),                           &
                               attribute_t("coordinates",   "lat lon")]
        end associate

        !>------------------------------------------------------------
        !!  Temperature on interface
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%temperature_interface))
            var%name        = "temperature_i"
            var%dimensions  = three_d_t_interface_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_temperature"),                    &
                               attribute_t("long_name",     "Temperature"),                        &
                               attribute_t("units",         "K"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        end associate



        !>------------------------------------------------------------
        !!  Downward Shortwave Radiation at the Surface (positive down)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%shortwave))
            var%name        = "rsds"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_downwelling_shortwave_flux_in_air"), &
                               attribute_t("units",         "W m-2"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Downward Direct Shortwave Radiation at the Surface (positive down)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%shortwave_direct))
            var%name        = "shortwave_direct"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_direct_downwelling_shortwave_flux_in_air"), &
                               attribute_t("units",         "W m-2"),                                            &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Downward Diffuse Shortwave Radiation at the Surface (positive down)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%shortwave_diffuse))
            var%name        = "shortwave_diffuse"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_diffuse_downwelling_shortwave_flux_in_air"), &
                               attribute_t("units",         "W m-2"),                                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Downward Longwave Radiation at the Surface (positive down)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%longwave))
            var%name        = "rlds"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_downwelling_longwave_flux_in_air"), &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Total Absorbed Solar Radiation
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%rad_absorbed_total))
            var%name        = "rad_absorbed_total"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "total_absorbed_radiation"),             &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Solar Radiation Absorbed by Vegetation
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%rad_absorbed_veg))
            var%name        = "rad_absorbed_veg"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "radiation_absorbed_by_vegetation"),     &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Solar Radiation Absorbed by Bare Ground
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%rad_absorbed_bare))
            var%name        = "rad_absorbed_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "radiation_absorbed_by_bare_ground"),    &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Net Longwave Radiation (positive up)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%rad_net_longwave))
            var%name        = "rad_net_longwave"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "net_upward_longwave_flux_in_air"),          &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Upward Longwave Radiation at the Surface (positive up)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%longwave_up))
            var%name        = "rlus"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_upwelling_longwave_flux_in_air"),   &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Ground Heat Flux (positive down)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ground_heat_flux))
            var%name        = "hfgs"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "upward_heat_flux_at_ground_level_in_soil"), &
                               attribute_t("units",         "W m-2"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Vegetation fraction
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%vegetation_fraction))
            var%name        = "vegetation_fraction"
            var%dimensions  = two_d_month_dimensions
            var%attributes  = [attribute_t("standard_name", "vegetation_area_fraction"),            &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Annual Maximum Vegetation Fraction
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%vegetation_fraction_max))
            var%name        = "vegetation_fraction_max"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "max_vegetation_area_fraction"),    &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Noah-MP Output Vegetation Fraction
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%vegetation_fraction_out))
            var%name        = "vegetation_fraction_out"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "vegetation_fraction_out"),         &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Land cover type
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%veg_type))
            var%name        = "veg_type"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "vegetation_type"),                 &
                               attribute_t("units",      "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Leaf Mass
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%mass_leaf))
            var%name        = "mass_leaf"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "leaf_mass"),                      &
                               attribute_t("units",      "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Root Mass
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%mass_root))
            var%name        = "mass_root"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "root_mass"),                      &
                               attribute_t("units",      "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Stem Mass
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%mass_stem))
            var%name        = "mass_stem"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "stem_mass"),                      &
                               attribute_t("units",      "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Wood Mass
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%mass_wood))
            var%name        = "mass_wood"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "wood_mass"),                      &
                               attribute_t("units",      "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Leaf Area Index
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%lai))
            var%name        = "lai"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "leaf_area_index"),                     &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Stem Area Index
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%sai))
            var%name        = "sai"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "stem_area_index"),                 &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Planting Date
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%date_planting))
            var%name        = "date_planting"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "planting_date"),                   &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Harvesting Date
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%date_harvest))
            var%name        = "date_harvest"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "harvest_date"),                    &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Crop Category
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%crop_category))
            var%name        = "crop_category"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "crop_category"),                   &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Crop Type (Noah-MP initialization variable)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%crop_type))
            var%name        = "crop_type"
            var%dimensions  = three_d_crop_dimensions
            var%attributes  = [attribute_t("non_standard_name", "crop_type"),                       &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Growing Season Growing Degree Days
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%growing_season_gdd))
            var%name        = "growing_season_gdd"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "growing_season_gdd"),              &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Event Number, Sprinkler
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_eventno_sprinkler))
            var%name        = "irr_eventno_sprinkler"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_eventno_sprinkler"),           &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Fraction
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_frac_total))
            var%name        = "irr_frac_total"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_frac_total"),                  &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Fraction, Sprinkler
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_frac_sprinkler))
            var%name        = "irr_frac_sprinkler"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_frac_sprinkler"),              &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Fraction, Micro
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_frac_micro))
            var%name        = "irr_frac_micro"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_frac_micro"),                  &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Fraction, Flood
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_frac_flood))
            var%name        = "irr_frac_flood"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_frac_flood"),                  &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Event Number, Micro
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_eventno_micro))
            var%name        = "irr_eventno_micro"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_eventno_micro"),               &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Event Number, Flood
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_eventno_flood))
            var%name        = "irr_eventno_flood"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_eventno_flood"),               &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Water Amount to be Applied, Sprinkler
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_alloc_sprinkler))
            var%name        = "irr_alloc_sprinkler"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_alloc_sprinkler"),             &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Water Amount to be Applied, Micro
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_alloc_micro))
            var%name        = "irr_alloc_micro"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_alloc_micro"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Irrigation Water Amount to be Applied, Flood
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_alloc_flood))
            var%name        = "irr_alloc_flood"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_alloc_flood"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Loss of Sprinkler Irrigation Water to Evaporation
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_evap_loss_sprinkler))
            var%name        = "irr_evap_loss_sprinkler"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_evap_loss_sprinkler"),         &
                               attribute_t("units",         "m timestep-1"),                        &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Amount of Irrigation, Sprinkler
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_amt_sprinkler))
            var%name        = "irr_amt_sprinkler"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_amt_sprinkler"),               &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Amount of Irrigation, Micro
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_amt_micro))
            var%name        = "irr_amt_micro"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_amt_micro"),                   &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Amount of Irrigation, Flood
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%irr_amt_flood))
            var%name        = "irr_amt_flood"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "irr_amt_flood"),                   &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Growing Season Growing Degree Days
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%growing_season_gdd))
            var%name        = "growing_season_gdd"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "growing_season_gdd"),              &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Mass of Agricultural Grain Produced
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%mass_ag_grain))
            var%name        = "mass_ag_grain"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "mass_agricultural_grain"),         &
                               attribute_t("units",         "g m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Growing Degree Days (based on 10C)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%growing_degree_days))
            var%name        = "growing_degree_days"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "growing_degree_days"),             &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Plant Growth Stage
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%plant_growth_stage))
            var%name        = "plant_growth_stage"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "plant_growth_stage"),              &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Net Ecosystem Exchange
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%net_ecosystem_exchange))
            var%name        = "net_ecosystem_exchange"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "net_ecosystem_exchange_expressed_as_carbon_dioxide"), &
                               attribute_t("units",         "g m-2 s-1"),                                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Gross Primary Productivity
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%gross_primary_prod))
            var%name        = "gross_primary_prod"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "gross_primary_productivity_of_biomass_expressed_as_carbon"), &
                               attribute_t("units",         "g m-2 s-1"),                                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Net Primary Productivity
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%net_primary_prod))
            var%name        = "net_primary_prod"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "net_primary_productivity_of_biomass_expressed_as_carbon"), &
                               attribute_t("units",         "g m-2 s-1"),                                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Absorbed Photosynthetically Active Radiation
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%apar))
            var%name        = "apar"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "absorbed_photosynthetically_active_radiation"), &
                               attribute_t("units",         "W m-2"),                                            &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Total Photosynthesis
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%photosynthesis_total))
            var%name        = "photosynthesis_total"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "total_photosynthesis_expressed_as_carbon_dioxide"), &
                               attribute_t("units",         "mmol m-2 s-1"),                                         &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Total Leaf Stomatal Resistance
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%stomatal_resist_total))
            var%name        = "stomatal_resist_total"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "total_leaf_stomatal_resistance"),                   &
                               attribute_t("units",         "s m-1"),                                                &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sunlit Leaf Stomatal Resistance
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%stomatal_resist_sun))
            var%name        = "stomatal_resist_sun"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "sunlif_leaf_stomatal_resistance"),                  &
                               attribute_t("units",         "s m-1"),                                                &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Shaded Leaf Stomatal Resistance
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%stomatal_resist_shade))
            var%name        = "stomatal_resist_shade"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "shaded_leaf_stomatal_resistance"),                 &
                               attribute_t("units",         "s m-1"),                                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  GECROS (Genotype-by-Envrionment interaction on CROp growth Simulator [Yin and van Laar, 2005]) crop model state
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%gecros_state))
            var%name        = "gecros_state"
            var%dimensions  = three_d_t_gecros_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "gecros_state"),                                    &
                               attribute_t("units",         "N/A"),                                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Canopy Total Water Content
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%canopy_water))
            var%name        = "canopy_water"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "canopy_water_amount"),                 &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Canopy Frozen Water Content
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%canopy_water_ice))
            var%name        = "canopy_ice"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "canopy_snow_amount"),                  &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Canopy Liquid Water Content (in snow)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%canopy_water_liquid))
            var%name        = "canopy_liquid"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "canopy_liquid_water_amount"),      &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Canopy Air Vapor Pressure
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%canopy_vapor_pressure))
            var%name        = "canopy_vapor_pressure"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "canopy_air_vapor_pressure"),       &
                               attribute_t("units",         "Pa"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Canopy Air Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%canopy_temperature))
            var%name        = "canopy_temperature"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "canopy_air_temperature"),          &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Canopy Wetted/Snowed Fraction
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%canopy_fwet))
            var%name        = "canopy_fwet"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "canopy_wetted_fraction"),          &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Vegetation Leaf Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%veg_leaf_temperature))
            var%name        = "veg_leaf_temperature"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "veg_leaf_temperature"),            &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Ground Surface Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ground_surf_temperature))
            var%name        = "ground_surf_temperature"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ground_surface_temperature"),      &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Between Gap Fraction
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%frac_between_gap))
            var%name        = "frac_between_gap"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "between_gap_fraction"),            &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Within Gap Fraction
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%frac_within_gap))
            var%name        = "frac_within_gap"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "within_gap_fraction"),             &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Under-Canopy Ground Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ground_temperature_canopy))
            var%name        = "ground_temperature_canopy"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "under_canopy_ground_temperature"), &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Bare Ground Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ground_temperature_bare))
            var%name        = "ground_temperature_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "bare_ground_temperature"),         &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snowfall on the Ground
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snowfall_ground))
            var%name        = "snowfall_ground"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ground_snow_rate"),                &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Rainfall on the Ground
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%rainfall_ground))
            var%name        = "rainfall_ground"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ground_rain_rate"),                &
                               attribute_t("units",         "mm s-1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow water equivalent on the surface
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_water_equivalent))
            var%name        = "swe"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_snow_amount"),                 &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow water equivalent from previous timestep
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_water_eq_prev))
            var%name        = "swe_0"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "surface_snow_amount_prev"),        &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow albedo from previous timestep
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_albedo_prev))
            var%name        = "snow_albedo_0"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "snowpack_albedo_prev"),            &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_temperature))
            var%name        = "snow_temperature"
            var%dimensions  = three_d_t_snow_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "temperature_in_surface_snow"),         &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow Layer Depth
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_layer_depth))
            var%name        = "snow_layer_depth"
            var%dimensions  = three_d_t_snowsoil_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "snow_layer_depth"),                &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow Layer Ice
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_layer_ice))
            var%name        = "snow_layer_ice"
            var%dimensions  = three_d_t_snow_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "snow_layer_ice_content"),          &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow Layer Liquid Water
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_layer_liquid_water))
            var%name        = "snow_layer_liquid_water"
            var%dimensions  = three_d_t_snow_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "snow_layer_liquid_water_content"), &
                               attribute_t("units",         "mm"),                                  &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow Age Factor
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_age_factor))
            var%name        = "tau_ss"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "snow_age_factor"),                 &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Snow height
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_height))
            var%name        = "snow_height"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_snow_height"),                 &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Number of Snowpack Layers
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_nlayers))
            var%name        = "snow_nlayers"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "snow_nlayers"),                    &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil water content
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_water_content))
            var%name        = "soil_water_content"
            var%dimensions  = three_d_t_soil_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "moisture_content_of_soil_layer"),      &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Equilibrium Volumetric Soil Moisture
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%eq_soil_moisture))
            var%name        = "eq_soil_moisture"
            var%dimensions  = three_d_t_soil_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "equilibrium_volumetric_soil_moisture"), &
                               attribute_t("units",         "m3 m-3"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil Moisture Content in the Layer Draining to Water Table when Deep
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%smc_watertable_deep))
            var%name        = "smc_watertable_deep"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "soil_moisture_content_in_layer_to_water_table_when_deep"), &
                               attribute_t("units",         "m3 m-3"),                                                      & !units not defined in noahmpdrv (guess)
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Groundwater Recharge
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%recharge))
            var%name        = "recharge"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "groundwater_recharge"),            &
                               attribute_t("units",         "mm"),                                  & !units not defined in noahmpdrv (guess)
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Groundwater Recharge when Water Table is Deep
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%recharge_deep))
            var%name        = "recharge_deep"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "groundwater_recharge_deep"),           &
                               attribute_t("units",         "mm"),                                  & !units not defined in noahmpdrv (guess)
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Canopy Evaporation Rate
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%evap_canopy))
            var%name        = "evap_canopy"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "water_evaporation_flux_from_canopy"),  &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil Surface Evaporation Rate
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%evap_soil_surface))
            var%name        = "evap_soil_surface"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "water_evaporation_flux_from_soil"),    &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Transpiration Rate
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%transpiration_rate))
            var%name        = "transpiration_rate"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "transpiration_rate"),              &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, Vegetated
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ch_veg))
            var%name        = "ch_veg"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ch_vegetated"),                    &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient at 2m, Vegetated
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ch_veg_2m))
            var%name        = "ch_veg_2m"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ch_vegetated_2m"),                 &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, Bare
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ch_bare))
            var%name        = "ch_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ch_bare_ground"),                  &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient at 2m, Bare
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ch_bare_2m))
            var%name        = "ch_bare_2m"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ch_bare_2m"),                      &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, Under Canopy
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ch_under_canopy))
            var%name        = "ch_under_canopy"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ch_under_canopy"),                 &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient, Leaf
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ch_leaf))
            var%name        = "ch_leaf"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ch_leaf"),                         &
                               attribute_t("units",         "1"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Flux, Vegetated Ground (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%sensible_heat_veg))
            var%name        = "sensible_heat_veg"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "sensible_heat_veg"),               &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Flux, Bare Ground (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%sensible_heat_bare))
            var%name        = "sensible_heat_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "sensible_heat_bare"),              &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Flux, Canopy (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%sensible_heat_canopy))
            var%name        = "sensible_heat_canopy"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "sensible_heat_canopy"),            &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Evaporation Heat Flux, Vegetated Ground (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%evap_heat_veg))
            var%name        = "evap_heat_veg"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "evap_heat_veg"),                   &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Evaporation Heat Flux, Bare Ground (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%evap_heat_bare))
            var%name        = "evap_heat_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "evap_heat_bare"),                  &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !! Evaporation Heat Flux, Canopy (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%evap_heat_canopy))
            var%name        = "evap_heat_canopy"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "evap_heat_canopy"),                &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Transpiration Heat Flux (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%transpiration_heat))
            var%name        = "transpiration_heat"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "transpiration_heat"),              &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Ground Heat Flux, Vegetated Ground (+ to soil)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ground_heat_veg))
            var%name        = "ground_heat_veg"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ground_heat_veg"),                 &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Ground Heat Flux, Bare Ground (+ to soil)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%ground_heat_bare))
            var%name        = "ground_heat_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "ground_heat_bare"),                &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Net Longwave Radiation Flux, Vegetated Ground (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%net_longwave_veg))
            var%name        = "net_longwave_veg"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "net_longwave_veg"),                &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Net Longwave Radiation Flux, Bare Ground (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%net_longwave_bare))
            var%name        = "net_longwave_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "net_longwave_bare"),               &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Net Longwave Radiation Flux, Canopy (+ to atmosphere)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%net_longwave_canopy))
            var%name        = "net_longwave_canopy"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "net_longwave_canopy"),             &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Surface Runoff Rate
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%runoff_surface))
            var%name        = "runoff_surface"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_runoff_flux"),                 &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Subsurface Runoff Rate
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%runoff_subsurface))
            var%name        = "runoff_subsurface"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "subsurface_runoff_flux"),          &
                               attribute_t("units",         "mm s-1"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Total Column Soil water content
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_totalmoisture))
            var%name        = "soil_column_total_water"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "soil_moisture_content"),               &
                               attribute_t("units",         "kg m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_temperature))
            var%name        = "soil_temperature"
            var%dimensions  = three_d_t_soil_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "soil_temperature"),                    &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Deep Soil Temperature (time constant)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_deep_temperature))
            var%name        = "soil_deep_temperature"
            var%dimensions  = two_d_t_dimensions
            var%attributes  = [attribute_t("non_standard_name", "deep_soil_temperature"),           &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Stable Carbon Mass in Deep Soil
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_carbon_stable))
            var%name        = "soil_carbon_stable"
            var%dimensions  = two_d_t_dimensions
            var%attributes  = [attribute_t("standard_name", "slow_soil_pool_mass_content_of_carbon"), &
                               attribute_t("units",         "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Short-lived Carbon Mass in Shallow Soil
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_carbon_fast))
            var%name        = "soil_carbon_fast"
            var%dimensions  = two_d_t_dimensions
            var%attributes  = [attribute_t("standard_name", "fast_soil_pool_mass_content_of_carbon"), &
                               attribute_t("units",         "g m-2"),                                 &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil Class, Layer 1
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_texture_1))
            var%name        = "soil_class_1"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "soil_class_layer1"),                 &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil Class, Layer 2
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_texture_2))
            var%name        = "soil_class_2"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "soil_class_layer2"),                 &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil Class, Layer 3
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_texture_3))
            var%name        = "soil_class_3"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "soil_class_layer3"),                 &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil Class, Layer 4
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_texture_4))
            var%name        = "soil_class_4"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "soil_class_layer4"),                 &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Soil Sand and Clay Composition by Layer
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%soil_sand_and_clay))
            var%name        = "soil_sand_and_clay_composition"
            var%dimensions  = three_d_soilcomp_dimensions
            var%attributes  = [attribute_t("non_standard_name", "soil_sand_and_clay_composition"),    &
                               attribute_t("units",         "1"),                                     &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Water Table Depth
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%water_table_depth))
            var%name        = "water_table_depth"
            var%dimensions  = two_d_t_dimensions
            var%attributes  = [attribute_t("standard_name", "water_table_depth"),                    &
                               attribute_t("units",         "m"),                                    &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Water in Aquifer
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%water_aquifer))
            var%name        = "water_aquifer"
            var%dimensions  = two_d_t_dimensions
            var%attributes  = [attribute_t("non_standard_name", "water_in_aquifer"),                 &
                               attribute_t("units",         "mm"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Groundwater Storage
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%storage_gw))
            var%name        = "storage_gw"
            var%dimensions  = two_d_t_dimensions
            var%attributes  = [attribute_t("non_standard_name", "groundwater_storage"),              &
                               attribute_t("units",         "mm"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake Storage
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%storage_lake))
            var%name        = "storage_lake"
            var%dimensions  = two_d_t_dimensions
            var%attributes  = [attribute_t("non_standard_name", "lake_storage"),                     &
                               attribute_t("units",         "mm"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Surface roughness length z0
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%roughness_z0))
            var%name        = "surface_roughness"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_roughness_length"),            &
                               attribute_t("long_name",     "Surface roughness length"),            &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Surface Radiative Temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%surface_rad_temperature))
            var%name        = "surface_rad_temperature"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_radiative_temperature"),       &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  2 meter air temperture
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%temperature_2m))
            var%name        = "ta2m"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_temperature"),                     &
                               attribute_t("long_name",     "Bulk air temperature at 2m"),          &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  2 meter Air Temperature over Vegetation
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%temperature_2m_veg))
            var%name        = "temperature_2m_veg"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_temperature"),                     &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  2 meter Air Temperature over Bare Ground
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%temperature_2m_bare))
            var%name        = "temperature_2m_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "air_temperature"),                     &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  2 meter Mixing Ratio over Vegetation
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%mixing_ratio_2m_veg))
            var%name        = "mixing_ratio_2m_veg"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "mixing_ratio"),                    &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  2 meter Mixing Ratio over Bare Ground
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%mixing_ratio_2m_bare))
            var%name        = "mixing_ratio_2m_bare"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "mixing_ratio"),                    &
                               attribute_t("units",         "kg kg-1"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  2 meter specific humidity
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%humidity_2m))
            var%name        = "hus2m"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "specific_humidity"),                   &
                               attribute_t("units",         "kg kg-2"),                             &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  10 meter height V component of wind field
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%v_10m))
            var%name        = "v10m"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "northward_10m_wind_speed"),            &
                               attribute_t("units",         "m s-1"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  10 meter height U component of the wind field
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%u_10m))
            var%name        = "u10m"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "eastward_10m_wind_speed"),             &
                               attribute_t("units",         "m s-1"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Momentum Drag Coefficient
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%coeff_momentum_drag))
            var%name        = "coeff_momentum_drag"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_drag_coefficient_for_momentum_in_air"), &
                               attribute_t("units",         "1"),                                            &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%coeff_heat_exchange))
            var%name        = "coeff_heat_exchange"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "sensible_heat_exchange_coefficient"), &
                               attribute_t("units",         "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible Heat Exchange Coefficient 3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%coeff_heat_exchange_3d))
            var%name        = "coeff_heat_exchange_3d"
            var%dimensions  = three_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "sensible_heat_exchange_coefficient_3d"), &
                               attribute_t("units",         "1"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  PBL height
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%hpbl))
            var%name        = "hpbl"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "height_of_planetary_boundary_layer"), &
                               attribute_t("units",         "m"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  PBL layer index
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%kpbl))
            var%name        = "kpbl"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "index_of_planetary_boundary_layer_height"), &
                               attribute_t("units",         "-"),                                      &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Land surface radiative skin temperature
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%skin_temperature))
            var%name        = "ts"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_temperature"),                 &
                               attribute_t("units",         "K"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Sensible heat flux from the surface (positive up)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%sensible_heat))
            var%name        = "hfss"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_upward_sensible_heat_flux"),   &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Latent heat flux from the surface (positive up)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%latent_heat))
            var%name        = "hfls"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_upward_latent_heat_flux"),     &
                               attribute_t("units",         "W m-2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake temperature 3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%t_lake3d))
            var%name        = "t_lake3d"
            var%dimensions  = three_d_t_lake_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "lake_water_temperature"),     &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake lake_icefraction_3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%lake_icefrac3d))
            var%name        = "lake_icefrac3d"
            var%dimensions  = three_d_t_lake_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "lake_icefraction_3d"),     &
                               attribute_t("units",         "-"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake z_lake3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%z_lake3d))
            var%name        = "z_lake3d"
            var%dimensions  = three_d_t_lake_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "lake_layer_depth"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake dz_lake3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%dz_lake3d))
            var%name        = "dz_lake3d"
            var%dimensions  = three_d_t_lake_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "lake_layer_thickness"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  lake snl2d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snl2d))
            var%name        = "snl2d"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "lake_snow_layer_2d"),           &
                               attribute_t("units",         "-"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate 
        !>------------------------------------------------------------
        !!  lake_t_grnd2d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%t_grnd2d))
            var%name        = "t_grnd2d"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "t_grnd2d"),           &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake t_soisno3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%t_soisno3d))
            var%name        = "t_soisno3d"
            var%dimensions  = three_d_t_lake_soisno_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "temperature_soil_snow_below_or_above_lake"),     &
                               attribute_t("units",         "K"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake h2osoi_ice3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%h2osoi_ice3d))
            var%name        = "h2osoi_ice3d"
            var%dimensions  = three_d_t_lake_soisno_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "h2osoi_ice3d"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake soil/snowliquid water (kg/m2)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%h2osoi_liq3d))
            var%name        = "h2osoi_liq3d"
            var%dimensions  = three_d_t_lake_soisno_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "lake_soil_or_snow_liquid water_content"),     &
                               attribute_t("units",         "kg/m2"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake h2osoi_vol3d volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%h2osoi_vol3d))
            var%name        = "h2osoi_vol3d"
            var%dimensions  = three_d_t_lake_soisno_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "volumetric_soil_water"),     &
                               attribute_t("units",         "m3/m3"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake z3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%z3d))
            var%name        = "z3d"
            var%dimensions  = three_d_t_lake_soisno_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "layer_depth_for_lake_snow&soil"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake layer_thickness_for_lake_snow&soil
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%dz3d))
            var%name        = "dz3d"
            var%dimensions  = three_d_t_lake_soisno_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "layer_thickness_for_lake_snow&soil"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake z3d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%zi3d))
            var%name        = "zi3d"
            var%dimensions  = three_d_t_lake_soisno_1_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "interface_layer_depth_for_lake_snow&soil"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake watsat3d: volumetric soil water at saturation (porosity)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%watsat3d))
            var%name        = "watsat3d"
            var%dimensions  = three_d_t_lake_soi_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "volumetric soil water at saturation (porosity)"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake csol3d: heat capacity, soil solids (J/m**3/Kelvin)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%csol3d))
            var%name        = "csol3d"
            var%dimensions  = three_d_t_lake_soi_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "heat capacity, soil solids "),     &
                               attribute_t("units",         "(J/m**3/Kelvin)"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake: thermal conductivity, soil minerals  [W/m-K]
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%tkmg3d))
            var%name        = "tkmg3d"
            var%dimensions  = three_d_t_lake_soi_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "thermal conductivity, soil minerals  [W/m-K]"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake lakemask
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%lakemask))
            var%name        = "lakemask"
            var%dimensions  = two_d_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "lakemask"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake lakedepth2d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%lakedepth2d))
            var%name        = "lakedepth2d"
            var%dimensions  = two_d_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "lake_depth"),     &
                               attribute_t("units",         "m"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  lake savedtke12d
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%savedtke12d))
            var%name        = "savedtke12d"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "savedtke12d"),           &
                               attribute_t("units",         "-?"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate 
        !>------------------------------------------------------------
        !!  Lake: thermal conductivity, saturated soil [W/m-K]
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%tksatu3d))
            var%name        = "tksatu3d"
            var%dimensions  = three_d_t_lake_soi_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "thermal conductivity, saturated soil [W/m-K]"),     &
                               attribute_t("units",         ""),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Lake tkdry3d: thermal conductivity, dry soil (W/m/Kelvin)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%tkdry3d))
            var%name        = "tkdry3d"
            var%dimensions  = three_d_t_lake_soi_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "thermal conductivity, dry soil (W/m/Kelvin)"),     &
                               attribute_t("units",         "?"),                               &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        




        
        ! type(variable_t) :: h2osoi_ice3d
        ! type(variable_t) :: h2osoi_liq3d! liquid water (kg/m2)
        ! type(variable_t) :: h2osoi_vol3d! volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
        ! type(variable_t) :: z3d ! layer depth for snow & soil (m)
        ! type(variable_t) :: dz3d

        ! type(variable_t) :: watsat3d
        ! type(variable_t) :: csol3d
        ! type(variable_t) :: tkmg3d
        ! type(variable_t) :: lakemask
        ! type(variable_t) :: tksatu3d
        ! type(variable_t) :: tkdry3d
        ! type(variable_t) :: zi3d
        !>------------------------------------------------------------
        !!  Binary land mask (water vs land)
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%land_mask))
            var%name        = "land_mask"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("non_standard_name", "land_water_mask"),                 &
                               attribute_t("coordinates",       "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Height of the terrain
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%terrain))
            var%name        = "terrain"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("standard_name", "height_above_reference_ellipsoid"),    &
                               attribute_t("units",         "m"),                                   &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Latitude y coordinate
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%latitude))
            var%name        = "lat"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("standard_name", "latitude"),                            &
                               attribute_t("units",         "degrees_north"),                       &
                               attribute_t("axis","Y")]
        end associate
        !>------------------------------------------------------------
        !!  Longitude x coordinate
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%longitude))
            var%name        = "lon"
            var%dimensions  = two_d_dimensions
            var%attributes  = [attribute_t("standard_name", "longitude"),                           &
                               attribute_t("units",         "degrees_east"),                        &
                               attribute_t("axis","X")]
        end associate
        !>------------------------------------------------------------
        !!  Latitude y coordinate on the U-grid
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%u_latitude))
            var%name        = "u_lat"
            var%dimensions  = two_d_u_dimensions
            var%attributes  = [attribute_t("non_standard_name", "latitude_on_u_grid"),              &
                               attribute_t("units",         "degrees_north")]
        end associate
        !>------------------------------------------------------------
        !!  Longitude x coordinate on the U-grid
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%u_longitude))
            var%name        = "u_lon"
            var%dimensions  = two_d_u_dimensions
            var%attributes  = [attribute_t("non_standard_name", "longitude_on_u_grid"),             &
                               attribute_t("units",         "degrees_east")]
        end associate
        !>------------------------------------------------------------
        !!  Latitude y coordinate on the V-grid
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%v_latitude))
            var%name        = "v_lat"
            var%dimensions  = two_d_v_dimensions
            var%attributes  = [attribute_t("non_standard_name", "latitude_on_v_grid"),              &
                               attribute_t("units",         "degrees_north")]
        end associate
        !>------------------------------------------------------------
        !!  Longitude x coordinate on the V-grid
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%v_longitude))
            var%name        = "v_lon"
            var%dimensions  = two_d_v_dimensions
            var%attributes  = [attribute_t("non_standard_name", "longitude_on_v_grid"),             &
                               attribute_t("units",         "degrees_east")]
        end associate

        ! loop through entire array setting n_dimensions and n_attrs based on the data that were supplied
        do i=1,size(var_meta)
            var_meta(i)%n_dimensions = size(var_meta(i)%dimensions)
            var_meta(i)%n_attrs      = size(var_meta(i)%attributes)
        enddo

    end subroutine init_var_meta

end module output_metadata
