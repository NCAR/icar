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
        character(len=16) :: three_d_soil_dimensions(3)         = [character(len=16) :: "lon_x","lat_y","nsoil"]
        character(len=16) :: three_d_t_soil_dimensions(4)       = [character(len=16) :: "lon_x","lat_y","nsoil","time"]

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
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "vegetation_fraction"),                 &
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
        !!  Leaf Area Index
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%lai))
            var%name        = "lai"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("non_standard_name", "leaf_area_index"),                 &
                               attribute_t("units",         "m2 m-2"),                              &
                               attribute_t("coordinates",   "lat lon")]
        end associate
        !>------------------------------------------------------------
        !!  Canopy Water Content
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
        !!  Snow height
        !!------------------------------------------------------------
        associate(var=>var_meta(kVARS%snow_height))
            var%name        = "snow_height"
            var%dimensions  = two_d_t_dimensions
            var%unlimited_dim=.True.
            var%attributes  = [attribute_t("standard_name", "surface_snow_height"),                 &
                               attribute_t("units",         "m"),                              &
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
