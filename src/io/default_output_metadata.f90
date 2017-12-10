module output_metadata

    use icar_constants
    implicit none

    type(variable_t), allocatable :: var_meta(:)

    !>------------------------------------------------------------
    !! Generic interface to the netcdf read routines
    !!------------------------------------------------------------
    interface get_metadata
        module procedure get_metadata_2d, get_metadata_3d
    end interface


contains

    function get_metadata_2d(var_idx, input_data) result(meta_data)
        implicit none
        integer, intent(in) :: var_idx
        real,    intent(in), pointer, optional :: input_data(:,:)
        type(variable_t) :: meta_data

        if (var_idx>kMAX_STORAGE_VARS) then
            stop "Invalid variable metadata requested"
        endif

        if (.not.allocated(var_meta)) call init_var_meta()

        meta_data = var_meta(var_idx)

        if (present(input_data)) then
            meta_data%data_2d => input_data
            meta_data%data_2d = .True.
            allocate(meta_data%dim_len, source=shape(input_data))
        endif

    end function get_metadata_2d

    function get_metadata_3d(var_idx, input_data) result(meta_data)
        implicit none
        integer, intent(in) :: var_idx
        real,    intent(in), pointer, optional :: input_data(:,:,:)
        type(variable_t) :: meta_data

        if (var_idx>kMAX_STORAGE_VARS) then
            stop "Invalid variable metadata requested"
        endif

        if (.not.allocated(var_meta)) call init_var_meta()

        meta_data = var_meta(var_idx)

        if (present(input_data)) then
            meta_data%data_3d => input_data
            meta_data%data_3d = .True.
            allocate(meta_data%dim_len, source=shape(input_data))
        endif

    end function get_metadata_3d


    subroutine init_var_meta()
        implicit none

        if (allocated(var_meta)) deallocate(var_meta)

        allocate(var_meta(kMAX_STORAGE_VARS))

        ! need to specify : name, n_attrs, attributes%name %value
        ! n_dimensions, dimensions

        associate(var=>var_meta(kVARS%u))
            var%name        = "u"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","east"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%v))
            var%name        = "v"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%w))
            var%name        = "w"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%pressure))
            var%name        = "pressure"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%pressure_interface))
            var%name        = "pressure_interface"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%potential_temperature))
            var%name        = "potential_temperature"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%temperature))
            var%name        = "temperature"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%water_vapor))
            var%name        = "water_vapor"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%cloud_water))
            var%name        = "cloud_water"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%cloud_number_concentration))
            var%name        = "cloud_number_concentration"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%cloud_ice))
            var%name        = "cloud_ice"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%ice_number_concentration))
            var%name        = "ice_number_concentration"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%rain_in_air))
            var%name        = "rain_in_air"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%rain_number_concentration))
            var%name        = "rain_number_concentration"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%snow_in_air))
            var%name        = "snow_in_air"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%snow_number_concentration))
            var%name        = "snow_number_concentration"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%graupel_in_air))
            var%name        = "graupel_in_air"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%graupel_number_concentration))
            var%name        = "graupel_number_concentration"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%precipitation))
            var%name        = "precipitation"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%snowfall))
            var%name        = "snowfall"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%graupel))
            var%name        = "graupel"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%exner))
            var%name        = "exner"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%density))
            var%name        = "density"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%z))
            var%name        = "z"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%z_interface))
            var%name        = "z_interface"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%dz))
            var%name        = "dz"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%dz_interface))
            var%name        = "dz_interface"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%shortwave))
            var%name        = "shortwave"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%longwave))
            var%name        = "longwave"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%vegetation_fraction))
            var%name        = "vegetation_fraction"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","vegetation_fraction"), &
                               attribute_t("units", "m2 m-2"),                     &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%land_cover))
            var%name        = "land_cover"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("non_standard_name","land_cover_type"), &
                               attribute_t("units", ""),                           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%lai))
            var%name        = "lai"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("non_standard_name","leaf_area_index"), &
                               attribute_t("units", "m2 m-2"),                     &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%canopy_water))
            var%name        = "canopy_water"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","canopy_water_amount"),  &
                               attribute_t("units", "kg m-2"),                      &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%snow_water_equivalent))
            var%name        = "snow_water_equivalent"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","surface_snow_amount"), &
                               attribute_t("units", "kg m-2"),                     &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%soil_water_content))
            var%name        = "soil_water_content"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","nsoil"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","soil_moisture_content"), &
                               attribute_t("units", "kg m-2"),                       &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%soil_temperature))
            var%name        = "soil_temperature"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","nsoil"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","soil_temperature"), &
                               attribute_t("units", "K"),                       &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%skin_temperature))
            var%name        = "skin_temperature"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%land_mask))
            var%name        = "land_mask"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","unassigned"), &
                               attribute_t("units", "unknown"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%terrain))
            var%name        = "terrain"
            var%n_dimensions= 3
            var%dimensions  = ["lon_x","lat_x","level"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","height_above_reference_ellipsoid"), &
                               attribute_t("units", "m"),           &
                               attribute_t("coordinates","lat lon")]
        end associate
        associate(var=>var_meta(kVARS%latitude))
            var%name        = "latitude"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","latitude"), &
                               attribute_t("units", "degrees_north"),           &
                               attribute_t("axis","Y")]
        end associate
        associate(var=>var_meta(kVARS%longitude))
            var%name        = "longitude"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 3
            var%attributes  = [attribute_t("standard_name","longitude"), &
                               attribute_t("units", "degrees_east"),     &
                               attribute_t("axis","X")]
        end associate
        associate(var=>var_meta(kVARS%u_latitude))
            var%name        = "u_latitude"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 2
            var%attributes  = [attribute_t("non_standard_name","latitude_on_u_grid"), &
                               attribute_t("units", "degrees_north")]
        end associate
        associate(var=>var_meta(kVARS%u_longitude))
            var%name        = "u_longitude"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 2
            var%attributes  = [attribute_t("non_standard_name","longitude_on_u_grid"), &
                               attribute_t("units", "degrees_east")]
        end associate
        associate(var=>var_meta(kVARS%v_latitude))
            var%name        = "v_latitude"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 2
            var%attributes  = [attribute_t("non_standard_name","latitude_on_v_grid"), &
                               attribute_t("units", "degrees_north")]
        end associate
        associate(var=>var_meta(kVARS%v_longitude))
            var%name        = "v_longitude"
            var%n_dimensions= 2
            var%dimensions  = ["lon_x","lat_x"]
            var%n_attrs     = 2
            var%attributes  = [attribute_t("non_standard_name","longitude_on_v_grid"), &
                               attribute_t("units", "degrees_east")]
        end associate

            ! need to specify : name, n_attrs, attributes%name %value
            ! n_dimensions, dimensions

    end subroutine init_var_meta

end module output_metadata
