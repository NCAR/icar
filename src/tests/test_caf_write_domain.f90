program test_caf_write_domain

    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use output_interface,   only : output_t
    use grid_interface,     only : grid_t
    use mod_atm_utilities,  only : exner_function, pressure_at_elevation, sat_mr
    use icar_constants,     only : kVARS

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    type(output_t)  :: dataset
    character(len=128) :: file_name
    integer :: i

    print*, "Reading Options Files"
    call options%init()
    sync all

    call options%alloc_vars([kVARS%temperature])


    print*, "Initializing Domain"
    call domain%init(options)
    sync all

    print*, "Adding domain to output dataset"
    call dataset%set_domain(domain)
    sync all

    print*, "Adding variables to output dataset"
    call dataset%add_variables(options%vars_for_restart, domain)
    sync all

    print*, "Initializing domain with ideal values"
    call initialize_ideal_domain(domain)
    sync all

    print*, "Writing sample output file"
    write(file_name, '("Initial_output_",I3.3,".nc")') this_image()
    call dataset%save_file(file_name, 1, domain%model_time)
    sync all

    print*, ""
    print*, "Test Completed"
    print*, ""

contains

    subroutine initialize_ideal_domain(domain)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(grid_t) :: grid

        integer :: nx,ny,nz, i

        nx=100
        ny=150
        nz=20

        call grid%set_grid_dimensions(nx, ny, nz)

        associate(                                              &
            u_test_val=>0.5, v_test_val=>0.5, w_test_val=>0.0,  &
            water_vapor_test_val            => 0.001,           &
            potential_temperature_test_val  => 300.0,           &
            cloud_water_mass_test_val       => 0.0,             &
            cloud_ice_mass_test_val         => 0.0,             &
            cloud_ice_number_test_val       => 0.0,             &
            rain_mass_test_val              => 0.0,             &
            rain_number_test_val            => 0.0,             &
            snow_mass_test_val              => 0.0,             &
            graupel_mass_test_val           => 0.0)

            ! domain%u%data_3d = u_test_val
            ! domain%v%data_3d = v_test_val
            ! domain%w%data_3d = w_test_val
            domain%water_vapor%data_3d = water_vapor_test_val
            domain%potential_temperature%data_3d = potential_temperature_test_val
            domain%cloud_water_mass%data_3d = cloud_water_mass_test_val
            ! domain%cloud_ice_mass%data_3d = cloud_ice_mass_test_val
            ! domain%cloud_ice_number%data_3d = cloud_ice_number_test_val
            domain%rain_mass%data_3d = rain_mass_test_val
            ! domain%rain_number%data_3d = rain_number_test_val
            domain%snow_mass%data_3d = snow_mass_test_val
            ! domain%graupel_mass%data_3d = graupel_mass_test_val

        end associate


        associate(                                                &
            kms=>domain%grid%kms, kme=>domain%grid%kme,           &
            sealevel_pressure     => 100000.0,                    &
            pressure              => domain%pressure%data_3d,     &
            exner                 => domain%exner%data_3d,        &
            temperature           => domain%temperature%data_3d,  &
            potential_temperature => domain%potential_temperature%data_3d,   &
            water_vapor           => domain%water_vapor%data_3d,  &
            z                     => domain%z%data_3d,            &
            dz                    => domain%dz_mass%data_3d )

            domain%accumulated_precipitation%data_2d = 0
            domain%accumulated_snowfall%data_2d      = 0

            z(:,kms,:) = domain%terrain%data_2d + dz(:,kms,:)/2
            do i=kms,kme
                if (i>kms) then
                    z(:,i,:) = z(:,i-1,:) + dz(:,i-1,:)/2 + dz(:,i,:)/2
                endif

                pressure(:,i,:)    = pressure_at_elevation(sealevel_pressure, z(:,i,:))
            enddo

            exner       = exner_function(pressure)
            temperature = exner * potential_temperature
            water_vapor = sat_mr(temperature,pressure)
        end associate

    end subroutine

end program
