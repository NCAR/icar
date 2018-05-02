!>----------------------------------------------------------
!! This module provides a wrapper to call various radiation models
!! It sets up variables specific to the physics package to be used
!!
!! The main entry point to the code is rad(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  radiation_init->[ external initialization routines]
!!  rad->[  external radiation routines]
!!
!! High level routine descriptions / purpose
!!   radiation_init     - initializes physics package
!!   rad                - sets up and calls main physics package
!!
!! Inputs: domain, options, dt
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module radiation
    use module_ra_simple, only: ra_simple, ra_simple_init
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use data_structures
    use icar_constants, only : kVARS

    implicit none

contains
    subroutine radiation_init(domain,options)
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options

        if (this_image()==1) write(*,*) "Initializing Radiation"

        if (options%physics%radiation==kRA_BASIC) then
            if (this_image()==1) write(*,*) "    Basic Radiation"
        endif
        if (options%physics%radiation==kRA_SIMPLE) then
            if (this_image()==1) write(*,*) "    Simple Radiation"
            call ra_simple_init(domain, options)
        endif

    end subroutine radiation_init


    subroutine ra_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        if (options%physics%radiation == kRA_SIMPLE) then
            call ra_simple_var_request(options)
        endif

    end subroutine ra_var_request


    !> ----------------------------------------------
    !! Communicate to the master process requesting the variables requred to be allocated, used for restart files, and advected
    !!
    !! ----------------------------------------------
    subroutine ra_simple_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options


        ! List the variables that are required to be allocated for the simple radiation code
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%cloud_fraction,   &
                      kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,  kVARS%snow_in_air,      &
                      kVARS%shortwave,   kVARS%longwave,                kVARS%cloud_ice,    kVARS%graupel_in_air])

        ! List the variables that are required to be advected for the simple radiation code
        call options%advect_vars( &
                      [kVARS%potential_temperature] )

        ! List the variables that are required when restarting for the simple radiation code
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature, kVARS%shortwave,   kVARS%longwave, kVARS%cloud_fraction] )

    end subroutine ra_simple_var_request


    subroutine rad(domain, options, dt, halo, subset)
        implicit none

        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt
        integer,        intent(in), optional :: halo, subset

        integer :: ims, ime, jms, jme, kms, kme
        integer :: its, ite, jts, jte, kts, kte

        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        kts = domain%grid%kts
        kte = domain%grid%kte

        if (options%physics%radiation==kRA_SIMPLE) then
            call ra_simple(theta = domain%potential_temperature%data_3d,         &
                           pii= domain%exner%data_3d,                            &
                           qv = domain%water_vapor%data_3d,                      &
                           qc = domain%cloud_water_mass%data_3d,                 &
                           qs = domain%snow_mass%data_3d                         &
                                + domain%cloud_ice_mass%data_3d                  &
                                + domain%graupel_mass%data_3d,                   &
                           qr = domain%rain_mass%data_3d,                        &
                           p =  domain%pressure%data_3d,                         &
                           swdown =  domain%shortwave%data_2d,                   &
                           lwdown =  domain%longwave%data_2d,                    &
                           cloud_cover =  domain%cloud_fraction%data_2d,         &
                           lat = domain%latitude%data_2d,                        &
                           lon = domain%longitude%data_2d,                       &
                           date = domain%model_time,                             &
                           options = options,                                    &
                           dt = dt,                                              &
                           ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                           its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
        endif

    end subroutine rad
end module radiation
