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
    use data_structures

    implicit none

contains
    subroutine radiation_init(domain,options)
        type(domain_type), intent(in) :: domain
        type(options_type),intent(in)    :: options

        write(*,*) "Initializing Radiation"

        if (options%physics%radiation==kRA_BASIC) then
            write(*,*) "    Basic Radiation"
        endif
        if (options%physics%radiation==kRA_SIMPLE) then
            write(*,*) "    Simple Radiation"
            call ra_simple_init(domain,options)
        endif

    end subroutine radiation_init

    subroutine rad(domain,options,date,dt)
        implicit none

        type(domain_type), intent(inout) :: domain
        type(options_type),intent(in)    :: options
        double precision, intent(in) :: date
        real, intent(in) :: dt

        if (options%physics%radiation==kRA_SIMPLE) then
            call ra_simple(domain%th,domain%pii,domain%qv,domain%cloud+domain%ice,domain%qsnow,&
                        domain%qrain,domain%p,domain%swdown,domain%lwdown,domain%cloudfrac,&
                        domain%lat,domain%lon,date,options,dt)
        endif

    end subroutine rad
end module radiation
