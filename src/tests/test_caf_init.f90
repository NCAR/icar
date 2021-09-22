program test_caf_init

    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    type(boundary_t):: boundary
    integer :: i

    call options%init()

    call domain%init(options)

    call boundary%init(options)

    call domain%get_initial_conditions(boundary, options)


end program test_caf_init
