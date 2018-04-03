program test_caf_linear_winds_setup

    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use linear_theory_winds,only : setup_linwinds

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    integer :: i

    call options%init()

    call domain%init(options)

    call setup_linwinds(domain, options, reverse=.False., useDensity=.False.)


end program test_caf_linear_winds_setup
