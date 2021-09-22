program test_caf_init

    use options_interface,  only : options_t
    use domain_interface,   only : domain_t

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    integer :: i

    call options%init()

    call domain%init(options)


end program test_caf_init
