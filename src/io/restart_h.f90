module restart_interface

    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use output_interface,   only : output_t

    implicit none

    interface
        module subroutine restart_model(domain, dataset, options)
            implicit none
            type(domain_t),  intent(inout) :: domain
            type(output_t),  intent(inout) :: dataset
            type(options_t), intent(inout) :: options

        end subroutine
    end interface

end module restart_interface
