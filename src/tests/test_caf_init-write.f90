program test_caf_init_write

    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use output_interface,   only : output_t

    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    type(boundary_t):: boundary
    type(output_t)  :: dataset

    character(len=1024) :: file_name
    integer :: i

    if (this_image()==1) print*, "Reading options"
    call options%init()

    if (this_image()==1) print*, "Initializing Domain"
    call domain%init(options)

    if (this_image()==1) print*, "Initializing boundary condition data structure"
    call boundary%init(options)

    if (this_image()==1) print*, "Reading Initial conditions from boundary dataset"
    call domain%get_initial_conditions(boundary, options)

    if (this_image()==1) print*, "Adding domain to output dataset"
    call dataset%set_domain(domain)

    if (this_image()==1) print*, "Adding variables to output dataset"
    call dataset%add_variables(options%vars_for_restart, domain)

    if (this_image()==1) print*, "Writing sample output file"
    write(file_name, '("icar_restart_output_",I3.3,".nc")') this_image()
    call dataset%save_file(file_name, 1, domain%model_time)

    if (this_image()==1) print*, "Model time = ", trim(domain%model_time%as_string())

end program
