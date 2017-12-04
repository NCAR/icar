program test_output

    use output_interface, only : output_t
    use variable_interface, only : variable_t

    implicit none

    type(output_t) :: dataset
    type(variable_t) :: vars(4)

    character(len=128) :: output_file
    real, pointer :: local(:,:,:)

    allocate(local(20,2,30))

    call setup(vars(1), "test_a", local)
    call setup(vars(2), "test_b", local)
    call setup(vars(3), "test_c", local)
    call setup(vars(4), "test_d", local)

    call dataset%add_to_output(vars(2))
    call dataset%add_to_output(vars(4))

    write(output_file,"(A,I0.3,A)") "test_", this_image(),".nc"
    call dataset%save_file(output_file)

contains
    subroutine setup(var, name, input_data)
        implicit none
        type(variable_t), intent(inout) :: var
        character(len=*), intent(in) :: name
        real, pointer :: input_data(:,:,:)

        var%name = name


        var%n_dimensions = 3
        var%dim_len =  [size(input_data,1), size(input_data,2), size(input_data,3)]
        var%dimensions=[character(len=1024) :: "x", "y", "z"]

        var%n_attrs=3
        allocate(var%attributes(var%n_attrs))
        var%attributes(:)%name  = ["a","b","c"]
        var%attributes(:)%value = ["a1","b1","c1"]

        var%three_d = .True.
        var%data_3d => input_data

    end subroutine setup


end program
