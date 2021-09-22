program test_output

    use icar_constants
    use output_interface,   only : output_t
    use variable_interface, only : variable_t
    use output_metadata,    only : get_metadata

    implicit none

    type(output_t) :: dataset
    type(variable_t) :: vars(4)

    character(len=128) :: output_file
    real, pointer :: local(:,:,:)
    real, pointer :: local_2d(:,:)
    integer :: i

    allocate(local(20,4,30))
    allocate(local_2d(20,30))

    do i=1,4
        local(:,i,:) = i
    enddo
    local_2d = 5

    call setup(vars(1), "test_data", local)
    vars(2) = get_metadata( kVARS%potential_temperature, input_data=local)
    vars(3) = get_metadata( kVARS%water_vapor,           input_data=local)
    vars(4) = get_metadata( kVARS%skin_temperature,      input_data=local_2d)

    do i=1,4
        call dataset%add_to_output(vars(i))
    enddo

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
        var%dim_len =  [size(input_data,1), size(input_data,3), size(input_data,2)]
        var%dimensions=[character(len=1024) :: "x", "y", "z"]

        var%n_attrs=3
        allocate(var%attributes(var%n_attrs))
        var%attributes(:)%name  = ["a","b","c"]
        var%attributes(:)%value = ["a1","b1","c1"]

        var%three_d = .True.
        var%data_3d => input_data

    end subroutine setup


end program
