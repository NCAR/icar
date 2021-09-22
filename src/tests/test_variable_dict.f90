program test_variable_dict

    use variable_dict_interface,    only : var_dict_t
    use variable_interface,         only : variable_t
    use grid_interface,             only : grid_t

    implicit none

    type(var_dict_t)    :: var_collection
    type(variable_t)    :: var1, var2, var3, output_var
    type(grid_t)        :: grid
    logical :: passed = .True.

    ! define a grid for a test variable
    call grid%set_grid_dimensions(nx=10,ny=10,nz=5)
    ! initialize a test variable with that grid
    call var1%initialize(grid)

    call grid%set_grid_dimensions(nx=2,ny=15,nz=5)
    ! initialize a test variable with that grid
    call var2%initialize(grid)

    ! initialize a test variable with a shape spec instead of a grid spec
    call var3%initialize( [2,5] )

    ! add the test variable to the dictionary for the key "data"
    call var_collection%add_var("data_first",var1)
    call var_collection%add_var("data",var2)
    call var_collection%add_var("data_third",var3)

    output_var = var_collection%get_var("data")

    if (.not.(output_var%three_d .eqv. var2%three_d))          passed = .False.
    if (.not.(output_var%n_dimensions == var2%n_dimensions))   passed = .False.
    if (.not.all(output_var%dim_len == var2%dim_len))          passed = .False.

    if (passed) then
        print*, "test_variable_dict Passed"
    else
        print*, "test_variable_dict Failed"
    endif
end program
