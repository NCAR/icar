module array_utilities

    implicit none

contains
    subroutine linear_space(input_array, min_value, max_value, n_values)
        implicit none
        real,   intent(inout), allocatable :: input_array(:)
        real,   intent(in)  :: min_value, max_value
        integer,intent(in)  :: n_values

        integer :: i

        if (allocated(input_array)) then
            if (size(input_array)/=n_values) then
                deallocate(input_array)
            endif
        endif
        ! note, this can't be an else statement because the above block might change the result
        if (.not.allocated(input_array)) then
            allocate(input_array(n_values))
        endif

        do i=1,n_values
            input_array(i) = (i-1.0)/real(n_values-1.0) * (max_value - min_value) + min_value
        enddo

    end subroutine linear_space

end module array_utilities
