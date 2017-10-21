module array_utilities

    implicit none

contains
    !>----------------------------------------------------------
    !! Generate an array with n_values elements spanning the range min to max
    !!
    !! This is similar to the matlab/numpy/etc linspace function
    !!
    !!----------------------------------------------------------
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

    pure function check_array_dims(input_array, d1, d2, d3, d4, d5) result(passed)
        implicit none
        real,    intent(in) :: input_array(:,:,:,:,:)
        integer, intent(in) :: d1, d2, d3, d4, d5
        logical :: passed

        passed = .True.

        if (size(input_array,1)/=d1) then
            passed = .False.
            return
        endif

        if (size(input_array,2)/=d2) then
            passed = .False.
            return
        endif

        if (size(input_array,3)/=d3) then
            passed = .False.
            return
        endif

        if (size(input_array,4)/=d4) then
            passed = .False.
            return
        endif

        if (size(input_array,5)/=d5) then
            passed = .False.
            return
        endif

    end function check_array_dims

    !>----------------------------------------------------------
    !! Calculate the weights between the positions bestpos and nextpos
    !! based on the distance between match and indata(nextpos) (normalized by nextpos - bestpos)
    !! assumes indata is monotonically increasing,
    !! bestpos must be set prior to entry
    !! nextpos is calculated internally (either 1, bestpos+1, or n)
    !!
    !!----------------------------------------------------------
    function calc_weight(indata, bestpos, nextpos, match) result(weight)
        implicit none
        real :: weight
        real, dimension(:), intent(in) :: indata
        integer, intent(in) :: bestpos
        integer, intent(inout) :: nextpos
        real, intent(in) :: match

        integer :: n

        n=size(indata)

        if (match<indata(1)) then
            nextpos=1
            weight=1
        else
            if (bestpos==n) then
                nextpos=n
                weight=1
            else
                nextpos=bestpos+1
                weight=(indata(nextpos)-match) / (indata(nextpos) - indata(bestpos))
            endif
        endif

    end function


end module array_utilities
