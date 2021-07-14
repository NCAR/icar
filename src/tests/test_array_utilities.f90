program test_array_utilities

    use array_utilities

    implicit none

    real, parameter :: SMALL_VALUE = 1e-6

    call test_linear_space

contains

    subroutine test_linear_space
        implicit none
        real, allocatable :: test_array(:)

        real    :: vmax, vmin, dv, current_dv
        integer :: i, n
        logical :: passed = .True.

        n    = 100
        vmin = 1e-4
        vmax = 3

        call linear_space(test_array,vmin,vmax,n)

        if (.not.allocated(test_array)) then
            stop "FAILED: linear_space failed to allocate input array"
        endif

        if (size(test_array)/=n) then
            print*, size(test_array), n
            stop "FAILED: linear_space allocated the wrong size array"
        endif

        if (abs(test_array(1)-vmin) > SMALL_VALUE) then
            print*, "FAILED: linear_space array start value /= minimum value"
            print*, test_array(1), vmin, test_array(1)-vmin
            passed = .False.
        endif

        if (abs(test_array(n)-vmax) > SMALL_VALUE) then
            print*, "FAILED: linear_space array end value /= maximum value"
            print*, test_array(n), vmax, test_array(n)-vmax
            passed = .False.
        endif

        dv = (vmax-vmin) / (n-1)

        do i=2,n
            current_dv = test_array(i)-test_array(i-1)

            if (abs(current_dv-dv) > SMALL_VALUE) then
                print*, "FAILED: linear_space wrong step size in side array"
                print*, test_array(i-1), test_array(i), current_dv, dv
                passed = .False.
            endif
        enddo

        if (passed) print*, "PASSED: linear_space"

    end subroutine test_linear_space
end program test_array_utilities
