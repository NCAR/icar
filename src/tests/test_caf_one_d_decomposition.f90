program caf_one_d_decomposition_test

    implicit none

    integer, parameter :: n_speeds = 31
    integer, parameter :: n_dir = 7
    integer, parameter :: n_stability = 10

    integer :: total = n_speeds * n_dir * n_stability

    integer :: start_pos[*], stop_pos[*]
    integer :: i
    logical :: passed[*] = .True.

    start_pos = nint((real(this_image()-1) / num_images()) * total) + 1

    if (this_image()==num_images()) then
        stop_pos = total
    else
        stop_pos  = nint((real(this_image()) / num_images()) * total)
    endif

    sync all

    if (this_image() < num_images()) then
        passed = (stop_pos+1 == start_pos[this_image()+1])
    endif

    sync all

    if (this_image() == 1) then

        passed = passed .and. ((start_pos[1]==1) .and. (stop_pos[num_images()]==total))

        do i = 1, num_images()
            passed = passed .and. passed[i]
        enddo

        if (passed) print*, "Test Passed"
    endif


end program caf_one_d_decomposition_test
