submodule (restart_interface) restart_implementation

    use io_routines, only : io_read
    use icar_constants, only : kMAX_FILE_LENGTH
    use time_object,        only : Time_type

! options%restart_step_in_file = restart_step
! options%restart_file         = restart_file
! options%restart_date         = restart_date
! options%restart_time         = restart_time

    implicit none


contains

module subroutine restart_model(domain, dataset, options)
    implicit none
    class(domain_t),  intent(inout) :: domain
    class(output_t),  intent(inout) :: dataset
    class(options_t), intent(inout) :: options

    character(len=kMAX_FILE_LENGTH) :: restart_file

    ! options%parameters%restart_file,
    restart_file = get_image_filename(this_image(), options%parameters%output_file, options%parameters%restart_time)

    call read_restart_data(domain, dataset, restart_file, options%parameters%restart_step_in_file)


end subroutine

subroutine read_restart_data(domain, dataset, filename, time_step)
    implicit none
    class(domain_t),  intent(inout) :: domain
    class(output_t),  intent(inout) :: dataset
    character(len=*), intent(in)    :: filename
    integer,          intent(in)    :: time_step

    integer :: i
    integer :: dim_3d(3)
    real, allocatable :: data_3d(:,:,:)
    real, allocatable :: data_2d(:,:)

    do i=1,dataset%n_variables

        associate(var => dataset%variables(i))

            if (var%three_d) then
                dim_3d = var%dim_len
                call io_read(filename, var%name, data_3d, extradim=time_step)
                if (associated(var%data_3d)) then
                    if (size(var%data_3d) /= size(data_3d)) then
                        call restart_domain_error(var%name)
                    endif
                    var%data_3d = reshape(data_3d, shape=dim_3d, order=[1,3,2])
                else
                    print*, "ERROR, variable not ready to be used:"//trim(var%name)
                endif

            elseif (var%two_d) then
                call io_read(filename, var%name, data_2d, extradim=time_step)
                if (associated(var%data_2d)) then
                    if (size(var%data_2d) /= size(data_2d)) then
                        call restart_domain_error(var%name)
                    endif
                    var%data_2d(:,:) = data_2d(:,:)
                else
                    print*, "ERROR, variable not ready to be used:"//trim(var%name)
                endif

            endif
        end associate
    end do

end subroutine read_restart_data

function get_image_filename(image_number, initial_filename, restart_time) result(file_name)
    implicit none
    integer,            intent(in) :: image_number
    character(len=*),   intent(in) :: initial_filename
    type(Time_type),    intent(in) :: restart_time

    character(len=kMAX_FILE_LENGTH) :: file_name
    integer :: n, i

    character(len=49)   :: file_date_format = '(I4,"-",I0.2,"-",I0.2,"_",I0.2,"-",I0.2,"-",I0.2)'

    write(file_name, '(A,I6.6,"_",A,".nc")') trim(initial_filename), image_number, trim(restart_time%as_string(file_date_format))

    n = len(trim(file_name))
    file_name(n-10:n-9) = "00"

end function get_image_filename

subroutine restart_domain_error(varname)
    implicit none
    character(len=*), intent(in) :: varname

    print*, "Error reading restart variable: ", trim(varname)
    print*, "The domain of the restart file does not match the current run"
    print*, "This can happen if you run a different number of parallel processes"

    error stop

end subroutine restart_domain_error

end submodule
