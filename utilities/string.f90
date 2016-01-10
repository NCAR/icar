!>------------------------------------------------------------
!!  Various functions to convert a number to a string and a string
!!  to a number
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module string

    implicit none
    interface str
        module procedure str_d
        module procedure str_r
        module procedure str_i
    end interface
    
    integer,parameter::MAXSTRINGLENGTH=100
contains
    function get_double(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        double precision :: get_double
        read(str_in,*) get_double
    end function get_double
    
    function get_real(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        real :: get_real
        read(str_in,*) get_real
    end function get_real

    function get_integer(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        integer :: get_integer
        read(str_in,*) get_integer
    end function get_integer
    
    
    
    function str_d(value,fmt) result(output_string)
        implicit none
        double precision :: value
        character(len=*), optional :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
    end function str_d

    function str_r(value,fmt) result(output_string)
        implicit none
        real :: value
        character(len=*), optional :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
    end function str_r

    function str_i(value,fmt) result(output_string)
        implicit none
        integer :: value
        character(len=*), optional :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
    end function str_i


end module string
