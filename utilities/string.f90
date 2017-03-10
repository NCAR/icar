!>------------------------------------------------------------
!!
!!  Various functions to convert a number to a string and a string
!!  to a number
!!
!!  Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module string

    use iso_fortran_env, only : real32, real64

    implicit none
    interface str
        module procedure str_d
        module procedure str_r
        module procedure str_i
    end interface
    
    integer,parameter::MAXSTRINGLENGTH=100
contains
    elemental function get_double(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        real(real64) :: get_double
        read(str_in,*) get_double
    end function get_double
    
    elemental function get_real(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        real(real32) :: get_real
        read(str_in,*) get_real
    end function get_real

    elemental function get_integer(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        integer :: get_integer
        read(str_in,*) get_integer
    end function get_integer
    
    elemental function str_d(value,fmt) result(output_string)
        implicit none
        real(real64), intent(in) :: value
        character(len=*), optional, intent(in) :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
    end function str_d

    elemental function str_r(value,fmt) result(output_string)
        implicit none
        real(real32), intent(in) :: value
        character(len=*), optional, intent(in) :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
    end function str_r

    elemental function str_i(value,fmt) result(output_string)
        implicit none
        integer, intent(in) :: value
        character(len=*), optional, intent(in) :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
    end function str_i

end module string
