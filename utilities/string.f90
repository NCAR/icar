!>------------------------------------------------------------
!!
!!  Various functions to convert a number to a string and a string
!!  to a number
!!
!!  Author: Ethan Gutmann (gutmann@ucar.edu)
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
        character(len=*), intent(in) :: str_in
        double precision :: get_double
        read(str_in,*) get_double
    end function get_double
    
    function get_real(str_in)
        character(len=*), intent(in) :: str_in
        real :: get_real
        read(str_in,*) get_real
    end function get_real

    function get_integer(str_in)
        character(len=*), intent(in) :: str_in
        integer :: get_integer
        read(str_in,*) get_integer
    end function get_integer
    
    
    
    function str_d(value)
        double precision :: value
        character(len=MAXSTRINGLENGTH) :: str_d
        write(str_d,*) value
        str_d=adjustl(str_d)
    end function str_d

    function str_r(value)
        real :: value
        character(len=MAXSTRINGLENGTH) :: str_r
        write(str_r,*) value
        str_r=adjustl(str_r)
    end function str_r
    
    function str_i(value)
        integer :: value
        character(len=MAXSTRINGLENGTH) :: str_i
        write(str_i,*) value
        str_i=adjustl(str_i)
    end function str_i

end module string
