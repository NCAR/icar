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
    !>-----------------------------
    !!  Generic interface to various types of string conversion functions
    !!
    !!-----------------------------
    interface str
        module procedure str_d
        module procedure str_r
        module procedure str_i
    end interface
    
    integer,parameter::MAXSTRINGLENGTH=100
contains
    !>------------------------------
    !! Convert a string to a double precision real number
    !!
    !!------------------------------
    function get_double(str_in)
        implicit none
        character(len=*), intent(in) :: str_in      ! input string
        double precision :: get_double              ! double precision return value
        
        read(str_in,*) get_double
        
    end function get_double
    
    !>------------------------------
    !! Convert a string to a single precision real number
    !!
    !!------------------------------
    function get_real(str_in)
        implicit none
        character(len=*), intent(in) :: str_in      ! input string
        real :: get_real                            ! return real value
        
        read(str_in,*) get_real
        
    end function get_real

    !>------------------------------
    !! Convert a string to a single precision integer
    !!
    !!------------------------------
    function get_integer(str_in)
        implicit none
        character(len=*), intent(in) :: str_in      ! input string
        integer :: get_integer                      ! return integer value
        
        read(str_in,*) get_integer
        
    end function get_integer
    
    
    !>------------------------------
    !! Convert a double precision real number to a string
    !!
    !!  Optionally specify a format string
    !!
    !!------------------------------
    function str_d(value,fmt) result(output_string)
        implicit none
        double precision :: value                       ! double precision value to be converted
        character(len=*), optional :: fmt               ! optional format string for conversion
        character(len=MAXSTRINGLENGTH) :: output_string ! return value
        
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        
        output_string=adjustl(output_string)
    end function str_d

    !>------------------------------
    !! Convert a single precision real number to a string
    !!
    !!  Optionally specify a format string
    !!
    !!------------------------------
    function str_r(value,fmt) result(output_string)
        implicit none
        real :: value                                   ! single precision value to be converted
        character(len=*), optional :: fmt               ! optional format string for conversion
        character(len=MAXSTRINGLENGTH) :: output_string ! return value
        
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        
        output_string=adjustl(output_string)
    end function str_r

    !>------------------------------
    !! Convert a single precision integer to a string
    !!
    !!  Optionally specify a format string
    !!
    !!------------------------------
    function str_i(value,fmt) result(output_string)
        implicit none
        integer :: value                                ! integer value to be converted
        character(len=*), optional :: fmt               ! optional format string for conversion
        character(len=MAXSTRINGLENGTH) :: output_string ! return value
        
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        
        output_string=adjustl(output_string)
    end function str_i


end module string
