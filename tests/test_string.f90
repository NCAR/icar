program main
  use string, only : get_double, get_real, get_integer, str
  use iso_fortran_env, only : real32, real64, int32, output_unit, error_unit
  implicit none

  ! Enumerate the test-result array indices:
  enum, bind(C)
    enumerator ::  get_double_, get_real_, get_integer_, str_d_, str_r_, str_i_
  end enum
  logical :: test_passed(get_double_:str_i_)=.false.

  ! Test data
  real(real64), parameter :: real64_datum=0.3_real64
  real(real32), parameter :: real32_datum=0.3_real32
  integer,      parameter :: integer_datum=1234567

  ! Construct expected results
  character(len=32) :: real64_string, real32_string, integer_string

  write(real64_string,*) real64_datum
  write(real32_string,*) real32_datum
  write(integer_string,*) integer_datum

  ! Test each string conversion function
  if ( get_double(real64_string)   == real64_datum  ) test_passed(get_double_) = .true.
  if ( get_real(real32_string)     == real32_datum  ) test_passed(get_real_) = .true.
  if ( get_integer(integer_string) == integer_datum ) test_passed(get_integer_) = .true.

  ! Left-justify the expected results for comparison
  if ( adjustl(real64_string)      == str(real64_datum)  ) test_passed(str_d_) = .true.
  if ( adjustl(real32_string)      == str(real32_datum)  ) test_passed(str_r_) = .true.
  if ( adjustl(integer_string)     == str(integer_datum) ) test_passed(str_i_) = .true.

  ! Test each string conversion function
  if (all(test_passed)) then
     write(output_unit,*) "Test passed." 
  else
     if ( .not. test_passed(get_double_ )) write(error_unit,*) "Test get_double failed."
     if ( .not. test_passed(get_real_   )) write(error_unit,*) "Test get_real failed."
     if ( .not. test_passed(get_integer_)) write(error_unit,*) "Test get_integer failed."
     if ( .not. test_passed(str_d_      )) write(error_unit,*) "Test str_d failed."
     if ( .not. test_passed(str_r_      )) write(error_unit,*) "Test str_r failed."
     if ( .not. test_passed(str_i_      )) write(error_unit,*) "Test str_i failed."
  end if

end program
