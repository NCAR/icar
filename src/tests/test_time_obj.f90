! Test if the precision of the time_delta_t object will lead to
!   precision errors
program test_time_object
  use iso_fortran_env, only: real64, real128, int64
  use time_object, only: time_type
  use time_delta_object, only: time_delta_t
  implicit none

  type(time_type) :: t_obj
  type(time_delta_t) :: dt
  real(real128) :: t128, err_r128
  real :: delta_r32, tmp1
  real(real64) :: delta_r64, tmp2
  integer(int64) :: i, n, i64 , err_i64
  character(len=44) :: c_i

  ! add delta_t n times to the time object
  n = 10000000
  write(c_i, '(i44)') n

  ! --- setup experiment ---
  ! tmp variables used so compiler doesn't optimize out
  tmp1 = 3600 / 86400.0
  delta_r32 = tmp1 * 86400.0

  tmp2 = 3600 / 86400.0
  delta_r64 = tmp2 * 86400.0
  call t_obj%set(1980, 1, 1, 0, 0, 0)
  call dt%set(seconds=delta_r32)
  t128 = t_obj%seconds()
  i64 = int(t_obj%seconds(), kind=int64)

  print *, "iterating ", trim(adjustl(c_i)), " times"
  print *, "delta_r32 = ", delta_r32
  print *, "delta_r64 = ", delta_r64
  print *, "--- running ---"

  ! do loop incrementing time object
  do i = 0, n
      ! increment time
      t_obj  = t_obj + dt
      t128  = t128 + delta_r32
      i64  = i64 + 3600

      ! calculate difference between time object and variables of type real128
      !   and type int64
      err_r128 = abs(t128 - t_obj%seconds())
      err_i64 = int(abs(i64 - t_obj%seconds()), kind=int64)

      ! report if error difference too large
      if (err_r128 .gt. 1.0) then
          print *, "err = ", err_r128
          stop "ERROR ERR_R128 > 1.0"
      end if
      if (err_i64 .gt. 0) then
          print *, "err_i64 = ", err_i64
          stop "ERROR ERR_I64 > 0"
      end if
  end do

  ! report or expirement
  write(c_i, '(i44)') i64
  print *, "SUCCESS"
  print *, "t_obj = ", t_obj%seconds(), "seconds"
  print *, "t128  = ", t128, "seconds"
  print *, "i64   =    ", trim(adjustl(c_i)), " seconds"

  print *, "final err_r128 = ", err_r128
  print *, "final err_i64  = ", err_i64

end program test_time_object
