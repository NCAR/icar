
! Small test for mp_simple code to be sure it computes rain and or snow
program main
  use module_mp_simple
  real, dimension(5) :: pressure, temperature, qv, qc, qr, qs, dz, rho
  real   :: rain, snow, dt
  integer:: nz, i

  nz = 5
  pressure = 80000.0
  temperature = 280.0
  rho = 1
  qv = 0.005
  qc = 0.
  qr = 0.
  qs = 0.
  dz = 200.0
  dt = 20.0
  rain=0
  snow=0

  cloud2snow = exp(-1.0*snow_formation_time_const*dt)
  cloud2rain = exp(-1.0*rain_formation_time_const*dt)

  print*, "dt=",dt, "   snow_formation_time_const=", snow_formation_time_const, "   cloud2snow=", cloud2snow
  print*, "dt=",dt, "   rain_formation_time_const=", rain_formation_time_const, "   cloud2rain=", cloud2rain

  ! slowly cool down and see if we keep getting precipitation out (including snow)
  do i=1,100
      call mp_simple(pressure, temperature, qv, qc, qr, qs, rain, snow, dt, dx2, dz, 1, nz, 1, nz)
      write(*,*) temperature(1),qv(1),qc(1),qr(1),qs(1),rain,snow
      temperature = temperature - 0.1
  end do
end program
