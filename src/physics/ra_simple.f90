!>----------------------------------------------------------
!! Very simple radiation code modeled after the description in Finch and Best(2004?)
!! of Reiff et al. 1984 shortwave and Idso and Jackson (1969) longwave.
!!
!! Clearsky Shortwave radiation is calculated as a function of day of year and time of day.
!! Cloudy Shortwave is calculated as clearsky SW * f(cloud cover) [0.25-1]
!!
!! Cloud cover is calculated as in Xu and Randal (1996) as f(surface_RH, qc+qs+qr)
!!
!! Clearsky Longwave radiation is calculated as f(Tair)
!! Cloudy longwave is scaled up by a f(cloud cover) [1-1.2]
!!
!! The entry point to the code is ra_simple.
!!
!! <pre>
!! Call tree graph :
!! ra_simple->
!!  [cloudfrac->],
!!  [shortwave->],
!!  [longwave->]
!!
!! High level routine descriptions / purpose
!!   ra_simple           - loops over X,Y grid cells, calls cloudfrac, shortwave,longwave on columns
!!   cloudfrac           - calculates the cloud fraction following Xu and Randall (1996)
!!   shortwave           - calculates shortwave at the surface following Reiff et al (1984)
!!   longwave            - calculates longwave at the surface following Idso and Jackson (1969)
!!
!! Driver inputs: p,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz
!!   p   = pressure                      - 3D - input  - Pa     - (nx,nz,ny)
!!   th  = potential temperature         - 3D - in/out - K      - (nx,nz,ny)
!!   pii = inverse exner function        - 3D - input  - []     - (nx,nz,ny)
!!   rho = air density                   - 3D - input  - kg/m^3 - (nx,nz,ny)
!!   qv  = specific humidity             - 3D - input  - kg/kg  - (nx,nz,ny)
!!   qc  = cloud water content           - 3D - input  - kg/kg  - (nx,nz,ny)
!!   qr  = rain water content            - 3D - input  - kg/kg  - (nx,nz,ny)
!!   qs  = snow water content            - 3D - input  - kg/kg  - (nx,nz,ny)
!!   swdown = shortwave down at surface  - 2D - output - W/m^2  - (nx,ny)
!!   lwdown = longwave down at surface   - 2D - output - W/m^2  - (nx,ny)
!!   dt = time step                      - 0D - input  - seconds    - scalar
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_ra_simple
    use time_object,        only : Time_type
    use mod_atm_utilities,  only : relative_humidity
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use data_structures
    ! use time

    implicit none

    real, allocatable, dimension(:,:) :: cos_lat_m,sin_lat_m
    integer :: nrad_layers
    real, parameter :: So = 1367.0  ! Solar "constant" W/m^2
    real, parameter :: qcmin = 1e-6 ! Minimum "cloud" water content to affect radiation
contains

    subroutine ra_simple_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in)    :: options

        associate(ims => domain%grid%ims,   &
                  ime => domain%grid%ime,   &
                  jms => domain%grid%jms,   &
                  jme => domain%grid%jme)

        allocate(cos_lat_m(ims:ime, jms:jme))
        allocate(sin_lat_m(ims:ime, jms:jme))

        cos_lat_m = cos(domain%latitude%data_2d / 360.0 * 2*pi)
        sin_lat_m = sin(domain%latitude%data_2d / 360.0 * 2*pi)

        nrad_layers = 5

        end associate
    end subroutine ra_simple_init


    function shortwave(day_frac, cloud_cover, solar_elevation, ims,ime, its,ite)
        ! compute shortwave down at the surface based on solar elevation, fractional day of the year, and cloud fraction
        ! based on Reiff et al. (1984)
        implicit none
        real             :: shortwave       (ims:ime)
        real, intent(in) :: day_frac        (ims:ime)
        real, intent(in) :: cloud_cover     (ims:ime)
        real, intent(in) :: solar_elevation (ims:ime)
        integer, intent(in) :: ims,ime, its,ite

        real :: sin_solar_elev(its:ite)

        sin_solar_elev = sin(solar_elevation(its:ite))

        shortwave(its:ite) = So * (1 + 0.035 * cos(day_frac(its:ite) * 2*pi)) * sin_solar_elev * (0.48 + 0.29 * sin_solar_elev)

        ! note it is cloud_cover**3.4 in Reiff, but this makes almost no difference and integer powers are fast so could use **3
        shortwave(its:ite) = shortwave(its:ite) * (1 - (0.75 * (cloud_cover(its:ite)**3.4)) )

    end function shortwave

    function longwave(T_air, cloud_cover, ims,ime, its,ite)
        ! compute longwave down at the surface based on air temperature and cloud fraction
        ! based on Idso and Jackson (1969)
        implicit none
        real                :: longwave(ims:ime)
        real,    intent(in) :: T_air(ims:ime), cloud_cover(ims:ime)
        integer, intent(in) :: ims,ime, its,ite
        real :: effective_emissivity(its:ite)

        effective_emissivity = 1 - 0.261 * exp((-7.77e-4) * (273.16-T_air(its:ite))**2)

        longwave(its:ite) = effective_emissivity * stefan_boltzmann * T_air(its:ite)**4

        longwave(its:ite) = min(longwave(its:ite) * (1 + 0.2 * cloud_cover(its:ite)), 600.0)

    end function longwave

    function cloudfrac(rh, qc, ims,ime, its,ite)
        ! Calculate the cloud fraction based on cloud water content (qc) and relative humidity
        ! based on equations from Xu and Randal (1996)
        implicit none
        real                :: cloudfrac(ims:ime)
        real,    intent(in) :: rh(ims:ime), qc(ims:ime)
        integer, intent(in) :: ims, ime, its, ite

        real :: temporary(its:ite)

        cloudfrac = 0

        temporary = ((1 - rh(its:ite)) * qc(its:ite))**0.25
        where(temporary > 1) temporary=1
        where(temporary < 0.0001) temporary=0.0001

        cloudfrac(its:ite) = qc(its:ite) - qcmin
        where(cloudfrac < 5e-8) cloudfrac = 5e-8

        cloudfrac(its:ite) = (rh(its:ite)**0.25) * (1 - exp((-2000*(cloudfrac(its:ite))) / temporary))

        where(cloudfrac < 0) cloudfrac = 0
        where(cloudfrac > 1) cloudfrac = 1

    end function

    function calc_solar_elevation(date, lon, j, ims,ime, jms,jme, its,ite, day_frac)
        implicit none
        real                       :: calc_solar_elevation(ims:ime)
        type(Time_type),intent(in) :: date
        real,           intent(in) :: lon(ims:ime, jms:jme)
        integer,        intent(in) :: j
        integer,        intent(in) :: ims, ime, jms, jme
        integer,        intent(in) :: its, ite
        real,           intent(out):: day_frac(ims:ime)

        integer :: i
        real, dimension(ims:ime) :: declination, day_of_year, hour_angle

        calc_solar_elevation = 0

        do i = its, ite
            day_of_year(i) = date%day_of_year(lon=lon(i,j))

            ! hour angle is 0 at noon
            hour_angle(i) = 2*pi* mod(day_of_year(i)+0.5, 1.0)

            day_frac(i) = date%year_fraction(lon=lon(i,j))
        end do

        ! fast approximation see : http://en.wikipedia.org/wiki/Position_of_the_Sun
        declination = (-0.4091) * cos(2.0*pi/365.0*(day_of_year+10))

        calc_solar_elevation(its:ite) = sin_lat_m(its:ite,j) * sin(declination(its:ite)) + &
                               cos_lat_m(its:ite,j) * cos(declination(its:ite)) * cos(hour_angle(its:ite))

        ! due to float precision errors, it is possible to exceed (-1 - 1) in which case asin will break
        where(calc_solar_elevation < -1)
            calc_solar_elevation = -1
        elsewhere(calc_solar_elevation > 1)
            calc_solar_elevation = 1
        endwhere

        calc_solar_elevation = asin(calc_solar_elevation)

        ! if the sun is below the horizon just set elevation to 0
        where(calc_solar_elevation<0) calc_solar_elevation=0
    end function calc_solar_elevation

    subroutine ra_simple(theta, pii, qv, qc, qs, qr, p, swdown, lwdown, cloud_cover, lat, lon, date, options, dt, &
                ims, ime, jms, jme, kms, kme, &
                its, ite, jts, jte, kts, kte, &
                F_runlw)
        implicit none
        real, intent(inout):: theta (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: pii   (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: qv    (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: qc    (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: qs    (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: qr    (ims:ime, kms:kme, jms:jme)
        real, intent(in)   :: p     (ims:ime, kms:kme, jms:jme)
        real, intent(out)  :: swdown     (ims:ime, jms:jme)
        real, intent(out)  :: lwdown     (ims:ime, jms:jme)
        real, intent(out)  :: cloud_cover(ims:ime, jms:jme)
        real, intent(in)   :: lat        (ims:ime, jms:jme)
        real, intent(in)   :: lon        (ims:ime, jms:jme)
        type(Time_type),    intent(in) :: date
        type(options_t), intent(in) :: options
        real,               intent(in) :: dt
        integer,            intent(in) :: ims, ime, jms, jme, kms, kme
        integer,            intent(in) :: its, ite, jts, jte, kts, kte
        logical,            intent(in), optional :: F_runlw

        logical :: runlw
        real :: coolingrate
        integer :: j, k
        real, allocatable, dimension(:) :: rh, T_air, solar_elevation, hydrometeors, day_frac


        runlw = .True.
        if (present(F_runlw)) runlw = F_runlw

        !$omp parallel private(j,k,rh,T_air,solar_elevation,hydrometeors,day_frac,coolingrate) &
        !$omp shared(theta,pii,qv,p,qc,qs,qr,date,lon,cloud_cover,swdown,lwdown)                      &
        !$omp firstprivate(runlw, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        allocate(rh             (ims:ime))
        allocate(T_air          (ims:ime))
        allocate(solar_elevation(ims:ime))
        allocate(hydrometeors   (ims:ime))
        allocate(day_frac       (ims:ime))

        ! 1.5K/day radiative cooling rate (300 = W/m^2 at 270K)
        coolingrate = 1.5 * (dt / 86400.0) * stefan_boltzmann / 300.0

        !$omp do
        do j = jts, jte

            T_air = 0
            rh = 0
            do k = kts, kts + nrad_layers - 1
                T_air = T_air + (theta(:,k,j)*pii(:,k,j))
                rh    = rh    + relative_humidity((theta(:,k,j)*pii(:,k,j)), qv(:,k,j), p(:,k,j))
            enddo
            T_air = T_air / nrad_layers
            rh    = rh    / nrad_layers
            where(rh > 1) rh = 1

            hydrometeors = qc(:,kts,j) + qs(:,kts,j) + qr(:,kts,j)
            do k = kts+1, kte
                hydrometeors = hydrometeors + qc(:,k,j) + qs(:,k,j) + qr(:,k,j)
            end do
            where(hydrometeors<0) hydrometeors = 0

            solar_elevation  = calc_solar_elevation(date, lon, j, ims,ime, jms,jme, its,ite, day_frac)

            cloud_cover(:,j) = cloudfrac(rh, hydrometeors, ims,ime, its,ite)
            swdown(:,j)      = shortwave(day_frac, cloud_cover(:,j), solar_elevation, ims,ime, its,ite)
            if (runlw) then
                lwdown(:,j)      = longwave(T_air, cloud_cover(:,j), ims,ime, its,ite)

                ! apply a simple radiative cooling to the atmosphere
                theta(its:ite, kts:kte, j) = theta(its:ite, kts:kte, j) - (((theta(its:ite, kts:kte, j) * pii(its:ite, kts:kte, j)) ** 4) * coolingrate)
            endif
        end do
        !$omp end do

        deallocate(rh,T_air, solar_elevation, hydrometeors, day_frac)
        !$omp end parallel

    end subroutine ra_simple
end module module_ra_simple
