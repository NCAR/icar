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
!!   ra_simple          - loops over X,Y grid cells, calls cloudfrac, shortwave,longwave on columns
!!   cloudfrac           - calculates the cloud fraction following Xu and Randall (1996)
!!   shortwave           - calculates shortwave at the surface following Reiff et al (1984)
!!   longwave               - calculates longwave at the surface following Idso and Jackson (1969)
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
!!   swdown = shortwave down at surface - 2D - output - W/m^2   - (nx,ny)
!!   lwdown = longwave down at surface  - 2D - output - W/m^2   - (nx,ny)
!!   dt = time step                      - 0D - input  - seconds    - scalar
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_ra_simple
    use data_structures
    use time
    implicit none
    real, allocatable, dimension(:,:) :: cos_lat_m,sin_lat_m
    real, parameter :: So=1367.0 ! Solar "constant" W/m^2
    real, parameter :: qcmin=1e-5 ! arbitrarily selected minimum "cloud" water content to affect radiation
    real, parameter :: MINIMUM_RH=1e-10 ! bounds on relative humidity = [min_rh to 1-min_rh]
contains
    subroutine ra_simple_init(domain,options)
        implicit none
        type(domain_type), intent(in) :: domain
        type(options_type),intent(in)    :: options
        integer :: nx,ny

        nx=size(domain%lat,1)
        ny=size(domain%lat,2)
        allocate(cos_lat_m(nx,ny))
        allocate(sin_lat_m(nx,ny))

        cos_lat_m=cos(domain%lat/360.0 * 2*pi)
        sin_lat_m=sin(domain%lat/360.0 * 2*pi)

    end subroutine ra_simple_init

!   day_of_year=calc_day_of_year(date)
!   time_of_day=calc_time_of_day(date)


    function relative_humidity(t,qv,p,j,nx)
        implicit none
        real,dimension(nx) :: relative_humidity
        real,intent(in),dimension(nx)::t
        real,intent(in),dimension(:,:,:)::qv,p
        integer,intent(in)::j,nx
        real,dimension(nx) :: mr, e, es

        ! convert sensible humidity to mixing ratio
        mr=qv(:,1,j)/(1-qv(:,1,j))
        ! convert mixing ratio to vapor pressure
        e=mr*p(:,1,j)/(0.62197+mr)
        ! convert temperature to saturated vapor pressure
        es=611.2*exp(17.67*(t-273.15)/(t-29.65))
        ! finally return relative humidity
        relative_humidity= e/es
        ! because it is an approximation things could go awry and rh outside or reasonable bounds could break something else.
        ! alternatively air could be supersaturated (esp. on boundary cells) but cloud fraction calculations will break.
        where(relative_humidity>(1-MINIMUM_RH)) relative_humidity=1-MINIMUM_RH
        where(relative_humidity<MINIMUM_RH) relative_humidity=MINIMUM_RH
    end function relative_humidity

    function shortwave(day_frac, cloud_cover, solar_elevation,nx)
!       compute shortwave down at the surface based on solar elevation, fractional day of the year, and cloud fraction
!       based on Reiff et al. (1984)
        implicit none
        real, dimension(nx) :: shortwave
        real, intent(in), dimension(nx) :: day_frac, cloud_cover, solar_elevation

        real, dimension(nx) :: sin_solar_elev
        integer, intent(in) :: nx

        sin_solar_elev = sin(solar_elevation)
        shortwave=So * (1+0.035*cos(day_frac * 2*pi)) * sin_solar_elev * (0.48+0.29*sin_solar_elev)

        ! note it is cloud_cover**3.4 in Reiff, but this makes almost no difference and integer powers are fast
        shortwave=shortwave * (1 - (0.75 * (cloud_cover**3)) )

    end function shortwave

    function longwave(T_air, cloud_cover,nx)
!       compute longwave down at the surface based on air temperature and cloud fraction
!       based on Idso and Jackson (1969)
        implicit none
        real, dimension(nx) :: longwave
        real, intent(in), dimension(nx) :: T_air,cloud_cover
        real, dimension(nx) :: effective_emissivity
        integer, intent(in) :: nx

        effective_emissivity = 1 - 0.261 * exp((-7.77e-4) * (273.16-T_air)**2)
        longwave = effective_emissivity * stefan_boltzmann * T_air**4

        longwave = longwave * (1+0.2*cloud_cover)
    end function longwave

    function cloudfrac(rh,qc,nx)
!       Calculate the cloud fraction based on cloud water content (qc) and relative humidity
!       based on equations from Xu and Randal (1996)
        implicit none
        real,dimension(nx)::cloudfrac
        real,intent(in),dimension(nx)::rh,qc
        integer, intent(in) :: nx
        real,dimension(nx)::temporary

        temporary=((1-rh)*qc)**0.25
        where(temporary >1) temporary=1
        where(temporary <0.0001) temporary=0.0001

        cloudfrac=qc-qcmin
        where(cloudfrac<0) cloudfrac=0

        cloudfrac=(rh**0.25) * (1-exp((-2000*(cloudfrac)) / temporary))
        where(cloudfrac<0) cloudfrac=0
        where(cloudfrac>1) cloudfrac=1

    end function

    function calc_solar_elevation(date,lon,j,nx,day_frac)
        implicit none
        real, dimension(nx) :: calc_solar_elevation
        double precision, intent(in) :: date
        real, dimension(:,:), intent(in) :: lon
        integer,intent(in) :: j, nx
        real, dimension(:), intent(out) :: day_frac

        integer :: i
        real, dimension(nx) :: declination,day_of_year,hour_angle

        do i=1,nx
            day_of_year(i) = floor(calc_day_of_year(date + lon(i,j)/360.0))
            hour_angle(i) = 2*pi* mod(date+0.5 + lon(i,j)/360.0,1.0)
        end do
        day_frac=day_of_year/365.25

        ! fast approximation see : http://en.wikipedia.org/wiki/Position_of_the_Sun
        declination = (-0.4091) * cos(2.0*pi/365.0*(day_of_year+10))
        calc_solar_elevation = sin_lat_m(:,j) * sin(declination) + &
                               cos_lat_m(:,j) * cos(declination) * cos(hour_angle)

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

    subroutine ra_simple(theta,pii,qv,qc,qs,qr,p,swdown,lwdown,cloud_cover,lat,lon,date,options,dt)
        implicit none
        real,dimension(:,:,:), intent(inout) :: theta
        real,dimension(:,:,:), intent(in) :: pii,qv,qc,qs,qr,p
        real,dimension(:,:), intent(out) :: swdown,lwdown,cloud_cover
        real,dimension(:,:), intent(in) :: lat,lon
        double precision, intent(in) :: date
        type(options_type),intent(in)    :: options
        real, intent(in) :: dt
        real :: coolingrate
        integer :: nx,ny,j,k,nz
        real, allocatable, dimension(:) :: rh,T_air,solar_elevation, hydrometeors,day_frac


        !$omp parallel private(nx,ny,nz,j,rh,T_air,solar_elevation,hydrometeors,day_frac,coolingrate) &
        !$omp shared(theta,pii,qv,p,qc,qs,qr,date,lon,cloud_cover,swdown,lwdown)
        nx=size(lat,1)
        ny=size(lat,2)
        nz=size(qv,2)

        allocate(rh(nx))
        allocate(T_air(nx))
        allocate(solar_elevation(nx))
        allocate(hydrometeors(nx))
        allocate(day_frac(nx))

        ! 1.5K/day radiative cooling rate (300 = W/m^2 at 270K)
        coolingrate=1.5*(dt/86400.0) *stefan_boltzmann / 300.0


        !$omp do
        do j=2,ny-1

            T_air=theta(:,1,j)*pii(:,1,j)
            rh=relative_humidity(T_air,qv,p,j,nx)

            hydrometeors=qc(:,1,j)+qs(:,1,j)+qr(:,1,j)
            do k=2,nz
                hydrometeors=hydrometeors+qc(:,k,j)+qs(:,k,j)+qr(:,k,j)
            end do
            where(hydrometeors<0) hydrometeors = 0

            solar_elevation=calc_solar_elevation(date,lon,j,nx,day_frac)

            cloud_cover(:,j) = cloudfrac(rh,hydrometeors,nx)
            swdown(:,j) = shortwave(day_frac,cloud_cover(:,j),solar_elevation,nx)
            lwdown(:,j) = longwave(T_air,cloud_cover(:,j),nx)
            ! apply a simple radiative cooling to the atmosphere
            theta(2:nx-1,:,j)=theta(2:nx-1,:,j) - (((theta(2:nx-1,:,j)*pii(2:nx-1,:,j))**4) * coolingrate)

        end do
        !$omp end do

        deallocate(rh,T_air,solar_elevation, hydrometeors, day_frac)
        !$omp end parallel

    end subroutine ra_simple
end module module_ra_simple
