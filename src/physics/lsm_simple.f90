!>----------------------------------------------------------
!! NOTE: CODE IS INCOMPLETE AND WILL NOT WORK
!!
!! Very simple land surface model code
!!
!! Rain is partitioned into infiltration and runoff
!! Snow is accumulated on the surface, then melts, runsoff, or sublimates
!! Soil moisture is permitted to be lost to ET or subsurface flow
!!
!! ET, Sensible Heat Flux, and Longwave are partitioned using Penman Monteith.
!!
!! The entry point to the code is lsm_simple.
!!
!! <pre>
!! Call tree graph :
!! lsm_simple->
!!  [->],
!!  [->],
!!  [->]
!!
!! High level routine descriptions / purpose
!!   lsm_simple         - loops over X,Y grid cells, calls a, b, c
!!
!! Driver inputs: p,th,pii,rho,qv,rain,snow,dt,dz
!!   p   = pressure                      - 3D - input  - Pa     - (nx,nz,ny)
!!   th  = potential temperature         - 3D - in/out - K      - (nx,nz,ny)
!!   pii = inverse exner function        - 3D - input  - []     - (nx,nz,ny)
!!   rho = air density                   - 3D - input  - kg/m^3 - (nx,nz,ny)
!!   qv  = specific humidity             - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   rain= rainfall                      - 2D - input  - mm     - (nx,ny)
!!   snow= snowfall                      - 2D - input  - mm     - (nx,ny)
!!   wind= wind speed                    - 2D - input  - m/s    - (nx,ny)
!!   swdown = shortwave down at surface - 2D - input  - W/m^2   - (nx,ny)
!!   lwdown = longwave down at surface  - 2D - input  - W/m^2   - (nx,ny)
!!   dt = time step                      - 0D - input  - seconds    - scalar
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_lsm_simple
    use data_structures
    implicit none
    private
    public lsm_simple_init, lsm_simple

    real :: albedo  ! Default surface albedo
    real :: zo      ! surface roughness
    real :: vegfrac ! vegetation fraction
    real :: emissivity ! land surface emissivity
    real :: Tdeep   ! Deep soil temperature
    real, allocatable,dimension(:,:) :: base_exchange_term,exchange_C,z_atm,lnz_atm_term,Ri,airt_m
    real :: soil_dz ! soil surface layer thickness
    real :: soil_thermal_conductivity
    real :: soil_density
    real :: soil_specific_heat
    real :: damping_time ! damping parameter to return soil temperature to deep soil temperature
    real :: vwc_min ! residual moisture content
    real :: vwc_max ! saturated moisture content
    real :: soil_hydro_conductivity ! soil hydraulic conductivity
    real :: soil_hydro_dz ! soil hydraulic thickness
    real, parameter :: kappa=0.4
contains
    subroutine lsm_simple_init(domain,options)
        implicit none
        type(domain_type), intent(in) :: domain
        type(options_type),intent(in)    :: options
        integer :: nx,ny

        nx=size(domain%terrain,1)
        ny=size(domain%terrain,2)

        albedo=0.2   ! arbitrary albedo
        zo=0.05      ! arbitrary roughness
!       vegfrac=0.75 ! aribtrary starting point
!       Ch=0.005     ! median value from Chen and Zhang 2009
        Tdeep=275.15 ! mean annual air temperature should be specified as an input field

        allocate(exchange_C(nx,ny))
        allocate(base_exchange_term(nx,ny))
        allocate(z_atm(nx,ny))
        allocate(lnz_atm_term(nx,ny))
        allocate(Ri(nx,ny))
        allocate(airt_m(nx,ny))
        z_atm=domain%z(:,1,:)
        lnz_atm_term = log((z_atm+zo)/zo)
        base_exchange_term=(75*kappa**2 * sqrt((z_atm+zo)/zo)) / (lnz_atm_term**2)
        lnz_atm_term=(kappa/lnz_atm_term)**2

        soil_dz=0.10 ! m
        soil_thermal_conductivity=0.28 ! W/m/K
        soil_density=1600 !kg/m^3
        soil_specific_heat=890 ! J/kg/K
        damping_time=300 ! days? sort of...

        vwc_min=0.03 ! residual moisture content
        vwc_max=0.45 ! saturated moisture content
        soil_hydro_conductivity=1e-7 ! m/s
        soil_hydro_dz=1.0     ! m
        emissivity=0.95

!       domain%soil_vwc=0.25
    end subroutine lsm_simple_init

    subroutine calc_ground_heat(tsoil,tskin,ground_heat,dt)
        implicit none
        real,dimension(:,:,:),intent(inout):: tsoil
        real,dimension(:,:),intent(in):: tskin
        real,dimension(:,:),intent(out)  :: ground_heat
        real,intent(in)::dt

        ground_heat = soil_thermal_conductivity * (tskin-tsoil(:,1,:)) / soil_dz

        tsoil(:,1,:)=tsoil(:,1,:) + ground_heat * dt * soil_dz*soil_density*soil_specific_heat &
                    + (Tdeep-tsoil(:,1,:))/damping_time
    end subroutine calc_ground_heat

    subroutine calc_exchange_coefficient(wind,tskin,airt)
        implicit none
        real, dimension(:,:),intent(in) :: wind,tskin,airt

        Ri = gravity/airt * (airt-tskin)*z_atm/wind**2

        where(Ri<0)  exchange_C=lnz_atm_term * (1.0-(15.0*Ri)/(1.0+(base_exchange_term * sqrt((-1.0)*Ri))))
        where(Ri>=0) exchange_C=lnz_atm_term * 1.0/((1.0+15.0*Ri)*sqrt(1.0+5.0*Ri))
    end subroutine calc_exchange_coefficient

    subroutine calc_sensible_heat(tskin,airt,wind,sensible_heat)
        implicit none
        real,dimension(:,:),intent(in)  :: tskin,airt,wind
        real,dimension(:,:),intent(out) :: sensible_heat

        sensible_heat = exchange_C * wind * (tskin-airt)

    end subroutine calc_sensible_heat

    subroutine calc_latent_heat(vwc,tskin,qv,wind,latent_heat)
        implicit none
        real,dimension(:,:,:),intent(in) :: vwc
        real,dimension(:,:),  intent(in) :: tskin,wind
        real,dimension(:,:,:),intent(in) :: qv
        real,dimension(:,:),  intent(out):: latent_heat

!       latent_heat = (vwc-vwc_min)/(vwc_max-vwc_min) * exchange_C * wind * (saturated(tskin) - qv(:,1,:))
        where(latent_heat<0) latent_heat=0

    end subroutine calc_latent_heat

    subroutine add_rain(vwc,rain,dt)
        implicit none
        real,dimension(:,:,:),intent(inout):: vwc
        real,dimension(:,:),  intent(in)   :: rain
        real,intent(in)::dt

        where((rain/1000.0/dt)>=soil_hydro_conductivity) &
            vwc(:,1,:)=vwc(:,1,:)+soil_hydro_conductivity*dt/soil_hydro_dz
        where((rain/1000.0/dt)<soil_hydro_conductivity) &
            vwc(:,1,:)=vwc(:,1,:)+rain/1000.0/dt/soil_hydro_dz
        where(vwc>vwc_max) vwc=vwc_max
    end subroutine add_rain

    subroutine lsm_simple(theta,pii,qv,rain,snow,p,swdown,lwdown, wind, &
                          sensible_heat, latent_heat, ground_heat,      &
                          tskin, tsoil, vwc, swe, options,dt)
        implicit none
        real,dimension(:,:,:), intent(inout) :: theta, qv, tsoil, vwc
        real,dimension(:,:,:), intent(in) :: pii,p
        real,dimension(:,:), intent(in) :: rain,snow,swdown,lwdown,wind
        real,dimension(:,:), intent(inout) :: sensible_heat, latent_heat, ground_heat, &
                                              tskin, swe
        type(options_type),intent(in)    :: options
        real, intent(in) :: dt
        integer :: nx,ny,j,k,nz

        airt_m=theta(:,1,:)*pii(:,1,:)

        call calc_ground_heat(tsoil,tskin,ground_heat,dt)
        call calc_exchange_coefficient(wind,tskin,airt_m)
        call calc_sensible_heat(tskin,airt_m,wind,sensible_heat)
        call calc_latent_heat(vwc,tskin,qv,wind,latent_heat)

        tskin = ((swdown*(1-albedo)+lwdown*emissivity    &
                  -ground_heat-sensible_heat-latent_heat) &
                 /stefan_boltzmann/emissivity)**0.25

!       call apply_fluxes(vwc,tsoil,qv,airt_m,latent_heat,sensible_heat,dt)
!       theta(:,1,:)=airt_m/pii(:,1,:)

        call add_rain(vwc,rain,dt)

    end subroutine lsm_simple
end module module_lsm_simple
