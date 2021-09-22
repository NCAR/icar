!>----------------------------------------------------------
!! This module provides a wrapper to call various land surface models
!!
!! It sets up variables specific to the LSM to be used including both
!! history variables not currently stored in the domain level data
!! structure, and runtime parameters
!!
!! The main entry point to the code is lsm(domain,options,dt,model_time)
!!
!! <pre>
!! Call tree graph :
!!  lsm_init->[ allocate_noah_data,
!!              external initialization routines]
!!  lsm->[  sat_mr,
!!          calc_exchange_coefficient,
!!          external LSM routines]
!!
!! High level routine descriptions / purpose
!!   lsm_init           - allocates module data and initializes physics package
!!   lsm                - sets up and calls main physics package
!!  calc_exchange_coefficient - calculates surface exchange coefficient (for Noah)
!!  allocate_noah_data  - allocate module level data for Noah LSM
!!  apply_fluxes        - apply LSM fluxes (e.g. sensible and latent heat fluxes) to atmosphere
!!
!! Inputs: domain, options, dt, model_time
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!!      model_time      = time since beginning date (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module land_surface
    use module_sf_noahdrv,   only : lsm_noah, lsm_noah_init
    ! use module_lsm_basic,    only : lsm_basic
    ! use module_lsm_simple,   only : lsm_simple, lsm_simple_init
    use module_water_simple, only : water_simple
    use mod_atm_utilities,   only : sat_mr
    use time_object,         only : Time_type
    use data_structures
    use options_interface,   only : options_t
    use domain_interface,    only : domain_t

    implicit none

    private
    public :: lsm_init, lsm, lsm_var_request

    ! Noah LSM required variables.  Some of these should be stored in domain, but tested here for now
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions

    ! LOTS of variables required by Noah, placed here "temporarily", this may be where some stay.
    ! Keeping them at the module level prevents having to allocate/deallocate every call
    ! also avoids adding LOTS of LSM variables to the main domain datastrucurt
    real,allocatable, dimension(:,:)    :: SMSTAV,SFCRUNOFF,UDRUNOFF,                               &
                                           SNOW,SNOWC,SNOWH, ACSNOW, ACSNOM, SNOALB, QFX,           &
                                           QGH, GSW, ALBEDO, ALBBCK, Z0, XICE, EMISS,               &
                                           EMBCK, QSFC, RAINBL, CHS, CHS2, CQS2, CPM, SR,           &
                                           CHKLOWQ, LAI, QZ0, VEGFRAC, SHDMIN,SHDMAX,SNOTIME,SNOPCX,&
                                           POTEVP,RIB, NOAHRES,FLX4_2D,FVB_2D,FBUR_2D,              &
                                           FGSN_2D, z_atm,lnz_atm_term,Ri,base_exchange_term,       &
                                           current_precipitation

    integer,allocatable, dimension(:,:) :: rain_bucket ! used to start the previous time step rain bucket

    logical :: MYJ, FRPCPN,ua_phys,RDLAI2D,USEMONALB
    real,allocatable, dimension(:,:,:)  :: SH2O,SMCREL
    real,allocatable, dimension(:,:)    :: dTemp,lhdQV, windspd
    real,allocatable, dimension(:)      :: Zs,DZs
    real :: XICE_THRESHOLD
    integer,allocatable, dimension(:,:) :: IVGTYP,ISLTYP
    integer :: ITIMESTEP, update_interval, cur_vegmonth

!     real, parameter :: kappa=0.4 ! this should be karman from data_structure
    real, parameter :: freezing_threshold=273.15
    real, parameter :: SMALL_PRESSURE=0.1 !note: 0.1Pa is very small 1e-10 wouldn't affect a single-precision float
    real, parameter :: SMALL_QV=1e-10
    real, parameter :: MAX_EXCHANGE_C = 0.5
    real, parameter :: MIN_EXCHANGE_C = 0.004

    character(len=MAXVARLENGTH) :: MMINLU
    logical :: FNDSOILW,FNDSNOWH,RDMAXALB
    integer :: num_soil_layers,ISURBAN,ISICE,ISWATER
    integer :: exchange_term
    real*8  :: last_model_time

contains


    subroutine lsm_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%landsurface > 1) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%snow_height,                                                                             &  ! BK 2020/10/26
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,                                      &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent,                               &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height,                                                                             &  ! BK 2020/10/26
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%veg_type])    ! BK uncommented 2021/03/20
                         ! kVARS%soil_type, kVARS%land_mask, kVARS%vegetation_fraction]
        endif
        if (options%physics%watersurface > 1) then
            call options%alloc_vars( &
                         [kVARS%sst, kVARS%ustar, kVARS%surface_pressure, kVARS%water_vapor,            &
                         kVARS%temperature, kVARS%sensible_heat, kVARS%latent_heat, kVARS%land_mask,    &
                         kVARS%humidity_2m, kVARS%temperature_2m, kVARS%skin_temperature, kVARS%u_10m, kVARS%v_10m])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%sst, kVARS%potential_temperature, kVARS%water_vapor, kVARS%skin_temperature,        &
                         kVARS%surface_pressure, kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,  &
                         kVARS%humidity_2m, kVARS%temperature_2m])
        endif


    end subroutine lsm_var_request

    subroutine calc_exchange_coefficient(wind,tskin,airt,exchange_C)
        implicit none
        real, dimension(:,:),intent(inout) :: wind,tskin
        real, dimension(:,:,:),intent(inout) :: airt
        real,dimension(:,:),intent(inout) :: exchange_C

        ! Richardson number
        where(wind==0) wind=1e-5
        exchange_C = 0

        Ri = gravity/airt(:,1,:) * (airt(:,1,:)-tskin)*z_atm/(wind**2)

        ! print*,"--------------------------------------------------"
        ! print*, "Surface Richardson number"
        where(Ri<0)  exchange_C = lnz_atm_term * (1.0-(15.0*Ri)/(1.0+(base_exchange_term * sqrt((-1.0)*Ri))))
        where(Ri>=0) exchange_C = lnz_atm_term * 1.0/((1.0+15.0*Ri)*sqrt(1.0+5.0*Ri))

        where(exchange_C > MAX_EXCHANGE_C) exchange_C=MAX_EXCHANGE_C
        where(exchange_C < MIN_EXCHANGE_C) exchange_C=MIN_EXCHANGE_C
    end subroutine calc_exchange_coefficient


! eqn A11 in Appendix A.2 of Chen et al 1997 (see below for reference)
    subroutine F2_formula(F2, z_atm, zo, Ri)
        real, dimension(:,:), intent(inout) :: F2, z_atm, zo, Ri

        ! for the stable case from Mahrt (1987)
        where(Ri>=0) F2=exp(-Ri)
        ! for the unstable case from Holtslag and Beljaars (1989)
        where(Ri<0)  F2=(1-(15*Ri)/(1+((70.5*karman**2 * sqrt(-Ri * z_atm/zo))/(lnz_atm_term**2))) )

    end subroutine F2_formula
!From Appendix A.2 in Chen et al 1997
! Impact of Atmospheric Surface-layer Parameterizations in the new Land-surface Scheme of the Ncep Mesoscale ETA Model
! Boundary-Layer Meteorology 85:391-421
    subroutine calc_mahrt_holtslag_exchange_coefficient(wind,tskin,airt,znt, exchange_C)
        implicit none
        real, dimension(:,:),intent(inout) :: wind,tskin
        real, dimension(:,:,:),intent(inout) :: airt
        real,dimension(:,:),intent(inout) :: exchange_C, znt

        ! Richardson number
        where(wind==0) wind=1e-10
        Ri = gravity/airt(:,1,:) * (airt(:,1,:)-tskin)*z_atm/(wind**2)

        call F2_formula(base_exchange_term, z_atm,znt,Ri)

        exchange_C = karman**2 * base_exchange_term / lnz_atm_term**2

        where(exchange_C > MAX_EXCHANGE_C) exchange_C=MAX_EXCHANGE_C
        where(exchange_C < MIN_EXCHANGE_C) exchange_C=MIN_EXCHANGE_C
    end subroutine calc_mahrt_holtslag_exchange_coefficient

    subroutine surface_diagnostics(HFX, QFX, TSK, QSFC, CHS2, CQS2,T2, Q2, PSFC)
        ! taken almost directly / exactly from WRF's module_sf_sfcdiags.F
        implicit none
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  HFX, QFX, TSK, QSFC
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(INOUT) ::  Q2, T2
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  PSFC, CHS2, CQS2
        integer :: i,j, nx,ny
        real :: rho

        !$omp parallel default(shared), private(nx,ny,i,j,rho)
        nx=size(HFX,1)
        ny=size(HFX,2)
        !$omp do
        do j=jts,jte
            do i=its,ite
                RHO = PSFC(I,J)/(Rd * TSK(I,J))
                if(CQS2(I,J).lt.1.E-3) then
                   Q2(I,J) = QSFC(I,J)
                else
                   Q2(I,J) = QSFC(I,J) - QFX(I,J)/(RHO*CQS2(I,J))
                endif
                if(CHS2(I,J).lt.1.E-3) then
                   T2(I,J) = TSK(I,J)
                else
                   T2(I,J) = TSK(I,J) - HFX(I,J)/(RHO*CP*CHS2(I,J))
                endif
                ! TH2(I,J) = T2(I,J)*(1.E5/PSFC(I,J))**ROVCP
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine surface_diagnostics

    subroutine apply_fluxes(domain,dt)
        ! add sensible and latent heat fluxes to the first atm level
        implicit none
        type(domain_t), intent(inout) :: domain
        real, intent(in) :: dt

        associate(density       => domain%density%data_3d,       &
                  sensible_heat => domain%sensible_heat%data_2d, &
                  latent_heat   => domain%latent_heat%data_2d,   &
                  dz            => domain%dz_interface%data_3d,  &
                  pii           => domain%exner%data_3d,         &
                  th            => domain%potential_temperature%data_3d, &
                  qv            => domain%water_vapor%data_3d    &
            )

        ! convert sensible heat flux to a temperature delta term
        ! (J/(s*m^2) * s / (J/(kg*K)) => kg*K/m^2) ... /((kg/m^3) * m) => K
        dTemp=(sensible_heat(its:ite,jts:jte) * dt/cp)  &
             / (density(its:ite,kts,jts:jte) * dz(its:ite,kts,jts:jte))
        ! add temperature delta and convert back to potential temperature
        th(its:ite,kts,jts:jte) = th(its:ite,kts,jts:jte) + (dTemp / pii(its:ite,kts,jts:jte))

        ! convert latent heat flux to a mixing ratio tendancy term
        ! (J/(s*m^2) * s / (J/kg) => kg/m^2) ... / (kg/m^3 * m) => kg/kg
        lhdQV=(latent_heat(its:ite,jts:jte) / LH_vaporization * dt) &
             / (density(its:ite,kts,jts:jte) * dz(its:ite,kts,jts:jte))
        ! add water vapor in kg/kg
        qv(its:ite,kts,jts:jte) = qv(its:ite,kts,jts:jte) + lhdQV

        ! enforce some minimum water vapor content... just in case
        where(qv(its:ite,kts,jts:jte) < SMALL_QV) qv(its:ite,kts,jts:jte) = SMALL_QV

        end associate

    end subroutine apply_fluxes

    subroutine allocate_noah_data(num_soil_layers)
        implicit none
        integer, intent(in) :: num_soil_layers
        integer :: i

        ITIMESTEP=1

        allocate(SMSTAV(ims:ime,jms:jme))
        SMSTAV = 0.5 !average soil moisture available for transp (between SMCWLT and SMCMAX)
        allocate(SFCRUNOFF(ims:ime,jms:jme))
        SFCRUNOFF = 0
        allocate(UDRUNOFF(ims:ime,jms:jme))
        UDRUNOFF = 0
        allocate(SNOW(ims:ime,jms:jme))
        SNOW = 0
        allocate(SNOWC(ims:ime,jms:jme))
        SNOWC = 0
        allocate(SNOWH(ims:ime,jms:jme))
        SNOWH = 0
        allocate(ACSNOW(ims:ime,jms:jme))
        ACSNOW = 0
        allocate(ACSNOM(ims:ime,jms:jme))
        ACSNOM = 0
        allocate(SNOALB(ims:ime,jms:jme))
        SNOALB = 0.8

        allocate(QGH(ims:ime,jms:jme))
        QGH = 0.02 ! saturated mixing ratio at ~20C
        allocate(GSW(ims:ime,jms:jme))
        GSW = 0

        allocate(ALBEDO(ims:ime,jms:jme))
        ALBEDO = 0.17
        allocate(ALBBCK(ims:ime,jms:jme))
        ALBBCK = 0.17 !?
        allocate(XICE(ims:ime,jms:jme))
        XICE = 0
        allocate(EMISS(ims:ime,jms:jme))
        EMISS = 0.95
        allocate(EMBCK(ims:ime,jms:jme))
        EMBCK = 0.95
        allocate(CPM(ims:ime,jms:jme))
        CPM = 0
        allocate(SR(ims:ime,jms:jme))
        SR = 0
        allocate(CHKLOWQ(ims:ime,jms:jme))
        CHKLOWQ = 0
        allocate(LAI(ims:ime,jms:jme))
        LAI = 3
        allocate(QZ0(ims:ime,jms:jme))
        QZ0 = 0 ! used to check for saturation? but only relevant if myj == True

        allocate(FLX4_2D(ims:ime,jms:jme))
        allocate(FVB_2D(ims:ime,jms:jme))
        allocate(FBUR_2D(ims:ime,jms:jme))
        allocate(FGSN_2D(ims:ime,jms:jme))

        allocate(SHDMIN(ims:ime,jms:jme))
        SHDMIN = 0
        allocate(SHDMAX(ims:ime,jms:jme))
        SHDMAX = 100
        allocate(SNOTIME(ims:ime,jms:jme))
        SNOTIME = 0
        allocate(SNOPCX(ims:ime,jms:jme))
        SNOPCX = 0
        allocate(POTEVP(ims:ime,jms:jme))
        POTEVP = 0
        allocate(SMCREL(ims:ime,num_soil_layers,jms:jme))
        SMCREL = 0
        allocate(RIB(ims:ime,jms:jme))
        RIB = 0
        allocate(NOAHRES(ims:ime,jms:jme))
        NOAHRES = 0
        allocate(VEGFRAC(ims:ime,jms:jme))
        VEGFRAC = 50


        XICE_THRESHOLD = 1
        RDLAI2D = .false.
        USEMONALB = .false.
        MYJ = .false.
        FRPCPN = .false. ! set this to true and calculate snow ratio to use microphysics based snow/rain partitioning
        ua_phys = .false.

        allocate(SH2O(ims:ime,num_soil_layers,jms:jme))
        SH2O = 0.25

        allocate(Zs(num_soil_layers))
        allocate(DZs(num_soil_layers))
        DZs = [0.1,0.3,0.6,1.0]
        Zs(1) = DZs(1)/2
        do i = 2,num_soil_layers
            Zs(i) = Zs(i-1) + DZs(i)/2 + DZs(i-1)/2
        end do

    end subroutine allocate_noah_data

    subroutine lsm_init(domain,options)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        integer :: i

        if (options%physics%landsurface == 0) return

        if (this_image()==1) write(*,*) "Initializing LSM"

        if (this_image()==1) write(*,*) "    max soil_deep_temperature on init: ", maxval(domain%soil_deep_temperature%data_2d)

        exchange_term = 1

        ! module level variables for easy access... need to think about tiling to permit halo processing separately.
        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        kds = domain%grid%kds
        kde = domain%grid%kde
        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        kts = domain%grid%kts
        kte = domain%grid%kte

        allocate(dTemp(its:ite,jts:jte))
        dTemp = 0
        allocate(lhdQV(its:ite,jts:jte))
        lhdQV = 0
        allocate(Z0(ims:ime,jms:jme))
        Z0 = domain%roughness_z0%data_2d ! this should get updated by the LSM(?)
        allocate(QSFC(ims:ime,jms:jme))
        QSFC = domain%water_vapor%data_3d(:,kms,:) ! this should get updated by the lsm
        allocate(Ri(ims:ime,jms:jme))
        Ri = 0
        allocate(z_atm(ims:ime,jms:jme))
        z_atm = 50

        allocate(lnz_atm_term(ims:ime,jms:jme))
        lnz_atm_term = 0.1

        allocate(base_exchange_term(ims:ime,jms:jme))
        base_exchange_term = 0.01

        allocate(QFX(ims:ime,jms:jme))
        QFX = 0

        allocate(current_precipitation(ims:ime,jms:jme))
        current_precipitation = 0

        allocate(windspd(ims:ime,jms:jme))
        windspd = 3

        ! NOTE, these fields have probably not been initialized yet...
        ! windspd = sqrt(domain%u10**2+domain%v10**2)

        allocate(CHS(ims:ime,jms:jme))
        CHS = 0.01

        allocate(CHS2(ims:ime,jms:jme))
        CHS2 = 0.01

        allocate(CQS2(ims:ime,jms:jme))
        CQS2 = 0.01


        allocate(RAINBL(ims:ime,jms:jme))
        RAINBL = domain%accumulated_precipitation%data_2d  ! used to store last time step accumulated precip so that it can be subtracted from the current step
                            ! set to domain%rain incase this is a restart run and rain is non-zero to start
        allocate(rain_bucket(ims:ime,jms:jme))
        rain_bucket = domain%precipitation_bucket  ! used to store last time step accumulated precip so that it can be subtracted from the current step


        ! initial guesses (not needed?)
        domain%temperature_2m%data_2d = domain%temperature%data_3d(:,kms,:)
        domain%humidity_2m%data_2d = domain%water_vapor%data_3d(:,kms,:)

        if (options%physics%landsurface==kLSM_SIMPLE) then
            write(*,*) "    Simple LSM (may not work?)"
            stop "Simple LSM not settup, choose a different LSM options"
            ! call lsm_simple_init(domain,options)
        endif
        ! Noah Land Surface Model
        if (options%physics%landsurface==kLSM_NOAH) then
            if (this_image()==1) write(*,*) "    Noah LSM"

            num_soil_layers=4

            ! if (this_image()==1) then
            !     write(*,*) "    options%parameters%external_files: ", trim(options%parameters%external_files)
            !     write(*,*) "    options%parameters%restart: ", options%parameters%restart
            !     write(*,*) "    options%parameters%rho_snow_ext ", trim(options%parameters%rho_snow_ext)
            !     write(*,*) "    options%parameters%swe_ext ", trim(options%parameters%swe_ext )
            ! endif

            if (options%parameters%rho_snow_ext /="" .AND. options%parameters%swe_ext /="") then ! calculate snowheight from external swe and density, but only if both are provided. (Swe alone will give FNDSNW = F)
                FNDSNOWH = .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            elseif(options%parameters%hsnow_ext /="" ) then  ! read in external snowheight if supplied
                FNDSNOWH= .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            elseif(options%parameters%restart) then   ! If restarting read in snow height, but only if this is in restart file?
                FNDSNOWH= .True.
                if (this_image()==1) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
            else
                FNDSNOWH=.False. ! calculate SNOWH from SNOW
            endif

            FNDSOILW=.False. ! calculate SOILW (this parameter is ignored in LSM_NOAH_INIT)
            RDMAXALB=.False.

            ISURBAN = options%lsm_options%urban_category
            ISICE   = options%lsm_options%ice_category
            ISWATER = options%lsm_options%water_category
            MMINLU  = options%lsm_options%LU_Categories !"MODIFIED_IGBP_MODIS_NOAH"

            call allocate_noah_data(num_soil_layers)

            if (options%lsm_options%monthly_vegfrac) then
                VEGFRAC = domain%vegetation_fraction%data_3d(:, domain%model_time%month, :)
            else
                VEGFRAC = domain%vegetation_fraction%data_3d(:, 1, :)
            endif
            cur_vegmonth = domain%model_time%month

            ! save the canopy water in a temporary variable in case this is a restart run because lsm_init resets it to 0
            CQS2 = domain%canopy_water%data_2d
            ! prevents init from failing when processing water points that may have "soil_t"=0
            where(domain%soil_temperature%data_3d<200) domain%soil_temperature%data_3d=200
            where(domain%soil_water_content%data_3d<0.0001) domain%soil_water_content%data_3d=0.0001

            call LSM_NOAH_INIT( VEGFRAC,                             &
                                domain%snow_water_equivalent%data_2d,         & !SNOW, &  BK 18/03/2021
                                SNOWC,                              &
                                domain%snow_height%data_2d,         & !SNOWH, &   BK 18/03/2021
                                domain%canopy_water%data_2d,        &
                                domain%soil_temperature%data_3d,    & !-- SMSTAV      Soil moisture availability for evapotranspiration ( fraction between SMCWLT and SMCMXA)
                                domain%soil_water_content%data_3d,  &
                                SFCRUNOFF,                          &
                                UDRUNOFF,                           &
                                ACSNOW,                             &
                                ACSNOM,                             &
                                domain%veg_type,                    &
                                domain%soil_type,                   &
                                domain%soil_temperature%data_3d,    &
                                domain%soil_water_content%data_3d,  &
                                SH2O,                               &
                                ZS,                                 &
                                DZS,                                &
                                MMINLU,                             &
                                SNOALB,                             &
                                FNDSOILW,                           &
                                FNDSNOWH,                           &
                                RDMAXALB,                           &
                                num_soil_layers,                    &
                                .False.,                            & ! nlayers, is_restart (can't yet)
                                .True. ,                            & ! allowed_to_read (e.g. soilparm.tbl)
                                ids,ide, jds,jde, kds,kde,          &
                                ims,ime, jms,jme, kms,kme,          &
                                its,ite, jts,jte, kts,kte  )

            domain%canopy_water%data_2d = CQS2
            CQS2=0.01
            where(domain%veg_type==ISWATER) domain%land_mask=kLC_WATER ! ensure VEGTYPE (land cover) and land-sea mask are consistent
        endif

        ! defines the height of the middle of the first model level
        z_atm = domain%z%data_3d(:,kts,:) - domain%terrain%data_2d
        lnz_atm_term = log((z_atm+Z0)/Z0)
        if (exchange_term==1) then
            base_exchange_term=(75*karman**2 * sqrt((z_atm+Z0)/Z0)) / (lnz_atm_term**2)
            lnz_atm_term=(karman/lnz_atm_term)**2
        endif

        update_interval=options%lsm_options%update_interval
        last_model_time=-999

    end subroutine lsm_init


    subroutine lsm(domain,options,dt)
        implicit none

        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real, intent(in) :: dt
        real ::lsm_dt
        integer :: nx,ny,i,j

        if (options%physics%landsurface == 0) return

        if (last_model_time==-999) then
            last_model_time = domain%model_time%seconds()-update_interval
        endif
        nx=size(domain%water_vapor%data_3d,1)
        ny=size(domain%water_vapor%data_3d,3)
        if ((domain%model_time%seconds() - last_model_time) >= update_interval) then
            lsm_dt = domain%model_time%seconds() - last_model_time
            last_model_time = domain%model_time%seconds()

            ! if (this_image()==1) write(*,*) "    lsm start: snow_water_equivalent max:", MAXVAL(domain%snow_water_equivalent%data_2d)

            ! exchange coefficients
            windspd = sqrt(domain%u_10m%data_2d**2 + domain%v_10m%data_2d**2)
            if (exchange_term==1) then
                call calc_exchange_coefficient(windspd,domain%skin_temperature%data_2d,domain%temperature%data_3d,CHS)
            elseif (exchange_term==2) then
                call calc_mahrt_holtslag_exchange_coefficient(windspd,domain%skin_temperature%data_2d,domain%temperature%data_3d,domain%roughness_z0%data_2d,CHS)
            endif

            ! --------------------------------------------------
            ! First handle the open water surface options
            ! --------------------------------------------------
            ! if (options%physics%watersurface==kWATER_BASIC) then
                ! Note, do nothing because QFX and QSFC are only used for to calculate diagnostic
                !    T2m and Q2m.  However, the fluxes and stability terms are not coordinated, so
                !    This leads to problems in the current formulation and this has been removed.
                ! do j=1,ny
                !     do i=1,nx
                !         if (domain%landmask(i,j)==kLC_WATER) then
                !             QFX(i,j) = domain%latent_heat(i,j) / LH_vaporization
                !             QSFC(i,j)=sat_mr(domain%T2m(i,j),domain%psfc(i,j))
                !         endif
                !     enddo
                ! enddo
            ! else
            if (options%physics%watersurface==kWATER_SIMPLE) then

                call water_simple(domain%sst%data_2d,                   &
                                  domain%surface_pressure%data_2d,      &
                                  windspd,                              &
                                  domain%ustar,                         &
                                  domain%water_vapor%data_3d,           &
                                  domain%temperature%data_3d,           &
                                  domain%sensible_heat%data_2d,         &
                                  domain%latent_heat%data_2d,           &
                                  z_atm, domain%roughness_z0%data_2d,   &
                                  domain%land_mask,                     &
                                  QSFC,                                 &
                                  QFX,                                  &
                                  domain%skin_temperature%data_2d)
            endif

            where(windspd<1) windspd=1 ! minimum wind speed to prevent the exchange coefficient from blowing up
            CHS = CHS * windspd
            CHS2 = CHS
            CQS2 = CHS

            ! --------------------------------------------------
            ! Now handle the land surface options
            ! --------------------------------------------------
            ! if (options%physics%landsurface==kLSM_BASIC) then
                ! call lsm_basic(domain,options,lsm_dt)
                ! Note, do nothing because QFX and QSFC are only used for to calculate diagnostic
                !    T2m and Q2m.  However, the fluxes and stability terms are not coordinated, so
                !    This leads to problems in the current formulation and this has been removed.
                ! do j=1,ny
                !     do i=1,nx
                !         if (domain%landmask(i,j)==kLC_LAND) then
                !             QFX(i,j) = domain%latent_heat(i,j) / LH_vaporization
                !             QSFC(i,j)=max(domain%water_vapor%data_3d(i,1,j),0.5*sat_mr(domain%T2m(i,j),domain%psfc(i,j)))
                !         endif
                !     enddo
                ! enddo


            ! else
            if (options%physics%landsurface == kLSM_SIMPLE) then
                write(*,*) "--------------------------"
                stop "Simple LSM not implemented yet"
                ! call lsm_simple(domain%th,domain%pii,domain%qv,domain%current_rain, domain%current_snow,domain%p_inter, &
                !                 domain%swdown,domain%lwdown, sqrt(domain%u10**2+domain%v10**2), &
                !                 domain%sensible_heat%data_2d, domain%latent_heat, domain%ground_heat_flux, &
                !                 domain%skin_t, domain%soil_t, domain%soil_vwc, domain%snow_swe, &
                !                 options,lsm_dt)

            else if (options%physics%landsurface == kLSM_NOAH) then
                ! Call the Noah Land Surface Model

                ! 2m saturated mixing ratio
                do j=jms,jme
                    do i=ims,ime
                        if (domain%land_mask(i,j) == kLC_LAND) then
                            QGH(i,j) = sat_mr(domain%temperature_2m%data_2d(i,j),domain%surface_pressure%data_2d(i,j))
                        endif
                    enddo
                enddo
                if (options%lsm_options%monthly_vegfrac) then
                    if (cur_vegmonth /= domain%model_time%month) then
                        VEGFRAC = domain%vegetation_fraction%data_3d(:, domain%model_time%month, :)
                        cur_vegmonth = domain%model_time%month
                    endif
                endif

                ! if (this_image()==1) write(*,*) "    lsm start: accumulated_precipitation max:", MAXVAL(domain%accumulated_precipitation%data_2d)
                ! if (this_image()==1) write(*,*) "    lsm start: RAINBL max:", MAXVAL(RAINBL)
                ! if (this_image()==1) write(*,*) "    lsm start: domain%precipitation_bucket max:", MAXVAL(domain%precipitation_bucket)
                ! if (this_image()==1) write(*,*) "    lsm start: rain_bucket max:", MAXVAL(rain_bucket)


                ! RAINBL(i,j) = [kg m-2]   RAINBL = domain%accumulated_precipitation%data_2d  ! used to store last time step accumulated precip so that it can be subtracted from the current step
                current_precipitation = (domain%accumulated_precipitation%data_2d-RAINBL)+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE
                call lsm_noah(domain%dz_interface%data_3d,                &
                            domain%water_vapor%data_3d,                   &
                            domain%pressure_interface%data_3d,            &
                            domain%temperature%data_3d,                   &
                            domain%skin_temperature%data_2d,              &
                            domain%sensible_heat%data_2d,                 &
                            QFX,                                          &
                            domain%latent_heat%data_2d,                   &
                            domain%ground_heat_flux%data_2d,              &
                            QGH,                                          &
                            GSW,                                          &
                            domain%shortwave%data_2d,                     &
                            domain%longwave%data_2d,                      &
                            SMSTAV,                                       &
                            domain%soil_totalmoisture%data_2d,            &  ! this is not defined on lsm_init (BK 2021/03/20)
                            SFCRUNOFF,                                    &
                            UDRUNOFF,                                     &
                            domain%veg_type,                              &
                            domain%soil_type,                             &
                            ISURBAN,                                      &
                            ISICE,                                        &
                            VEGFRAC,                                      &
                            ALBEDO,                                       &
                            ALBBCK,                                       &
                            domain%roughness_z0%data_2d,                  &
                            Z0,                                           &
                            domain%soil_deep_temperature%data_2d,         &
                            real(domain%land_mask),                       &
                            XICE,                                         &
                            EMISS,                                        &
                            EMBCK,                                        &
                            SNOWC,                                        &
                            QSFC,                                         &
                            current_precipitation,                        &  ! RAINBL
                            MMINLU,                                       &
                            num_soil_layers,                              &
                            lsm_dt,                                       &
                            DZS,                                          &
                            ITIMESTEP,                                    &
                            domain%soil_water_content%data_3d,            &
                            domain%soil_temperature%data_3d,              &
                            domain%snow_water_equivalent%data_2d,         &
                            domain%canopy_water%data_2d,                  &
                            CHS,                                          &
                            CHS2,                                         &
                            CQS2,                                         &
                            CPM,                                          &
                            ROVCP,                                        &
                            SR,                                           &
                            chklowq,                                      &
                            lai,                                          &
                            qz0,                                          & !H
                            myj,frpcpn,                                   &
                            SH2O,                                         &
                            domain%snow_height%data_2d, &     !SNOWH,                                   & !H
                            SNOALB,SHDMIN,SHDMAX,                         & !I
                            SNOTIME,                                      & !?
                            ACSNOM,ACSNOW,                                & !O
                            SNOPCX,                                       & !O
                            POTEVP,                                       & !O
                            SMCREL,                                       & !O
                            XICE_THRESHOLD,                               &
                            RDLAI2D,USEMONALB,                            &
                            Ri,                                           & !I
                            NOAHRES,                                      &
                            ua_phys,flx4_2d,fvb_2d,fbur_2d,fgsn_2d,       & ! Noah UA changes
                            ids,ide, jds,jde, kds,kde,                    &
                            ims,ime, jms,jme, kms,kme,                    &
                            its,ite, jts,jte, kts,kte)

                ! now that znt (roughness_z0) has been updated, we need to recalculate terms
                lnz_atm_term = log((z_atm+domain%roughness_z0%data_2d)/domain%roughness_z0%data_2d)
                if (exchange_term==1) then
                    base_exchange_term=(75*karman**2 * sqrt((z_atm+domain%roughness_z0%data_2d)/domain%roughness_z0%data_2d)) / (lnz_atm_term**2)
                    lnz_atm_term=(karman/lnz_atm_term)**2
                endif

                ! note this is more or less just diagnostic and could be removed
                domain%longwave_up%data_2d = stefan_boltzmann * EMISS * domain%skin_temperature%data_2d**4
                RAINBL = domain%accumulated_precipitation%data_2d
                rain_bucket = domain%precipitation_bucket

            endif


            if (options%physics%landsurface > kLSM_BASIC) then
                ! accumulate soil moisture over the entire column
                domain%soil_totalmoisture%data_2d = domain%soil_water_content%data_3d(:,1,:) * DZS(1) * 1000
                do i = 2,num_soil_layers
                    domain%soil_totalmoisture%data_2d = domain%soil_totalmoisture%data_2d + domain%soil_water_content%data_3d(:,i,:) * DZS(i)
                enddo

                ! 2m Air T and Q are not well defined if Tskin is not coupled with the surface fluxes
                call surface_diagnostics(domain%sensible_heat%data_2d,    &
                                         QFX,                             &
                                         domain%skin_temperature%data_2d, &
                                         QSFC,                            &
                                         CHS2,                            &
                                         CQS2,                            &
                                         domain%temperature_2m%data_2d,   &
                                         domain%humidity_2m%data_2d,      &
                                         domain%surface_pressure%data_2d)
            endif
        endif
        if (options%physics%landsurface>0) then
            call apply_fluxes(domain, dt)
        endif

    end subroutine lsm
end module land_surface
