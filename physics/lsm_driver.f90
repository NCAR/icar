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
!!  sat_mr              - calculate saturated mixing ratio (should be moved to )
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
    use module_lsm_basic,    only : lsm_basic
    use module_lsm_simple,   only : lsm_simple, lsm_simple_init
    use module_water_simple, only : water_simple
    use io_routines,         only : io_write3d, io_write2d
    use data_structures

    implicit none

    private
    public :: lsm_init, lsm

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
                                           FGSN_2D, z_atm,lnz_atm_term,Ri,base_exchange_term

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

    real function sat_mr(t,p)
    ! Calculate the saturated mixing ratio at a temperature (K), pressure (Pa)
        implicit none
        real,intent(in) :: t,p
        real :: e_s,mr_s,a,b

        ! from http://www.dtic.mil/dtic/tr/fulltext/u2/778316.pdf
        !     Lowe, P.R. and J.M. Ficke., 1974: THE COMPUTATION OF SATURATION VAPOR PRESSURE
        !         Environmental Prediction Research Facility, Technical Paper No. 4-74
        ! which references:
        !     Murray, F. W., 1967: On the computation of saturation vapor pressure.
        !         Journal of Applied Meteorology, Vol. 6, pp. 203-204.
        ! Also notes a 6th order polynomial and look up table as viable options.
        if (t<freezing_threshold) then
            a=21.8745584
            b=7.66
        else
            a=17.2693882
            b=35.86
        endif
        e_s = 610.78* exp(a*(t-273.16)/(t-b)) !(Pa)

        ! alternate formulations
        ! Polynomial:
        ! e_s = ao + t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t))))) a0-6 defined separately for water and ice
        ! e_s = 611.2*exp(17.67*(t-273.15)/(t-29.65)) ! (Pa)
        ! from : http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf
        ! e_s = 611.0*10.0**(7.5*(t-273.15)/(t-35.45))

        ! enforce e_s < air pressure incase we are out on one edge of a polynomial
        ! I'm not sure this should ever be encounted anymore, but left in for now...
        if ((p-e_s)<=0) then
            e_s=p*0.99999
        endif
        ! e_s=min(e_s,p-SMALL_PRESSURE) ! this is harder to cover a reasonable range of pressure in single precision
        !from : http://www.srh.noaa.gov/images/epz/wxcalc/mixingRatio.pdf
        sat_mr=0.6219907*e_s/(p-e_s) !(kg/kg)
    end function sat_mr


    subroutine calc_exchange_coefficient(wind,tskin,airt,exchange_C)
        implicit none
        real, dimension(:,:),intent(inout) :: wind,tskin
        real, dimension(:,:,:),intent(inout) :: airt
        real,dimension(:,:),intent(inout) :: exchange_C

        ! Richardson number
        where(wind==0) wind=1e-10

        Ri = gravity/airt(:,1,:) * (airt(:,1,:)-tskin)*z_atm/(wind**2)

!         print*,"--------------------------------------------------"
!         print*, "Surface Richardson number"
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
        REAL, DIMENSION( :,: ), INTENT(IN)    ::  HFX, QFX, TSK, QSFC
        REAL, DIMENSION( :,: ), INTENT(INOUT) ::  Q2, T2
        REAL, DIMENSION( :,: ), INTENT(IN)    ::  PSFC, CHS2, CQS2
        integer :: i,j, nx,ny
        real :: rho

        !$omp parallel default(shared), private(nx,ny,i,j,rho)
        nx=size(HFX,1)
        ny=size(HFX,2)
        !$omp do
        do j=1,ny
            do i=1,nx
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
        type(domain_type), intent(inout) :: domain
        real, intent(in) :: dt
        integer :: nx,ny
        nx = size(domain%sensible_heat,1)
        ny = size(domain%sensible_heat,2)

        ! convert sensible heat flux to a temperature delta term
        ! (J/(s*m^2) * s / (J/(kg*K)) => kg*K/m^2) ... /((kg/m^3) * m) => K
        dTemp=(domain%sensible_heat(2:nx-1,2:ny-1)*dt/cp)  &
             / (domain%rho(2:nx-1,1,2:ny-1)*domain%dz_inter(2:nx-1,1,2:ny-1))
        ! add temperature delta and convert back to potential temperature
        domain%th(2:nx-1,1,2:ny-1)=domain%th(2:nx-1,1,2:ny-1) + (dTemp / domain%pii(2:nx-1,1,2:ny-1))

        ! convert latent heat flux to a mixing ratio tendancy term
        ! (J/(s*m^2) * s / (J/kg) => kg/m^2) ... / (kg/m^3 * m) => kg/kg
        lhdQV=(domain%latent_heat(2:nx-1,2:ny-1) / LH_vaporization * dt) &
             / (domain%rho(2:nx-1,1,2:ny-1) * domain%dz_inter(2:nx-1,1,2:ny-1))
        ! add water vapor in kg/kg
        domain%qv(2:nx-1,1,2:ny-1)=domain%qv(2:nx-1,1,2:ny-1) + lhdQV

        ! enforce some minimum water vapor content... just in case
        where(domain%qv<SMALL_QV) domain%qv=SMALL_QV

    end subroutine apply_fluxes

    subroutine allocate_noah_data(ime,jme,kme,num_soil_layers)
        implicit none
        integer, intent(in) :: ime,jme,kme,num_soil_layers
        integer :: i

        ITIMESTEP=1

        allocate(SMSTAV(ime,jme))
        SMSTAV=0.5 !average soil moisture available for transp (between SMCWLT and SMCMAX)
        allocate(SFCRUNOFF(ime,jme))
        SFCRUNOFF=0
        allocate(UDRUNOFF(ime,jme))
        UDRUNOFF=0
        allocate(SNOW(ime,jme))
        SNOW=0
        allocate(SNOWC(ime,jme))
        SNOWC=0
        allocate(SNOWH(ime,jme))
        SNOWH=0
        allocate(ACSNOW(ime,jme))
        ACSNOW=0
        allocate(ACSNOM(ime,jme))
        ACSNOM=0
        allocate(SNOALB(ime,jme))
        SNOALB=0.8

        allocate(QGH(ime,jme))
        QGH=0.02 ! saturated mixing ratio at ~20C
        allocate(GSW(ime,jme))
        GSW=0

        allocate(ALBEDO(ime,jme))
        ALBEDO=0.17
        allocate(ALBBCK(ime,jme))
        ALBBCK=0.17 !?
        allocate(XICE(ime,jme))
        XICE=0
        allocate(EMISS(ime,jme))
        EMISS=0.95
        allocate(EMBCK(ime,jme))
        EMBCK=0.95
        allocate(CPM(ime,jme))
        CPM=0
        allocate(SR(ime,jme))
        SR=0
        allocate(CHKLOWQ(ime,jme))
        CHKLOWQ=0
        allocate(LAI(ime,jme))
        LAI=3
        allocate(QZ0(ime,jme))
        QZ0=0 ! used to check for saturation? but only relevant if myj==True

        allocate(FLX4_2D(ime,jme))
        allocate(FVB_2D(ime,jme))
        allocate(FBUR_2D(ime,jme))
        allocate(FGSN_2D(ime,jme))

        allocate(SHDMIN(ime,jme))
        SHDMIN=0
        allocate(SHDMAX(ime,jme))
        SHDMAX=100
        allocate(SNOTIME(ime,jme))
        SNOTIME=0
        allocate(SNOPCX(ime,jme))
        SNOPCX=0
        allocate(POTEVP(ime,jme))
        POTEVP=0
        allocate(SMCREL(ime,num_soil_layers,jme))
        SMCREL=0
        allocate(RIB(ime,jme))
        RIB=0
        allocate(NOAHRES(ime,jme))
        NOAHRES=0
        allocate(VEGFRAC(ime,jme))
        VEGFRAC=50


        XICE_THRESHOLD=1
        RDLAI2D=.false.
        USEMONALB=.false.
        MYJ=.false.
        FRPCPN=.false. ! set this to true and calculate snow ratio to use microphysics based snow/rain partitioning
        ua_phys=.false.

        allocate(SH2O(ime,num_soil_layers,jme))
        SH2O=0.25

        allocate(Zs(num_soil_layers))
        allocate(DZs(num_soil_layers))
        DZs=[0.1,0.3,0.6,1.0]
        Zs(1)=DZs(1)/2
        do i=2,num_soil_layers
            Zs(i)=Zs(i-1) + DZs(i)/2 + DZs(i-1)/2
        end do

    end subroutine allocate_noah_data

    subroutine lsm_init(domain,options)
        implicit none
        type(domain_type), intent(inout) :: domain
        type(options_type),intent(in)    :: options
        integer :: i

        write(*,*) "Initializing LSM"

        exchange_term = 1

        ime=size(domain%th,1)
        jme=size(domain%th,3)

        allocate(dTemp(ime-2,jme-2))
        dTemp=0
        allocate(lhdQV(ime-2,jme-2))
        lhdQV=0
        allocate(Z0(ime,jme))
        Z0=0.01 ! this should get updated by the LSM
        allocate(QSFC(ime,jme))
        QSFC=domain%qv(:,1,:) ! this should get updated by the lsm
        allocate(Ri(ime,jme))
        Ri=0
        allocate(z_atm(ime,jme))
        z_atm=50
        allocate(lnz_atm_term(ime,jme))
        lnz_atm_term=0.1
        allocate(base_exchange_term(ime,jme))
        base_exchange_term=0.01
        allocate(QFX(ime,jme))
        QFX=0
        allocate(windspd(ime,jme))
        windspd=0
        ! NOTE, these fields have probably not been initialized yet...
        ! windspd=sqrt(domain%u10**2+domain%v10**2)
        allocate(CHS(ime,jme))
        CHS=0.01
        allocate(CHS2(ime,jme))
        CHS2=0.01
        allocate(CQS2(ime,jme))
        CQS2=0.01

        allocate(RAINBL(ime,jme))
        RAINBL=domain%rain  ! used to store last time step accumulated precip so that it can be subtracted from the current step
                            ! set to domain%rain incase this is a restart run and rain is non-zero to start
        allocate(rain_bucket(ime,jme))
        rain_bucket=domain%rain_bucket  ! used to store last time step accumulated precip so that it can be subtracted from the current step


        domain%T2m = domain%th(:,1,:) * domain%pii(:,1,:)
        domain%Q2m = domain%qv(:,1,:)

        if (options%physics%landsurface==kLSM_SIMPLE) then
            write(*,*) "    Simple LSM (may not work?)"
            stop "Simple LSM not settup, choose a different LSM options"
            call lsm_simple_init(domain,options)
        endif
        ! Noah Land Surface Model
        if (options%physics%landsurface==kLSM_NOAH) then
            write(*,*) "    Noah LSM"
            ids=1;jds=1;kds=1
            ims=1;jms=1;kms=1
            its=2;jts=2;kts=1

            ide=size(domain%p,1);jde=size(domain%p,3);kde=size(domain%p,2)
            ime=ide;jme=jde;kme=kde
            ite=ide-1;jte=jde-1;kte=kde

            num_soil_layers=4
            FNDSNOWH=.False. ! calculate SNOWH from SNOW
            FNDSOILW=.False. ! calculate SOILW (this parameter is ignored in LSM_NOAH_INIT)
            RDMAXALB=.False.

            ISURBAN = options%lsm_options%urban_category
            ISICE   = options%lsm_options%ice_category
            ISWATER = options%lsm_options%water_category
            MMINLU  = options%lsm_options%LU_Categories !"MODIFIED_IGBP_MODIS_NOAH"

            call allocate_noah_data(ime,jme,kme,num_soil_layers)

            if (options%lsm_options%monthly_vegfrac) then
                VEGFRAC=domain%vegfrac(:,:,domain%current_month)
            else
                VEGFRAC=domain%vegfrac(:,:,1)
            endif
            cur_vegmonth=domain%current_month

            ! save the canopy water in a temporary variable in case this is a restart run because lsm_init resets it to 0
            CQS2=domain%canopy_water
            ! prevents init from failing when processing water poitns that may have "soil_t"=0
            where(domain%soil_t<200) domain%soil_t=200
            where(domain%soil_vwc<0.0001) domain%soil_vwc=0.0001
            call LSM_NOAH_INIT(VEGFRAC,SNOW,SNOWC,SNOWH,domain%canopy_water,domain%soil_t,    &
                                domain%soil_vwc, SFCRUNOFF,UDRUNOFF,ACSNOW,  &
                                ACSNOM,domain%veg_type,domain%soil_type,     &
                                domain%soil_t,                               &
                                domain%soil_vwc,SH2O,ZS,DZS,                 &
                                MMINLU,                                      &
                                SNOALB, FNDSOILW, FNDSNOWH, RDMAXALB,        &
                                num_soil_layers, .False.,                    & ! nlayers, is_restart (can't yet)
                                .True. ,                                     & ! allowed_to_read (e.g. soilparm.tbl)
                                ids,ide, jds,jde, kds,kde,                   &
                                ims,ime, jms,jme, kms,kme,                   &
                                its,ite, jts,jte, kts,kte  )

            domain%canopy_water=CQS2
            CQS2=0.01
            where(domain%veg_type==ISWATER) domain%landmask=kLC_WATER ! ensure VEGTYPE (land cover) and land-sea mask are consistent
        endif

        ! defines the height of the middle of the first model level
        z_atm=domain%z(:,1,:)-domain%terrain
        lnz_atm_term = log((z_atm+Z0)/Z0)
        if (exchange_term==1) then
            base_exchange_term=(75*karman**2 * sqrt((z_atm+Z0)/Z0)) / (lnz_atm_term**2)
            lnz_atm_term=(karman/lnz_atm_term)**2
        endif

        update_interval=options%lsm_options%update_interval
        last_model_time=-999

    end subroutine lsm_init


    subroutine lsm(domain,options,dt,model_time)
        implicit none

        type(domain_type), intent(inout) :: domain
        type(options_type),intent(in)    :: options
        real, intent(in) :: dt
        real ::lsm_dt
        double precision, intent(in) :: model_time
        integer :: nx,ny,i,j

        if (last_model_time==-999) then
            last_model_time=model_time-update_interval
        endif
        nx=size(domain%qv,1)
        ny=size(domain%qv,3)
        if ((model_time-last_model_time)>=update_interval) then
            lsm_dt=model_time-last_model_time
            last_model_time=model_time

            ! exchange coefficients
            windspd=sqrt(domain%u10**2+domain%v10**2)
            if (exchange_term==1) then
                call calc_exchange_coefficient(windspd,domain%skin_t,domain%T,CHS)
            elseif (exchange_term==2) then
                call calc_mahrt_holtslag_exchange_coefficient(windspd,domain%skin_t,domain%T,domain%znt,CHS)
            endif
!             print*, CHS(128,103)

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
                call water_simple(domain%sst, domain%psfc, windspd, domain%ustar,  &
                                  domain%qv, domain%t,                             &
                                  domain%sensible_heat, domain%latent_heat,        &
                                  domain%z(:,1,:)-domain%terrain, Z0,              &
                                  domain%landmask, QSFC, QFX, domain%skin_t)
            endif

            where(windspd<1) windspd=1
            CHS = CHS * windspd
            CHS2=CHS
            CQS2=CHS

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
                !             QSFC(i,j)=max(domain%qv(i,1,j),0.5*sat_mr(domain%T2m(i,j),domain%psfc(i,j)))
                !         endif
                !     enddo
                ! enddo


            ! else
            if (options%physics%landsurface==kLSM_SIMPLE) then
                write(*,*) "--------------------------"
                stop "Simple LSM not implemented yet"
                call lsm_simple(domain%th,domain%pii,domain%qv,domain%current_rain, domain%current_snow,domain%p_inter, &
                                domain%swdown,domain%lwdown, sqrt(domain%u10**2+domain%v10**2), &
                                domain%sensible_heat, domain%latent_heat, domain%ground_heat, &
                                domain%skin_t, domain%soil_t, domain%soil_vwc, domain%snow_swe, &
                                options,lsm_dt)

            else if (options%physics%landsurface==kLSM_NOAH) then
                ! Call the Noah Land Surface Model

                ! 2m saturated mixing ratio
                do j=1,ny
                    do i=1,nx
                        if (domain%landmask(i,j)==kLC_LAND) then
                            QGH(i,j)=sat_mr(domain%T2m(i,j),domain%psfc(i,j))
                        endif
                    enddo
                enddo
                if (options%lsm_options%monthly_vegfrac) then
                    if (cur_vegmonth/=domain%current_month) then
                        VEGFRAC=domain%vegfrac(:,:,domain%current_month)
                        cur_vegmonth=domain%current_month
                    endif
                endif

                call lsm_noah(domain%dz_inter,domain%qv,domain%p_inter,domain%th*domain%pii,domain%skin_t,  &
                            domain%sensible_heat,QFX,domain%latent_heat,domain%ground_heat, &
                            QGH,GSW,domain%swdown,domain%lwdown,SMSTAV,domain%soil_totalmoisture, &
                            SFCRUNOFF, UDRUNOFF, &
                            domain%veg_type,domain%soil_type, &
                            ISURBAN,ISICE, &
                            VEGFRAC, &
                            ALBEDO,ALBBCK,domain%znt,Z0,domain%soil_tdeep,domain%landmask,XICE,EMISS,EMBCK,     &
                            SNOWC,QSFC,                                   &
                            (domain%rain-RAINBL)+(domain%rain_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE, &
                            MMINLU, &
                            num_soil_layers,lsm_dt,DZS,ITIMESTEP,         &
                            domain%soil_vwc,domain%soil_t,domain%snow_swe,&
                            domain%canopy_water,            &
                            CHS,CHS2,CQS2,CPM,ROVCP,SR,chklowq,lai,qz0,   & !H
                            myj,frpcpn,                                   &
                            SH2O,SNOWH,                                   & !H
                            domain%Um, domain%Vm,                         & !I
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

                ! now that znt has been updated, we need to recalculate terms
                lnz_atm_term = log((z_atm+domain%znt)/domain%znt)
                if (exchange_term==1) then
                    base_exchange_term=(75*karman**2 * sqrt((z_atm+domain%znt)/domain%znt)) / (lnz_atm_term**2)
                    lnz_atm_term=(karman/lnz_atm_term)**2
                endif

                ! note this is more or less just diagnostic and could be removed
                domain%lwup=stefan_boltzmann*EMISS*domain%skin_t**4
                RAINBL=domain%rain
                rain_bucket=domain%rain_bucket

                ! i=32
                ! j=82
                ! print*,"            ---------------------------"
                ! print*, "   soil_t ", "       skin_t ", "         T2m"
                ! print*, domain%soil_t(i,1,j), domain%skin_t(i,j), domain%T2m(i,j)
                ! print*, "   Tair ", "        sensible ", "        latent   ", "        CHS"
                ! print*, domain%T(i,1,j), domain%sensible_heat(i,j), domain%latent_heat(i,j), CHS(i,j)
                ! print*, "    SWD  ", "          LWD ", "         LWU ", "          G  "
                ! print*, domain%swdown(i,j), domain%lwdown(i,j), domain%lwup(i,j), domain%ground_heat(i,j)
                ! print*, "ln z term,        base exchange term        wind"
                ! print*, lnz_atm_term(i,j), base_exchange_term(i,j), windspd(i,j)

            endif
            ! 2m Air T and Q are not well defined if Tskin is not coupled with the surface fluxes
            if (options%physics%landsurface > kLSM_BASIC) then
                call surface_diagnostics(domain%sensible_heat, QFX, domain%skin_t, QSFC,  &
                                         CHS2, CQS2,domain%T2m, domain%Q2m, domain%psfc)
            endif
        endif
        if (options%physics%landsurface>0) then
            call apply_fluxes(domain,dt)
        endif

    end subroutine lsm
end module land_surface
