module land_surface
    use module_sf_noahdrv, only : lsm_noah, lsm_noah_init
    use module_lsm_basic,  only : lsm_basic
    use module_lsm_simple, only : lsm_simple, lsm_simple_init
    use io_routines,       only : io_write3d, io_write2d
	use output, 		   only : write_domain
    use data_structures
    
    implicit none
    ! Noah LSM required variables.  Some of these should be stored in domain, but tested here for now
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions
    
    ! LOTS of variables required by Noah, placed here temporarily to get it to compile, this may be where some stay.
    real,allocatable, dimension(:,:)    :: SMSTAV,SFCRUNOFF,UDRUNOFF,   &
                                           SNOW,SNOWC,SNOWH, ACSNOW, ACSNOM, SNOALB, QFX,   &
                                           QGH, GSW, ALBEDO, ALBBCK, ZNT, Z0, XICE, EMISS, &
                                           EMBCK, QSFC, RAINBL, CHS, CHS2, CQS2, CPM, SR,       &
                                           CHKLOWQ, LAI, QZ0, SHDMIN,SHDMAX,SNOTIME,SNOPCX,     &
                                           POTEVP,RIB, NOAHRES,FLX4_2D,FVB_2D,FBUR_2D,   &
                                           FGSN_2D, z_atm,lnz_atm_term,Ri,base_exchange_term, T2m
                                           
    logical :: MYJ, FRPCPN,ua_phys,RDLAI2D,USEMONALB
    real,allocatable, dimension(:,:,:)  :: SH2O,SMCREL
	real,allocatable, dimension(:,:)    :: dTemp,lhdQV
    real,allocatable, dimension(:)      :: Zs,DZs
    real :: ROVCP,XICE_THRESHOLD
    integer,allocatable, dimension(:,:) :: IVGTYP,ISLTYP
    integer :: ITIMESTEP
    
    real, parameter :: kappa=0.4
    real, parameter :: freezing_threshold=273.15
    real, parameter :: SMALL_PRESSURE=1e-10
    real, parameter :: MAX_EXCHANGE_C = 0.5
    
    character(len=MAXVARLENGTH) :: MMINLU
    logical :: FNDSOILW,FNDSNOWH,RDMAXALB
    integer :: num_soil_layers,ISURBAN,ISICE,counter,steps_per_lsm_step
    real*8  :: last_model_time
    
contains
    real function sat_mr(t,p)
    ! Calculate the saturated mixing ratio at a temperature (K), pressure (Pa)
        implicit none
        real,intent(in) :: t,p
        real :: e_s,mr_s,a,b

!       from http://www.dtic.mil/dtic/tr/fulltext/u2/778316.pdf
!           Lowe, P.R. and J.M. Ficke., 1974: THE COMPUTATION OF SATURATION VAPOR PRESSURE 
!               Environmental Prediction Research Facility, Technical Paper No. 4-74
!       which references:
!           Murray, F. W., 1967: On the computation of saturation vapor pressure. 
!               Journal of Applied Meteorology, Vol. 6, pp. 203-204.
!       Also notes a 6th order polynomial and look up table as viable options. 
        if (t<freezing_threshold) then
            a=21.8745584
            b=7.66
        else
            a=17.2693882
            b=35.86
        endif
        e_s = 610.78* exp(a*(t-273.16)/(t-b)) !(Pa)

!       alternate formulations
!       Polynomial:
!       e_s = ao + t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t))))) a0-6 defined separately for water and ice
!       e_s = 611.2*exp(17.67*(t-273.15)/(t-29.65)) ! (Pa)
        !from : http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf
!       e_s = 611.0*10.0**(7.5*(t-273.15)/(t-35.45))
        
        
        e_s=min(e_s,p-SMALL_PRESSURE)
        !from : http://www.srh.noaa.gov/images/epz/wxcalc/mixingRatio.pdf
        sat_mr=0.6219907*e_s/(p-e_s) !(kg/kg)
    end function sat_mr
    
    subroutine calc_exchange_coefficient(wind,tskin,airt,exchange_C)
        implicit none
        real, dimension(:,:),intent(in) :: wind,tskin,airt
        real,dimension(:,:),intent(inout) :: exchange_C
        
        Ri = g/airt * (airt-tskin)*z_atm/wind**2
        
        where(Ri<0)  exchange_C=lnz_atm_term * (1.0-(15.0*Ri)/(1.0+(base_exchange_term * sqrt((-1.0)*Ri))))
        where(Ri>=0) exchange_C=lnz_atm_term * 1.0/((1.0+15.0*Ri)*sqrt(1.0+5.0*Ri))
        where(exchange_C>MAX_EXCHANGE_C) exchange_C=MAX_EXCHANGE_C
    end subroutine calc_exchange_coefficient
    
    
	subroutine apply_fluxes(domain,dt)
		implicit none
		type(domain_type), intent(inout) :: domain
		real, intent(in) :: dt
		integer :: nx,ny
		nx=ime
		ny=jme
		! convert sensible heat flux to a temperature delta term
        ! J/s/m^2 * s / J/(kg*K) => kg*K/m^2 ... /((kg/m^3) * m) => K
        dTemp=(domain%sensible_heat(2:nx-1,2:ny-1)*dt/cp)  &
			 /(domain%rho(2:nx-1,1,2:ny-1)*domain%dz(2:nx-1,1,2:ny-1))
		! add temperature delta and convert back to potential temperature
        domain%th(2:nx-1,1,2:ny-1)=domain%th(2:nx-1,1,2:ny-1)+dTemp/domain%pii(2:nx-1,1,2:ny-1)
		
		! convert latent heat flux to a mixing ratio tendancy term
        ! J/s/m^2 * s / J/kg => kg/m^2 ... / (kg/m^3 * m) => kg/kg
		lhdQV=(domain%latent_heat(2:nx-1,2:ny-1)/LH_vaporization*dt) / (domain%rho(2:nx-1,1,2:ny-1)*domain%dz(2:nx-1,1,2:ny-1))
		domain%qv(2:nx-1,1,2:ny-1)=domain%qv(2:nx-1,1,2:ny-1)+lhdQV
		where(domain%qv<SMALL_PRESSURE) domain%qv=SMALL_PRESSURE
		
	end subroutine apply_fluxes
	
    subroutine allocate_noah_data(ime,jme,kme,num_soil_layers)
        implicit none
        integer, intent(in) :: ime,jme,kme,num_soil_layers
        integer :: i
        
        ITIMESTEP=1
        
        allocate(T2m(ime,jme))
        T2m=270
        allocate(Ri(ime,jme))
        Ri=0
        allocate(z_atm(ime,jme))
        z_atm=50
        allocate(lnz_atm_term(ime,jme))
        lnz_atm_term=0.1
        allocate(base_exchange_term(ime,jme))
        base_exchange_term=0.01
        
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

        allocate(QFX(ime,jme))
        QFX=0
        allocate(QGH(ime,jme))
        QGH=0.02 ! saturated mixing ratio at ~20C
        allocate(GSW(ime,jme))
        GSW=0

        allocate(ALBEDO(ime,jme))
        ALBEDO=0.17
        allocate(ALBBCK(ime,jme))
        ALBBCK=0.17 !?
        allocate(ZNT(ime,jme))
        ZNT=0.01 !?
        allocate(Z0(ime,jme))
        Z0=0.01 !?
        allocate(XICE(ime,jme))
        XICE=0
        allocate(EMISS(ime,jme))
        EMISS=0.95
        allocate(EMBCK(ime,jme))
        EMBCK=0.95
        allocate(QSFC(ime,jme))
        QSFC=0
        allocate(RAINBL(ime,jme))
        RAINBL=0 ! used to store last time step accumulated precip so that it can be subtracted from the current step
        allocate(CHS(ime,jme))
        CHS=0.01
        allocate(CHS2(ime,jme))
        CHS2=0.01
        allocate(CQS2(ime,jme))
        CQS2=0.01
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
        SHDMAX=0
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
        
        
        
        ROVCP=R/cp
        XICE_THRESHOLD=1
        RDLAI2D=.false.
        USEMONALB=.false.
        MYJ=.false.
        FRPCPN=.false. ! set this to true and calculate snow ratio to use microphysics based snow/rain partitioning
        ua_phys=.false.
        

        MMINLU="USGS"
		
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
        
        write(*,*) "Initializing land surface model"
        
		ime=size(domain%th,1)
		jme=size(domain%th,3)
		allocate(dTemp(ime-2,jme-2))
		dTemp=0
		allocate(lhdQV(ime-2,jme-2))
		lhdQV=0
        ! Noah Land Surface Model
        if (options%physics%landsurface==3) then
            
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
            ISURBAN=0
            ISICE=0
            call allocate_noah_data(ime,jme,kme,num_soil_layers)
            
            call LSM_NOAH_INIT(domain%vegfrac,SNOW,SNOWC,SNOWH,domain%canopy_water,domain%soil_t,    &
	                            domain%soil_vwc, SFCRUNOFF,UDRUNOFF,ACSNOW,  &
	                            ACSNOM,domain%veg_type,domain%soil_type,	 &
								domain%soil_t, 								 &
								domain%soil_vwc,SH2O,ZS,DZS, 				 &
	                            MMINLU,                                      &
	                            SNOALB, FNDSOILW, FNDSNOWH, RDMAXALB,        &
	                            num_soil_layers, .False.,                    & ! nlayers, is_restart (can't yet)
	                            .True. ,                                     & ! allowed_to_read (e.g. soilparm.tbl)
	                            ids,ide, jds,jde, kds,kde,                   &
	                            ims,ime, jms,jme, kms,kme,                   &
	                            its,ite, jts,jte, kts,kte  )
							
            z_atm=domain%z(:,1,:)
            lnz_atm_term = log((z_atm+Z0)/Z0)
            base_exchange_term=(75*kappa**2 * sqrt((z_atm+Z0)/Z0)) / (lnz_atm_term**2)
            lnz_atm_term=(kappa/lnz_atm_term)**2
			where(domain%veg_type==16) domain%landmask=2 ! ensure VEGTYPE (land cover) and land-sea mask are consistent
        endif
        if (options%physics%landsurface==2) then
            call lsm_simple_init(domain,options)
        endif
        counter=0
        steps_per_lsm_step=10
        last_model_time=-999
        
    end subroutine lsm_init
    
	
    subroutine lsm(domain,options,dt,model_time)
        implicit none
        
        type(domain_type), intent(inout) :: domain
        type(options_type),intent(in)    :: options
        real, intent(in) :: dt
        real ::lsm_dt
        real*8, intent(in) :: model_time
        integer :: nx,ny,i,j
        
        if (last_model_time==-999) then
            last_model_time=model_time-dt
        endif
        nx=size(domain%qv,1)
        ny=size(domain%qv,3)
        counter=counter+1
        if (counter>steps_per_lsm_step) then
            counter=0
            lsm_dt=model_time-last_model_time
            last_model_time=model_time
            
            if (options%physics%landsurface==1) then
                call lsm_basic(domain,options,lsm_dt)
            else if (options%physics%landsurface==2) then
                call lsm_simple(domain%th,domain%pii,domain%qv,domain%current_rain, domain%current_snow,domain%p, &
                                domain%swdown,domain%lwdown, sqrt(domain%u(1:ny,1,1:nx)**2+domain%v(1:ny,1,1:nx)**2), &
                                domain%sensible_heat, domain%latent_heat, domain%ground_heat, &
                                domain%skin_t, domain%soil_t, domain%soil_vwc, domain%snow_swe, &
                                options,lsm_dt)
                                
            else if (options%physics%landsurface==3) then
                ! Call the Noah Land Surface Model
                
				! it would be better to call the MYJ SFC scheme here...
                T2m=domain%th(:,1,:)*domain%pii(:,1,:)
                ! 2m saturated mixing ratio
                do j=1,ny
                    do i=1,nx
                        QGH(i,j)=sat_mr(T2m(i,j),domain%p(i,1,j))
                    enddo
                enddo
                ! shortwave down
                ! GSW=domain%swdown   ! This does not actually get used in Noah
                
				! exchange coefficients
                call calc_exchange_coefficient(sqrt(domain%u(1:nx,1,1:ny)**2+domain%v(1:nx,1,1:ny)**2),domain%skin_t,T2m,CHS)
                CHS2=CHS
                CQS2=CHS
				
                call lsm_noah(domain%dz,domain%qv,domain%p,domain%th*domain%pii,domain%skin_t,  &
							domain%sensible_heat,QFX,domain%latent_heat,domain%ground_heat, &
							QGH,GSW,domain%swdown,domain%lwdown,SMSTAV,domain%soil_totalmoisture, &
							SFCRUNOFF, UDRUNOFF, &
							domain%veg_type,domain%soil_type, &
							ISURBAN,ISICE, &
							domain%vegfrac, &
							ALBEDO,ALBBCK,ZNT,Z0,domain%soil_tdeep,domain%landmask,XICE,EMISS,EMBCK,     &
							SNOWC,QSFC,domain%rain-RAINBL,MMINLU,         &
							num_soil_layers,lsm_dt,DZS,ITIMESTEP,         &
							domain%soil_vwc,domain%soil_t,domain%snow_swe,&
							domain%canopy_water,            &
							CHS,CHS2,CQS2,CPM,ROVCP,SR,chklowq,lai,qz0,   & !H
							myj,frpcpn,                                   &
							SH2O,SNOWH,                                   & !H
							domain%u(1:nx,:,1:ny), domain%v(1:nx,:,1:ny), & !I
							SNOALB,SHDMIN,SHDMAX,                         & !I
							SNOTIME,                                      & !?
							ACSNOM,ACSNOW,                                & !O
							SNOPCX,                                       & !O
							POTEVP,                                       & !O
							SMCREL,                                       & !O
							XICE_THRESHOLD,                               &
							RDLAI2D,USEMONALB,                            &
							RIB,                                          & !?
							NOAHRES,                                      &
							ua_phys,flx4_2d,fvb_2d,fbur_2d,fgsn_2d,       & ! Noah UA changes
							ids,ide, jds,jde, kds,kde,                    &
							ims,ime, jms,jme, kms,kme,                    &
							its,ite, jts,jte, kts,kte)
                
				! note this is more or less just diagnostic and could be removed
				domain%lwup=stefan_boltzmann*EMISS*domain%skin_t**4
                RAINBL=domain%rain
				call apply_fluxes(domain,lsm_dt)
            endif
        endif
        
    end subroutine lsm
end module land_surface
