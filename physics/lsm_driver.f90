module land_surface
	use module_sf_noahdrv, only : lsm_noah, lsm_noah_init
	use module_lsm_basic,  only : lsm_basic
	use module_lsm_simple, only : lsm_simple, lsm_simple_init
	use data_structures
	
	implicit none
	! Noah LSM required variables.  Some of these should be stored in domain, but tested here for now
	integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
	integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
	integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions
	
	! LOTS of variables required by Noah, placed here temporarily to get it to compile, this may be where they stay...
	real,allocatable, dimension(:,:)    :: VEGFRA, CANWAT, SMSTAV, SMSTOT,SFCRUNOFF,UDRUNOFF,   &
										   SNOW,SNOWC,SNOWH, ACSNOW, ACSNOM, SNOALB,TSK, QFX,   &
										   QGH, GSW, ALBEDO, ALBBCK, ZNT, Z0, TMN, XICE, EMISS, &
										   EMBCK, QSFC, RAINBL, CHS, CHS2, CQS2, CPM, SR,       &
										   CHKLOWQ, LAI, QZ0, SHDMIN,SHDMAX,SNOTIME,SNOPCX,     &
										   POTEVP,SMCREL,RIB, NOAHRES,FLX4_2D,FVB_2D,FBUR_2D,   &
										   FGSN_2D
										   
	logical :: MYJ, FRPCPN,ua_phys,RDLAI2D,USEMONALB
	real,allocatable, dimension(:,:,:)  :: TSLB,SMOIS,SH2O
	real,allocatable, dimension(:)      :: Zs,DZs
	real :: ROVCP,XICE_THRESHOLD
	integer,allocatable, dimension(:,:) :: IVGTYP,ISLTYP
	integer :: ITIMESTEP
	
	character(len=MAXVARLENGTH) :: MMINLU
	logical :: FNDSOILW,FNDSNOWH,RDMAXALB
	integer :: num_soil_layers,ISURBAN,ISICE
	
contains
	subroutine allocate_noah_data(ime,jme,kme,num_soil_layers)
		implicit none
		integer, intent(in) :: ime,jme,kme,num_soil_layers
		integer :: i
		
		ITIMESTEP=1
		
		allocate(VEGFRA(ime,jme))
		VEGFRA=0.5
		allocate(CANWAT(ime,jme))
		CANWAT=0
		allocate(SMSTAV(ime,jme))
		SMSTAV=280
		allocate(SMSTOT(ime,jme))
		SMSTOT=0
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
		allocate(TSK(ime,jme))
		TSK=280.0

		allocate(QFX(ime,jme))
		QFX=0
		allocate(QGH(ime,jme))
		QGH=0
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
		allocate(TMN(ime,jme))
		TMN=0 !?
		allocate(XICE(ime,jme))
		XICE=0 !?
		allocate(EMISS(ime,jme))
		EMISS=0.95
		allocate(EMBCK(ime,jme))
		EMBCK=0.95
		allocate(QSFC(ime,jme))
		QSFC=0
		allocate(RAINBL(ime,jme))
		RAINBL=0
		allocate(CHS(ime,jme))
		CHS=0.1
		allocate(CHS2(ime,jme))
		CHS2=0.1
		allocate(CQS2(ime,jme))
		CQS2=0.1
		allocate(CPM(ime,jme))
		CPM=0
		allocate(SR(ime,jme))
		SR=0
		allocate(CHKLOWQ(ime,jme))
		CHKLOWQ=0
		allocate(LAI(ime,jme))
		LAI=3
		allocate(QZ0(ime,jme))
		QZ0=0
		
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
		allocate(SMCREL(ime,jme))
		SMCREL=0
		allocate(RIB(ime,jme))
		RIB=0
		allocate(NOAHRES(ime,jme))
		NOAHRES=0
		
		
		
		ROVCP=0.01
		XICE_THRESHOLD=0
		RDLAI2D=.false.
		USEMONALB=.false.
		MYJ=.false.
		FRPCPN=.false.
		ua_phys=.false.
		

		MMINLU="USGS"
		allocate(IVGTYP(ime,jme))
		IVGTYP=7 ! Grassland
		allocate(ISLTYP(ime,jme))
		ISLTYP=6 ! Loam

		
		allocate(TSLB(ime,num_soil_layers,jde))
		TSLB=280.0
		allocate(SMOIS(ime,num_soil_layers,jde))
		SMOIS=0.2
		allocate(SH2O(ime,num_soil_layers,jde))
		SH2O=SMOIS
		
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
		type(domain_type), intent(in) :: domain
		type(options_type),intent(in)    :: options
		
		
		write(*,*) "Initializing land surface model"
		
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

		    call LSM_NOAH_INIT(VEGFRA,SNOW,SNOWC,SNOWH,CANWAT,SMSTAV,    &
							SMSTOT, SFCRUNOFF,UDRUNOFF,ACSNOW,           &
							ACSNOM,IVGTYP,ISLTYP,TSLB,SMOIS,SH2O,ZS,DZS, &
							MMINLU,                                      &
							SNOALB, FNDSOILW, FNDSNOWH, RDMAXALB,        &
							num_soil_layers, options%restart,            &
							.True. ,                                     & ! allowed_to_read (e.g. soilparm.tbl)
							ids,ide, jds,jde, kds,kde,                   &
							ims,ime, jms,jme, kms,kme,                   &
							its,ite, jts,jte, kts,kte  )
		endif
		if (options%physics%landsurface==2) then
			call lsm_simple_init(domain,options)
		endif
		
	end subroutine lsm_init
	
	subroutine lsm(domain,options,dt)
		implicit none
		
		type(domain_type), intent(inout) :: domain
		type(options_type),intent(in)    :: options
		real, intent(in) :: dt

		if (options%physics%landsurface==1) then
			call lsm_basic(domain,options,dt)
		else if (options%physics%landsurface==2) then
			call lsm_simple(domain%th,domain%pii,domain%qv,domain%current_rain, domain%current_snow,domain%p, &
							domain%swdown,domain%lwdown, sqrt(domain%u(:,1,:)**2+domain%v(:,1,:)**2), &
							domain%sensible_heat, domain%latent_heat, domain%ground_heat, &
							domain%skin_t, domain%soil_t, domain%soil_vwc, domain%snow_swe, &
							options,dt)
		else if (options%physics%landsurface==3) then
! Variables I still haven't set (and may not know what they are)
! QGH,GSW,ALBEDO,ALBBCK,ZNT,Z0,TMN, XICE,EMISS,EMBCK, QSFC,RAINBL
! ITIMESTEP, CHS,CHS2,CQS2,CPM,ROVCP,SR,chklowq,lai,qz0,
! myj,frpcpn, SNOALB,SHDMIN,SHDMAX, SNOTIME

			call lsm_noah(domain%dz,domain%qv,domain%p,domain%th,TSK,  &
                  domain%sensible_heat,QFX,domain%latent_heat,domain%ground_heat, &
				  QGH,GSW,domain%swdown,domain%lwdown,SMSTAV,SMSTOT, &
                  SFCRUNOFF, UDRUNOFF,IVGTYP,ISLTYP,ISURBAN,ISICE,VEGFRA, &
                  ALBEDO,ALBBCK,ZNT,Z0,TMN,domain%landmask,XICE,EMISS,EMBCK,     &
                  SNOWC,QSFC,RAINBL,MMINLU,                     &
                  num_soil_layers,dt,DZS,ITIMESTEP,             &
                  SMOIS,TSLB,domain%snow_swe,CANWAT,            &
                  CHS,CHS2,CQS2,CPM,ROVCP,SR,chklowq,lai,qz0,   & !H
                  myj,frpcpn,                                   &
                  SH2O,SNOWH,                                   & !H
                  domain%u, domain%v,                           & !I
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
		endif
		
	end subroutine lsm
end module land_surface
