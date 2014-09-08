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
	real,allocatable, dimension(:,:)    :: VEGFRA, CANWAT, SMSTAV, SMSTOT,SFCRUNOFF,UDRUNOFF, &
										   SNOW,SNOWC,SNOWH, ACSNOM, SNOALB
	real,allocatable, dimension(:,:,:)  :: TSLB,SMOIS,SH2O
	real,allocatable, dimension(:)      :: Zs,DZs
	integer,allocatable, dimension(:,:) :: IVGTYP,ISLTYP
	
	character(len=MAXVARLENGTH) :: MMINLU
	logical :: FNDSOILW,FNDSNOWH,RDMAXALB
	integer :: num_soil_layers,ISURBAN,ISICE
	
contains
	subroutine allocate_noah_data(ime,jme,kme,num_soil_layers)
		implicit none
		integer, intent(in) :: ime,jme,kme,num_soil_layers
		integer :: i
		
		allocate(VEGFRA(ime,jme))
		VEGFRAC=0.5
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
		allocate(ACSNOM(ime,jme))
		ACSNOM=0
		allocate(SNOWALB(ime,jme))
		SNOWALB=0.8

		MMINLU="USGS"
		allocate(IVGTYP(ime,jme))
		IVGTYP=7 ! Grassland
		allocate(ISLTYP(ime,jme))
		ISLTYP=6 ! Loam

		
		allocate(TSLB(ime,num_soil_layers,jde))
		TSLB=280
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
							SMSTOT, SFCRUNOFF,UDRUNOFF,ACSNOW,        &
							ACSNOM,IVGTYP,ISLTYP,TSLB,SMOIS,SH2O,ZS,DZS, &
							MMINLU,                                   &
							SNOALB, FNDSOILW, FNDSNOWH, RDMAXALB,     &
							num_soil_layers, options%restart,         &
							.True ,                                   & ! allowed_to_read (e.g. soilparm.tbl)
							ids,ide, jds,jde, kds,kde,                &
							ims,ime, jms,jme, kms,kme,                &
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
                  SMOIS,TSLB,domain%swe,CANWAT,                 &
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
                  its,ite, jts,jte, kts,kte,                    &
                  sf_urban_physics,                             &
                  CMR_SFCDIF,CHR_SFCDIF,CMC_SFCDIF,CHC_SFCDIF)

		endif
		
	end subroutine lsm
end module land_surface
