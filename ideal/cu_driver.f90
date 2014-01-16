module convection
	use data_structures
	use module_cu_tiedtke
	implicit none
	real,allocatable,dimension(:,:,:)::pii,p8,U3D,V3D,T3D,zed1,zed2,rho,                 &
									   RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,    &
					                   RUCUTEN,RVCUTEN
	logical,allocatable,dimension(:,:)::CU_ACT_FLAG
    REAL,allocatable,DIMENSION(:,:) ::XLAND, RAINCV, PRATEC
	real,dimension(44) :: ZNU = (/0.9975, 0.9915, 0.984, 0.975, 0.965, 0.9525, 0.9375, &
											 0.92, 0.9, 0.88, 0.86, 0.835, 0.805, 0.775, 0.745,   &
											 0.71, 0.67, 0.63, 0.59, 0.55, 0.51, 0.47, 0.43, 0.39,&
											 0.355, 0.325, 0.295, 0.27, 0.25, 0.23, 0.21, 0.19,   &
											 0.17, 0.15, 0.13, 0.11, 0.09100001, 0.074, 0.059,    &
											 0.046, 0.035, 0.025, 0.015, 0.005/)

contains
	subroutine init_convection(domain,options)
		type(domain_type), intent(in) :: domain
		type(options_type), intent(in) :: options
		integer::ids,ide,jds,jde,kds,kde

		ids=1
		ide=size(domain%qv,1)
		kds=1
		kde=size(domain%qv,2)
		jds=1
		jde=size(domain%qv,3)
		
		if (options%physics%convection==1) then
			write(*,*) "Initializing Tiedtke scheme"
			allocate(RTHCUTEN(ids:ide,kds:kde,jds:jde))
			RTHCUTEN=0
			allocate(RQVCUTEN(ids:ide,kds:kde,jds:jde))
			RQVCUTEN=0
			allocate(RQCCUTEN(ids:ide,kds:kde,jds:jde))
			RQCCUTEN=0
			allocate(RQICUTEN(ids:ide,kds:kde,jds:jde))
			RQICUTEN=0
			allocate(RUCUTEN(ids:ide,kds:kde,jds:jde))
			RUCUTEN=0
			allocate(RVCUTEN(ids:ide,kds:kde,jds:jde))
			RVCUTEN=0
			call tiedtkeinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,   &
			                     RUCUTEN,RVCUTEN,                   &
			                     .false.,1,1,0,                     &
			                     .true.,                            &
			                     ids, ide, jds, jde, kds, kde-1,      &
			                     ids, ide, jds, jde, kds, kde,      &
			                     ids, ide, jds, jde, kds, kde-1)
		 endif
		
	end subroutine init_convection
	
subroutine convect(domain,options,dt_in)
	type(domain_type), intent(inout) :: domain
	type(options_type), intent(in)   :: options
	real, intent(in) ::dt_in
	integer ::ids,ide,jds,jde,kds,kde,i,itimestep=1,STEPCU=1
	
	ids=1
	ide=size(domain%qv,1)
	kds=1
	kde=size(domain%qv,2)
	jds=1
	jde=size(domain%qv,3)
	
	if (.not.allocated(pii)) then
		allocate(pii(ids:ide,kds:kde,jds:jde))
		pii=1
		allocate(rho(ids:ide,kds:kde,jds:jde))
		rho=1
		allocate(zed1(ids:ide,kds:kde,jds:jde))
		zed1=0
		allocate(zed2(ids:ide,kds:kde,jds:jde))
		zed2=0
		allocate(p8(ids:ide,kds:kde,jds:jde))
		p8=0
		allocate(U3D(ids:ide,kds:kde,jds:jde))
		U3D=0
		allocate(V3D(ids:ide,kds:kde,jds:jde))
		V3D=0
		allocate(T3D(ids:ide,kds:kde,jds:jde))
		T3D=0
		allocate(CU_ACT_FLAG(ids:ide,jds:jde))
		CU_ACT_FLAG=.True.
		
		allocate(XLAND(ids:ide,jds:jde))
		allocate(RAINCV(ids:ide,jds:jde))
		allocate(PRATEC(ids:ide,jds:jde))
		
		where(domain%landmask==1) XLAND=1 ! land
		where(domain%landmask==0) XLAND=2 ! ocean (using LANDMASK as input)
		where(domain%landmask==2) XLAND=2 ! ocean (using XLAND as input)
		RAINCV=0
		PRATEC=0
	endif
	if (options%physics%convection==1) then
! 		used to convert potential temperature to sensible temperature

		!$omp parallel private(i) &
		!$omp firstprivate(ids,ide,jds,jde,kds,kde) &
		!$omp default(shared)
		!$omp do schedule(static)
		do i=ids,ide
			pii(i,:,:)=1.0/((100000.0/domain%p(i,:,:))**(R/cp))
			RAINCV(i,:)=0
			p8(i,:kde-1,:)=(domain%p(i,:kde-1,:)+domain%p(i,kds+1:,:))/2
			p8(i,kde,:)=p8(i,kde-1,:)+(p8(i,kde-1,:)-p8(i,kde-2,:))
			if (i>ids) then
				U3D(i,:,:)=(domain%u(i-1,:,:)+domain%u(i,:,:))/2
			endif
			V3D(i,:,:)=(domain%v(i,:,jds:jde)+domain%v(i,:,jds+1:jde+1))/2
			T3D(i,:,:)=domain%th(i,:,:)*pii(i,:,:)
			rho(i,:,:)=domain%p(i,:,:)/(R*T3D(i,:,:))
		enddo
		!$omp end do
		!$omp end parallel	
		
! 		write(*,*) "Entering Tiedtke"
		call CU_TIEDTKE(                                    &
                 dt_in,itimestep,STEPCU                            &
                ,RAINCV,PRATEC,domain%latent_heat,domain%sensible_heat,ZNU(kds:kde) &
                ,U3D,V3D,domain%w,T3D,domain%qv,domain%cloud,domain%ice,pii,rho &
                ,zed1,zed2 &
                ,domain%dz,p8,domain%p,XLAND,CU_ACT_FLAG                &
                ,ids,ide, jds,jde, kds,kde-1                      &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ids+1,ide-1, jds+1,jde-1, kds,kde-1					&
                ,RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN            &
                ,RUCUTEN, RVCUTEN                               &
                ,.True.    ,.True.    ,.True.    ,.True.    ,.True.       &
				)
		domain%qv=domain%qv+RQVCUTEN*dt_in
		domain%cloud=domain%cloud+RQCCUTEN*dt_in
		domain%th=domain%th+RTHCUTEN*dt_in
		domain%ice=domain%ice+RQICUTEN*dt_in
		domain%rain=domain%rain+RAINCV
		domain%crain=domain%crain+RAINCV
	endif
	
	
end subroutine convect	
end module convection