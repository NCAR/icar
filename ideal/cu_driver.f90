module convection
	use data_structures
	use module_cu_tiedtke
	implicit none
	real,allocatable,dimension(:,:,:)::pii,p8,U3D,V3D,zed,rho
	logical,allocatable,dimension(:,:)::CU_ACT_FLAG
    REAL,allocatable,DIMENSION(:,:) ::XLAND, RAINCV, PRATEC

contains
subroutine convect(domain,options,dt_in)
	type(domain_type), intent(inout) :: domain
	type(options_type), intent(in)   :: options
	real, intent(in) ::dt_in
	integer ::ids,ide,jds,jde,kds,kde,itimestep=1,STEPCU=1
	
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
		allocate(zed(ids:ide,kds:kde,jds:jde))
		zed=0
		allocate(p8(ids:ide,kds:kde,jds:jde))
		p8=0
		allocate(U3D(ids:ide,kds:kde,jds:jde))
		U3D=0
		allocate(V3D(ids:ide,kds:kde,jds:jde))
		V3D=0
		allocate(CU_ACT_FLAG(ids:ide,jds:jde))
		CU_ACT_FLAG=.True.
		
		allocate(XLAND(ids:ids,jds:jde))
		allocate(RAINCV(ids:ids,jds:jde))
		allocate(PRATEC(ids:ids,jds:jde))
		XLAND=1
		RAINCV=0
		PRATEC=0
	endif
	if (options%physics%convection==1) then
! 		used to convert potential temperature to sensible temperature
		pii=1.0/((100000.0/domain%p)**(R/cp))
		RAINCV=0
		p8(:,:kde-1,:)=(domain%p(:,:kde-1,:)+domain%p(:,kds+1:,:))/2
		p8(:,kde,:)=p8(:,kde-1,:)+(p8(:,kde-1,:)-p8(:,kde-2,:))
		U3D=(domain%u(ids:ide,:,:)+domain%u(ids+1:ide+1,:,:))/2
		V3D=(domain%v(:,:,jds:jde)+domain%v(:,:,jds+1:jde+1))/2
		rho=domain%p/(R*domain%th*pii)
		call CU_TIEDTKE(                                    &
                 dt_in,itimestep,STEPCU                            &
                ,RAINCV,PRATEC,domain%latent_heat,domain%sensible_heat,domain%z(1,:,1) &
                ,U3D,V3D,domain%w,domain%th,domain%qv,domain%cloud,domain%ice,pii,rho &
                ,zed,zed &
                ,domain%dz,p8,domain%p,XLAND,CU_ACT_FLAG                &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ids,ide, jds,jde, kds,kde)

		domain%rain=domain%rain+RAINCV
	endif
	
	
end subroutine convect	
end module convection