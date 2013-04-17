module microphysics
	use data_structures
	use module_mp_thompson
	implicit none
	real, parameter :: LH_vaporization=2260000.0 ! J/kg
	real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
	real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
	real, parameter :: g=9.81 ! gravity m/s^2

contains
	subroutine mp_init(physics_level)
		integer, intent(in)::physics_level
		call thompson_init()
	end subroutine mp_init
	
	subroutine mp(domain,options,dt_in)
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		real,intent(in)::dt_in
		real,allocatable,dimension(:,:,:)::pii,SR! ,qv,qc,qr,qi,qs,qg, &
! 					ni,nr,th,p,dz
! 		real,allocatable,dimension(:,:)::rain,snow,rainc,snowc,graupel,graupelc
		integer ::ids,ide,jds,jde,kds,kde,itimestep=1
		
		ids=1
		ide=size(domain%qv,1)
		kds=1
		kde=size(domain%qv,2)
		jds=1
		jde=size(domain%qv,3)
		
		allocate(pii(ids:ide,kds:kde,jds:jde))
		allocate(SR(ids:ide,kds:kde,jds:jde))
		SR=0
		pii=1.0/((100000.0/domain%p)**(R/cp))
! 		write(*,"(20F6.3)") pii(4,0,20:40)
		call mp_gt_driver(domain%qv, domain%cloud, domain%qrain, domain%ice, &
		                domain%qsnow, domain%qgrau, domain%nice, domain%nrain, &
						domain%th, pii, domain%p, domain%dz, dt_in, itimestep, &
						domain%rain, domain%rain, &
						domain%snow, domain%snow, &
						domain%graupel, domain%graupel, &
						SR, &
						ids,ide, jds,jde, kds,kde, &    ! domain dims
						ids,ide, jds,jde, kds,kde, &    ! memory dims
						ids+1,ide-1, jds+1,jde-1, kds,kde)      ! tile dims
		deallocate(pii,SR)

! 		allocate(qv(ids:ide,jds:jde,1))
! 		allocate(qc(ids:ide,jds:jde,1))
! 		allocate(qr(ids:ide,jds:jde,1))
! 		allocate(qi(ids:ide,jds:jde,1))
! 		allocate(qs(ids:ide,jds:jde,1))
! 		allocate(qg(ids:ide,jds:jde,1))
! 		allocate(ni(ids:ide,jds:jde,1))
! 		allocate(nr(ids:ide,jds:jde,1))
! 		allocate(th(ids:ide,jds:jde,1))
! 		allocate(p(ids:ide,jds:jde,1))
! 		allocate(dz(ids:ide,jds:jde,1))
! 		allocate(rain(ids:ide,1))
! 		allocate(rainc(ids:ide,1))
! 		allocate(snow(ids:ide,1))
! 		allocate(snowc(ids:ide,1))
! 		allocate(graupel(ids:ide,1))
! 		allocate(graupelc(ids:ide,1))
! 		
! 		qv(:,:,1)=domain%qv
! 		qc(:,:,1)=domain%cloud
! 		qr(:,:,1)=domain%qrain
! 		qi(:,:,1)=domain%ice
! 		qs(:,:,1)=domain%qsnow
! 		qg(:,:,1)=domain%qgrau
! 		ni(:,:,1)=domain%nice
! 		nr(:,:,1)=domain%nrain
! 		th(:,:,1)=domain%th
! 		p(:,:,1)=domain%p
! 		dz(:,:,1)=domain%dz
! 		rain(:,1)=domain%rain
! 		rainc=0
! 		snow(:,1)=domain%snow
! 		snowc=0
! 		graupel=0
! 		graupelc=0
! 		call mp_gt_driver(qv, qc, qr, qi, qs,qg,ni,nr, &
! 						th, pii, p, dz, dt_in, itimestep, &
! 						rain, rainc, &
! 						snow, snowc, &
! 						graupel, graupelc, &
! 						SR, &
! 						ids,ide, jds,jde, kds,kde, &    ! domain dims
! 						ids,ide, jds,jde, kds,kde, &    ! memory dims
! 						ids,ide, jds,jde, kds,kde)      ! tile dims
! 
! 
! 		domain%qv=qv(:,:,1)
! 		domain%cloud=qc(:,:,1)
! 		domain%qrain=qr(:,:,1)
! 		domain%ice=qi(:,:,1)
! 		domain%qsnow=qs(:,:,1)
! 		domain%qgrau=qg(:,:,1)
! 		domain%nice=ni(:,:,1)
! 		domain%nrain=nr(:,:,1)
! 		domain%th=th(:,:,1)
! 		domain%p=p(:,:,1)
! 		domain%dz=dz(:,:,1)
! 		domain%rain=rain(:,1)
! 		domain%snow=snow(:,1)
! 		domain%graupel=graupel(:,1)
! 
! 		deallocate(qv,qc,qr,qi,qs,qg,ni,nr,th,p,dz,rain,rainc,snow,snowc,graupel,graupelc)
						
	end subroutine mp
end module microphysics