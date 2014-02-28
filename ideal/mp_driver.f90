module microphysics
	use data_structures
	use module_mp_thompson, only: mp_gt_driver,thompson_init
	implicit none
! 	these are now defined in data_structures.f90
! 	real, parameter :: LH_vaporization=2260000.0 ! J/kg
! 	real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
! 	real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
! 	real, parameter :: g=9.81 ! gravity m/s^2

	real,allocatable,dimension(:,:,:)::SR

contains
	subroutine mp_init(options)
		implicit none
		type(options_type), intent(in)::options
		if (options%physics%microphysics==1) then
			call thompson_init()
		endif
	end subroutine mp_init
	
	subroutine mp(domain,options,dt_in)
		implicit none
		type(domain_type),intent(inout)::domain
		type(options_type),intent(in)::options
		real,intent(in)::dt_in
		integer ::ids,ide,jds,jde,kds,kde,itimestep=1
		
		ids=1
		ide=size(domain%qv,1)
		kds=1
		kde=size(domain%qv,2)
		jds=1
		jde=size(domain%qv,3)

		if (.not.allocated(SR)) then
			allocate(SR(ids:ide,kds:kde,jds:jde))
	! 		snow rain ratio
			SR=0
		endif
		if (options%physics%microphysics==1) then
	! 		used to convert potential temperature to sensible temperature
			
			call mp_gt_driver(domain%qv, domain%cloud, domain%qrain, domain%ice, &
			                domain%qsnow, domain%qgrau, domain%nice, domain%nrain, &
							domain%th, domain%pii, domain%p, domain%dz, dt_in, itimestep, &
							domain%rain, domain%rain, &
							domain%snow, domain%snow, &
							domain%graupel, domain%graupel, &
							SR, &
							ids,ide, jds,jde, kds,kde, &    ! domain dims
							ids,ide, jds,jde, kds,kde, &    ! memory dims
							ids+1,ide-1, jds+1,jde-1, kds,kde)      ! tile dims
		endif
						
	end subroutine mp
	
	subroutine mp_finish()
		implicit none
		if (allocated(SR)) then
			deallocate(SR)
		endif
	end subroutine mp_finish
end module microphysics