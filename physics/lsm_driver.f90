module land_surface
	use module_lsm_basic, only : lsm_basic
	use module_lsm_simple, only: lsm_simple, lsm_simple_init
	use data_structures
	
	implicit none
	
contains
	subroutine lsm_init(domain,options)
		implicit none
		type(domain_type), intent(in) :: domain
		type(options_type),intent(in)    :: options
		
		write(*,*) "Initializing land surface model"
		
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
							domain%sensible_heat, domain%latent_heat, &
							domain%skin_t, domain%soil_t, domain%soil_vwc, domain%snow_swe, &
							options,dt)
		endif
		
	end subroutine lsm
end module land_surface
