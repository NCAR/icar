module radiation
	use module_ra_simple, only: ra_simple, ra_simple_init
	use data_structures
	
	implicit none
	
contains
	subroutine radiation_init(domain,options)
		type(domain_type), intent(in) :: domain
		type(options_type),intent(in)    :: options
		
		write(*,*) "Initializing radiation"
		
		if (options%physics%radiation==2) then
			call ra_simple_init(domain,options)
		endif
		
	end subroutine radiation_init
	
	subroutine rad(domain,options,date,dt)
		implicit none
		
		type(domain_type), intent(inout) :: domain
		type(options_type),intent(in)    :: options
		double precision, intent(in) :: date
		real, intent(in) :: dt

		if (options%physics%radiation==2) then
			call ra_simple(domain%th,domain%pii,domain%qv,domain%cloud+domain%ice,domain%qsnow,&
						domain%qrain,domain%p,domain%swdown,domain%lwdown,domain%cloudfrac,&
						domain%lat,domain%lon,date,options,dt)
		endif
		
	end subroutine rad
end module radiation
