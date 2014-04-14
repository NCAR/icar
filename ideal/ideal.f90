program ideal
	use init                ! init_*  routines
! 	use boundary_conditions ! bc_*    routines, constant in an ideal simulation? 
	use data_structures     ! *_type  types
	use microphysics        ! mp*     routines
	use wind                ! update_winds
	use advection           ! advect
	use output              ! write_* routines
	
	implicit none
	type(options_type) :: options
	type(domain_type)  :: domain
	type(bc_type)  :: boundary
	real :: dt
	integer::i
	
	call init_model("ideal_options.namelist",options,domain,boundary)
	call mp_init(options%physics%microphysics)
	
	! linear theory winds if specified else just balance dudx+dvdy+dwdz = 0
	call update_winds(domain,options)
	
	! courant condition for 3D advection... should make 3 x 1D to maximize dt?
	dt=options%dx/maxval(domain%u)/3.0
	write(*,*) "dt = ",dt
	
	do i=0,options%ntimesteps
		call advect(domain,options,dt,options%dx)
		call mp(domain,options,dt)
! 		call lsm(domain)
! 		call pbl(domain)
! 		call radiation(domain)
		call write_domain(domain,options,i)
	end do
	
end program ideal
