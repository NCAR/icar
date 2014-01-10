program real
	use init                ! init_model
	use boundary_conditions ! bc_*    routines, constant in an ideal simulation? 
	use data_structures     ! *_type  types
	use microphysics        ! mp_init
	use output              ! write_domain
	use time_step			! step
	use convection
	
	implicit none
	type(options_type) :: options
	type(domain_type)  :: domain
	type(bc_type)      :: boundary
	integer::i,nx,ny,start_point
	real*8::model_time,next_output
	
! 	initialize model including options, terrain, lat, lon data. 
	call init_model("real_options.namelist",options,domain,boundary)
! 	initialize microphysics code (e.g. compute look up tables in Thompson et al)
	write(*,*) "Initializing microphysics"
! 	write(*,*) "WARNING: NOT Initializing microphysics"
	call mp_init(options%physics%microphysics) !this could easily be moved to init_model...
	call init_convection(domain,options)
! 	read initial conditions from the boundary file
	write(*,*) "Initializing Boundary conditions"
	call bc_init(domain,boundary,options)
	if (options%restart) then
		start_point=options%restart_step
	else
		start_point=1
	endif
	model_time=(start_point-1)*options%in_dt
	next_output=model_time+options%out_dt
	
! 	note that a timestep here is an IO timestep O(1hr), not a physics timestep O(20s)
	do i=start_point,options%ntimesteps
		write(*,*) "Timestep:", i, "  of ", options%ntimesteps
		write(*,*) "  Model time=",model_time/3600.0,"hrs"
! 		update boundary conditions (dXdt variables)
		call bc_update(domain,boundary,options)
! 		this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
		call step(domain,options,boundary,model_time,next_output)
! 		finally write the output for this timestep
!       this is now handled internal to step to make sub-Input timesteps in output easy. 
! 		call write_domain(domain,options,i)
	end do
	
end program real

