!-----------------------------------------
!
! Main Program
!
! Initialize options and memory in init_model
! Read initial conditions in bc_init (from a restart file if requested)
! initialize physics packages in init_physics (e.g. tiedke and thompson if used)
! If this run is a restart run, then set start to the restart timestep
!      in otherwords, ntimesteps is the number of BC updates from the beginning of the entire model 
!      run, not just from the begining of this restart run
! calculate model time in seconds based on the time between BC updates (in_dt)
! Calculate the next model output time from current model time + output time delta (out_dt)
!
! Finally, loop until ntimesteps are reached updating boundary conditions and stepping the model forward
!
!-----------------------------------------
program real
	use init                ! init_model, init_physics             Initialize model (not initial conditions)
	use boundary_conditions ! bc_init,bc_update                    Boundary and initial conditions
	use data_structures     ! *_type datatypes                     Data-types and physical "constants"
	use output              ! write_domain                         Used to output initial model state
	use time_step			! step                                 Advance the model forward in time
	
	implicit none
	type(options_type) :: options
	type(domain_type)  :: domain
	type(bc_type)      :: boundary
	integer            :: i,nx,ny,start_point
	double precision   :: model_time,next_output
		
!-----------------------------------------
!  Model Initialization
!
! 	initialize model including options, terrain, lat, lon data. 
	call init_model(options,domain,boundary)
! 	read initial conditions from the boundary file
	write(*,*) "Initializing Boundary conditions"
	call bc_init(domain,boundary,options)
	write(*,*) maxval(domain%u),minval(domain%u)
	
	write(*,*) "Initializing Physics packages"
	call init_physics(options,domain)
	
	if (options%restart) then
		start_point=options%restart_step
	else
		start_point=1
	endif
	!note, startpoint at the beginning is 1, but we want model time at this point to be 0, thus start_point-1
	model_time=(start_point-1)*options%in_dt
	next_output=model_time+options%out_dt
	call bc_update(domain,boundary,options)
	
	! write the initial state of the model (primarily useful for debugging)
	if (.not.options%restart) then
		call write_domain(domain,options,nint(model_time/options%out_dt))
	endif
!
!-----------------------------------------
!-----------------------------------------
!  Time Loop
!
! 	note that a timestep here is a forcing input timestep O(1-3hr), not a physics timestep O(20-100s)
	do i=start_point,options%ntimesteps
		write(*,*) ""
		write(*,*) " ----------------------------------------------------------------------"
		write(*,*) "Timestep:", i, "  of ", options%ntimesteps
		write(*,*) "  Model time=",dnint(100*model_time/3600.0)/100.0,"hrs"
		
		if (.not.options%ideal) then
	! 		update boundary conditions (dXdt variables)
			call bc_update(domain,boundary,options)
		endif
! 		this is the meat of the model physics, run all the physics for the current time step looping over internal timesteps
		call step(domain,options,boundary,model_time,next_output)

	end do
!
!-----------------------------------------
	
end program real

