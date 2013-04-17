module init
	use data_structures
	use io_routines
	use geo
	implicit none
	contains
	
	subroutine init_model(options_filename,options,domain,boundary)
		implicit none
		character(len=*), intent(in) :: options_filename
		type(options_type), intent(out) :: options
		type(domain_type), intent(out):: domain
		type(bc_type), intent(out):: boundary
		
		call init_options(options_filename,options)
		call init_domain(options,domain)
		call init_bc(options,domain,boundary)
		
	end subroutine init_model
	
	subroutine init_options(options_filename,options)
		implicit none
		character(len=*), intent(in) :: options_filename
		type(options_type), intent(out) :: options
		
		character(len=100) :: init_conditions_file, output_file, boundary_file
		character(len=100) :: latvar,lonvar
		real :: dx
		integer :: name_unit,ntimesteps,outputinterval,timestep
		integer :: pbl,lsm,mp,rad,conv,adv,wind
		
		namelist /files_list/ init_conditions_file,output_file,boundary_file
		namelist /var_list/ latvar,lonvar
		namelist /parameters/ ntimesteps,outputinterval,timestep,dx
		namelist /physics/ pbl,lsm,mp,rad,conv,adv,wind
		
		
		open(io_newunit(name_unit), file=options_filename)
		read(name_unit,nml=files_list)
		read(name_unit,nml=var_list)
		read(name_unit,nml=parameters)
		read(name_unit,nml=physics)
		close(name_unit)
		
		options%init_conditions_file=init_conditions_file
		options%boundary_file=boundary_file
		options%output_file=output_file
		options%latvar=latvar
		options%lonvar=lonvar
		options%ntimesteps=ntimesteps
		options%outputinterval=outputinterval
		options%timestep=timestep
		options%dx=dx
		options%physics%boundarylayer=pbl
		options%physics%convection=conv
		options%physics%advection=adv
		options%physics%landsurface=lsm
		options%physics%microphysics=mp
		options%physics%radiation=rad
		options%physics%windtype=wind
		
	end subroutine init_options
			
	subroutine init_domain(options, domain)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(out):: domain
		integer:: ny,nz,nx
		
! 		allocate(domain%reader)
! 		call domain%reader%init(domain)
! 		eventually may want to make these read from a low res BC or init file? 

! 		these are the only fields we read initial conditions for
		call io_read3d(options%init_conditions_file,"p",domain%p)
		call io_read3d(options%init_conditions_file,"u",domain%u)
		call io_read3d(options%init_conditions_file,"v",domain%v)
! 		w is actually calculated form the divergence in u and v
! 		call io_read3d(options%init_conditions_file,"w",domain%w)
		call io_read3d(options%init_conditions_file,"th",domain%th)
		call io_read3d(options%init_conditions_file,"qv",domain%qv)
		call io_read3d(options%init_conditions_file,"qc",domain%cloud)
		
		call io_read3d(options%init_conditions_file,"z",domain%z)
		call io_read3d(options%init_conditions_file,"dz",domain%dz)
		call io_read2d(options%init_conditions_file,"hgt",domain%terrain)

		call io_read2d(options%init_conditions_file,options%latvar,domain%lat)
		call io_read2d(options%init_conditions_file,options%lonvar,domain%lon)
				
! 		any variables that aren't read from an initial conditions file should be allocated and initialized to 0
		ny=size(domain%p,1)
		nz=size(domain%p,2)
		nx=size(domain%p,3)
		allocate(domain%w(ny,nz,nx))
		domain%w=0
		allocate(domain%ice(ny,nz,nx))
		domain%ice=0
		allocate(domain%nice(ny,nz,nx))
		domain%nice=0
		allocate(domain%qrain(ny,nz,nx))
		domain%qrain=0
		allocate(domain%nrain(ny,nz,nx))
		domain%nrain=0
		allocate(domain%qsnow(ny,nz,nx))
		domain%qsnow=0
		allocate(domain%qgrau(ny,nz,nx))
		domain%qgrau=0
		allocate(domain%rain(ny,nx))
		domain%rain=0
		allocate(domain%snow(ny,nx))
		domain%snow=0
		allocate(domain%graupel(ny,nx))
		domain%graupel=0
		
		domain%dx=options%dx
		domain%dt=options%dt
		
	end subroutine init_domain
	
	subroutine init_bc_data(options,boundary)
		implicit none
		type(options_type), intent(in) :: options
		type(bc_type), intent(out):: boundary
		
		call io_read2d(options%boundary_file,options%latvar,boundary%lat)
		call io_read2d(options%boundary_file,options%lonvar,boundary%lon)
		
	end subroutine init_bc_data
	
	subroutine init_bc(options,domain,boundary)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(in):: domain
		type(bc_type), intent(out):: boundary
		
		call init_bc_data(options,boundary)
		call geo_LUT(domain,boundary)
		
		
	end subroutine init_bc
end module
