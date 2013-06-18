module init
	use data_structures
	use io_routines
	use geo
	
	implicit none
	private
	public::init_model
	
contains
	subroutine init_model(options_filename,options,domain,boundary)
		implicit none
		character(len=*), intent(in) :: options_filename
		type(options_type), intent(out) :: options
		type(domain_type), intent(out):: domain
		type(bc_type), intent(out):: boundary
		
! 		read in options file
		write(*,*) "Init Options"
		call init_options(options_filename,options)
! 		allocate and initialize the domain
		write(*,*) "Init Domain"
		call init_domain(options,domain)
! 		allocate and initialize the boundary conditions structure (includes 3D grids too...)
!		this might be more apropriately though of as a forcing data structure (for low res model)
		write(*,*) "Init Boundaries"
		call init_bc(options,domain,boundary)
		write(*,*) "Finished Initialization"
		
	end subroutine init_model
	
	subroutine init_options(options_filename,options)
! 		reads a series of options from a namelist file and stores them in the 
! 		options data structure
		implicit none
		character(len=*), intent(in) :: options_filename
		type(options_type), intent(out) :: options
		
		character(len=MAXFILELENGTH) :: init_conditions_file, output_file
		character(len=MAXFILELENGTH),allocatable:: boundary_files(:),ext_wind_files(:)
		character(len=MAXVARLENGTH) :: latvar,lonvar
		real :: dx,outputinterval,dz
		integer :: name_unit,ntimesteps,nfiles
		integer :: pbl,lsm,mp,rad,conv,adv,wind,nz,n_ext_winds
		logical :: readz,debug,external_winds
		
! 		set up namelist structures
		namelist /var_list/ latvar,lonvar
		namelist /parameters/ ntimesteps,outputinterval,dx,readz,nz,debug,dz,nfiles,external_winds,n_ext_winds
		namelist /files_list/ init_conditions_file,output_file,boundary_files,ext_wind_files
		namelist /physics/ pbl,lsm,mp,rad,conv,adv,wind
		
! 		read namelists
		open(io_newunit(name_unit), file=options_filename)
		read(name_unit,nml=var_list)
		read(name_unit,nml=parameters)
		read(name_unit,nml=physics)
		allocate(boundary_files(nfiles))
		if(external_winds)then
			allocate(ext_wind_files(n_ext_winds))
		endif
		read(name_unit,nml=files_list)
		close(name_unit)
		
! 		could probably simplify and read these all right from the namelist file, 
! 		but this way we can change the names in the file independant of the internal variable names
		options%init_conditions_file=init_conditions_file
		options%nfiles=nfiles
		allocate(options%boundary_files(nfiles))
		options%boundary_files=boundary_files
		options%output_file=output_file
		options%latvar=latvar
		options%lonvar=lonvar
		options%ntimesteps=ntimesteps
		options%io_dt=outputinterval
		options%dx=dx
		options%dz=dz
		options%readz=readz
		options%external_winds=external_winds
		options%n_ext_winds=n_ext_winds
		if(external_winds)then
			allocate(options%ext_wind_files(n_ext_winds))
		endif
		options%nz=nz
		options%debug=debug
		options%physics%boundarylayer=pbl
		options%physics%convection=conv
		options%physics%advection=adv
		options%physics%landsurface=lsm
		options%physics%microphysics=mp
		options%physics%radiation=rad
		options%physics%windtype=wind
		
	end subroutine init_options
	
	subroutine remove_edges(domain,edgesize)
		type(domain_type), intent(inout) :: domain
		integer, intent(in)::edgesize
		
		write(*,*) "reduce edges here..."
		
	end subroutine remove_edges
			
	subroutine init_domain(options, domain)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(out):: domain
		integer:: ny,nz,nx,i
		
! 		these are the only required variables on a high-res grid, lat, lon, and terrain elevation
		call io_read2d(options%init_conditions_file,"HGT",domain%terrain,1)
		call io_read2d(options%init_conditions_file,options%latvar,domain%lat,1)
		call io_read2d(options%init_conditions_file,options%lonvar,domain%lon,1)
		
		if(options%buffer>0) then
			call remove_edges(domain,options%buffer)
		endif
! 		use the lat variable to define the x and y dimensions for all other variables
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
! 		assumes nz is defined in the options
		nz=options%nz
		
! 		if a 3d grid was also specified, then read those data in
		if (options%readz) then
			if (options%debug) then
				write(*,*) "Reading 3D Z data"
			endif
			call io_read3d(options%init_conditions_file,"z", domain%z)
! 			dz also has to be calculated from the 3d z file
			allocate(domain%dz(nx,nz,ny))
			domain%dz(:,1:nz-1,:)=domain%z(:,2:nz,:)-domain%z(:,1:nz-1,:)
			domain%dz(:,nz,:)=domain%dz(:,nz-1,:)
		else
! 			otherwise, set up the z grid to be evenly spaced in z using the terrain +dz/2 for the base
! 			and z[1]+i*dz for the res
			allocate(domain%z(nx,nz,ny))
			do i=1,nz
				domain%z(:,i,:)=domain%terrain+(i*options%dz)-(options%dz/2)
			enddo
! 			here dz is just constant, but must be on a 3d grid for microphysics code
			allocate(domain%dz(nx,nz,ny))
			domain%dz=options%dz
		endif
		if (options%debug) then
			write(*,*) "allocating domain wide memory"
		endif
! 		all other variables should be allocated and initialized to 0
		allocate(domain%p(nx,nz,ny))
		domain%p=0
		allocate(domain%u(nx,nz,ny))
		domain%u=0
		allocate(domain%v(nx,nz,ny))
		domain%v=0
		allocate(domain%th(nx,nz,ny))
		domain%th=0
		allocate(domain%qv(nx,nz,ny))
		domain%qv=0
		allocate(domain%cloud(nx,nz,ny))
		domain%cloud=0
		allocate(domain%w(nx,nz,ny))
		domain%w=0
		allocate(domain%ice(nx,nz,ny))
		domain%ice=0
		allocate(domain%nice(nx,nz,ny))
		domain%nice=0
		allocate(domain%qrain(nx,nz,ny))
		domain%qrain=0
		allocate(domain%nrain(nx,nz,ny))
		domain%nrain=0
		allocate(domain%qsnow(nx,nz,ny))
		domain%qsnow=0
		allocate(domain%qgrau(nx,nz,ny))
		domain%qgrau=0
		allocate(domain%rain(nx,ny))
		domain%rain=0
		allocate(domain%snow(nx,ny))
		domain%snow=0
		allocate(domain%graupel(nx,ny))
		domain%graupel=0
		
! 		store dx in domain as well as options, read as an option, but it is more appropriate in domain
		domain%dx=options%dx
		
	end subroutine init_domain
	
	subroutine init_bc_data(options,boundary,domain)
		implicit none
		type(options_type), intent(in) :: options
		type(bc_type), intent(out):: boundary
		type(domain_type), intent(in):: domain
		integer::nx,ny,nz
		
! 		these variables are required for any boundary/forcing file type
		call io_read2d(options%boundary_files(1),options%latvar,boundary%lat)
		call io_read2d(options%boundary_files(1),options%lonvar,boundary%lon)
		call io_read2d(options%boundary_files(1),"HGT",boundary%terrain)
		
		nx=size(domain%lat,1)
		nz=options%nz
		ny=size(domain%lat,2)
! 		all other structures must be allocated and initialized, but will be set on a forcing timestep
! 		this also makes it easier to change how these variables are read from various forcing model file structures
		allocate(boundary%dudt(nx,nz,ny))
		boundary%dudt=0
		allocate(boundary%dvdt(nx,nz,ny))
		boundary%dvdt=0
		allocate(boundary%dwdt(nx,nz,ny))
		boundary%dwdt=0
		allocate(boundary%dpdt(nx,nz,ny))
		boundary%dpdt=0
		allocate(boundary%dthdt(nz,max(nx,ny),4))
		boundary%dthdt=0
		allocate(boundary%dqvdt(nz,max(nx,ny),4))
		boundary%dqvdt=0
		allocate(boundary%dqcdt(nz,max(nx,ny),4))
		boundary%dqcdt=0
		
	end subroutine init_bc_data
	
	subroutine init_ext_winds(options,bc)
		type(options_type), intent(in) :: options
		type(bc_type),intent(inout) :: bc
			
		real, allocatable, dimension(:,:,:) :: u,v
		real, allocatable, dimension(:,:) :: lat,lon
		
		call io_read2d(options%ext_wind_files(1),options%latvar,bc%ext_winds%lat)
		call io_read2d(options%ext_wind_files(1),options%lonvar,bc%ext_winds%lon)
		call geo_LUT(bc%next_domain, bc%ext_winds)
	end subroutine init_ext_winds
	
	
	subroutine init_bc(options,domain,boundary)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(in):: domain
		type(bc_type), intent(out):: boundary
			
! 		set up base data
		call init_bc_data(options,boundary,domain)
		call init_domain(options,boundary%next_domain) !set up a domain to hold the forcing for the next time step
! 		create the geographic look up table used to calculate boundary forcing data
		call geo_LUT(domain,boundary)
		
		if (options%external_winds) then
			call init_ext_winds(options,boundary)
		endif
		
		call geo_interp2d(boundary%next_domain%terrain,boundary%terrain,boundary%geolut)
		
	end subroutine init_bc
end module
