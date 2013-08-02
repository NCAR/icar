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
		type(options_type), intent(inout) :: options
		type(domain_type), intent(inout):: domain
		type(bc_type), intent(inout):: boundary
		
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
		type(options_type), intent(inout) :: options
		
		character(len=MAXFILELENGTH) :: init_conditions_file, output_file,restart_file
		character(len=MAXFILELENGTH),allocatable:: boundary_files(:),ext_wind_files(:)
		character(len=MAXVARLENGTH) :: latvar,lonvar
		real :: dx,outputinterval,dz
		integer :: name_unit,ntimesteps,nfiles
		integer :: pbl,lsm,mp,rad,conv,adv,wind,nz,n_ext_winds,buffer,restart_step
		logical :: readz,debug,external_winds,remove_lowres_linear,mean_winds,mean_fields,restart,add_low_topo
		n_ext_winds=200
		
! 		set up namelist structures
		namelist /var_list/ latvar,lonvar
		namelist /parameters/ ntimesteps,outputinterval,dx,readz,nz,debug,dz,nfiles, &
							  external_winds,buffer,n_ext_winds,add_low_topo,&
							  remove_lowres_linear,mean_winds,mean_fields,restart
		namelist /files_list/ init_conditions_file,output_file,boundary_files
		namelist /restart_info/ restart_step,restart_file
		namelist /ext_winds_info/ n_ext_winds,ext_wind_files
		namelist /physics/ pbl,lsm,mp,rad,conv,adv,wind
		
! 		read namelists
		open(io_newunit(name_unit), file=options_filename)
		read(name_unit,nml=var_list)
		read(name_unit,nml=parameters)
		
		if (restart) then
			read(name_unit,nml=restart_info)
			options%restart=restart
			options%restart_step=restart_step
			options%restart_file=restart_file
		endif

		read(name_unit,nml=physics)
		
		allocate(boundary_files(nfiles))
		read(name_unit,nml=files_list)
		
		if (external_winds) then
			write(*,*) n_ext_winds
			allocate(ext_wind_files(n_ext_winds))
			write(*,*) n_ext_winds
			read(name_unit,nml=ext_winds_info)
			write(*,*) n_ext_winds
			options%external_winds=external_winds
			options%ext_winds_nfiles=n_ext_winds
			allocate(options%ext_wind_files(n_ext_winds))
			options%ext_wind_files=ext_wind_files
			deallocate(ext_wind_files)
		endif
		
		close(name_unit)
		
! 		could probably simplify and read these all right from the namelist file, 
! 		but this way we can change the names in the file independant of the internal variable names
		options%init_conditions_file=init_conditions_file
		options%nfiles=nfiles
		allocate(options%boundary_files(nfiles))
		options%boundary_files=boundary_files
		deallocate(boundary_files)
		options%output_file=output_file
		options%latvar=latvar
		options%lonvar=lonvar
		options%ntimesteps=ntimesteps
		options%io_dt=outputinterval
		options%dx=dx
		options%dz=dz
		options%readz=readz
		options%buffer=buffer
		options%remove_lowres_linear=remove_lowres_linear
		options%add_low_topo=add_low_topo
		options%mean_winds=mean_winds
		options%mean_fields=mean_fields
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
		integer::nx1,ny1,nx2,ny2
		real,allocatable,dimension(:,:)::temp_data
		
		nx1=size(domain%lat,1)
		ny1=size(domain%lat,2)
		nx2=nx1-(edgesize*2)
		ny2=ny1-(edgesize*2)
		allocate(temp_data(nx1,ny1))
		
		temp_data=domain%lat
		deallocate(domain%lat)
		allocate(domain%lat(nx2,ny2))
		domain%lat=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		temp_data=domain%lon
		deallocate(domain%lon)
		allocate(domain%lon(nx2,ny2))
		domain%lon=temp_data(edgesize:nx1-edgesize,edgesize:ny1-edgesize)

		temp_data=domain%terrain
		deallocate(domain%terrain)
		allocate(domain%terrain(nx2,ny2))
		domain%terrain=temp_data(edgesize:nx1-edgesize,edgesize:ny1-edgesize)

		deallocate(temp_data)
	end subroutine remove_edges
	
! 	allocate all arrays in domain
	subroutine domain_allocation(domain,nx,nz,ny)
		type(domain_type), intent(inout) :: domain
		integer,intent(in)::nx,nz,ny
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
		
	end subroutine domain_allocation
	
! 	initialize the domain e.g. lat,lon,terrain, 3D z coordinate
	subroutine init_domain(options, domain)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(inout):: domain
		real,dimension(45)::fulldz
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
			allocate(domain%dz(nx,nz,ny))
! 			domain%dz=options%dz
! 			mean layer thicknesses from a 36km WRF run over the "CO-headwaters" domain
			fulldz=[36.,   51.,   58.,   73.,   74.,  111.,  113.,  152.,  155.,  157.,  160.,  245., &
				   251.,  258.,  265.,  365.,  379.,  395.,  413.,  432.,  453.,  476.,  503.,  533., &
				   422.,  443.,  467.,  326.,  339.,  353.,  369.,  386.,  405.,  426.,  450.,  477., &
				   455.,  429.,  396.,  357.,  311.,  325.,  340.,  356.,  356.]
			domain%dz(:,1,:)=fulldz(1)
			domain%z(:,1,:)=domain%terrain+fulldz(1)/2
			do i=2,nz
				domain%z(:,i,:)=domain%z(:,i-1,:)+(fulldz(i)+fulldz(i-1))/2
				domain%dz(:,i,:)=fulldz(i)
			enddo
! 			here dz is just constant, but must be on a 3d grid for microphysics code
		endif
		if (options%debug) then
			write(*,*) "allocating domain wide memory"
		endif
! 		all other variables should be allocated and initialized to 0
		call domain_allocation(domain,nx,nz,ny)
		
! 		store dx in domain as well as options, read as an option, but it is more appropriate in domain
		domain%dx=options%dx
		
	end subroutine init_domain
	
! 	allocate arrays in boundary condition data structure
	subroutine boundary_allocate(boundary,nx,nz,ny)
		type(bc_type), intent(inout) :: boundary
		integer,intent(in)::nx,nz,ny
		
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
	end subroutine boundary_allocate
	
! 	initialize the boundary condition data structure e.g. lat,lon,terrain,3D Z coord
	subroutine init_bc_data(options,boundary,domain)
		implicit none
		type(options_type), intent(in) :: options
		type(bc_type), intent(inout):: boundary
		type(domain_type), intent(in):: domain
		real,dimension(45)::fulldz
		integer::nx,ny,nz,i
		
! 		these variables are required for any boundary/forcing file type
		call io_read2d(options%boundary_files(1),options%latvar,boundary%lat)
		call io_read2d(options%boundary_files(1),options%lonvar,boundary%lon)
		call io_read2d(options%boundary_files(1),"HGT",boundary%terrain)
		
		nx=size(boundary%lat,1)
		nz=options%nz
		ny=size(boundary%lat,2)
		
		allocate(boundary%z(nx,nz,ny))
		allocate(boundary%dz(nx,nz,ny))
! 			domain%dz=options%dz
! 			mean layer thicknesses from a 36km WRF run over the "CO-headwaters" domain
		fulldz=[36.,   51.,   58.,   73.,   74.,  111.,  113.,  152.,  155.,  157.,  160.,  245., &
			   251.,  258.,  265.,  365.,  379.,  395.,  413.,  432.,  453.,  476.,  503.,  533., &
			   422.,  443.,  467.,  326.,  339.,  353.,  369.,  386.,  405.,  426.,  450.,  477., &
			   455.,  429.,  396.,  357.,  311.,  325.,  340.,  356.,  356.]
		boundary%dz(:,1,:)=fulldz(1)
		boundary%z(:,1,:)=boundary%terrain+fulldz(1)/2
		do i=2,nz
			boundary%z(:,i,:)=boundary%z(:,i-1,:)+(fulldz(i)+fulldz(i-1))/2
			boundary%dz(:,i,:)=fulldz(i)
		enddo
	
		
! 		all other structures must be allocated and initialized, but will be set on a high-res grid
! 		u/v are seperate so we can read them on the low res grid and adjust/rm-linearwinds before interpolating
! 		this also makes it easier to change how these variables are read from various forcing model file structures
		allocate(boundary%u(nx,nz,ny))
		boundary%u=0
		allocate(boundary%v(nx,nz,ny))
		boundary%v=0
		
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
		
		call boundary_allocate(boundary,nx,nz,ny)
	end subroutine init_bc_data
	
! 	initialize the external wind system (GEOLUT)
	subroutine init_ext_winds(options,bc)
		type(options_type), intent(in) :: options
		type(bc_type),intent(inout) :: bc
			
		real, allocatable, dimension(:,:,:) :: u,v
		real, allocatable, dimension(:,:) :: lat,lon
		
		call io_read2d(options%ext_wind_files(1),options%latvar,bc%ext_winds%lat)
		call io_read2d(options%ext_wind_files(1),options%lonvar,bc%ext_winds%lon)
		write(*,*) maxval(bc%ext_winds%lat),minval(bc%ext_winds%lat)
		write(*,*) maxval(bc%next_domain%lat),minval(bc%next_domain%lat)
		write(*,*) "Setting up ext wind geoLUT"
		call geo_LUT(bc%next_domain, bc%ext_winds)
! 		if (options%debug) then
! 			call io_write3di("geolut_x.nc","data",bc%ext_winds%geolut%x)
! 			call io_write3di("geolut_y.nc","data",bc%ext_winds%geolut%y)
! 			call io_write3d("geolut_w.nc","data",bc%ext_winds%geolut%w)
! 		endif
	end subroutine init_ext_winds
	
! 	initialize the boundary condiditions (init data structures and GEOLUT)
	subroutine init_bc(options,domain,boundary)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(inout):: domain
		type(bc_type), intent(inout):: boundary
			
		write(*,*) "WARNING hardcoded low-res dx=36km"
		boundary%dx=36000.0
! 		set up base data
		call init_bc_data(options,boundary,domain)
		call init_domain(options,boundary%next_domain) !set up a domain to hold the forcing for the next time step
! 		create the geographic look up table used to calculate boundary forcing data
		call geo_LUT(domain,boundary)
! 		if (options%debug) then
! 			call io_write3di("bcgeolut_x.nc","data",boundary%geolut%x)
! 			call io_write3di("bcgeolut_y.nc","data",boundary%geolut%y)
! 			call io_write3d("bcgeolut_w.nc","data",boundary%geolut%w)
! 		endif
		
		if (options%external_winds) then
			call init_ext_winds(options,boundary)
		endif
		
! 		interpolate the low-res terrain to the high-res grid for pressure adjustments. 
! 		the correct way would probably be to adjust all low-res pressures to Sea level before interpolating
! 		then pressure adjustments all occur from SLP. 
		call geo_interp2d(boundary%next_domain%terrain,boundary%terrain,boundary%geolut)
		if (options%add_low_topo) then
			domain%terrain=domain%terrain+(boundary%next_domain%terrain-sum(boundary%next_domain%terrain)/size(boundary%next_domain%terrain))
		endif
		
	end subroutine init_bc
end module
