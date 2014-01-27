module init
! ----------------------------------------------------------------------------
! 	NOTE: This module was initially written to read WRF output files as input.
! 		This should serve as a basis for any additional file types, and can be 
! 		readily modified.  At some point this could be modified to check the 
! 		type if input file being specified, and call an appropriate routine. 
! 		e.g. if options%inputtype=="WRF" then call wrf_init/update()
! ----------------------------------------------------------------------------
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
		character(len=MAXVARLENGTH) :: landvar="",latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar,&
										hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,     &
										pvar,tvar,qvvar,qcvar,qivar,hgtvar,shvar,lhvar,pblhvar,version
		real :: dx,dxlow,outputinterval,inputinterval
		real, allocatable, dimension(:) :: dz_levels
		integer :: name_unit,ntimesteps,nfiles,xmin,xmax,ymin,ymax
		integer :: pbl,lsm,mp,rad,conv,adv,wind,nz,n_ext_winds,buffer,restart_step
		logical :: ideal, readz,readdz,debug,external_winds,remove_lowres_linear,&
		           mean_winds,mean_fields,restart,add_low_topo
   		real,dimension(45)::fulldz
		
! 		set up namelist structures
		namelist /model_version/ version
		namelist /var_list/ pvar,tvar,qvvar,qcvar,qivar,hgtvar,shvar,lhvar,pblhvar,&
							landvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,zvar, &
							hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi 
		namelist /parameters/ ntimesteps,outputinterval,inputinterval,dx,dxlow,ideal,readz,readdz,nz,debug,nfiles, &
							  external_winds,buffer,n_ext_winds,add_low_topo,&
							  remove_lowres_linear,mean_winds,mean_fields,restart,xmin,xmax,ymin,ymax
		namelist /z_info/ dz_levels
		namelist /files_list/ init_conditions_file,output_file,boundary_files
		namelist /restart_info/ restart_step,restart_file
		namelist /ext_winds_info/ n_ext_winds,ext_wind_files
		namelist /physics/ pbl,lsm,mp,rad,conv,adv,wind
		n_ext_winds=200
		
! 		read namelists
		open(io_newunit(name_unit), file=options_filename)
		read(name_unit,nml=model_version)
		if (version.ne."0.7.2") then
			write(*,*) "Model version does not match namelist version"
			write(*,*) "  Model version: 0.7.2"
			write(*,*) "  Namelist version:",version
			stop
		endif
		write(*,*) "Model version: ",version, "uses landmask, working on better vertical interpolation"
		read(name_unit,nml=var_list)
		read(name_unit,nml=parameters)
		
		if (readdz) then
			allocate(dz_levels(nz),options%dz_levels(nz))
			read(name_unit,nml=z_info)
			options%dz_levels=dz_levels
			deallocate(dz_levels)
		else
! 			mean layer thicknesses from a 36km WRF run over the "CO-headwaters" domain
			fulldz=[36.,   51.,   58.,   73.,   74.,  111.,  113.,  152.,  155.,  157.,  160.,  245., &
				   251.,  258.,  265.,  365.,  379.,  395.,  413.,  432.,  453.,  476.,  503.,  533., &
				   422.,  443.,  467.,  326.,  339.,  353.,  369.,  386.,  405.,  426.,  450.,  477., &
				   455.,  429.,  396.,  357.,  311.,  325.,  340.,  356.,  356.]
! 			mean layer thicknesses from ERAi domain
! 		   fulldz=[   24.8,  36.5,  51.8,  70.1,  90.8, 113.5, 137.9, 163.7, 190.5, 218.1, 246.4, &
! 		   			 275.1, 304.3, 333.6, 363.0, 392.4, 421.7, 450.8, 479.6, 508.0, 535.9, 563.2, &
! 					 589.8, 615.7, 640.9, 665.5, 689.8, 714.1, 739.4, 767.2, 796.8, 826.6, 856.2, &
! 					 885.1, 912.5, 937.9, 961.4, 979.4, 990.1, 976.6, 937.6, 900.1, 864.2, 829.6, 796.5]
 			allocate(options%dz_levels(nz))
			options%dz_levels=fulldz(1:nz)
		endif
		
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
			allocate(ext_wind_files(n_ext_winds))
			read(name_unit,nml=ext_winds_info)
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
		options%latvar=latvar
		options%lonvar=lonvar
		options%hgtvar=hgtvar
		options%uvar=uvar
		options%shvar=shvar
		options%lhvar=lhvar
		options%pblhvar=pblhvar
		options%ulat=ulat
		options%ulon=ulon
		options%vvar=vvar
		options%vlat=vlat
		options%vlon=vlon
		options%pvar=pvar
		options%tvar=tvar
		options%qvvar=qvvar
		options%qcvar=qcvar
		options%qivar=qivar
		
		options%landvar=landvar
		options%zvar=zvar
		options%hgt_hi=hgt_hi
		options%lat_hi=lat_hi
		options%lon_hi=lon_hi
		options%ulat_hi=ulat_hi
		options%ulon_hi=ulon_hi
		options%vlat_hi=vlat_hi
		options%vlon_hi=vlon_hi
		
		options%ntimesteps=ntimesteps
		options%in_dt=inputinterval
		options%out_dt=outputinterval
		options%dx=dx
		options%dxlow=dxlow
		options%ideal=ideal
		if (ideal) then
			write(*,*) "Running Idealized simulation (time step does not advance)"
		endif
		options%readz=readz
		options%buffer=buffer
		options%remove_lowres_linear=remove_lowres_linear
		options%add_low_topo=add_low_topo
		options%mean_winds=mean_winds
		options%mean_fields=mean_fields
		options%nz=nz
		options%xmin=xmin
		options%xmax=xmax
		options%ymin=ymin
		options%ymax=ymax
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
		integer::nx1,ny1,nx2,ny2,nz
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

		domain%lon=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		temp_data=domain%terrain
		deallocate(domain%terrain)
		allocate(domain%terrain(nx2,ny2))
		domain%terrain=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		temp_data=domain%landmask
		deallocate(domain%landmask)
		allocate(domain%landmask(nx2,ny2))
		domain%landmask=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		
		deallocate(temp_data)
		
		nx1=size(domain%u_geo%lat,1)
		ny1=size(domain%u_geo%lat,2)
		nx2=nx1-(edgesize*2)
		ny2=ny1-(edgesize*2)
		
		allocate(temp_data(nx1,ny1))
		temp_data=domain%u_geo%lat
		deallocate(domain%u_geo%lat)
		allocate(domain%u_geo%lat(nx2,ny2))
		domain%u_geo%lat=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
		
		temp_data=domain%u_geo%lon
		deallocate(domain%u_geo%lon)
		allocate(domain%u_geo%lon(nx2,ny2))
		domain%u_geo%lon=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

		deallocate(temp_data)
		
		nx1=size(domain%v_geo%lat,1)
		ny1=size(domain%v_geo%lat,2)
		nx2=nx1-(edgesize*2)
		ny2=ny1-(edgesize*2)
		
		allocate(temp_data(nx1,ny1))
		temp_data=domain%v_geo%lat
		deallocate(domain%v_geo%lat)
		allocate(domain%v_geo%lat(nx2,ny2))
		domain%v_geo%lat=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
		
		temp_data=domain%v_geo%lon
		deallocate(domain%v_geo%lon)
		allocate(domain%v_geo%lon(nx2,ny2))
		domain%v_geo%lon=temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
		
		deallocate(temp_data)
	end subroutine remove_edges
	
! 	allocate all arrays in domain
	subroutine domain_allocation(domain,nx,nz,ny)
		type(domain_type), intent(inout) :: domain
		integer,intent(in)::nx,nz,ny
		allocate(domain%p(nx,nz,ny))
		domain%p=0
		allocate(domain%u(nx+1,nz,ny))
		domain%u=0
		allocate(domain%v(nx,nz,ny+1))
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
		allocate(domain%crain(nx,ny))
		domain%crain=0
		allocate(domain%snow(nx,ny))
		domain%snow=0
		allocate(domain%graupel(nx,ny))
		domain%graupel=0
		allocate(domain%sensible_heat(nx,ny))
		domain%sensible_heat=0
		allocate(domain%latent_heat(nx,ny))
		domain%latent_heat=0
		allocate(domain%pbl_height(nx,ny))
		domain%pbl_height=0
		
	end subroutine domain_allocation
	
! 	initialize the domain e.g. lat,lon,terrain, 3D z coordinate
	subroutine init_domain(options, domain)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(inout):: domain
		real,dimension(:,:,:),allocatable::temporary_z
		integer:: ny,nz,nx,i,buf
		
! 		these are the only required variables on a high-res grid, lat, lon, and terrain elevation
		call io_read2d(options%init_conditions_file,options%hgt_hi,domain%terrain,1)
		call io_read2d(options%init_conditions_file,options%lat_hi,domain%lat,1)
		call io_read2d(options%init_conditions_file,options%lon_hi,domain%lon,1)
		call io_read2d(options%init_conditions_file,options%ulat_hi,domain%u_geo%lat,1)
		call io_read2d(options%init_conditions_file,options%ulon_hi,domain%u_geo%lon,1)
		call io_read2d(options%init_conditions_file,options%vlat_hi,domain%v_geo%lat,1)
		call io_read2d(options%init_conditions_file,options%vlon_hi,domain%v_geo%lon,1)
		if (options%landvar/="") then
			call io_read2d(options%init_conditions_file,options%landvar,domain%landmask,1)
		else
			nx=size(domain%lat,1)
			ny=size(domain%lat,2)
			allocate(domain%landmask(nx,ny))
			domain%landmask=1 !if we weren't supplied a landmask field, assume all is land (what we care about anyway)
		endif
		
		if(options%buffer>0) then
			call remove_edges(domain,options%buffer)
		endif
! 		use the lat variable to define the x and y dimensions for all other variables
		nx=size(domain%lat,1)
		ny=size(domain%lat,2)
! 		assumes nz is defined in the options
		nz=options%nz
		
! 		if a 3d grid was also specified, then read those data in
		if ((options%readz).and.(options%ideal)) then
			if (.not.options%ideal) then
				write(*,*) "Reading Z only recommended for ideal runs at the moment"
			endif
			
			call io_read3d(options%init_conditions_file,options%zvar, domain%z)
! 			dz also has to be calculated from the 3d z file
			buf=options%buffer
			allocate(domain%dz(nx,nz,ny))
			allocate(temporary_z(nx,nz,ny))
			do i=1,nz-1
				domain%dz(:,i,:)=domain%z(buf+1:nx+buf,buf+1:ny+buf,i+1)-domain%z(buf+1:nx+buf,buf+1:ny+buf,i)
				temporary_z(:,i,:)=domain%z(buf+1:nx+buf,buf+1:ny+buf,i)
			enddo
			temporary_z(:,nz,:)=domain%z(buf+1:nx+buf,buf+1:ny+buf,nz)
			domain%dz(:,nz,:)=domain%dz(:,nz-1,:)
			deallocate(domain%z)
			allocate(domain%z(nx,nz,ny))
			domain%z=temporary_z
			deallocate(temporary_z)
		else
! 			otherwise, set up the z grid to be evenly spaced in z using the terrain +dz/2 for the base
! 			and z[1]+i*dz for the res
			allocate(domain%z(nx,nz,ny))
			allocate(domain%dz(nx,nz,ny))
			domain%dz(:,1,:)=options%dz_levels(1)
			domain%z(:,1,:)=domain%terrain+options%dz_levels(1)/2
			do i=2,nz
				domain%z(:,i,:)=domain%z(:,i-1,:)+(options%dz_levels(i)+options%dz_levels(i-1))/2
				domain%dz(:,i,:)=options%dz_levels(i)
			enddo

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
		
		allocate(boundary%dudt(nx+1,nz,ny))
		boundary%dudt=0
		allocate(boundary%dvdt(nx,nz,ny+1))
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
		allocate(boundary%dlhdt(nx,ny))
		boundary%dlhdt=0
		allocate(boundary%dshdt(nx,ny))
		boundary%dshdt=0
		allocate(boundary%dpblhdt(nx,ny))
		boundary%dpblhdt=0
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
		call io_read2d(options%boundary_files(1),options%ulat,boundary%u_geo%lat)
		call io_read2d(options%boundary_files(1),options%ulon,boundary%u_geo%lon)
		call io_read2d(options%boundary_files(1),options%vlat,boundary%v_geo%lat)
		call io_read2d(options%boundary_files(1),options%vlon,boundary%v_geo%lon)
		call io_read2d(options%boundary_files(1),options%hgtvar,boundary%terrain)
		
		nx=size(boundary%lat,1)
		nz=options%nz
		ny=size(boundary%lat,2)
		
		allocate(boundary%z(nx,nz,ny))
		allocate(boundary%dz(nx,nz,ny))
!		mean layer thicknesses from a 36km WRF run over the "CO-headwaters" domain
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
	
    subroutine setup_extwinds(domain)
        implicit none
        type(wind_type),intent(inout)::domain
        integer::nx,ny
        
        nx=size(domain%terrain,1)
        ny=size(domain%terrain,2)
		
		! dzdx/y used in rotating windfield back to terrain following grid in a simple fashion
		allocate(domain%dzdx(nx-1,ny))
		allocate(domain%dzdy(nx,ny-1))
		domain%dzdx=sqrt((domain%terrain(2:nx,:)-domain%terrain(1:nx-1,:))**2+domain%dx**2)/domain%dx
		domain%dzdy=sqrt((domain%terrain(:,2:ny)-domain%terrain(:,1:ny-1))**2+domain%dx**2)/domain%dx
    end subroutine
	
	
! 	initialize the external wind system (GEOLUT)
	subroutine init_ext_winds(options,bc)
		type(options_type), intent(in) :: options
		type(bc_type),intent(inout) :: bc
			
		real, allocatable, dimension(:,:,:) :: u,v
		real, allocatable, dimension(:,:) :: lat,lon
		
		call io_read2d(options%ext_wind_files(1),options%hgt_hi,bc%ext_winds%terrain,1)
		call io_read2d(options%ext_wind_files(1),options%latvar,bc%ext_winds%lat)
		call io_read2d(options%ext_wind_files(1),options%lonvar,bc%ext_winds%lon)
		call io_read2d(options%ext_wind_files(1),options%ulat_hi,bc%ext_winds%u_geo%lat)
		call io_read2d(options%ext_wind_files(1),options%ulon_hi,bc%ext_winds%u_geo%lon)
		call io_read2d(options%ext_wind_files(1),options%vlat_hi,bc%ext_winds%v_geo%lat)
		call io_read2d(options%ext_wind_files(1),options%vlon_hi,bc%ext_winds%v_geo%lon)
		write(*,*) "Setting up ext wind geoLUTs"
		call geo_LUT(bc%next_domain, bc%ext_winds%u_geo)
		call geo_LUT(bc%next_domain, bc%ext_winds%v_geo)
		
! 		force all weight to be on the first x,y pair...
! 		this assumes the "external winds" file is on the exact same grid as the high res model grid
		bc%ext_winds%u_geo%geolut%w(2:,:,:)=0
		bc%ext_winds%u_geo%geolut%w(1,:,:)=1
		bc%ext_winds%v_geo%geolut%w(2:,:,:)=0
		bc%ext_winds%v_geo%geolut%w(1,:,:)=1
		bc%ext_winds%dx=bc%next_domain%dx
		call setup_extwinds(bc%ext_winds)
	end subroutine init_ext_winds
	
! 	initialize the boundary condiditions (init data structures and GEOLUT)
	subroutine init_bc(options,domain,boundary)
		implicit none
		type(options_type), intent(in) :: options
		type(domain_type), intent(inout):: domain
		type(bc_type), intent(inout):: boundary
		integer::i
			
		boundary%dx=options%dxlow
! 		set up base data
		call init_bc_data(options,boundary,domain)
		call init_domain(options,boundary%next_domain) !set up a domain to hold the forcing for the next time step
! 		create the geographic look up table used to calculate boundary forcing data
		write(*,*) "Setting up domain geographic Look Up Tables"
		call geo_LUT(domain,boundary)
		call geo_LUT(domain%u_geo,boundary%u_geo)
		call geo_LUT(domain%v_geo,boundary%v_geo)
		
		if (options%external_winds) then
			call init_ext_winds(options,boundary)
		endif
		
! 		interpolate the low-res terrain to the high-res grid for pressure adjustments. 
! 		the correct way would probably be to adjust all low-res pressures to Sea level before interpolating
! 		then pressure adjustments all occur from SLP. 
! 		This should be done on a separate lowres terrain grid so the embedded high res terrain grid can also be used in pressure adjustments on each time step...
! 		allocate(boundary%lowres_terrain(nx,ny))
		call geo_interp2d(boundary%next_domain%terrain,boundary%terrain,boundary%geolut)
		call geo_interp(boundary%next_domain%z,boundary%z,boundary%geolut,.false.)
		if (options%add_low_topo) then
			domain%terrain=domain%terrain+(boundary%next_domain%terrain-sum(boundary%next_domain%terrain) &
											 /size(boundary%next_domain%terrain))/2.0
			do i=1,size(domain%z,2)
				domain%z(:,i,:)=domain%z(:,i,:)+(boundary%next_domain%terrain-sum(boundary%next_domain%terrain) &
												 /size(boundary%next_domain%terrain))/2.0
			enddo
		endif
		
	end subroutine init_bc
end module
