module boundary_conditions
! ----------------------------------------------------------------------------
! 	NOTE: This module was initially written to read WRF output files as input.
! 		This should serve as a basis for any additional file types, and can be 
! 		readily modified.  At some point this could be modified to check the 
! 		type if input file being specified, and call an appropriate routine. 
! 		e.g. if options%inputtype=="WRF" then call wrf_init/update()
! ----------------------------------------------------------------------------
	use io_routines
	use data_structures
	use wind
	use linear_theory_winds
	use geo
	use output
	
	implicit none
	private
! 	these could be better stored in bc_type and initialized in init_bc?
	character (len=255),dimension(:),allocatable :: file_list
! 	manage file pointer and position in file for boundary conditions
	integer::curfile,curstep
	integer::steps_in_file,nfiles
! 	manage file pointer and position in file for external winds
	character (len=255),dimension(:),allocatable :: ext_winds_file_list
	integer::ext_winds_curfile,ext_winds_curstep
	integer::ext_winds_steps_in_file,ext_winds_nfiles

	integer,parameter::smoothing_window=14
!   these are now specified in data_structures.f90
!     real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
!     real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
	
	public::bc_init
	public::bc_update
contains
	
! 	smooth an array (written for wind but will work for anything)
	subroutine smooth_wind(wind,windowsize,ydim)
		real, intent(inout), dimension(:,:,:):: wind
		integer,intent(in)::windowsize
		integer,intent(in)::ydim
		real,allocatable,dimension(:,:,:)::inputwind
		integer::i,j,k,nx,ny,nz
		nx=size(wind,1)
		ny=size(wind,2)
		nz=size(wind,3)
		allocate(inputwind(nx,ny,nz))
		inputwind=wind
		
		!$omp parallel firstprivate(windowsize,nx,ny,nz),private(i,j,k),shared(wind,inputwind)
		!$omp do
		do k=1,nz
			do j=1,ny
				do i=1,nx
					if (ydim==2) then
						wind(i,j,k)=sum(inputwind(max(i-windowsize,1):min(nx,i+windowsize),max(1,j-windowsize):min(ny,j+windowsize),k)) &
							/ ((min(nx,i+windowsize)-max(i-windowsize,1)+1) * (min(ny,j+windowsize)-max(1,j-windowsize)+1))
					else
						wind(i,j,k)=sum(inputwind(max(i-windowsize,1):min(nx,i+windowsize),j,max(1,k-windowsize):min(nz,k+windowsize))) &
							/ ((min(nx,i+windowsize)-max(i-windowsize,1)+1) * (min(nz,k+windowsize)-max(1,k-windowsize)+1))
					endif
				enddo
			enddo
		enddo
		!$omp end do
		!$omp end parallel
		
		deallocate(inputwind)
	end subroutine smooth_wind
	
! 	generic routine to read a low res variable (varname) from a netcdf file (filename) at the current time step (curstep)
!   then interpolate it to the high res grid either in 3D or at the boundaries only (boundary_only)
! 	applies modifications specificaly for U,V,T,P and PH variables
	subroutine read_var(highres,filename,varname,geolut,curstep,boundary_only,low_res_terrain)
		implicit none
		real,dimension(:,:,:),intent(inout)::highres
		character(len=*),intent(in) :: filename,varname
		type(geo_look_up_table),intent(in) :: geolut
		integer,intent(in)::curstep
		logical, intent(in) :: boundary_only
		real,optional,dimension(:,:,:) :: low_res_terrain
		
		real,dimension(:,:,:),allocatable :: inputdata,extra_data
		integer,dimension(io_maxDims)::dims
		integer::nx,ny,nz,i
		
! 		Read the data in
		call io_getdims(filename,varname, dims)
		call io_read3d(filename,varname,inputdata,curstep)
		
! 		note dims(1)=ndims
		nx=dims(2)
		ny=dims(3)
		nz=dims(4)
		
! 		perform unit conversion on specific variables
		if (varname=="U") then
			call smooth_wind(inputdata,1,2)
		else if (varname=="V") then
			call smooth_wind(inputdata,1,2)
		else if (varname=="T") then
			inputdata=inputdata+300
		else if (varname=="P") then
			call io_read3d(filename,"PB",extra_data,curstep)
			inputdata=inputdata+extra_data
			deallocate(extra_data)
		else if (varname=="PH") then
			call io_read3d(filename,"PHB",extra_data,curstep)
			inputdata=(inputdata+extra_data)/9.8
			nz=nz-1
			inputdata(:,:,1:nz)=(inputdata(:,:,1:nz)+inputdata(:,:,2:nz+1))/2
			deallocate(extra_data)
		endif
		
! 		interpolate data onto the high resolution grid after re-arranging the dimensions. 
		nz=min(nz,size(highres,2))
		allocate(extra_data(nx,nz,ny))
		do i=1,nz
			extra_data(:,i,:)=inputdata(:,:,i)
		enddo
! 		extra_data=reshape(inputdata,[nx,nz,ny],order=[1,3,2])
		call geo_interp(highres, &
						extra_data, &
						geolut,boundary_only)
		deallocate(extra_data)
		deallocate(inputdata)
						
	end subroutine read_var
	
! 	rotate winds from real space back to terrain following grid (approximately)
!   assumes a simple slope transform in u and v independantly
	subroutine rotate_ext_wind_field(domain,ext_winds)
        implicit none
        type(domain_type),intent(inout)::domain
        type(wind_type),intent(inout)::ext_winds
		integer :: nx,ny,nz,i
	
		ny=size(domain%z,1)
		nz=size(domain%z,2)
		nx=size(domain%z,3)
		do i=1,nz
			domain%u(:,i,1:nx-1)=domain%u(:,i,1:nx-1)*ext_winds%dzdx
			domain%v(1:ny-1,i,:)=domain%v(1:ny-1,i,:)*ext_winds%dzdy
		end do
	
	end subroutine rotate_ext_wind_field

	
! 	initialize the eternal winds information (filenames, nfiles, etc) and read the initial conditions
	subroutine ext_winds_init(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		! MODULE variables : ext_winds_ curstep, curfile, nfiles, steps_in_file, file_list
		ext_winds_curfile=1
		if (options%restart) then
			ext_winds_curstep=options%restart_step
		else
			ext_winds_curstep=1
		endif
		ext_winds_nfiles=options%ext_winds_nfiles
		allocate(ext_winds_file_list(ext_winds_nfiles))
		ext_winds_file_list=options%ext_wind_files
		call io_getdims(ext_winds_file_list(ext_winds_curfile),"U", dims)
		if (dims(1)==3) then
			ext_winds_steps_in_file=1
		else
			ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims
		endif
		
! 		ext_winds_curstep=ext_winds_curstep+1
		if (.not.options%ideal) then
			do while (ext_winds_curstep>ext_winds_steps_in_file)
				ext_winds_curfile=ext_winds_curfile+1
				if (ext_winds_curfile>ext_winds_nfiles) then
					stop "Ran out of files to process!"
				endif
				ext_winds_curstep=ext_winds_curstep-ext_winds_steps_in_file 
				!instead of setting=1, this way we can set an arbitrary starting point multiple files in
				call io_getdims(ext_winds_file_list(ext_winds_curfile),"U", dims)
				if (dims(1)==3) then
					ext_winds_steps_in_file=1
				else
					ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
				endif
			enddo
		endif
		write(*,*) "Initial ext wind file:step=",ext_winds_curfile," : ",ext_winds_curstep
		call read_var(domain%u,    ext_winds_file_list(ext_winds_curfile),"U",      bc%ext_winds%u_geo%geolut,ext_winds_curstep,.FALSE.)
		call read_var(domain%v,    ext_winds_file_list(ext_winds_curfile),"V",      bc%ext_winds%v_geo%geolut,ext_winds_curstep,.FALSE.)
		call rotate_ext_wind_field(domain,bc%ext_winds)
	end subroutine ext_winds_init
	
! 	remove linear theory topographic winds perturbations from the low resolution wind field. 
	subroutine remove_linear_winds(domain,bc,options,filename,curstep)
        implicit none
		type(domain_type), intent(inout) :: domain
		type(bc_type), intent(inout) :: bc
		type(options_type),intent(in) :: options
		character(len=*),intent(in) :: filename
		integer,intent(in)::curstep
		integer,dimension(io_maxDims)::dims !note, io_maxDims is included from io_routines.
		real,allocatable,dimension(:,:,:)::inputdata,extra_data
		logical :: reverse_winds=.TRUE.
		integer :: nx,ny,nz,nz_output
		character(len=255) :: outputfilename
		
		call io_getdims(filename,"U", dims)
		nx=dims(2)
		ny=dims(3)
		nz=dims(4)
		nz_output=size(bc%u,2)
		
		allocate(inputdata(nx,ny,nz_output))
		
! 		first read in the low-res U and V data directly
! 		load low-res U data
		call io_read3d(filename,"U",extra_data,curstep)
		inputdata=extra_data(1:nx,1:ny,:nz_output)
		call smooth_wind(inputdata,2,2)
		bc%u=reshape(inputdata,[nx,nz_output,ny],order=[1,3,2])
		deallocate(extra_data,inputdata)
		
		allocate(inputdata(nx-1,ny+1,nz_output))
! 		load low-res V data
		call io_read3d(filename,"V",extra_data,curstep)
		inputdata=extra_data(1:nx-1,1:ny+1,:nz_output)
		call smooth_wind(inputdata,2,2)
		bc%v=reshape(inputdata,[nx,nz_output,ny],order=[1,3,2])
		deallocate(extra_data,inputdata)
! 		if (options%debug) then
! 			write(outputfilename,"(A,I5.5)") "U_pre",curstep
! 			call io_write3d(outputfilename,"data",bc%u)
! 			write(outputfilename,"(A,I5.5)") "V_pre",curstep
! 			call io_write3d(outputfilename,"data",bc%v)
! 		endif
! 		remove the low-res linear wind contribution effect
		call linear_perturb(bc,reverse_winds)
! 		if (options%debug) then
! 			write(outputfilename,"(A,I5.5)") "U_post",curstep
! 			call io_write3d(outputfilename,"data",bc%u)
! 			write(outputfilename,"(A,I5.5)") "V_post",curstep
! 			call io_write3d(outputfilename,"data",bc%v)
! 		endif
		
! 		finally interpolate low res winds to the high resolutions grid
		call geo_interp(domain%u, bc%u,bc%u_geo%geolut,.FALSE.)
		call geo_interp(domain%v, bc%v,bc%v_geo%geolut,.FALSE.)
	end subroutine remove_linear_winds
	
! for test cases compute the mean winds and make them constant everywhere...
	subroutine mean_winds(domain,filename,curstep)
		type(domain_type), intent(inout) :: domain
		character(len=*),intent(in)::filename
		integer,intent(in)::curstep
		
		real,allocatable,dimension(:,:,:)::extra_data
		integer,dimension(io_maxDims)::dims
		integer::nx,ny,nz
		
		nz=size(domain%u,2)
		
! 		load low-res U data
		call io_read3d(filename,"U",extra_data,curstep)
		domain%u=sum(extra_data(:,:,:nz))/size(extra_data(:,:,:nz))
		deallocate(extra_data)

! 		load low-res V data
		call io_read3d(filename,"V",extra_data,curstep)
		domain%v=sum(extra_data(:,:,:nz))/size(extra_data(:,:,:nz))
		deallocate(extra_data)
				
	end subroutine mean_winds
	
! 	if we are restarting from a given point, initialize the domain from the given restart file
	subroutine load_restart_file(domain,restart_file)
		type(domain_type), intent(inout) :: domain
		character(len=*),intent(in)::restart_file
		real,allocatable,dimension(:,:,:)::inputdata
		
		write(*,*) "WARNING Restart file may not be correct since setting U/V on staggered grids" 
		call io_read3d(restart_file,"u",inputdata)
		domain%u=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"v",inputdata)
		domain%v=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qv",inputdata)
		domain%qv=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qc",inputdata)
		domain%cloud=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qr",inputdata)
		domain%qrain=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qi",inputdata)
		domain%ice=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qs",inputdata)
		domain%qsnow=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qg",inputdata)
		domain%qgrau=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"nr",inputdata)
		domain%nrain=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"ni",inputdata)
		domain%nice=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"p",inputdata)
		domain%p=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"th",inputdata)
		domain%th=inputdata
		deallocate(inputdata)
	end subroutine load_restart_file
	
! 	initialize the boundary conditions (read inital conditions, etc.)
	subroutine bc_init(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		real,dimension(:,:,:),allocatable::inputdata
		logical :: boundary_value
		integer::nx,ny,nz,i
		real::domainsize
		! MODULE variables : curstep, curfile, nfiles, steps_in_file, file_list
		
! 		in case we are using a restart file we have some trickery to do here to find the proper file to be reading from
!       and set the current time step appropriately... should probably be moved to a subroutine. 
		curfile=1
		if (options%restart) then
			curstep=options%restart_step
		else
			curstep=1
		endif
		nfiles=options%nfiles
		allocate(file_list(nfiles))
		file_list=options%boundary_files
		call io_getdims(file_list(curfile),"P", dims)
		if (dims(1)==3) then
			steps_in_file=1
		else
			steps_in_file=dims(dims(1)+1) !dims(1) = ndims
		endif
		
		if (.not.options%ideal) then
			do while (curstep>steps_in_file)
				curfile=curfile+1
				if (curfile>nfiles) then
					stop "Ran out of files to process!"
				endif
				curstep=curstep-steps_in_file !instead of setting=1, this way we can set an arbitrary starting point multiple files in
				call io_getdims(file_list(curfile),"P", dims)
				if (dims(1)==3) then
					steps_in_file=1
				else
					steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
				endif
			enddo
		endif
! 		load the restart file
		if (options%restart) then
			call load_restart_file(domain,options%restart_file)
			if (options%external_winds) then
				call ext_winds_init(domain,bc,options)
			endif
			call update_winds(domain,options)
			call write_domain(domain,options,-1)
		else
! 			else load data from the first Boundary conditions file
			boundary_value=.False.
			nx=size(domain%u,1)
			ny=size(domain%u,3)
			if (options%external_winds) then
				call ext_winds_init(domain,bc,options)
! 				call smooth_wind(domain%u,1,3)
! 				call smooth_wind(domain%v,1,3)
			elseif (options%remove_lowres_linear) then
				call remove_linear_winds(domain,bc,options,file_list(curfile),curstep)
			elseif (options%mean_winds) then
				call mean_winds(domain,file_list(curfile),curstep)
			else
				call read_var(domain%u,    file_list(curfile),"U",      bc%u_geo%geolut,curstep,boundary_value)
				call read_var(domain%v,    file_list(curfile),"V",      bc%v_geo%geolut,curstep,boundary_value)
				if (.not.options%ideal)then
					call smooth_wind(domain%u,smoothing_window,3)
					call smooth_wind(domain%v,smoothing_window,3)
				endif
			endif
			call read_var(domain%p,    file_list(curfile),"P",      bc%geolut,curstep,boundary_value)
			call read_var(domain%th,   file_list(curfile),"T",      bc%geolut,curstep,boundary_value)
			call read_var(domain%qv,   file_list(curfile),"QVAPOR", bc%geolut,curstep,boundary_value)
			call read_var(domain%cloud,file_list(curfile),"QCLOUD", bc%geolut,curstep,boundary_value)
			call read_var(domain%ice,  file_list(curfile),"QICE",   bc%geolut,curstep,boundary_value)
		
			call update_pressure(domain%p,domain%th/((100000.0/domain%p)**(R/cp)), &
								 bc%next_domain%terrain,domain%terrain)
							 
	 		nz=size(domain%th,2)
			domainsize=size(domain%th,1)*size(domain%th,3)
	 		if (options%mean_fields) then
	 			do i=1,nz
	 				domain%th(:,i,:)=sum(domain%th(:,i,:))/domainsize
	 				domain%qv(:,i,:)=sum(domain%qv(:,i,:))/domainsize
	 				domain%cloud(:,i,:)=sum(domain%cloud(:,i,:))/domainsize
	 				domain%ice(:,i,:)=sum(domain%ice(:,i,:))/domainsize
	 			enddo
	 		endif
			
			call update_winds(domain,options)
		endif

	end subroutine bc_init


	subroutine update_edges(dxdt,d1,d2)
		implicit none
		real,dimension(:,:,:), intent(inout) :: dxdt
		real,dimension(:,:,:), intent(in) :: d1,d2
		integer :: nx,nz,ny,i

		nx=size(d1,1)
		nz=size(d1,2)
		ny=size(d1,3)
		do i=1,nz
			dxdt(i,:ny,1)=d1(1,i,:) -d2(1,i,:)
			dxdt(i,:ny,2)=d1(nx,i,:)-d2(nx,i,:)
			dxdt(i,:nx,3)=d1(:,i,1) -d2(:,i,1)
			dxdt(i,:nx,4)=d1(:,i,ny)-d2(:,i,ny)
		enddo
		dxdt(:,1,3:4)=0
		dxdt(:,nx,3:4)=0
	end subroutine update_edges
	
	
	subroutine update_dxdt(bc,domain)
		implicit none
		type(bc_type), intent(inout) :: bc
		type(domain_type), intent(in) :: domain
		
		bc%dudt=bc%next_domain%u-domain%u
		bc%dvdt=bc%next_domain%v-domain%v
		bc%dwdt=bc%next_domain%w-domain%w
		bc%dpdt=bc%next_domain%p-domain%p

		call update_edges(bc%dthdt,bc%next_domain%th,domain%th)
		call update_edges(bc%dqvdt,bc%next_domain%qv,domain%qv)
		call update_edges(bc%dqcdt,bc%next_domain%cloud,domain%cloud)
	end subroutine update_dxdt
	
	subroutine update_pressure(pressure,temperature,hgt_lo,hgt_hi)
		implicit none
		real,dimension(:,:,:), intent(inout) :: pressure
		real,dimension(:,:,:), intent(in) :: temperature
		real,dimension(:,:), intent(in) :: hgt_lo,hgt_hi
		real,dimension(:,:),allocatable::slp
		integer :: nx,ny,nz,i
		nx=size(pressure,1)
		nz=size(pressure,2)
		ny=size(pressure,3)
		allocate(slp(nx,ny))
		do i=1,nz
		    slp = pressure(:,i,:)/(1 - 2.25577E-5*hgt_lo)**5.25588
			pressure(:,i,:) = slp * (1 - 2.25577e-5 * hgt_hi)**5.25588
		enddo
		deallocate(slp)
	end subroutine update_pressure
	
	subroutine update_ext_winds(bc,options)
		implicit none
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		logical :: use_boundary,use_interior
		! MODULE variables : ext_winds_ curstep, curfile, nfiles, steps_in_file, file_list
		if (.not.options%ideal) then
			ext_winds_curstep=ext_winds_curstep+1
			if (ext_winds_curstep>ext_winds_steps_in_file) then
				ext_winds_curfile=ext_winds_curfile+1
				ext_winds_curstep=1
				call io_getdims(ext_winds_file_list(ext_winds_curfile),"U", dims)
				if (dims(1)==3) then
					ext_winds_steps_in_file=1
				else
					ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims
				endif
			endif
			if (ext_winds_curfile>ext_winds_nfiles) then
				stop "Ran out of files to process!"
			endif
		endif
		
		use_interior=.False.
		use_boundary=.True.
		call read_var(bc%next_domain%u,    ext_winds_file_list(ext_winds_curfile),"U", &
		              bc%ext_winds%u_geo%geolut,ext_winds_curstep,use_interior)
		call read_var(bc%next_domain%v,    ext_winds_file_list(ext_winds_curfile),"V", &
		              bc%ext_winds%v_geo%geolut,ext_winds_curstep,use_interior)
		call rotate_ext_wind_field(bc%next_domain,bc%ext_winds)
	
	end subroutine update_ext_winds
	
	
	subroutine bc_update(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		logical :: use_boundary,use_interior
		integer::i,nz,nx,ny
		! MODULE variables : curstep, curfile, nfiles, steps_in_file, file_list
		
		if (.not.options%ideal) then
			curstep=curstep+1
			do while (curstep>steps_in_file)
				curfile=curfile+1
				if (curfile>nfiles) then
					stop "Ran out of files to process!"
				endif
				curstep=curstep-steps_in_file !instead of setting=1, this way we can set an arbitrary starting point multiple files in
				call io_getdims(file_list(curfile),"P", dims)
				if (dims(1)==3) then
					steps_in_file=1
				else
					steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
				endif
			enddo
		endif
		use_interior=.False.
		use_boundary=.True.
		if (options%external_winds) then
			call update_ext_winds(bc,options)
! 			call smooth_wind(bc%next_domain%u,1,3)
! 			call smooth_wind(bc%next_domain%v,1,3)
		elseif (options%remove_lowres_linear) then
			call remove_linear_winds(bc%next_domain,bc,options,file_list(curfile),curstep)
		elseif (options%mean_winds) then
			call mean_winds(bc%next_domain,file_list(curfile),curstep)
		else
			call read_var(bc%next_domain%u,    file_list(curfile),"U",      bc%u_geo%geolut,curstep,use_interior)
			call read_var(bc%next_domain%v,    file_list(curfile),"V",      bc%v_geo%geolut,curstep,use_interior)
			if (.not.options%ideal)then
				call smooth_wind(bc%next_domain%u,smoothing_window,3)
				call smooth_wind(bc%next_domain%v,smoothing_window,3)
			endif
		endif
		call read_var(bc%next_domain%p,    file_list(curfile),"P",      bc%geolut,curstep,use_interior)
		call read_var(bc%next_domain%th,   file_list(curfile),"T",      bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%qv,   file_list(curfile),"QVAPOR", bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%cloud,file_list(curfile),"QCLOUD", bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%ice,  file_list(curfile),"QICE",   bc%geolut,curstep,use_boundary)
		
		call update_pressure(bc%next_domain%p,domain%th/((100000.0/domain%p)**(R/cp)), &
							 bc%next_domain%terrain,domain%terrain)
		nx=size(bc%next_domain%th,1)
		nz=size(bc%next_domain%th,2)
		ny=size(bc%next_domain%th,3)
		if (options%mean_fields) then
			do i=1,nz
				bc%next_domain%th(1,i,:)=sum(bc%next_domain%th(1,i,:))/ny
				bc%next_domain%th(nx,i,:)=sum(bc%next_domain%th(nx,i,:))/ny
				bc%next_domain%th(:,i,1)=sum(bc%next_domain%th(:,i,1))/nx
				bc%next_domain%th(:,i,ny)=sum(bc%next_domain%th(:,i,ny))/nx
				
				bc%next_domain%qv(1,i,:)=sum(bc%next_domain%qv(1,i,:))/ny
				bc%next_domain%qv(nx,i,:)=sum(bc%next_domain%qv(nx,i,:))/ny
				bc%next_domain%qv(:,i,1)=sum(bc%next_domain%qv(:,i,1))/nx
				bc%next_domain%qv(:,i,ny)=sum(bc%next_domain%qv(:,i,ny))/nx

				bc%next_domain%cloud(1,i,:)=sum(bc%next_domain%cloud(1,i,:))/ny
				bc%next_domain%cloud(nx,i,:)=sum(bc%next_domain%cloud(nx,i,:))/ny
				bc%next_domain%cloud(:,i,1)=sum(bc%next_domain%cloud(:,i,1))/nx
				bc%next_domain%cloud(:,i,ny)=sum(bc%next_domain%cloud(:,i,ny))/nx

				bc%next_domain%ice(1,i,:)=sum(bc%next_domain%ice(1,i,:))/ny
				bc%next_domain%ice(nx,i,:)=sum(bc%next_domain%ice(nx,i,:))/ny
				bc%next_domain%ice(:,i,1)=sum(bc%next_domain%ice(:,i,1))/nx
				bc%next_domain%ice(:,i,ny)=sum(bc%next_domain%ice(:,i,ny))/nx
			enddo
		endif
		
		call update_winds(bc%next_domain,options)
		call update_dxdt(bc,domain)
	end subroutine bc_update
end module boundary_conditions